#!/usr/bin/env python3
"""Inspect one NWB recording from the KK dataset (DANDI 000692, Kimura lab).

Mirrors scripts/inspect_ascent.py and scripts/inspect_ey.py section by section
so the three inspections are comparable.  Loading and annotation extraction are
adapted to this file's layout (ndx-multichannel-volume extension):

  acquisition/CalciumImageSeries/data       (T, X, Y, Z, C) uint16, gzip chunks; C = GCaMP, TagRFP-T
  processing/CalciumActivity/CalciumSeriesSegmentation/data
                                            (T, X, Y, Z) uint16 LABEL VOLUMES: voxel value i = ROI i, 0 = background
  processing/CalciumActivity/SegmentationVol1/Seg_tpoint_1
                                            one (x, y, z, weight) voxel per ROI for a single time point
  processing/CalciumActivity/NeuronIDs/labels   biological name per ROI ('unknown' allowed)
  processing/CalciumActivity/SignalRawFluor/... stored activity traces (T, N)
  acquisition/NeuroPALImageRaw/data         (X, Y, Z, C) structural stack with its own segmentation

Unlike EY there are no per-frame point tables: ROI positions are derived here
as the centroid of each label in each frame (the same reduction WormND's
3DeeCellTracker loader applies).  Deriving them streams the whole label series
once (~1.6k gzip-compressed volumes); the result can be cached with
--tracks-csv and reused by make_previews_kk.py.

Usage:
  python -P scripts/inspect_kk.py --nwb data/kk/<file>.nwb \
      --out-json previews/kk/inspection_summary.json \
      --tracks-csv previews/kk/calcium_tracks.csv \
      [--per-frame-means previews/kk/per_frame_means.csv]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))  # sibling scripts, not the data dir
from inspect_ascent import banner, jsonable, local_mean  # noqa: E402
from inspect_ey import decode, ragged_labels  # noqa: E402

CALCIUM = "acquisition/CalciumImageSeries"
CAL_VOL = "general/optophysiology/CalciumImVol"
SEG = "processing/CalciumActivity/CalciumSeriesSegmentation"
SEG_POINTS = "processing/CalciumActivity/SegmentationVol1"
LABELS = "processing/CalciumActivity/NeuronIDs/labels"
TRACE = "processing/CalciumActivity/SignalRawFluor/SignalCalciumImResponseSeries"
NEUROPAL_IMG = "acquisition/NeuroPALImageRaw"
NEUROPAL_VOL = "general/optophysiology/NeuroPALImVol"
NEUROPAL_SEG = "processing/NeuroPAL/NeuroPALSegmentation/NeuroPALNeurons"


# ----------------------------------------------------------------------------
# loading helpers (shared with make_previews_kk.py)
# ----------------------------------------------------------------------------
def channel_names(f: h5py.File) -> list[str]:
    """Optical channel names in data-axis order, from order_optical_channels if present."""
    g = f[CAL_VOL]
    if "order_optical_channels" in g and "channels" in g["order_optical_channels"]:
        return [decode(c) for c in g["order_optical_channels/channels"][()]]
    return [k for k in g.keys() if isinstance(g[k], h5py.Group) and "emission_lambda" in g[k]]


def read_frame(f: h5py.File, t: int, c: int) -> np.ndarray:
    """One channel of one calcium volume as (Z, Y, X), the axis order of the ASCENT helpers.

    Stored (T, X, Y, Z, C); reading one (t, c) touches 1/C of the frame's gzip chunks."""
    return f[f"{CALCIUM}/data"][t, :, :, :, c].transpose(2, 1, 0)


def read_frame_all(f: h5py.File, t: int) -> np.ndarray:
    """All channels of one volume as (C, Z, Y, X)."""
    return f[f"{CALCIUM}/data"][t].transpose(3, 2, 1, 0)


def read_seg(f: h5py.File, t: int) -> np.ndarray:
    """One label volume as (Z, Y, X) uint16; 0 = background, i = ROI i."""
    return f[f"{SEG}/data"][t].transpose(2, 1, 0)


def labels_list(f: h5py.File) -> list[str] | None:
    return [decode(v) for v in f[LABELS][()]] if LABELS in f else None


def seg_frame_table(seg: np.ndarray, t: int, n_labels: int) -> pd.DataFrame:
    """Centroid, voxel count and bounding box of every label present in one label volume.

    Columns: t, roi (0-based = label value - 1, the index into NeuronIDs/labels and the
    stored trace columns), label_value, x, y, z (centroid, float), n_voxels, x0..z1.
    Works on the foreground voxels only (about 1 % of the volume)."""
    Z, Y, X = seg.shape
    flat = seg.ravel()
    nz = np.flatnonzero(flat)
    lab = flat[nz].astype(np.int64)
    zz, yy, xx = np.unravel_index(nz, seg.shape)
    n = max(n_labels, int(lab.max()) if len(lab) else 0) + 1
    cnt = np.bincount(lab, minlength=n)
    present = np.nonzero(cnt[1:])[0] + 1
    rows = dict(t=t, roi=present - 1, label_value=present)
    for name, c in (("x", xx), ("y", yy), ("z", zz)):
        rows[name] = np.bincount(lab, weights=c, minlength=n)[present] / cnt[present]
    rows["n_voxels"] = cnt[present]
    for name, c, size in (("z", zz, Z), ("y", yy, Y), ("x", xx, X)):
        lo = np.full(n, size, dtype=np.int64); hi = np.full(n, -1, dtype=np.int64)
        np.minimum.at(lo, lab, c); np.maximum.at(hi, lab, c)
        rows[f"{name}0"] = lo[present]; rows[f"{name}1"] = hi[present]
    return pd.DataFrame(rows)


def load_calcium_tracks(f: h5py.File, cache: Path | None = None, verbose: bool = True) -> pd.DataFrame:
    """Per-frame ROI centroids from the label-volume series as one long table.

    Streams every label volume once (gzip chunks, ~5 MB each decompressed).  The
    label value is used as the track id, as the authors' description implies
    ("each value i indicates that voxel belongs to ROI i")."""
    if cache is not None and cache.exists():
        df = pd.read_csv(cache)
        if verbose:
            print(f"loaded cached tracks from {cache}: {len(df)} rows, {df.t.nunique()} frames")
        return df
    labels = labels_list(f)
    n_lab = len(labels) if labels else 0
    T = f[f"{SEG}/data"].shape[0]
    parts = []
    t0 = time.time()
    for t in range(T):
        parts.append(seg_frame_table(read_seg(f, t), t, n_lab))
        if verbose and t % 200 == 0:
            print(f"  label volume t={t}: {len(parts[-1])} ROIs ({time.time() - t0:.0f} s)")
    df = pd.concat(parts, ignore_index=True)
    if labels is not None:
        lab = np.array(labels + [""] * (int(df.roi.max()) + 1 - n_lab), dtype=object)
        df["label"] = lab[df.roi.to_numpy()]
    if verbose:
        print(f"  derived centroids for {T} label volumes in {time.time() - t0:.0f} s: {len(df)} rows")
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True); df.to_csv(cache, index=False); print(f"wrote {cache}")
    return df


def rounded(pts: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (pts.z.round().astype(int).to_numpy(), pts.y.round().astype(int).to_numpy(), pts.x.round().astype(int).to_numpy())


def mask_means(vol: np.ndarray, seg: np.ndarray, n: int) -> np.ndarray:
    """Mean of `vol` over each label 1..n of `seg` (NaN where the label is absent)."""
    lab = seg.ravel().astype(np.int64)
    cnt = np.bincount(lab, minlength=n + 1)[1:n + 1].astype(float)
    s = np.bincount(lab, weights=vol.ravel().astype(np.float64), minlength=n + 1)[1:n + 1]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, s / cnt, np.nan)


# ----------------------------------------------------------------------------
# 1. structure
# ----------------------------------------------------------------------------
def inspect_structure(f: h5py.File) -> dict:
    banner("1. NWB / HDF5 STRUCTURE (verified from file)")
    out = dict(root_keys=list(f.keys()), nwb_version=decode(f.attrs.get("nwb_version", "")))
    print(f"root keys: {out['root_keys']}; nwb_version {out['nwb_version']}")
    out["spec_namespaces"] = {k: list(f["specifications"][k].keys()) for k in f["specifications"].keys()} if "specifications" in f else {}
    print(f"specification namespaces embedded: {out['spec_namespaces']}")

    def dsinfo(d, extra=None):
        layout = d.id.get_create_plist().get_layout()
        r = dict(shape=list(d.shape), dtype=str(d.dtype), chunks=d.chunks, compression=d.compression,
                 compression_opts=d.compression_opts, layout={0: "compact", 1: "contiguous", 2: "chunked"}.get(layout, str(layout)),
                 nbytes_uncompressed=int(d.nbytes), storage_size=int(d.id.get_storage_size()),
                 data_attrs={k: jsonable(v) for k, v in d.attrs.items()})
        r.update(extra or {}); return r

    d = f[f"{CALCIUM}/data"]
    out["calcium"] = dsinfo(d, dict(bytes_per_frame=int(d.nbytes // d.shape[0]), dimension_dataset=f[f"{CALCIUM}/dimension"][()].tolist(),
                                    rate_hz=float(f[f"{CALCIUM}/starting_time"].attrs["rate"]), starting_time_s=float(f[f"{CALCIUM}/starting_time"][()]),
                                    attrs={k: jsonable(v) for k, v in f[CALCIUM].attrs.items() if k != "object_id"}, channels=channel_names(f)))
    c = out["calcium"]
    print(f"CalciumImageSeries/data: shape {c['shape']} dtype {c['dtype']} layout {c['layout']} chunks {c['chunks']} "
          f"compression {c['compression']} level {c['compression_opts']}; {c['nbytes_uncompressed']/1e9:.2f} GB uncompressed, "
          f"{c['storage_size']/1e9:.2f} GB on disk (ratio {c['nbytes_uncompressed']/c['storage_size']:.2f}); {c['bytes_per_frame']/1e6:.1f} MB/frame uncompressed")
    print(f"  dimension dataset = {c['dimension_dataset']} (per-frame axes X, Y, Z); channel axis last: {c['channels']}")
    print(f"  starting_time {c['starting_time_s']} s, rate {c['rate_hz']} Hz -> {d.shape[0]} frames = {d.shape[0]/c['rate_hz']:.1f} s; timestamps dataset present: {'timestamps' in f[CALCIUM]}")
    print(f"  attrs: {c['attrs']}\n  data attrs: {c['data_attrs']}")
    n_frames = d.shape[0]
    out["n_frames"] = n_frames
    out["axis_order"] = "T,X,Y,Z,C as documented by `dimension` and the channel groups; tested in sections 5-6"

    ni = f[f"{NEUROPAL_IMG}/data"]
    out["neuropal_image"] = dsinfo(ni, dict(RGBW_channels=f[f"{NEUROPAL_IMG}/RGBW_channels"][()].tolist(), description=decode(f[f"{NEUROPAL_IMG}/description"][()])))
    print(f"NeuroPALImageRaw/data: shape {ni.shape} dtype {ni.dtype} chunks {ni.chunks} {ni.compression} ({ni.nbytes/1e6:.0f} MB); "
          f"RGBW_channels {out['neuropal_image']['RGBW_channels']}; '{out['neuropal_image']['description']}'")

    s = f[f"{SEG}/data"]
    out["segmentation"] = dsinfo(s, dict(rate_hz=float(f[f"{SEG}/starting_time"].attrs["rate"]) if "starting_time" in f[SEG] else None,
                                         attrs={k: jsonable(v) for k, v in f[SEG].attrs.items() if k != "object_id"},
                                         members=list(f[SEG].keys())))
    sg = out["segmentation"]
    print(f"CalciumSeriesSegmentation/data: shape {sg['shape']} dtype {sg['dtype']} chunks {sg['chunks']} {sg['compression']}; "
          f"{sg['storage_size']/1e9:.2f} GB on disk; rate {sg['rate_hz']} Hz; members {sg['members']}")
    print(f"  description: {sg['attrs'].get('description')}")
    if SEG_POINTS in f:
        g = f[SEG_POINTS]
        keys = list(g.keys()); p = g[keys[0]]
        out["segmentation_points"] = dict(groups=keys, n_rois=int(p["id"].shape[0]), voxel_mask_dtype=str(p["voxel_mask"].dtype),
                                          description=decode(p.attrs.get("description", "")),
                                          voxel_mask_index_is_1_to_n=bool(np.array_equal(p["voxel_mask_index"][()], np.arange(1, p["id"].shape[0] + 1))))
        sp = out["segmentation_points"]
        print(f"SegmentationVol1: groups {keys}; {sp['n_rois']} ROIs, one voxel each: {sp['voxel_mask_index_is_1_to_n']}; '{sp['description']}'")
    nid = f["processing/CalciumActivity/NeuronIDs"]
    out["neuron_ids"] = dict(members=list(nid.keys()), n_labels=int(f[LABELS].shape[0]), description=decode(nid["description"][()]) if "description" in nid else None,
                             attrs={k: jsonable(v) for k, v in nid.attrs.items() if k != "object_id"})
    if "MCVSeriesSegmentation" in nid:
        out["neuron_ids"]["MCVSeriesSegmentation_links_to"] = nid["MCVSeriesSegmentation"].name if isinstance(nid["MCVSeriesSegmentation"], h5py.Group) else None
        try:
            out["neuron_ids"]["MCVSeriesSegmentation_same_object_as_seg"] = bool(nid["MCVSeriesSegmentation"].attrs.get("object_id") == f[SEG].attrs.get("object_id"))
        except Exception:
            pass
    print(f"NeuronIDs: {out['neuron_ids']}")
    td = f[f"{TRACE}/data"]
    out["SignalRawFluor"] = dict(shape=list(td.shape), dtype=str(td.dtype), rate_hz=float(f[f"{TRACE}/starting_time"].attrs["rate"]),
                                 starting_time_s=float(f[f"{TRACE}/starting_time"][()]), unit=decode(td.attrs.get("unit", "")),
                                 description=decode(f[TRACE].attrs.get("description", "")),
                                 rois_table=f[f[f"{TRACE}/rois"].attrs["table"]].name if "rois" in f[TRACE] and "table" in f[f"{TRACE}/rois"].attrs else None)
    tr = out["SignalRawFluor"]
    print(f"SignalRawFluor: shape {tr['shape']} {tr['dtype']} unit '{tr['unit']}' rate {tr['rate_hz']} Hz; '{tr['description']}'; rois -> {tr['rois_table']}")

    # pixel-value spot checks on three frames, both channels
    spot = {}
    t0 = time.time()
    for t in sorted({0, n_frames // 2, n_frames - 1}):
        v = read_frame_all(f, t)
        for ci, cn in enumerate(c["channels"]):
            a = v[ci]
            spot[f"t{t}_{cn}"] = dict(min=int(a.min()), max=int(a.max()), mean=float(a.mean()), frac_zero=float((a == 0).mean()),
                                      p50=float(np.percentile(a, 50)), p99=float(np.percentile(a, 99)), frac_saturated=float((a == 65535).mean()))
            # fully zero planes along each axis (a zero strip would indicate padding after registration)
            zero = (a == 0)
            sp = spot[f"t{t}_{cn}"]
            sp["all_zero_z_slices"] = int(zero.all(axis=(1, 2)).sum()); sp["all_zero_y_columns"] = int(zero.all(axis=(0, 2)).sum()); sp["all_zero_x_rows"] = int(zero.all(axis=(0, 1)).sum())
            sp["all_zero_y_columns_which"] = np.flatnonzero(zero.all(axis=(0, 2))).tolist()[:12]
            print(f"  t{t} {cn}: min={sp['min']} max={sp['max']} mean={sp['mean']:.1f} p50={sp['p50']:.0f} p99={sp['p99']:.0f} zero-fraction={sp['frac_zero']:.4f}; "
                  f"all-zero planes: z {sp['all_zero_z_slices']}, y columns {sp['all_zero_y_columns']} {sp['all_zero_y_columns_which']}, x rows {sp['all_zero_x_rows']}")
    out["spot_checks"] = spot
    out["seconds_to_read_3_frames"] = round(time.time() - t0, 1)
    print(f"  (reading 3 full frames took {out['seconds_to_read_3_frames']} s)")
    if all(s["max"] <= 255 for s in spot.values()):
        print("  NOTE: every sampled voxel is <= 255 although the dtype is uint16 -> 8-bit values inside a 16-bit container (section 7 gives the global maximum)")
    return out


# ----------------------------------------------------------------------------
# 2. metadata stored in the file
# ----------------------------------------------------------------------------
def inspect_metadata(f: h5py.File) -> dict:
    banner("2. METADATA STORED IN THE FILE (NWB fields)")
    out = {}

    def scalar(path):
        if path not in f:
            return None
        v = f[path][()]
        if isinstance(v, np.ndarray):
            return [decode(x) if isinstance(x, bytes) else jsonable(x) for x in v]
        return decode(v) if isinstance(v, bytes) else jsonable(v)

    for k in ["identifier", "session_description", "session_start_time", "timestamps_reference_time", "file_create_date"]:
        out[k] = scalar(k)
    for k in f["general"].keys():
        if isinstance(f["general"][k], h5py.Dataset):
            out[f"general/{k}"] = scalar(f"general/{k}")
    out["subject"] = {k: scalar(f"general/subject/{k}") for k in f["general/subject"].keys()}
    out["devices"] = {k: {a: jsonable(b) for a, b in f["general/devices"][k].attrs.items() if a != "object_id"} for k in f["general/devices"].keys()}
    vols = {}
    for name in f["general/optophysiology"].keys():
        g = f["general/optophysiology"][name]
        v = {k: scalar(f"general/optophysiology/{name}/{k}") for k in g.keys() if isinstance(g[k], h5py.Dataset)}
        v["grid_spacing_unit"] = decode(g["grid_spacing"].attrs.get("unit", "")) if "grid_spacing" in g else None
        v["channels"] = {}
        for k in g.keys():
            if isinstance(g[k], h5py.Group) and "emission_lambda" in g[k]:
                v["channels"][k] = {kk: scalar(f"general/optophysiology/{name}/{k}/{kk}") for kk in g[k].keys()}
        if "order_optical_channels" in g:
            v["order_optical_channels"] = scalar(f"general/optophysiology/{name}/order_optical_channels/channels")
        vols[name] = v
    out["imaging_volumes"] = vols
    out["stimuli"] = []  # no intervals/chemical_stimuli in this dataset
    out["has_intervals"] = "intervals" in f
    out["stimulus_presentation_members"] = list(f["stimulus/presentation"].keys()) if "stimulus/presentation" in f else None
    out["neuron_ids_description"] = decode(f["processing/CalciumActivity/NeuronIDs/description"][()]) if "processing/CalciumActivity/NeuronIDs/description" in f else None
    out["neuropal_seg_description"] = decode(f[NEUROPAL_SEG].attrs.get("description", "")) if NEUROPAL_SEG in f else None

    for k in ["identifier", "session_description", "session_start_time", "file_create_date", "general/experiment_description",
              "general/experimenter", "general/institution", "general/lab", "general/related_publications"]:
        print(f"  {k}: {out.get(k)}")
    print(f"  subject: {out['subject']}")
    print(f"  devices: {out['devices']}")
    for name, v in vols.items():
        print(f"  imaging volume {name}: grid_spacing {v.get('grid_spacing')} {v.get('grid_spacing_unit')}, location {v.get('location')}, "
              f"channels {list(v['channels'].keys())}, order {v.get('order_optical_channels')}, description '{v.get('description')}'")
        for ck, cv in v["channels"].items():
            print(f"      {ck}: ex {cv.get('excitation_lambda')} nm, em {cv.get('emission_lambda')} nm, '{cv.get('description')}'")
    print(f"  intervals group present: {out['has_intervals']}; stimulus/presentation members: {out['stimulus_presentation_members']} (no stimulus in this dataset)")
    print(f"  NeuronIDs description: {out['neuron_ids_description']}")
    print(f"  NeuroPAL segmentation description: {out['neuropal_seg_description']}")
    return out


def validate_with_pynwb(path: Path) -> dict:
    banner("2b. PYNWB READ CHECK")
    try:
        from pynwb import NWBHDF5IO
    except ImportError as e:
        print(f"pynwb not installed ({e}); skipped")
        return dict(ok=False, error=str(e))
    t0 = time.time()
    try:
        with NWBHDF5IO(str(path), mode="r", load_namespaces=True) as io:
            nwb = io.read()
            seg = nwb.processing["CalciumActivity"]["CalciumSeriesSegmentation"]
            info = dict(ok=True, identifier=nwb.identifier, session_description=nwb.session_description,
                        acquisition=list(nwb.acquisition.keys()), processing=list(nwb.processing.keys()),
                        imaging_planes=list(nwb.imaging_planes.keys()), seg_type=type(seg).__name__,
                        seg_has_plane_segmentations=hasattr(seg, "plane_segmentations"), seg_has_data=hasattr(seg, "data"),
                        seconds=round(time.time() - t0, 1))
        print(f"pynwb read OK in {info['seconds']} s: acquisition {info['acquisition']}, processing {info['processing']}; "
              f"CalciumSeriesSegmentation is a {info['seg_type']} (plane_segmentations: {info['seg_has_plane_segmentations']}, data: {info['seg_has_data']})")
        return info
    except Exception as e:
        print(f"pynwb read FAILED: {type(e).__name__}: {e}")
        return dict(ok=False, error=f"{type(e).__name__}: {e}")


# ----------------------------------------------------------------------------
# 3. channel behaviour (two channels: the ASCENT test applies)
# ----------------------------------------------------------------------------
def inspect_channels(f: h5py.File, tracks: pd.DataFrame, n_frames: int, n_sample: int, n_neurons: int, seed: int) -> dict:
    banner("3. CHANNEL BEHAVIOUR: which channel is calcium (GCaMP) vs reference (TagRFP-T)?")
    chans = channel_names(f)
    rng = np.random.default_rng(seed)
    ts = np.unique(np.linspace(0, n_frames - 1, n_sample).round().astype(int))
    n_full = tracks.groupby("roi")["t"].nunique()
    full = n_full[n_full == tracks.t.nunique()].index.to_numpy()
    rois = np.sort(rng.choice(full, size=min(n_neurons, len(full)), replace=False))
    sel = tracks[tracks.roi.isin(rois)].set_index(["t", "roi"]).sort_index()
    per_chan = {}
    neuron_ts = {c: np.zeros((len(ts), len(rois))) for c in chans}
    mask_ts = {c: np.zeros((len(ts), len(rois))) for c in chans}
    stats = {c: dict(vol_mean=[], p50=[], p99=[], vmax=[]) for c in chans}
    n_lab = int(tracks.roi.max()) + 1
    for i, t in enumerate(ts):
        v = read_frame_all(f, int(t)); seg = read_seg(f, int(t))
        pts = sel.loc[int(t)].loc[rois]
        z, y, x = rounded(pts)
        for ci, cn in enumerate(chans):
            a = v[ci]
            stats[cn]["vol_mean"].append(a.mean()); stats[cn]["p50"].append(np.percentile(a, 50)); stats[cn]["p99"].append(np.percentile(a, 99)); stats[cn]["vmax"].append(a.max())
            neuron_ts[cn][i] = local_mean(a, z, y, x, r=2)
            mask_ts[cn][i] = mask_means(a, seg, n_lab)[rois]
    for cn in chans:
        vol_mean = np.array(stats[cn]["vol_mean"])
        cv = []
        for j in range(len(rois)):
            yv = neuron_ts[cn][:, j]
            trend = np.polyval(np.polyfit(ts, yv, 1), ts)
            cv.append(np.std(yv - trend) / max(yv.mean(), 1e-9))
        cv = np.array(cv)
        per_chan[cn] = dict(volume_mean=vol_mean.tolist(), p50=np.array(stats[cn]["p50"]).tolist(), p99=np.array(stats[cn]["p99"]).tolist(),
                            max=np.array(stats[cn]["vmax"]).tolist(), volume_mean_temporal_cv=float(vol_mean.std() / vol_mean.mean()),
                            neuron_detrended_cv_median=float(np.median(cv)), neuron_detrended_cv_p90=float(np.percentile(cv, 90)),
                            neuron_mean_intensity_median=float(np.median(neuron_ts[cn])),
                            mask_mean_intensity_median=float(np.nanmedian(mask_ts[cn])))
        p = per_chan[cn]
        print(f"{cn}: volume mean over {len(ts)} sampled frames: {vol_mean.min():.1f}..{vol_mean.max():.1f} (temporal CV {p['volume_mean_temporal_cv']:.3f}); "
              f"p99 {np.min(p['p99']):.0f}..{np.max(p['p99']):.0f}; max {np.min(p['max'])}..{np.max(p['max'])}")
        print(f"{cn}: per-neuron detrended temporal CV over {len(rois)} full-length tracks (5x5 at rounded centroid): median {p['neuron_detrended_cv_median']:.3f}, "
              f"p90 {p['neuron_detrended_cv_p90']:.3f}; median neuron brightness {p['neuron_mean_intensity_median']:.0f} (mask mean {p['mask_mean_intensity_median']:.0f})")
    result = dict(sampled_frames=ts.tolist(), rois=rois.tolist(), per_channel=per_chan)
    if len(chans) == 2:
        a, b = chans
        cva, cvb = per_chan[a]["neuron_detrended_cv_median"], per_chan[b]["neuron_detrended_cv_median"]
        calcium = a if cva > cvb else b; ref = b if calcium == a else a
        corr = float(np.corrcoef(neuron_ts[a].mean(0), neuron_ts[b].mean(0))[0, 1])
        frac = float(np.mean(neuron_ts[a].std(0) / neuron_ts[a].mean(0) > neuron_ts[b].std(0) / neuron_ts[b].mean(0)))
        result.update(calcium_channel_guess=calcium, reference_channel_guess=ref, cv_ratio=float(max(cva, cvb) / max(min(cva, cvb), 1e-9)),
                      cross_channel_brightness_corr_over_neurons=corr, frac_neurons_with_more_variability_in=dict(channel=a, fraction=frac))
        print(f"-> more temporally variable channel (calcium-like): {calcium}; more stable channel (reference-like): {ref}; CV ratio {result['cv_ratio']:.2f}")
        print(f"-> correlation of mean neuron brightness across channels (same spots bright in both?): {corr:.2f}")
        print(f"-> fraction of sampled neurons more variable in {a} than {b}: {frac:.2f}")
        print(f"   (documented: channel 0 = {a}, channel 1 = {b}; description says tracking was done on the RFP series)")
    return result


# ----------------------------------------------------------------------------
# 4. annotations
# ----------------------------------------------------------------------------
def inspect_annotations(f: h5py.File, tracks: pd.DataFrame, n_frames: int) -> dict:
    banner("4. ANNOTATIONS: label volumes -> per-frame ROI centroids, track structure, labels")
    out = {}
    per_t = tracks.groupby("t").size()
    n_t = tracks.t.nunique()
    out["positions"] = dict(
        source="centroids of the label values in CalciumSeriesSegmentation/data (derived here, not stored as points)",
        columns=list(tracks.columns), n_rows=len(tracks), t_min=int(tracks.t.min()), t_max=int(tracks.t.max()), n_time_points=int(n_t),
        unannotated_frames=sorted(set(range(n_frames)) - set(tracks.t.unique())),
        rois_per_frame=dict(min=int(per_t.min()), median=float(per_t.median()), max=int(per_t.max())),
        ranges={c: [float(tracks[c].min()), float(tracks[c].max())] for c in "xyz"},
        fraction_of_movie_frames_annotated=float(n_t / n_frames),
    )
    p = out["positions"]
    print(f"{len(tracks)} ROI rows over {n_t} of {n_frames} frames (t {p['t_min']}..{p['t_max']}); frames with no label at all: {p['unannotated_frames'][:10]}{'...' if len(p['unannotated_frames']) > 10 else ''}")
    print(f"  ROIs per frame: {p['rois_per_frame']}; centroid ranges x{np.round(p['ranges']['x'], 1).tolist()} y{np.round(p['ranges']['y'], 1).tolist()} z{np.round(p['ranges']['z'], 1).tolist()}")

    # ROI size statistics (this is what masks add over point annotations)
    nv = tracks.n_voxels
    ext = pd.DataFrame({"dx": tracks.x1 - tracks.x0 + 1, "dy": tracks.y1 - tracks.y0 + 1, "dz": tracks.z1 - tracks.z0 + 1})
    out["roi_size"] = dict(voxels=dict(min=int(nv.min()), p10=float(nv.quantile(0.1)), median=float(nv.median()), p90=float(nv.quantile(0.9)), max=int(nv.max())),
                           bbox_extent_median=ext.median().to_dict(), bbox_extent_max=ext.max().to_dict(),
                           frac_rows_single_voxel=float((nv == 1).mean()), frac_rows_single_slice=float((ext.dz == 1).mean()))
    rs = out["roi_size"]
    print(f"  ROI size (voxels per label per frame): min {rs['voxels']['min']}, p10 {rs['voxels']['p10']:.0f}, median {rs['voxels']['median']:.0f}, "
          f"p90 {rs['voxels']['p90']:.0f}, max {rs['voxels']['max']}; bbox extent median {rs['bbox_extent_median']}, max {rs['bbox_extent_max']}; "
          f"{100*rs['frac_rows_single_voxel']:.2f}% single-voxel rows, {100*rs['frac_rows_single_slice']:.1f}% confined to one slice")

    # label values vs the NeuronIDs table
    n_labels = int(f[LABELS].shape[0])
    out["label_values"] = dict(min=int(tracks.label_value.min()), max=int(tracks.label_value.max()), n_distinct=int(tracks.label_value.nunique()),
                               n_neuron_id_labels=n_labels, max_equals_n_labels=bool(tracks.label_value.max() == n_labels))
    lv = out["label_values"]
    print(f"  label values {lv['min']}..{lv['max']} ({lv['n_distinct']} distinct) vs {n_labels} NeuronIDs labels -> max label value == number of labels: {lv['max_equals_n_labels']}")

    # track structure: label value as track id
    tl = tracks.groupby("roi")["t"].nunique()
    gaps = []
    for roi, g in tracks.groupby("roi"):
        tt = np.sort(g.t.to_numpy()); gaps.append(int(((np.diff(tt) - 1) > 0).sum()))
    out["tracks"] = dict(n_tracks=int(tracks.roi.nunique()), track_length=dict(min=int(tl.min()), median=float(tl.median()), max=int(tl.max())),
                         n_tracks_spanning_all_frames=int((tl == n_frames).sum()), n_tracks_present_in_all_annotated_frames=int((tl == n_t).sum()),
                         n_missing_roi_frames=int(tracks.roi.nunique() * n_frames - len(tracks)),
                         frac_missing=float(1 - len(tracks) / (tracks.roi.nunique() * n_frames)),
                         n_tracks_with_gaps=int(sum(1 for g in gaps if g > 0)),
                         shortest_tracks=[(int(r), tracks.loc[tracks.roi == r, "label"].iloc[0] if "label" in tracks else "", int(n)) for r, n in tl.sort_values().head(8).items()])
    tr = out["tracks"]
    print(f"  tracks (label value): {tr['n_tracks']}; length min/median/max {tr['track_length']['min']}/{tr['track_length']['median']}/{tr['track_length']['max']}; "
          f"{tr['n_tracks_spanning_all_frames']} present in all {n_frames} frames; {tr['n_missing_roi_frames']} (roi, frame) pairs missing "
          f"({100*tr['frac_missing']:.2f}%); {tr['n_tracks_with_gaps']} tracks have gaps; shortest (roi, label, n frames): {tr['shortest_tracks']}")

    # motion
    ts_ = tracks.sort_values(["roi", "t"])
    g = ts_.groupby("roi")
    step = g[["z", "y", "x"]].diff(); dt = g["t"].diff()
    step = step[dt == 1]  # consecutive frames only
    disp = np.sqrt((step[["y", "x"]] ** 2).sum(axis=1))
    out["per_step_motion_voxels"] = dict(n_steps=int(len(disp)), xy_median=float(disp.median()), xy_p90=float(disp.quantile(0.9)), xy_p99=float(disp.quantile(0.99)),
                                         xy_max=float(disp.max()), z_median=float(step.z.abs().median()), z_p99=float(step.z.abs().quantile(0.99)), z_max=float(step.z.abs().max()),
                                         frac_steps_below_half_voxel=float((disp < 0.5).mean()))
    m = out["per_step_motion_voxels"]
    print(f"  motion between consecutive frames (centroids, voxels): in-plane median {m['xy_median']:.2f}, p90 {m['xy_p90']:.2f}, p99 {m['xy_p99']:.1f}, max {m['xy_max']:.1f}; "
          f"|dz| median {m['z_median']:.2f}, p99 {m['z_p99']:.2f}, max {m['z_max']:.1f}; {100*m['frac_steps_below_half_voxel']:.0f}% of steps move < 0.5 voxel in-plane")
    # per-frame population motion: median displacement of all ROIs between t-1 and t (shows when the worm moves)
    fm = pd.DataFrame({"t": ts_.t[dt == 1], "d": disp}).groupby("t")["d"].median()
    out["population_motion_per_frame"] = dict(median_over_frames=float(fm.median()), p99_over_frames=float(fm.quantile(0.99)), max=float(fm.max()),
                                              n_frames_median_step_above_2vox=int((fm > 2).sum()), frames_with_largest_step=fm.sort_values(ascending=False).head(5).index.astype(int).tolist())
    pm = out["population_motion_per_frame"]
    print(f"  population median step per frame: median {pm['median_over_frames']:.2f}, p99 {pm['p99_over_frames']:.2f}, max {pm['max']:.1f} voxels; "
          f"{pm['n_frames_median_step_above_2vox']} frames where the median ROI moves > 2 voxels (largest at t={pm['frames_with_largest_step']})")
    first = g.first(); last = g.last()
    drift = np.sqrt(((last[["x", "y"]] - first[["x", "y"]]) ** 2).sum(axis=1))
    out["first_to_last_drift_voxels"] = dict(median=float(drift.median()), max=float(drift.max()))
    # range of motion of each track over the whole recording (bounding box of its centroid path)
    rng_ = g[["x", "y", "z"]].agg(lambda s: s.max() - s.min())
    out["track_path_extent_voxels"] = dict(xy_median=float(np.sqrt(rng_.x ** 2 + rng_.y ** 2).median()), xy_max=float(np.sqrt(rng_.x ** 2 + rng_.y ** 2).max()),
                                           z_median=float(rng_.z.median()), z_max=float(rng_.z.max()))
    te = out["track_path_extent_voxels"]
    print(f"  in-plane displacement first->last frame per track: median {drift.median():.1f}, max {drift.max():.1f} voxels; "
          f"extent of each track's path over the recording: xy median {te['xy_median']:.1f}, max {te['xy_max']:.1f}; z median {te['z_median']:.1f}, max {te['z_max']:.1f}")

    # z distribution of centroids
    zc = tracks.z.round().astype(int).value_counts().sort_index()
    out["z_histogram"] = {int(k): int(v) for k, v in zc.items()}
    print(f"  rows per (rounded) z slice: {out['z_histogram']}")

    # label volume sanity: overlapping labels are impossible by construction (one value per voxel); connectivity not checked
    # labels
    if "label" in tracks:
        labels = pd.Series(labels_list(f))
        nonempty = labels[(labels.str.strip() != "")]
        named = nonempty[~nonempty.str.lower().isin(["unknown", "none", "nan"])]
        lr = sum(1 for l in named if l.endswith("L") and (l[:-1] + "R") in set(named))
        out["labels"] = dict(n=int(len(labels)), n_nonempty=int(len(nonempty)), n_named=int(len(named)), n_unique_named=int(named.nunique()),
                             n_unknown=int((nonempty.str.lower() == "unknown").sum()), n_LR_pairs=int(lr), examples=labels.head(12).tolist(),
                             duplicated_names=sorted(named[named.duplicated()].unique().tolist()))
        lb = out["labels"]
        print(f"  biological labels (NeuronIDs/labels): {lb['n']} total, {lb['n_named']} real names ({lb['n_unique_named']} unique, duplicates {lb['duplicated_names']}), "
              f"{lb['n_unknown']} 'unknown', {lb['n_LR_pairs']} left/right pairs; e.g. {lb['examples']}")

    # the single stored point table (SegmentationVol1/Seg_tpoint_1): which frame does it match, and does ROI i sit inside label i+1?
    if SEG_POINTS in f:
        g1 = f[SEG_POINTS]; key = list(g1.keys())[0]; p1 = g1[key]
        vm = p1["voxel_mask"][()]
        pts = pd.DataFrame({"roi": p1["id"][()], "x": vm["x"].astype(int), "y": vm["y"].astype(int), "z": vm["z"].astype(int), "w": vm["weight"]})
        res = {}
        for t in (0, 1, 2):
            seg = read_seg(f, t)
            inside = (pts.x < seg.shape[2]) & (pts.y < seg.shape[1]) & (pts.z < seg.shape[0])
            val = np.full(len(pts), -1); val[inside] = seg[pts.z[inside], pts.y[inside], pts.x[inside]]
            c = tracks[tracks.t == t].set_index("roi")
            d = np.sqrt(((c.loc[pts.roi[pts.roi.isin(c.index)], ["x", "y", "z"]].to_numpy() - pts.set_index("roi").loc[pts.roi[pts.roi.isin(c.index)], ["x", "y", "z"]].to_numpy()) ** 2).sum(1))
            res[f"t{t}"] = dict(frac_point_inside_label_roi_plus_1=float((val == pts.roi + 1).mean()), frac_point_on_background=float((val == 0).mean()),
                                frac_point_inside_other_label=float(((val > 0) & (val != pts.roi + 1)).mean()),
                                centroid_distance_median=float(np.median(d)), centroid_distance_max=float(np.max(d)))
        out["segmentation_points_vs_label_volumes"] = dict(group=key, n_points=len(pts), weights_unique=np.unique(pts.w).tolist(),
                                                           ranges={c: [int(pts[c].min()), int(pts[c].max())] for c in "xyz"}, per_frame=res)
        print(f"  {SEG_POINTS}/{key}: {len(pts)} points, weights {np.unique(pts.w).tolist()}, ranges {out['segmentation_points_vs_label_volumes']['ranges']}")
        for t, r in res.items():
            print(f"      vs label volume {t}: point lies inside label (roi+1) for {100*r['frac_point_inside_label_roi_plus_1']:.1f}% of ROIs, on background {100*r['frac_point_on_background']:.1f}%, "
                  f"inside another label {100*r['frac_point_inside_other_label']:.1f}%; distance to our centroid median {r['centroid_distance_median']:.2f}, max {r['centroid_distance_max']:.1f} voxels")

    # NeuroPAL segmentation and its relation to calcium frame 0 (matched by name)
    if NEUROPAL_SEG in f:
        gg = f[NEUROPAL_SEG]
        vm = gg["voxel_mask"][()]
        npal = pd.DataFrame({"x": vm["x"].astype(int), "y": vm["y"].astype(int), "z": vm["z"].astype(int), "label": ragged_labels(gg)})
        ni_shape = list(f[f"{NEUROPAL_IMG}/data"].shape)
        cal_labels = set(tracks.label.unique()) if "label" in tracks else set()
        named_npal = npal[~npal.label.str.lower().isin(["", "unknown"])]
        out["neuropal_segmentation"] = dict(n_rois=len(npal), ranges={c: [int(npal[c].min()), int(npal[c].max())] for c in "xyz"}, image_shape=ni_shape,
                                            n_named=int(len(named_npal)), n_unique_named=int(named_npal.label.nunique()),
                                            n_names_shared_with_calcium=int(len(set(named_npal.label) & cal_labels)),
                                            n_calcium_names_not_in_neuropal=int(len({l for l in cal_labels if l.lower() not in ("", "unknown")} - set(named_npal.label))),
                                            voxels_per_roi_all_one=bool(np.array_equal(gg["voxel_mask_index"][()], np.arange(1, len(npal) + 1))))
        ns = out["neuropal_segmentation"]
        print(f"  NeuroPAL segmentation: {ns['n_rois']} ROIs (one voxel each: {ns['voxels_per_roi_all_one']}) in NeuroPAL image coords (image shape {ns['image_shape']}), ranges {ns['ranges']}; "
              f"{ns['n_named']} named ({ns['n_unique_named']} unique); names shared with the calcium labels: {ns['n_names_shared_with_calcium']}; "
              f"calcium names absent from NeuroPAL: {ns['n_calcium_names_not_in_neuropal']}")
        if "label" in tracks:
            t0 = tracks[tracks.t == tracks.t.min()]
            t0 = t0[~t0.label.str.lower().isin(["", "unknown"])].drop_duplicates("label").set_index("label")
            both = named_npal.drop_duplicates("label").set_index("label")
            common = both.index.intersection(t0.index)
            if len(common) >= 5:
                A = t0.loc[common, ["x", "y", "z"]].to_numpy(float); B = both.loc[common, ["x", "y", "z"]].to_numpy(float)
                cm = np.corrcoef(A.T, B.T)[:3, 3:]  # rows: calcium xyz, cols: neuropal xyz
                fit = {}
                for i, c in enumerate("xyz"):
                    j = int(np.argmax(np.abs(cm[i])))  # best-matching neuropal axis
                    a, b = np.polyfit(A[:, i], B[:, j], 1)
                    res = B[:, j] - (a * A[:, i] + b)
                    fit[c] = dict(neuropal_axis="xyz"[j], slope=float(a), intercept=float(b), residual_rms=float(np.sqrt(np.mean(res ** 2))), r=float(cm[i, j]))
                out["neuropal_vs_calcium_t0_fit"] = dict(n_matched_labels=int(len(common)), correlation_matrix_calcium_rows_neuropal_cols=cm.round(3).tolist(), per_axis=fit)
                print(f"  NeuroPAL position vs calcium t0 centroid, matched by name ({len(common)} neurons); correlation matrix (calcium xyz rows, NeuroPAL xyz cols):\n      {np.round(cm, 2).tolist()}")
                for c, v in fit.items():
                    print(f"      calcium {c} -> NeuroPAL {v['neuropal_axis']}: slope {v['slope']:.3f}, intercept {v['intercept']:.1f}, residual rms {v['residual_rms']:.1f} voxels, r {v['r']:.3f}")
    return out


# ----------------------------------------------------------------------------
# 5. bounds
# ----------------------------------------------------------------------------
def inspect_bounds(tracks: pd.DataFrame, data_shape: tuple) -> dict:
    banner("5. ARE ROI COORDINATES IN BOUNDS? (and do they identify the axes?)")
    T, X, Y, Z = data_shape[:4]
    print("Centroids are derived from label volumes with the same (X, Y, Z) shape as the images, so they are in bounds by construction;")
    print("the axis-order question is whether the label volume's axes are the image's axes, tested with image content in section 6.")
    print(f"data shape {data_shape}; centroid ranges x[{tracks.x.min():.1f},{tracks.x.max():.1f}] y[{tracks.y.min():.1f},{tracks.y.max():.1f}] z[{tracks.z.min():.1f},{tracks.z.max():.1f}]")
    inb = ((tracks.x >= 0) & (tracks.x < X) & (tracks.y >= 0) & (tracks.y < Y) & (tracks.z >= 0) & (tracks.z < Z))
    edge = dict(frac_rows_touching_z0=float((tracks.z0 == 0).mean()), frac_rows_touching_zmax=float((tracks.z1 == Z - 1).mean()),
                frac_rows_touching_xy_border=float(((tracks.x0 == 0) | (tracks.y0 == 0) | (tracks.x1 == X - 1) | (tracks.y1 == Y - 1)).mean()))
    print(f"  all centroids in bounds: {bool(inb.all())}; masks touching the first z slice: {100*edge['frac_rows_touching_z0']:.2f}% of rows, the last: "
          f"{100*edge['frac_rows_touching_zmax']:.2f}%, the xy border: {100*edge['frac_rows_touching_xy_border']:.2f}%")
    return dict(shape=list(data_shape), all_in_bounds=bool(inb.all()), mask_edge_contact=edge)


# ----------------------------------------------------------------------------
# 6. alignment
# ----------------------------------------------------------------------------
def inspect_alignment(f: h5py.File, tracks: pd.DataFrame, seed: int, n_time_test: int = 24) -> dict:
    banner("6. DO ROI MASKS / CENTROIDS LAND ON BRIGHT SPOTS? (alignment / axis-order test, both channels)")
    chans = channel_names(f)
    rng = np.random.default_rng(seed)
    ts = [int(t) for t in sorted({0, tracks.t.max() // 4, tracks.t.max() // 2, (3 * tracks.t.max()) // 4, int(tracks.t.max())})]
    out = dict(frames=ts, per_channel={})
    n_lab = int(tracks.roi.max()) + 1
    for ci, cn in enumerate(chans):
        ratios, ratios_sw, frac_local_max, mask_ratios = [], [], [], []
        for t in ts:
            v = read_frame(f, t, ci).astype(np.float64); seg = read_seg(f, t)
            Z, Y, X = v.shape
            pts = tracks[tracks.t == t]
            z, y, x = rounded(pts)
            rz = rng.integers(z.min(), z.max() + 1, len(z)); ry = rng.integers(y.min(), y.max() + 1, len(z)); rx = rng.integers(x.min(), x.max() + 1, len(z))
            gt = local_mean(v, z, y, x, 2); rnd = local_mean(v, rz, ry, rx, 2)
            ok = (x < Y) & (y < X)
            sw = local_mean(v, z[ok], x[ok], y[ok], 2) if ok.any() else np.array([np.nan])
            # mean ratios (as for EY): ASCENT's background-subtracted median ratio degenerates when the
            # median random window equals the background, which happens in these dim 8-bit images
            ratios.append(float(gt.mean() / max(rnd.mean(), 1e-9)))
            ratios_sw.append(float(np.nanmean(sw) / max(rnd.mean(), 1e-9)))
            # mask-based: mean inside each label vs the mean of the background (label 0)
            mm = mask_means(v, seg, n_lab); bgm = v[seg == 0].mean()
            mask_ratios.append(float(np.nanmedian(mm) / bgm))
            hit = 0
            for zi, yi, xi in zip(z, y, x):
                win = v[zi, max(0, yi - 5):yi + 6, max(0, xi - 5):xi + 6]
                py, px = np.unravel_index(np.argmax(win), win.shape)
                py += max(0, yi - 5); px += max(0, xi - 5)
                hit += (abs(py - yi) <= 3) and (abs(px - xi) <= 3)
            frac_local_max.append(hit / len(z))
        aligned = bool(np.median(ratios) > 3); order_ok = bool(np.median(ratios) > np.nanmedian(ratios_sw))
        out["per_channel"][cn] = dict(gt_over_random_contrast=ratios, swapped_xy_over_random_contrast=ratios_sw, mask_mean_over_background=mask_ratios,
                                      frac_gt_within_3px_of_local_max=frac_local_max, aligned=aligned, xyz_order_confirmed=order_ok)
        print(f"[{cn}] mean 5x5 brightness at centroids / at random points in the same box: {np.round(ratios, 1).tolist()} "
              f"(swapped x<->y: {np.round(ratios_sw, 2).tolist()})")
        print(f"[{cn}] median mask-mean / background mean: {np.round(mask_ratios, 1).tolist()}; fraction of centroids within 3 px of the local in-plane maximum: {np.round(frac_local_max, 2).tolist()}")
        print(f"[{cn}] -> ROIs much brighter than random points: {aligned}; labelled (x,y,z) <-> data axes (1,2,3) beats swapped: {order_ok}")

    # time-index test: on the frames with the most motion, which label volume (t-1, t, t+1) best fits image t?
    # (EY could not answer this because its worm did not move; here the centroids move several voxels on some frames)
    ts_ = tracks.sort_values(["roi", "t"]); g = ts_.groupby("roi")
    step = g[["x", "y"]].diff(); dt = g["t"].diff()
    disp = pd.DataFrame({"t": ts_.t[dt == 1], "d": np.sqrt((step[dt == 1] ** 2).sum(axis=1))}).groupby("t")["d"].median()
    cand = [int(t) for t in disp.sort_values(ascending=False).index if 1 <= t <= tracks.t.max() - 1][:n_time_test]
    ref_ci = 1 if len(chans) == 2 else 0  # the reference channel the masks were drawn on
    wins = {"t-1": 0, "t": 0, "t+1": 0}; margins = []
    for t in cand:
        v = read_frame(f, t, ref_ci).astype(np.float64)
        bg = np.median(v)
        scores = {}
        for name, tt in (("t-1", t - 1), ("t", t), ("t+1", t + 1)):
            seg = read_seg(f, tt)
            scores[name] = float(np.nanmean(mask_means(v, seg, n_lab)) - bg)
        best = max(scores, key=scores.get); wins[best] += 1
        others = sorted(scores.values())[-2]
        margins.append((scores[best] - others) / max(abs(scores["t"]), 1e-9))
    out["time_index_test"] = dict(frames=cand, reference_channel=chans[ref_ci], wins=wins, median_margin_of_winner=float(np.median(margins)),
                                  n_test_frames_with_pop_motion_above=float(disp.loc[cand].min()))
    print(f"time-index test on the {len(cand)} frames with the largest population motion (median step >= {disp.loc[cand].min():.1f} voxels), "
          f"{chans[ref_ci]} mask mean minus background using label volume t-1 / t / t+1: wins {wins}; median winning margin {100*np.median(margins):.1f}%")
    return out


# ----------------------------------------------------------------------------
# 7. stream through every frame
# ----------------------------------------------------------------------------
def per_frame_means(f: h5py.File, n_frames: int, out_csv: Path, tracks: pd.DataFrame) -> dict:
    banner(f"7. PER-FRAME STATISTICS over all {n_frames} frames (streaming: images + label volume, one frame at a time)")
    chans = channel_names(f)
    n_lab = int(tracks.roi.max()) + 1
    rows, track_rows = [], []
    t0 = time.time()
    for t in range(n_frames):
        v = read_frame_all(f, t); seg = read_seg(f, t)
        row = {"t": t}
        for ci, cn in enumerate(chans):
            a = v[ci]
            row[f"{cn}_mean"] = float(a.mean()); row[f"{cn}_p99"] = float(np.percentile(a, 99)); row[f"{cn}_max"] = int(a.max())
            mm = mask_means(a, seg, n_lab)
            row[f"{cn}_mask_mean"] = float(np.nanmean(mm)); row[f"{cn}_bg_mean"] = float(a[seg == 0].mean())
            if ci == 0:
                tr = pd.DataFrame({"t": t, "roi": np.arange(n_lab)})
            tr[f"{cn}_mask"] = mm
            pts = tracks[tracks.t == t].set_index("roi")
            lm = np.full(n_lab, np.nan)
            if len(pts):
                z, y, x = rounded(pts); lm[pts.index.to_numpy()] = local_mean(a, z, y, x, 2)
            tr[f"{cn}_local5"] = lm
        row["n_rois"] = int(len(tracks[tracks.t == t]))
        track_rows.append(tr); rows.append(row)
        if t % 100 == 0:
            print(f"  t={t}: " + ", ".join(f"{k}={v:.4g}" for k, v in row.items() if k != "t") + f"  ({time.time() - t0:.0f} s)")
    print(f"  streamed {n_frames} frames in {time.time() - t0:.0f} s")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pf = pd.DataFrame(rows); pf.to_csv(out_csv, index=False); print(f"wrote {out_csv}")
    tp = out_csv.with_name(out_csv.stem + "_tracks.csv")
    mine = pd.concat(track_rows, ignore_index=True)
    labels = labels_list(f)
    if labels:
        mine["label"] = np.array(labels + [""] * (n_lab - len(labels)), dtype=object)[mine.roi.to_numpy()]
    mine.to_csv(tp, index=False); print(f"wrote {tp}")
    out = {f"{cn}_global_max": int(pf[f"{cn}_max"].max()) for cn in chans}
    out["brightness_trend"] = {cn: dict(first_over_median=float(pf[f"{cn}_mask_mean"].iloc[0] / pf[f"{cn}_mask_mean"].median()),
                                        last_over_median=float(pf[f"{cn}_mask_mean"].iloc[-1] / pf[f"{cn}_mask_mean"].median())) for cn in chans}
    print(f"  mask-mean brightness first/last frame relative to the median: {out['brightness_trend']}")

    # compare with the traces stored in the file
    stored = f[f"{TRACE}/data"][()]  # (T, N)
    sp = out_csv.with_name("stored_traces.csv")
    pd.DataFrame(stored, columns=[f"roi{i}" for i in range(stored.shape[1])]).rename_axis("t").to_csv(sp); print(f"wrote {sp}")
    n = min(n_frames, stored.shape[0])
    comp = {}
    cands = {f"{chans[0]}_mask": mine.pivot(index="t", columns="roi", values=f"{chans[0]}_mask"),
             f"{chans[0]}_local5": mine.pivot(index="t", columns="roi", values=f"{chans[0]}_local5")}
    if len(chans) == 2:
        cands[f"{chans[0]}_mask/{chans[1]}_mask"] = cands[f"{chans[0]}_mask"] / mine.pivot(index="t", columns="roi", values=f"{chans[1]}_mask")
        cands[f"{chans[1]}_mask"] = mine.pivot(index="t", columns="roi", values=f"{chans[1]}_mask")
    for name, wide in cands.items():
        cors, rel = [], []
        for j in range(min(wide.shape[1], stored.shape[1])):
            a = wide.iloc[:n, j].to_numpy(float); b = stored[:n, j]; ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() > 10:
                cors.append(np.corrcoef(a[ok], b[ok])[0, 1]); rel.append(np.median(np.abs(a[ok] - b[ok]) / np.maximum(np.abs(b[ok]), 1e-9)))
            else:
                cors.append(np.nan); rel.append(np.nan)
        comp[name] = dict(corr_median=float(np.nanmedian(cors)), corr_p10=float(np.nanpercentile(cors, 10)), corr_min=float(np.nanmin(cors)),
                          frac_corr_above_0_9=float(np.nanmean(np.array(cors) > 0.9)), frac_corr_above_0_99=float(np.nanmean(np.array(cors) > 0.99)),
                          median_relative_abs_difference=float(np.nanmedian(rel)))
        c = comp[name]
        print(f"  stored SignalRawFluor vs our {name}: per-neuron corr median {c['corr_median']:.3f}, p10 {c['corr_p10']:.3f}, min {c['corr_min']:.3f}; "
              f"{100*c['frac_corr_above_0_9']:.0f}% > 0.9, {100*c['frac_corr_above_0_99']:.0f}% > 0.99; median |ours - stored|/stored {c['median_relative_abs_difference']:.3f}")
    out["stored_trace_comparison"] = dict(n_frames_compared=int(n), stored_value_range=[float(np.nanmin(stored)), float(np.nanmax(stored))],
                                          stored_has_nan=bool(np.isnan(stored).any()), n_stored_nan=int(np.isnan(stored).sum()), candidates=comp)
    print(f"  stored values {out['stored_trace_comparison']['stored_value_range']}, NaN entries {out['stored_trace_comparison']['n_stored_nan']}")
    # stored-trace jump statistic (the DANDI description warns of transient tracking failures)
    d = np.abs(np.diff(stored, axis=0)) / np.maximum(np.abs(stored[:-1]), 1e-9)
    out["stored_trace_jumps"] = dict(frac_steps_change_above_50pct=float(np.nanmean(d > 0.5)), frac_steps_change_above_100pct=float(np.nanmean(d > 1.0)),
                                     n_neurons_with_any_step_above_100pct=int((np.nanmax(d, axis=0) > 1.0).sum()))
    j = out["stored_trace_jumps"]
    print(f"  stored-trace frame-to-frame changes > 50%: {100*j['frac_steps_change_above_50pct']:.2f}% of steps; > 100%: {100*j['frac_steps_change_above_100pct']:.2f}%; "
          f"neurons with at least one > 100% step: {j['n_neurons_with_any_step_above_100pct']}")
    return out


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nwb", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--tracks-csv", type=Path, default=None, help="cache for the derived per-frame centroids (written if absent, read if present)")
    ap.add_argument("--n-sample-frames", type=int, default=40)
    ap.add_argument("--n-neurons", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--per-frame-means", type=Path, default=None, help="CSV path; streams through ALL frames (images and label volumes)")
    ap.add_argument("--skip-pynwb", action="store_true")
    args = ap.parse_args(argv)

    summary = dict(nwb=str(args.nwb))
    with h5py.File(args.nwb, "r") as f:
        s = inspect_structure(f); summary["structure"] = s
        summary["metadata"] = inspect_metadata(f)
        if not args.skip_pynwb:
            summary["pynwb"] = validate_with_pynwb(args.nwb)
        banner("4a. DERIVING PER-FRAME CENTROIDS FROM THE LABEL VOLUMES")
        tracks = load_calcium_tracks(f, args.tracks_csv)
        n = s["n_frames"]
        summary["annotations"] = inspect_annotations(f, tracks, n)
        summary["bounds"] = inspect_bounds(tracks, tuple(s["calcium"]["shape"]))
        summary["channels"] = inspect_channels(f, tracks, n, args.n_sample_frames, args.n_neurons, args.seed)
        summary["alignment"] = inspect_alignment(f, tracks, args.seed)
        if args.per_frame_means:
            summary["per_frame"] = per_frame_means(f, n, args.per_frame_means, tracks)

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(jsonable(summary), indent=2))
        print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
