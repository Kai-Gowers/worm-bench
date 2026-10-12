#!/usr/bin/env python3
"""Inspect one NWB recording from the EY dataset (DANDI 000541, Yemini et al.).

Mirrors scripts/inspect_ascent.py section by section so that the two
inspections are comparable, with the loading and annotation extraction adapted
to the NWB layout written by the ndx-multichannel-volume extension:

  acquisition/CalciumImageSeries/data            (T, X, Y, Z) uint16, one channel
  processing/CalciumActivity/CalciumSeriesSegmentation/Seg_tpoint_{t}/voxel_mask
                                                 one (x, y, z, weight) row per ROI
  processing/CalciumActivity/NeuronIDs/labels    biological name per ROI
  processing/CalciumActivity/SignalRawFluor/...  stored activity traces (T-1, N)
  acquisition/NeuroPALImageRaw/data              (X, Y, Z, C) structural stack
  intervals/chemical_stimuli                      stimulus timing

Only a few dozen frames are read unless --per-frame-means is given, in which
case every frame is streamed once (one frame in memory at a time).

Usage:
  python -P scripts/inspect_ey.py --nwb data/ey/<file>.nwb \
      --out-json previews/ey/inspection_summary.json \
      [--tracks-csv previews/ey/calcium_tracks.csv] \
      [--per-frame-means previews/ey/per_frame_means.csv]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))  # sibling scripts, not the data dir
from inspect_ascent import banner, jsonable, local_mean  # noqa: E402

CALCIUM = "acquisition/CalciumImageSeries"
SEG = "processing/CalciumActivity/CalciumSeriesSegmentation"
SEG_DNMF = "processing/CalciumActivity/CalciumSeriesSegmentationdNMF"
LABELS = "processing/CalciumActivity/NeuronIDs/labels"
TRACE = "processing/CalciumActivity/SignalRawFluor/SignalCalciumImResponseSeries"
TRACE_DNMF = "processing/CalciumActivity/SignalRawFluordNMF/dNMFCalciumImResponseSeries"
NEUROPAL_IMG = "acquisition/NeuroPALImageRaw"
NEUROPAL_SEG = "processing/NeuroPAL/NeuroPALSegmentation/NeuroPALNeurons"
STIM = "intervals/chemical_stimuli"
TPOINT_RE = re.compile(r"^Seg_tpoint_(\d+)$")


# ----------------------------------------------------------------------------
# loading helpers (shared with make_previews_ey.py)
# ----------------------------------------------------------------------------
def read_frame(f: h5py.File, t: int) -> np.ndarray:
    """One calcium volume as (Z, Y, X), the axis order used by the ASCENT helpers.

    The dataset is stored (T, X, Y, Z); the transpose is a view, no copy."""
    return f[f"{CALCIUM}/data"][t].transpose(2, 1, 0)


def decode(v):
    return v.decode() if isinstance(v, bytes) else str(v)


def ragged_labels(group: h5py.Group, name: str = "ID_labels") -> list[str]:
    """NeuroPAL labels are stored one character per element with a VectorIndex."""
    chars = [decode(c) for c in group[name][()]]
    idx = group[f"{name}_index"][()]
    starts = np.r_[0, idx[:-1]]
    return ["".join(chars[a:b]) for a, b in zip(starts, idx)]


def load_calcium_tracks(f: h5py.File, group: str = SEG) -> pd.DataFrame:
    """All per-frame ROI positions as one long table: t, roi, x, y, z, weight, label.

    Each Seg_tpoint_{t} PlaneSegmentation holds one voxel per ROI; ROI order is
    the same in every frame (checked in inspect_annotations), so `roi` doubles
    as a track id."""
    g = f[group]
    labels = [decode(v) for v in f[LABELS][()]] if LABELS in f else None
    parts = []
    for k in g.keys():
        m = TPOINT_RE.match(k)
        if not m:
            continue
        t = int(m.group(1))
        vm = g[k]["voxel_mask"][()]
        ids = g[k]["id"][()]
        parts.append(pd.DataFrame({"t": t, "roi": ids, "x": vm["x"].astype(np.int64), "y": vm["y"].astype(np.int64),
                                   "z": vm["z"].astype(np.int64), "weight": vm["weight"].astype(np.float64)}))
    df = pd.concat(parts, ignore_index=True).sort_values(["t", "roi"]).reset_index(drop=True)
    if labels is not None and len(labels) == df.roi.max() + 1:
        df["label"] = np.array(labels, dtype=object)[df.roi.to_numpy()]
    return df


def in_bounds(pts: pd.DataFrame, zyx_shape: tuple[int, int, int]) -> np.ndarray:
    """True for rows whose (z, y, x) index lies inside a (Z, Y, X) volume.  Some ROI
    rows in this dataset have z >= Z (see section 4/5); image-based tests skip them."""
    Z, Y, X = zyx_shape
    return ((pts.z >= 0) & (pts.z < Z) & (pts.y >= 0) & (pts.y < Y) & (pts.x >= 0) & (pts.x < X)).to_numpy()


def load_stimuli(f: h5py.File) -> list[dict]:
    if STIM not in f:
        return []
    g = f[STIM]
    return [dict(stimulus=decode(s), start_s=float(a), stop_s=float(b))
            for s, a, b in zip(g["stimulus"][()], g["start_time"][()], g["stop_time"][()])]


# ----------------------------------------------------------------------------
# 1. structure
# ----------------------------------------------------------------------------
def inspect_structure(f: h5py.File) -> dict:
    banner("1. NWB / HDF5 STRUCTURE (verified from file)")
    out = dict(root_keys=list(f.keys()), nwb_version=decode(f.attrs.get("nwb_version", "")))
    print(f"root keys: {out['root_keys']}; nwb_version {out['nwb_version']}")
    print(f"specification namespaces embedded: {list(f['specifications'].keys()) if 'specifications' in f else None}")
    out["spec_namespaces"] = {k: list(f["specifications"][k].keys()) for k in f["specifications"].keys()} if "specifications" in f else {}

    d = f[f"{CALCIUM}/data"]
    layout = d.id.get_create_plist().get_layout()  # 0 compact, 1 contiguous, 2 chunked
    out["calcium"] = dict(shape=list(d.shape), dtype=str(d.dtype), chunks=d.chunks, compression=d.compression,
                          layout={0: "compact", 1: "contiguous", 2: "chunked"}.get(layout, str(layout)),
                          nbytes=int(d.nbytes), bytes_per_frame=int(d.nbytes // d.shape[0]),
                          dimension_dataset=f[f"{CALCIUM}/dimension"][()].tolist(),
                          rate_hz=float(f[f"{CALCIUM}/starting_time"].attrs["rate"]),
                          starting_time_s=float(f[f"{CALCIUM}/starting_time"][()]),
                          attrs={k: jsonable(v) for k, v in f[CALCIUM].attrs.items() if k != "object_id"},
                          data_attrs={k: jsonable(v) for k, v in d.attrs.items()})
    c = out["calcium"]
    print(f"CalciumImageSeries/data: shape {c['shape']} dtype {c['dtype']} layout {c['layout']} chunks {c['chunks']} "
          f"compression {c['compression']} ({c['nbytes']/1e9:.2f} GB, {c['bytes_per_frame']/1e6:.2f} MB/frame)")
    print(f"  dimension dataset = {c['dimension_dataset']}  (documents the per-frame axes as X, Y, Z)")
    print(f"  starting_time {c['starting_time_s']} s, rate {c['rate_hz']} Hz -> {d.shape[0]} frames = {d.shape[0]/c['rate_hz']:.1f} s")
    print(f"  attrs: {c['attrs']}")
    print(f"  data attrs: {c['data_attrs']}")
    n_frames = d.shape[0]
    out["n_frames"] = n_frames
    out["axis_order"] = "T,X,Y,Z as documented by `dimension`; tested in section 6"

    if NEUROPAL_IMG in f:
        ni = f[f"{NEUROPAL_IMG}/data"]
        out["neuropal_image"] = dict(shape=list(ni.shape), dtype=str(ni.dtype), nbytes=int(ni.nbytes),
                                     RGBW_channels=f[f"{NEUROPAL_IMG}/RGBW_channels"][()].tolist(),
                                     description=decode(f[f"{NEUROPAL_IMG}/description"][()]))
        print(f"NeuroPALImageRaw/data: shape {ni.shape} dtype {ni.dtype} ({ni.nbytes/1e6:.0f} MB); RGBW_channels {out['neuropal_image']['RGBW_channels']}")

    seg_keys = [k for k in f[SEG].keys() if TPOINT_RE.match(k)]
    tps = sorted(int(TPOINT_RE.match(k).group(1)) for k in seg_keys)
    out["segmentation"] = dict(n_tpoints=len(tps), t_min=tps[0], t_max=tps[-1], contiguous=tps == list(range(tps[0], tps[-1] + 1)),
                               members_t0=list(f[f"{SEG}/Seg_tpoint_0"].keys()),
                               voxel_mask_dtype=str(f[f"{SEG}/Seg_tpoint_0/voxel_mask"].dtype),
                               description_t0=decode(f[f"{SEG}/Seg_tpoint_0"].attrs.get("description", "")),
                               has_dnmf_copy=SEG_DNMF in f)
    s = out["segmentation"]
    print(f"CalciumSeriesSegmentation: {s['n_tpoints']} PlaneSegmentation groups Seg_tpoint_{s['t_min']}..{s['t_max']} "
          f"(contiguous {s['contiguous']}); movie has {n_frames} frames -> {n_frames - s['n_tpoints']} frame(s) without annotation")
    print(f"  members: {s['members_t0']}; voxel_mask dtype {s['voxel_mask_dtype']}; dNMF copy present: {s['has_dnmf_copy']}")

    for name, path in [("SignalRawFluor", TRACE), ("SignalRawFluordNMF", TRACE_DNMF)]:
        if path in f:
            td = f[f"{path}/data"]
            out[name] = dict(shape=list(td.shape), dtype=str(td.dtype), rate_hz=float(f[f"{path}/starting_time"].attrs["rate"]),
                             unit=decode(td.attrs.get("unit", "")), description=decode(f[path].attrs.get("description", "")),
                             rois_table=f[f"{path}/rois"].attrs["table"] and f[f[f"{path}/rois"].attrs["table"]].name if "rois" in f[path] else None)
            print(f"{name}: shape {td.shape} {td.dtype} unit '{out[name]['unit']}' rate {out[name]['rate_hz']} Hz; rois -> {out[name]['rois_table']}")

    # pixel-value spot checks on three frames (reads 3 frames)
    spot = {}
    for t in sorted({0, n_frames // 2, n_frames - 1}):
        v = read_frame(f, t)
        spot[f"t{t}"] = dict(min=int(v.min()), max=int(v.max()), mean=float(v.mean()), frac_zero=float((v == 0).mean()),
                             frac_at_255=float((v == 255).mean()), p50=float(np.percentile(v, 50)), p99=float(np.percentile(v, 99)))
        print(f"  t{t}: min={v.min()} max={v.max()} mean={v.mean():.2f} p50={spot[f't{t}']['p50']:.0f} p99={spot[f't{t}']['p99']:.0f} "
              f"zero-fraction={spot[f't{t}']['frac_zero']:.3f} at-255-fraction={spot[f't{t}']['frac_at_255']:.2e}")
    out["spot_checks"] = spot
    if all(s["max"] <= 255 for s in spot.values()):
        print("  NOTE: every sampled voxel is <= 255 although the dtype is uint16 -> values look 8-bit (rescaled?), see section 7 for the full-movie maximum")
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
    out["stimuli"] = load_stimuli(f)
    out["stimuli_description"] = decode(f[STIM].attrs.get("description", "")) if STIM in f else None
    out["neuron_ids_description"] = decode(f["processing/CalciumActivity/NeuronIDs/description"][()]) if "processing/CalciumActivity/NeuronIDs" in f else None
    out["neuropal_seg_description"] = decode(f[NEUROPAL_SEG].attrs.get("description", "")) if NEUROPAL_SEG in f else None

    for k in ["identifier", "session_description", "session_start_time", "file_create_date", "general/experiment_description",
              "general/experimenter", "general/institution", "general/lab", "general/related_publications"]:
        print(f"  {k}: {out.get(k)}")
    print(f"  subject: {out['subject']}")
    print(f"  devices: {out['devices']}")
    for name, v in vols.items():
        print(f"  imaging volume {name}: grid_spacing {v.get('grid_spacing')} {v.get('grid_spacing_unit')}, location {v.get('location')}, "
              f"channels {list(v['channels'].keys())}, order {v.get('order_optical_channels')}")
        for ck, cv in v["channels"].items():
            print(f"      {ck}: ex {cv.get('excitation_lambda')} nm, em {cv.get('emission_lambda')} nm, filter '{cv.get('description')}'")
    print(f"  chemical stimuli: {out['stimuli']}")
    print(f"    description: {out['stimuli_description']}")
    print(f"  NeuronIDs description: {out['neuron_ids_description']}")
    print(f"  NeuroPAL segmentation description: {out['neuropal_seg_description']}")
    return out


def validate_with_pynwb(path: Path) -> dict:
    """Optional: confirm the file opens with pynwb and its embedded extension spec."""
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
            info = dict(ok=True, identifier=nwb.identifier, session_description=nwb.session_description,
                        acquisition=list(nwb.acquisition.keys()), processing=list(nwb.processing.keys()),
                        imaging_planes=list(nwb.imaging_planes.keys()), intervals=list(nwb.intervals.keys()),
                        seconds=round(time.time() - t0, 1))
            ps = nwb.processing["CalciumActivity"]["CalciumSeriesSegmentation"].plane_segmentations
            info["n_plane_segmentations"] = len(ps)
        print(f"pynwb read OK in {info['seconds']} s: acquisition {info['acquisition']}, processing {info['processing']}, "
              f"{info['n_plane_segmentations']} plane segmentations")
        return info
    except Exception as e:  # report, do not abort
        print(f"pynwb read FAILED: {type(e).__name__}: {e}")
        return dict(ok=False, error=f"{type(e).__name__}: {e}")


# ----------------------------------------------------------------------------
# 3. channel behaviour (single channel)
# ----------------------------------------------------------------------------
def inspect_channels(f: h5py.File, tracks: pd.DataFrame, n_frames: int, n_sample: int, n_neurons: int, seed: int) -> dict:
    banner("3. CHANNEL BEHAVIOUR (single GCaMP channel; no reference channel to compare against)")
    print("ASCENT's calcium-vs-reference comparison is inapplicable here: the series has one channel.")
    print("Instead: per-neuron temporal variability of that channel on sampled frames (same statistic as ASCENT),")
    print("and, in section 7, stimulus-locked responses against the stimulus times stored in the file.")
    rng = np.random.default_rng(seed)
    annotated = np.sort(tracks.t.unique())
    ts = np.unique(np.linspace(0, len(annotated) - 1, n_sample).round().astype(int))
    ts = annotated[ts]
    n_full = tracks.groupby("roi")["t"].nunique()
    full = n_full[n_full == len(annotated)].index.to_numpy()
    rois = rng.choice(full, size=min(n_neurons, len(full)), replace=False)
    sel = tracks[tracks.roi.isin(rois)].set_index(["t", "roi"]).sort_index()
    vol_mean, p50, p99, vmax = [], [], [], []
    neuron_ts = np.zeros((len(ts), len(rois)))
    for i, t in enumerate(ts):
        v = read_frame(f, int(t))
        vol_mean.append(v.mean()); p50.append(np.percentile(v, 50)); p99.append(np.percentile(v, 99)); vmax.append(v.max())
        pts = sel.loc[int(t)].loc[rois]
        ok = in_bounds(pts, v.shape)
        vals = np.full(len(rois), np.nan)
        vals[ok] = local_mean(v, pts["z"].to_numpy()[ok], pts["y"].to_numpy()[ok], pts["x"].to_numpy()[ok], r=2)
        neuron_ts[i] = vals
    vol_mean = np.array(vol_mean)
    cv = []
    for j in range(len(rois)):
        y = neuron_ts[:, j]; ok = np.isfinite(y)
        trend = np.polyval(np.polyfit(ts[ok], y[ok], 1), ts[ok])
        cv.append(np.std(y[ok] - trend) / max(y[ok].mean(), 1e-9))
    cv = np.array(cv)
    n_oob = int(np.isnan(neuron_ts).sum())
    out = dict(sampled_frames=ts.tolist(), rois=rois.tolist(), volume_mean=vol_mean.tolist(), p50=np.array(p50).tolist(),
               p99=np.array(p99).tolist(), max=np.array(vmax).tolist(),
               volume_mean_temporal_cv=float(vol_mean.std() / vol_mean.mean()),
               neuron_detrended_cv_median=float(np.median(cv)), neuron_detrended_cv_p90=float(np.percentile(cv, 90)),
               neuron_mean_intensity_median=float(np.nanmedian(neuron_ts)),
               neuron_mean_intensity_range=[float(np.nanmean(neuron_ts, 0).min()), float(np.nanmean(neuron_ts, 0).max())],
               n_sampled_points_out_of_bounds=n_oob)
    print(f"volume mean over {len(ts)} sampled frames: {vol_mean.min():.2f}..{vol_mean.max():.2f} (temporal CV {out['volume_mean_temporal_cv']:.3f}); "
          f"p99 {np.min(p99):.0f}..{np.max(p99):.0f}; max {np.min(vmax)}..{np.max(vmax)}")
    print(f"per-neuron detrended temporal CV over {len(rois)} full-length tracks: median {out['neuron_detrended_cv_median']:.3f}, "
          f"p90 {out['neuron_detrended_cv_p90']:.3f}; median neuron brightness {out['neuron_mean_intensity_median']:.1f} "
          f"(per-neuron means {out['neuron_mean_intensity_range'][0]:.1f}..{out['neuron_mean_intensity_range'][1]:.1f})")
    print(f"  ({n_oob} of {neuron_ts.size} sampled (frame, neuron) points skipped because the ROI lies outside the stack)")
    print("  (ASCENT, same statistic: c1/GCaMP median CV was the larger of its two channels; here there is no second channel to rank against)")
    return out


# ----------------------------------------------------------------------------
# 4. annotations
# ----------------------------------------------------------------------------
def inspect_annotations(f: h5py.File, tracks: pd.DataFrame, n_frames: int) -> dict:
    banner("4. ANNOTATIONS: per-frame ROI positions, track structure, labels")
    out = {}
    per_t = tracks.groupby("t").size()
    n_t = tracks.t.nunique()
    out["positions"] = dict(
        columns=list(tracks.columns), n_rows=len(tracks), t_min=int(tracks.t.min()), t_max=int(tracks.t.max()), n_time_points=int(n_t),
        missing_t_in_range=int(tracks.t.max() - tracks.t.min() + 1 - n_t),
        unannotated_frames=sorted(set(range(n_frames)) - set(tracks.t.unique())),
        points_per_frame=dict(min=int(per_t.min()), median=float(per_t.median()), max=int(per_t.max())),
        integer_coords=True,  # voxel_mask stores unsigned ints
        ranges={c: [int(tracks[c].min()), int(tracks[c].max())] for c in "xyz"},
        weights_unique=np.unique(tracks.weight).tolist(),
        fraction_of_movie_frames_annotated=float(n_t / n_frames),
    )
    p = out["positions"]
    print(f"{len(tracks)} ROI rows over {n_t} time points (t {p['t_min']}..{p['t_max']}, {p['missing_t_in_range']} missing in range); "
          f"unannotated movie frames: {p['unannotated_frames']}")
    print(f"  points per frame: {p['points_per_frame']}; coords stored as unsigned ints; ranges x{p['ranges']['x']} y{p['ranges']['y']} z{p['ranges']['z']}; "
          f"weights {p['weights_unique']}")

    # one voxel per ROI? (voxel_mask_index must be 1..n in every frame)
    one_voxel = True
    same_ids = True
    ids0 = f[f"{SEG}/Seg_tpoint_0/id"][()]
    for k in f[SEG].keys():
        if TPOINT_RE.match(k):
            vmi = f[f"{SEG}/{k}/voxel_mask_index"][()]
            one_voxel &= bool(np.array_equal(vmi, np.arange(1, len(vmi) + 1)))
            same_ids &= bool(np.array_equal(f[f"{SEG}/{k}/id"][()], ids0))
    out["one_voxel_per_roi_every_frame"] = one_voxel
    out["same_roi_ids_every_frame"] = same_ids
    print(f"  exactly one voxel per ROI in every frame: {one_voxel}; identical ROI id list in every frame: {same_ids}")

    # track structure: roi index as track id
    tl = tracks.groupby("roi")["t"].nunique()
    out["tracks"] = dict(n_tracks=int(tracks.roi.nunique()), roi_min=int(tracks.roi.min()), roi_max=int(tracks.roi.max()),
                         track_length=dict(min=int(tl.min()), median=float(tl.median()), max=int(tl.max())),
                         n_tracks_spanning_all_annotated_frames=int((tl == n_t).sum()),
                         duplicate_roi_time_pairs=int(tracks.duplicated(["roi", "t"]).sum()),
                         coincident_points=int(tracks.duplicated(["t", "x", "y", "z"]).sum()))
    tr = out["tracks"]
    print(f"  tracks (ROI index): {tr['n_tracks']} (ids {tr['roi_min']}..{tr['roi_max']}), length min/median/max "
          f"{tr['track_length']['min']}/{tr['track_length']['median']}/{tr['track_length']['max']}, "
          f"{tr['n_tracks_spanning_all_annotated_frames']} span all {n_t} annotated frames; duplicate (roi,t) {tr['duplicate_roi_time_pairs']}; "
          f"two ROIs at the same voxel in the same frame: {tr['coincident_points']} rows")

    # motion
    ts = tracks.sort_values(["roi", "t"])
    step = ts.groupby("roi")[["z", "y", "x"]].diff().dropna()
    disp = np.sqrt((step[["y", "x"]] ** 2).sum(axis=1))
    out["per_step_motion_voxels"] = dict(xy_median=float(disp.median()), xy_p99=float(disp.quantile(0.99)), xy_max=float(disp.max()),
                                         z_max=int(step.z.abs().max()), frac_steps_zero=float((disp == 0).mean()),
                                         frac_steps_z_nonzero=float((step.z != 0).mean()))
    m = out["per_step_motion_voxels"]
    print(f"  motion between consecutive annotated frames (voxels): in-plane median {m['xy_median']:.2f}, p99 {m['xy_p99']:.1f}, max {m['xy_max']:.1f}; "
          f"max |dz| {m['z_max']}; {100*m['frac_steps_zero']:.0f}% of steps have zero in-plane motion; {100*m['frac_steps_z_nonzero']:.1f}% change z")
    # total drift: displacement of each track between first and last frame
    first = ts.groupby("roi").first(); last = ts.groupby("roi").last()
    drift = np.sqrt(((last[["x", "y"]] - first[["x", "y"]]) ** 2).sum(axis=1))
    out["first_to_last_drift_voxels"] = dict(median=float(drift.median()), max=float(drift.max()))
    print(f"  in-plane displacement first->last annotated frame per track: median {drift.median():.1f}, max {drift.max():.1f} voxels")

    # rows outside the image stack (z >= Z etc.)
    T, X, Y, Z = f[f"{CALCIUM}/data"].shape
    oob = tracks[~in_bounds(tracks, (Z, Y, X))]
    per_roi = oob.groupby("roi").size().sort_values(ascending=False)
    out["out_of_bounds"] = dict(n_rows=int(len(oob)), frac_rows=float(len(oob) / len(tracks)), n_rois_affected=int(per_roi.size),
                                n_frames_affected=int(oob.t.nunique()), z_values=sorted(oob.z.unique().tolist()),
                                axes_affected={c: int(((tracks[c] < 0) | (tracks[c] >= dim)).sum()) for c, dim in zip("xyz", (X, Y, Z))},
                                worst_rois=[(int(r), tracks.loc[tracks.roi == r, "label"].iloc[0] if "label" in tracks else "", int(n)) for r, n in per_roi.head(8).items()])
    ob = out["out_of_bounds"]
    print(f"  rows outside the (X,Y,Z)=({X},{Y},{Z}) stack: {ob['n_rows']} ({100*ob['frac_rows']:.2f}%), all on the z axis: {ob['axes_affected']}, "
          f"z values {ob['z_values']}; {ob['n_rois_affected']} ROIs and {ob['n_frames_affected']} frames affected; "
          f"worst (roi, label, n frames): {ob['worst_rois']}")

    # z distribution
    zc = tracks.z.value_counts().sort_index()
    out["z_histogram"] = {int(k): int(v) for k, v in zc.items()}
    print(f"  rows per z slice: {out['z_histogram']}")

    # labels
    if "label" in tracks:
        labels = tracks.drop_duplicates("roi").set_index("roi")["label"].sort_index()
        nonempty = labels[labels.str.strip() != ""]
        lr = sum(1 for l in nonempty if l.endswith("L") and (l[:-1] + "R") in set(nonempty))
        out["labels"] = dict(n=int(len(labels)), n_nonempty=int(len(nonempty)), n_unique_nonempty=int(nonempty.nunique()),
                             n_LR_pairs=int(lr), examples=labels.head(12).tolist(),
                             duplicated_labels=nonempty[nonempty.duplicated()].tolist())
        print(f"  biological labels (NeuronIDs/labels, one per ROI): {len(labels)} total, {len(nonempty)} non-empty, "
              f"{nonempty.nunique()} unique, {lr} left/right pairs; duplicates {out['labels']['duplicated_labels']}; e.g. {labels.head(12).tolist()}")

    # dNMF copy of the segmentation
    if SEG_DNMF in f:
        d2 = load_calcium_tracks(f, SEG_DNMF).drop(columns=["label"], errors="ignore")
        mm = tracks.merge(d2, on=["t", "roi"], suffixes=("", "_dnmf"))
        diff = (mm[["x", "y", "z"]].to_numpy() - mm[["x_dnmf", "y_dnmf", "z_dnmf"]].to_numpy())
        ad = np.abs(diff)
        out["dnmf_segmentation_vs_main"] = dict(n_rows=len(mm), n_tpoints=int(d2.t.nunique()), n_rois=int(d2.roi.nunique()),
                                                frac_rows_identical=float((ad.sum(1) == 0).mean()),
                                                max_abs_diff_xyz=ad.max(0).tolist(), mean_abs_diff_xyz=ad.mean(0).tolist())
        dn = out["dnmf_segmentation_vs_main"]
        print(f"  dNMF segmentation copy: {dn['n_tpoints']} tpoints x {dn['n_rois']} ROIs; identical to main in {100*dn['frac_rows_identical']:.1f}% of rows; "
              f"max |diff| xyz {dn['max_abs_diff_xyz']}, mean |diff| {np.round(dn['mean_abs_diff_xyz'], 2).tolist()}")

    # NeuroPAL segmentation and its relation to calcium frame 0
    if NEUROPAL_SEG in f:
        g = f[NEUROPAL_SEG]
        vm = g["voxel_mask"][()]
        npal = pd.DataFrame({"x": vm["x"].astype(int), "y": vm["y"].astype(int), "z": vm["z"].astype(int), "label": ragged_labels(g)})
        ni_shape = list(f[f"{NEUROPAL_IMG}/data"].shape) if NEUROPAL_IMG in f else None
        out["neuropal_segmentation"] = dict(n_rois=len(npal), ranges={c: [int(npal[c].min()), int(npal[c].max())] for c in "xyz"},
                                            image_shape=ni_shape, n_labels_nonempty=int((npal.label.str.strip() != "").sum()),
                                            label_set_equals_calcium=bool(set(npal.label) == set(tracks.label.unique())) if "label" in tracks else None,
                                            same_order_as_calcium=bool(npal.label.tolist() == tracks[tracks.t == tracks.t.min()].sort_values("roi").label.tolist()) if "label" in tracks else None)
        ns = out["neuropal_segmentation"]
        print(f"  NeuroPAL segmentation: {ns['n_rois']} ROIs in NeuroPAL image coords (image shape {ns['image_shape']}), ranges {ns['ranges']}; "
              f"label set equals calcium label set: {ns['label_set_equals_calcium']}; same ROI order: {ns['same_order_as_calcium']}")
        if "label" in tracks:
            t0 = tracks[tracks.t == tracks.t.min()].set_index("label")
            both = npal[npal.label.isin(t0.index) & (npal.label.str.strip() != "")].drop_duplicates("label").set_index("label")
            fit = {}
            for c in "xyz":
                a, b = np.polyfit(t0.loc[both.index, c].to_numpy(float), both[c].to_numpy(float), 1)
                res = both[c].to_numpy(float) - (a * t0.loc[both.index, c].to_numpy(float) + b)
                fit[c] = dict(slope=float(a), intercept=float(b), residual_rms=float(np.sqrt(np.mean(res ** 2))), r=float(np.corrcoef(t0.loc[both.index, c], both[c])[0, 1]))
            out["neuropal_vs_calcium_t0_fit"] = dict(n_matched_labels=int(len(both)), per_axis=fit)
            print(f"  NeuroPAL position vs calcium t0 position, matched by label ({len(both)} neurons), per-axis linear fit neuropal = slope*calcium + intercept:")
            for c, v in fit.items():
                print(f"      {c}: slope {v['slope']:.3f}, intercept {v['intercept']:.1f}, residual rms {v['residual_rms']:.1f} voxels, r {v['r']:.3f}")
    return out


# ----------------------------------------------------------------------------
# 5. bounds
# ----------------------------------------------------------------------------
def inspect_bounds(tracks: pd.DataFrame, data_shape: tuple) -> dict:
    banner("5. ARE ROI COORDINATES IN BOUNDS? (and do they identify the axes?)")
    T, X, Y, Z = data_shape
    checks = {}
    for name, (xc, yc, zc) in {"as_labelled_xyz": ("x", "y", "z"), "swapped_xy": ("y", "x", "z")}.items():
        inb = ((tracks[xc] >= 0) & (tracks[xc] < X) & (tracks[yc] >= 0) & (tracks[yc] < Y) & (tracks[zc] >= 0) & (tracks[zc] < Z))
        checks[name] = dict(n_out_of_bounds=int((~inb).sum()), all_in_bounds=bool(inb.all()))
    print(f"data shape (T,X,Y,Z)={data_shape}; coords ranges x[{tracks.x.min()},{tracks.x.max()}] y[{tracks.y.min()},{tracks.y.max()}] z[{tracks.z.min()},{tracks.z.max()}]")
    for k, v in checks.items():
        print(f"  {k}: out-of-bounds rows = {v['n_out_of_bounds']}")
    print("  note: unlike ASCENT's square field, X=256 and Y=128 differ, so bounds alone already discriminate x from y;"
          " section 6 confirms with image content")
    return dict(shape=list(data_shape), checks=checks)


# ----------------------------------------------------------------------------
# 6. alignment
# ----------------------------------------------------------------------------
def inspect_alignment(f: h5py.File, tracks: pd.DataFrame, seed: int) -> dict:
    banner("6. DO ROI POINTS LAND ON BRIGHT SPOTS? (alignment / axis-order test)")
    rng = np.random.default_rng(seed)
    annotated = np.sort(tracks.t.unique())
    ts = [int(annotated[i]) for i in sorted({0, len(annotated) // 4, len(annotated) // 2, (3 * len(annotated)) // 4, len(annotated) - 1})]
    ratios, ratios_sw, frac_local_max, nonzero = [], [], [], []
    for t in ts:
        v = read_frame(f, t).astype(np.float64)
        Z, Y, X = v.shape
        pts = tracks[tracks.t == t]
        pts = pts[in_bounds(pts, v.shape)]
        z, y, x = pts.z.to_numpy(), pts.y.to_numpy(), pts.x.to_numpy()
        rz = rng.integers(z.min(), z.max() + 1, len(z)); ry = rng.integers(y.min(), y.max() + 1, len(z)); rx = rng.integers(x.min(), x.max() + 1, len(z))
        gt = local_mean(v, z, y, x, 2); rnd = local_mean(v, rz, ry, rx, 2)
        # swapped: treat stored x as the row index and y as the column index (only where in bounds)
        ok = (x < Y) & (y < X)
        sw = local_mean(v, z[ok], x[ok], y[ok], 2) if ok.any() else np.array([np.nan])
        # the median voxel is 0 in this dataset (background was zeroed, section 1), so a
        # background-subtracted median ratio degenerates; use mean ratios and nonzero fractions instead
        ratios.append(float(gt.mean() / max(rnd.mean(), 1e-9)))
        ratios_sw.append(float(np.nanmean(sw) / max(rnd.mean(), 1e-9)))
        nonzero.append([float((gt > 0).mean()), float((rnd > 0).mean()), float((sw > 0).mean())])
        hit = 0
        for zi, yi, xi in zip(z, y, x):
            win = v[zi, max(0, yi - 5):yi + 6, max(0, xi - 5):xi + 6]
            py, px = np.unravel_index(np.argmax(win), win.shape)
            py += max(0, yi - 5); px += max(0, xi - 5)
            hit += (abs(py - yi) <= 3) and (abs(px - xi) <= 3)
        frac_local_max.append(hit / len(z))
    print(f"mean 5x5 brightness at ROI points / at random points in the same bounding box: {np.round(ratios, 1).tolist()}  "
          f"(swapped x<->y, in-bounds subset only: {np.round(ratios_sw, 2).tolist()})")
    print(f"fraction of 5x5 windows with any nonzero voxel [ROI, random, swapped]: {np.round(nonzero, 2).tolist()}")
    print(f"fraction of ROI points within 3 px of the local in-plane maximum: {np.round(frac_local_max, 2).tolist()}")
    aligned = bool(np.median(ratios) > 3)
    order_ok = bool(np.median(ratios) > np.nanmedian(ratios_sw))
    print(f"-> ROI points are much brighter than random points: {aligned}")
    print(f"-> labelled (x,y,z) <-> data axes (1,2,3) beats swapped order: {order_ok}")
    return dict(frames=ts, gt_over_random_contrast=ratios, swapped_xy_over_random_contrast=ratios_sw,
                frac_nonzero_window_gt_random_swapped=nonzero,
                frac_gt_within_3px_of_local_max=frac_local_max, aligned=aligned, xyz_order_confirmed=order_ok)


# ----------------------------------------------------------------------------
# 7. stream through every frame
# ----------------------------------------------------------------------------
def per_frame_means(f: h5py.File, n_frames: int, out_csv: Path, tracks: pd.DataFrame, stimuli: list[dict], rate_hz: float) -> dict:
    banner(f"7. PER-FRAME STATISTICS over all {n_frames} frames (streaming, one frame at a time)")
    rows, track_rows = [], []
    gmax = 0
    t0 = time.time()
    for t in range(n_frames):
        v = read_frame(f, t)
        gmax = max(gmax, int(v.max()))
        row = {"t": t, "mean": float(v.mean()), "p99": float(np.percentile(v, 99)), "max": int(v.max()), "frac_zero": float((v == 0).mean())}
        pts = tracks[tracks.t == t]
        if len(pts):
            ok = in_bounds(pts, v.shape)
            vals = np.full(len(pts), np.nan)
            vals[ok] = local_mean(v, pts.z.to_numpy()[ok], pts.y.to_numpy()[ok], pts.x.to_numpy()[ok], 2)
            row["gt_mean"] = float(np.nanmean(vals)); row["n_gt_out_of_bounds"] = int((~ok).sum())
            track_rows.append(pd.DataFrame({"t": t, "roi": pts.roi.to_numpy(), "label": pts.label.to_numpy() if "label" in pts else "", "gcamp": vals}))
        rows.append(row)
        if t % 100 == 0:
            print(f"  t={t}: " + ", ".join(f"{k}={v:.3g}" for k, v in row.items() if k != "t"))
    print(f"  streamed {n_frames} frames in {time.time() - t0:.0f} s; global max voxel value = {gmax}")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pf = pd.DataFrame(rows); pf.to_csv(out_csv, index=False); print(f"wrote {out_csv}")
    tp = out_csv.with_name(out_csv.stem + "_tracks.csv")
    mine = pd.concat(track_rows, ignore_index=True); mine.to_csv(tp, index=False); print(f"wrote {tp}")
    out = dict(global_max=gmax, mean_frac_zero=float(pf.frac_zero.mean()))

    # compare with the traces stored in the file
    if TRACE in f:
        stored = f[f"{TRACE}/data"][()]  # (T-1, N)
        sp = out_csv.with_name("stored_traces.csv")
        pd.DataFrame(stored, columns=[f"roi{i}" for i in range(stored.shape[1])]).rename_axis("t").to_csv(sp); print(f"wrote {sp}")
        wide = mine.pivot(index="t", columns="roi", values="gcamp").sort_index()
        n = min(len(wide), stored.shape[0])
        cors = []
        for j in range(min(wide.shape[1], stored.shape[1])):
            a = wide.iloc[:n, j].to_numpy(); ok = np.isfinite(a)
            cors.append(np.corrcoef(a[ok], stored[:n, j][ok])[0, 1] if ok.sum() > 10 else np.nan)
        # does the stored trace look like a simple rescaling of local brightness? ratio stored/mine per neuron
        out["stored_trace_comparison"] = dict(n_frames_compared=int(n), corr_median=float(np.nanmedian(cors)), corr_p10=float(np.nanpercentile(cors, 10)),
                                              corr_min=float(np.nanmin(cors)), frac_corr_above_0_8=float(np.mean(np.array(cors) > 0.8)),
                                              stored_value_range=[float(stored.min()), float(stored.max())], stored_mean_over_all=float(stored.mean()))
        sc = out["stored_trace_comparison"]
        print(f"  stored SignalRawFluor traces vs 5x5x1 local mean at ROI positions: per-neuron correlation median {sc['corr_median']:.2f}, "
              f"p10 {sc['corr_p10']:.2f}, min {sc['corr_min']:.2f}; {100*sc['frac_corr_above_0_8']:.0f}% of neurons > 0.8; stored values {sc['stored_value_range']}")

    # stimulus-locked responses
    if stimuli:
        labels = tracks.drop_duplicates("roi").set_index("roi")["label"] if "label" in tracks else None
        sources = {"local_mean": mine.pivot(index="t", columns="roi", values="gcamp").sort_index()}
        if TRACE in f:
            sources["stored_trace"] = pd.DataFrame(stored, columns=range(stored.shape[1]))
        resp = {}
        for src, wide in sources.items():
            # baseline floor: a 5x5 window that is all zero has no defined fractional change, so the
            # denominator is floored at a tenth of the typical neuron brightness
            floor = 0.1 * float(np.nanmedian(wide.to_numpy()))
            resp[src] = {}
            for s in stimuli:
                a = int(round(s["start_s"] * rate_hz)); b = int(round(s["stop_s"] * rate_hz)); post = int(round(10 * rate_hz)); pre = int(round(20 * rate_hz))
                base = wide.iloc[max(0, a - pre):a].mean(); during = wide.iloc[a:b + post].max()
                ratio = (during - base) / base.clip(lower=floor)
                top = ratio.sort_values(ascending=False).head(6)
                resp[src][s["stimulus"]] = dict(frames=[a, b], baseline_floor=floor,
                                                n_neurons_peak_above_50pct=int((ratio > 0.5).sum()), n_neurons_peak_above_100pct=int((ratio > 1.0).sum()),
                                                top_responders=[(int(r), labels.get(r, "") if labels is not None else "", round(float(v), 2)) for r, v in top.items()])
                rr = resp[src][s["stimulus"]]
                print(f"  [{src}] {s['stimulus']} ({s['start_s']}-{s['stop_s']} s = frames {a}-{b}): neurons whose peak in [onset, offset+10 s] exceeds the 20 s "
                      f"pre-stimulus mean by >50%: {rr['n_neurons_peak_above_50pct']}, by >100%: {rr['n_neurons_peak_above_100pct']}; "
                      f"top (roi, label, (peak-base)/base): {rr['top_responders']}")
        out["stimulus_responses"] = resp
    return out


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nwb", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--tracks-csv", type=Path, default=None, help="write the per-frame ROI positions as a long CSV here")
    ap.add_argument("--n-sample-frames", type=int, default=40)
    ap.add_argument("--n-neurons", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--per-frame-means", type=Path, default=None, help="CSV path; streams through ALL frames")
    ap.add_argument("--skip-pynwb", action="store_true")
    args = ap.parse_args(argv)

    summary = dict(nwb=str(args.nwb))
    with h5py.File(args.nwb, "r") as f:
        s = inspect_structure(f); summary["structure"] = s
        summary["metadata"] = inspect_metadata(f)
        if not args.skip_pynwb:
            summary["pynwb"] = validate_with_pynwb(args.nwb)
        tracks = load_calcium_tracks(f)
        if args.tracks_csv:
            args.tracks_csv.parent.mkdir(parents=True, exist_ok=True); tracks.to_csv(args.tracks_csv, index=False); print(f"\nwrote {args.tracks_csv}")
        n = s["n_frames"]
        summary["annotations"] = inspect_annotations(f, tracks, n)
        summary["bounds"] = inspect_bounds(tracks, tuple(s["calcium"]["shape"]))
        summary["channels"] = inspect_channels(f, tracks, n, args.n_sample_frames, args.n_neurons, args.seed)
        summary["alignment"] = inspect_alignment(f, tracks, args.seed)
        if args.per_frame_means:
            summary["per_frame"] = per_frame_means(f, n, args.per_frame_means, tracks, summary["metadata"]["stimuli"], s["calcium"]["rate_hz"])

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(jsonable(summary), indent=2))
        print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
