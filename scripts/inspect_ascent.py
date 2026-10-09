#!/usr/bin/env python3
"""Inspect the ASCENT in-house C. elegans recording (Zenodo 17561700).

Reads only selected frames from the HDF5 movie (never the whole file) and
reports:
  * HDF5 structure: frame count, channels, shape, dtype, chunking, attributes
  * which channel behaves like calcium (GCaMP) vs stable reference (mNeptune)
  * ground-truth CSV coverage: neurons, time points, track lengths
  * whether GT coordinates are in bounds and align with bright image spots

Usage:
  python -P scripts/inspect_ascent.py --h5 data/ascent/<movie>.h5 \
      --coords data/ascent/InHouse_Opterra_coords_GT.csv \
      --tracks data/ascent/InHouse_Opterra_tracks_GT.csv \
      --out-json previews/ascent/inspection_summary.json

Add --per-frame-means out.csv to stream through every frame once (reads the
whole movie sequentially, one frame at a time) and save per-frame mean
intensities, which reveal the stimulus timing.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

FRAME_RE = re.compile(r"^t(\d+)$")
CHAN_RE = re.compile(r"^c(\d+)$")


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def jsonable(o):
    """Convert numpy scalars/arrays to plain Python for json.dump."""
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def read_frame(f: h5py.File, t: int, c: int) -> np.ndarray:
    """Read one channel of one time point: a single Z x Y x X volume."""
    return f[f"t{t}/c{c}"][()]


def local_mean(vol: np.ndarray, z: np.ndarray, y: np.ndarray, x: np.ndarray, r: int = 2) -> np.ndarray:
    """Mean of a (2r+1)x(2r+1) in-plane window around each (z,y,x) point."""
    Z, Y, X = vol.shape
    out = np.empty(len(z), dtype=np.float64)
    for i, (zi, yi, xi) in enumerate(zip(z, y, x)):
        y0, y1 = max(0, yi - r), min(Y, yi + r + 1)
        x0, x1 = max(0, xi - r), min(X, xi + r + 1)
        out[i] = vol[zi, y0:y1, x0:x1].mean()
    return out


# ----------------------------------------------------------------------------
# 1. HDF5 structure
# ----------------------------------------------------------------------------
def inspect_structure(f: h5py.File) -> dict:
    banner("1. HDF5 STRUCTURE (verified from file)")
    root_keys = list(f.keys())
    frame_ids = sorted(int(m.group(1)) for k in root_keys if (m := FRAME_RE.match(k)))
    other = [k for k in root_keys if not FRAME_RE.match(k)]
    n = len(frame_ids)
    contiguous = frame_ids == list(range(n))
    print(f"root keys: {len(root_keys)}  (t-groups: {n}, other: {other})")
    print(f"frame ids contiguous 0..{n-1}: {contiguous}")

    # channels inside the first frame group
    t0 = f[f"t{frame_ids[0]}"]
    chan_keys = sorted(t0.keys())
    print(f"t{frame_ids[0]} members: {chan_keys}")
    chan_info = {}
    for k in chan_keys:
        d = t0[k]
        if isinstance(d, h5py.Dataset):
            chan_info[k] = dict(
                shape=d.shape,
                dtype=str(d.dtype),
                chunks=d.chunks,
                compression=d.compression,
                compression_opts=d.compression_opts,
                nbytes=int(d.nbytes),
            )
            print(f"  {k}: shape={d.shape} dtype={d.dtype} chunks={d.chunks} "
                  f"compression={d.compression} ({d.nbytes/1e6:.1f} MB raw)")

    # Full scan of metadata (shape/dtype) for every frame.  This touches only
    # the HDF5 object headers, not the pixel data, so it is cheap.
    shapes, dtypes, members = set(), set(), set()
    for t in frame_ids:
        g = f[f"t{t}"]
        members.add(tuple(sorted(g.keys())))
        for k in g.keys():
            d = g[k]
            shapes.add(d.shape)
            dtypes.add(str(d.dtype))
    print(f"all frames: member sets={members}")
    print(f"all frames: shapes={shapes} dtypes={dtypes}")

    # Axis-order inference: documented layout is (channels, z, y, x) =
    # (2, 9, 512, 512) per frame, with z stored first inside each channel.
    shape = next(iter(shapes))
    axis_guess = "unknown"
    if len(shape) == 3:
        small = int(np.argmin(shape))
        axis_guess = "ZYX" if small == 0 else ("YXZ" if small == 2 else "ambiguous")
    print(f"axis-order guess from shape {shape}: {axis_guess} "
          f"(smallest axis = depth; 'Z' is the one with few slices)")

    # Actual pixel-value spot checks on a few frames
    spot = {}
    for t in sorted({frame_ids[0], frame_ids[n // 2], frame_ids[-1]}):
        for k in chan_keys:
            v = f[f"t{t}/{k}"][()]
            spot[f"t{t}/{k}"] = dict(min=int(v.min()), max=int(v.max()), mean=float(v.mean()))
            print(f"  t{t}/{k}: min={v.min()} max={v.max()} mean={v.mean():.1f}")

    return dict(
        n_frames=n, frame_ids_contiguous=contiguous, other_root_keys=other,
        channel_keys=chan_keys, channel_info=chan_info, all_shapes=[list(s) for s in shapes],
        all_dtypes=sorted(dtypes), member_sets=[list(m) for m in members],
        axis_order_guess=axis_guess, spot_checks=spot,
    )


def inspect_attrs(f: h5py.File) -> dict:
    banner("2. HDF5 ATTRIBUTES / METADATA stored in the file")
    found = {}

    def dump(name, obj):
        if len(obj.attrs):
            found[name or "/"] = {k: jsonable(v if not isinstance(v, np.ndarray) else v) for k, v in obj.attrs.items()}

    dump("/", f)
    f["t0"].visititems(dump)  # one frame group is representative
    # also any non-frame root member
    for k in f.keys():
        if not FRAME_RE.match(k):
            f[k].visititems(dump) if isinstance(f[k], h5py.Group) else dump(k, f[k])
    if not found:
        print("no HDF5 attributes found on root, t0, or t0/c* "
              "-> voxel size, frame rate, channel names are NOT stored in the file")
    else:
        for k, v in found.items():
            print(f"  {k}: {v}")
    return found


# ----------------------------------------------------------------------------
# 3. channels: which one is calcium?
# ----------------------------------------------------------------------------
def inspect_channels(f: h5py.File, tracks: pd.DataFrame, n_frames: int, chan_keys: list[str],
                     n_sample: int, n_neurons: int, seed: int) -> dict:
    banner("3. CHANNEL BEHAVIOUR: which channel is calcium (GCaMP) vs reference (mNeptune)?")
    rng = np.random.default_rng(seed)
    ts = np.unique(np.linspace(0, n_frames - 1, n_sample).round().astype(int))
    full_tracks = tracks.groupby("TrackID").filter(lambda g: g["t"].nunique() == n_frames)
    tids = rng.choice(full_tracks["TrackID"].unique(), size=min(n_neurons, full_tracks["TrackID"].nunique()), replace=False)
    sel = full_tracks[full_tracks["TrackID"].isin(tids)].set_index(["t", "TrackID"]).sort_index()

    per_chan = {}
    for ck in chan_keys:
        c = int(CHAN_RE.match(ck).group(1))
        vol_mean, p50, p99, vmax = [], [], [], []
        neuron_ts = np.zeros((len(ts), len(tids)))
        for i, t in enumerate(ts):
            v = read_frame(f, int(t), c)
            vol_mean.append(v.mean()); p50.append(np.percentile(v, 50)); p99.append(np.percentile(v, 99)); vmax.append(v.max())
            pts = sel.loc[int(t)].loc[tids]
            neuron_ts[i] = local_mean(v, pts["z"].to_numpy(), pts["y"].to_numpy(), pts["x"].to_numpy(), r=2)
        vol_mean = np.array(vol_mean)
        # per-neuron temporal variability (std / mean) after removing a linear
        # trend (bleaching), then the median over neurons
        cv = []
        for j in range(len(tids)):
            y = neuron_ts[:, j]
            trend = np.polyval(np.polyfit(ts, y, 1), ts)
            cv.append(np.std(y - trend) / max(y.mean(), 1e-9))
        cv = np.array(cv)
        per_chan[ck] = dict(
            sampled_frames=ts.tolist(),
            volume_mean=vol_mean.tolist(), p50=np.array(p50).tolist(), p99=np.array(p99).tolist(), max=np.array(vmax).tolist(),
            volume_mean_temporal_cv=float(vol_mean.std() / vol_mean.mean()),
            neuron_detrended_cv_median=float(np.median(cv)),
            neuron_detrended_cv_p90=float(np.percentile(cv, 90)),
            neuron_mean_intensity_median=float(np.median(neuron_ts)),
            neuron_ts=neuron_ts,
        )
        print(f"{ck}: volume mean over {len(ts)} sampled frames: {vol_mean.min():.1f}..{vol_mean.max():.1f} "
              f"(temporal CV {per_chan[ck]['volume_mean_temporal_cv']:.3f}); p99 {np.min(p99):.0f}..{np.max(p99):.0f}; "
              f"max {np.min(vmax)}..{np.max(vmax)}")
        print(f"{ck}: per-neuron detrended temporal CV over {len(tids)} full-length tracks: "
              f"median {per_chan[ck]['neuron_detrended_cv_median']:.3f}, p90 {per_chan[ck]['neuron_detrended_cv_p90']:.3f}; "
              f"median neuron brightness {per_chan[ck]['neuron_mean_intensity_median']:.0f}")

    # cross-channel comparison
    result = dict(sampled_frames=ts.tolist(), track_ids=tids.tolist(), per_channel={})
    for ck, d in per_chan.items():
        result["per_channel"][ck] = {k: v for k, v in d.items() if k != "neuron_ts"}
    if len(per_chan) == 2:
        a, b = chan_keys
        cva, cvb = per_chan[a]["neuron_detrended_cv_median"], per_chan[b]["neuron_detrended_cv_median"]
        calcium = a if cva > cvb else b
        ref = b if calcium == a else a
        corr = np.corrcoef(per_chan[a]["neuron_ts"].mean(0), per_chan[b]["neuron_ts"].mean(0))[0, 1]
        # fraction of neurons whose c-calcium trace varies more than their reference trace
        nts_a, nts_b = per_chan[a]["neuron_ts"], per_chan[b]["neuron_ts"]
        frac = float(np.mean(nts_a.std(0) / nts_a.mean(0) > nts_b.std(0) / nts_b.mean(0)))
        result.update(calcium_channel_guess=calcium, reference_channel_guess=ref,
                      cv_ratio=float(max(cva, cvb) / max(min(cva, cvb), 1e-9)),
                      cross_channel_brightness_corr_over_neurons=float(corr),
                      frac_neurons_with_more_variability_in=dict(channel=a, fraction=frac))
        print(f"-> more temporally variable channel (calcium-like): {calcium}; "
              f"more stable channel (reference-like): {ref}; CV ratio {result['cv_ratio']:.2f}")
        print(f"-> correlation of mean neuron brightness across channels (same spots bright in both?): {corr:.2f}")
        print(f"-> fraction of sampled neurons more variable in {a} than {b}: {frac:.2f}")
        print("   (documented: c0 = mNeptune reference, c1 = GCaMP6 calcium)")
    return result


# ----------------------------------------------------------------------------
# 4. annotations
# ----------------------------------------------------------------------------
def inspect_annotations(coords: pd.DataFrame, tracks: pd.DataFrame, n_frames: int | None) -> dict:
    banner("4. GROUND-TRUTH ANNOTATIONS (CSV)")
    out = {}
    for name, df in [("coords_GT", coords), ("tracks_GT", tracks)]:
        per_t = df.groupby("t").size()
        d = dict(
            columns=list(df.columns), dtypes={c: str(t) for c, t in df.dtypes.items()}, n_rows=len(df),
            t_min=int(df.t.min()), t_max=int(df.t.max()), n_time_points=int(df.t.nunique()),
            missing_t_in_range=int(df.t.max() - df.t.min() + 1 - df.t.nunique()),
            points_per_frame=dict(min=int(per_t.min()), median=float(per_t.median()), max=int(per_t.max())),
            integer_coords={c: bool(np.array_equal(df[c], np.round(df[c]))) for c in "zyx"},
            ranges={c: [int(df[c].min()), int(df[c].max())] for c in "zyx"},
        )
        print(f"{name}: {len(df)} rows, columns {list(df.columns)}, dtypes {set(d['dtypes'].values())}")
        print(f"  t: {d['t_min']}..{d['t_max']} ({d['n_time_points']} distinct, {d['missing_t_in_range']} missing in range)")
        print(f"  points per frame: {d['points_per_frame']}")
        print(f"  coords integer-valued: {d['integer_coords']}; ranges z{d['ranges']['z']} y{d['ranges']['y']} x{d['ranges']['x']}")
        out[name] = d

    oid = coords["object_id"]
    out["coords_GT"].update(n_unique_object_id=int(oid.nunique()), object_id_max=int(oid.max()),
                            duplicate_object_ids=int(oid.duplicated().sum()))
    print(f"coords_GT: {oid.nunique()} unique object_id (max id {oid.max()}, so ids are sparse / not row numbers)")

    tl = tracks.groupby("TrackID")["t"].nunique()
    n_t = tracks.t.nunique()
    out["tracks_GT"].update(
        n_tracks=int(tracks.TrackID.nunique()), track_id_min=int(tracks.TrackID.min()), track_id_max=int(tracks.TrackID.max()),
        track_length=dict(min=int(tl.min()), median=float(tl.median()), max=int(tl.max())),
        n_tracks_spanning_all_annotated_frames=int((tl == n_t).sum()),
        duplicate_track_time_pairs=int(tracks.duplicated(["TrackID", "t"]).sum()),
        tracks_shorter_than_full={int(k): int(v) for k, v in tl[tl < n_t].items()},
    )
    print(f"tracks_GT: {tracks.TrackID.nunique()} tracks (ids {tracks.TrackID.min()}..{tracks.TrackID.max()}), "
          f"length min/median/max {tl.min()}/{tl.median()}/{tl.max()}, "
          f"{(tl == n_t).sum()} span all {n_t} annotated frames; shorter tracks: {out['tracks_GT']['tracks_shorter_than_full']}")
    print(f"  duplicate (TrackID, t) pairs: {out['tracks_GT']['duplicate_track_time_pairs']}")

    m = tracks.merge(coords, left_on="ObjectID", right_on="object_id", suffixes=("_tr", "_co"))
    agree = bool(len(m) == len(tracks) == len(coords)
                 and np.array_equal(m[["t_tr", "z_tr", "y_tr", "x_tr"]].to_numpy(), m[["t_co", "z_co", "y_co", "x_co"]].to_numpy()))
    out["coords_tracks_consistent"] = agree
    print(f"coords_GT and tracks_GT describe the same points (ObjectID == object_id, identical t,z,y,x): {agree}")

    if n_frames is not None:
        out["fraction_of_movie_frames_annotated"] = float(n_t / n_frames)
        print(f"annotated time points / movie frames: {n_t}/{n_frames} = {n_t / n_frames:.3f}")

    # motion summary: per-step displacement of each track (how much do neurons move?)
    ts = tracks.sort_values(["TrackID", "t"])
    step = ts.groupby("TrackID")[["z", "y", "x"]].diff().dropna()
    disp = np.sqrt((step[["y", "x"]] ** 2).sum(axis=1))
    out["per_step_motion_voxels"] = dict(
        xy_median=float(disp.median()), xy_p99=float(disp.quantile(0.99)), xy_max=float(disp.max()),
        z_max=int(step.z.abs().max()), frac_steps_zero=float((disp == 0).mean()),
    )
    print(f"motion between consecutive frames (voxels): in-plane median {disp.median():.2f}, p99 {disp.quantile(0.99):.1f}, "
          f"max {disp.max():.1f}; max |dz| {int(step.z.abs().max())}; {100*(disp == 0).mean():.0f}% of steps have zero in-plane motion")
    return out


# ----------------------------------------------------------------------------
# 5. bounds + alignment
# ----------------------------------------------------------------------------
def inspect_bounds(coords: pd.DataFrame, shape: tuple[int, int, int]) -> dict:
    banner("5. ARE GT COORDINATES IN BOUNDS?")
    Z, Y, X = shape
    checks = {}
    for name, (zc, yc, xc) in {"as_labelled_zyx": ("z", "y", "x"), "swapped_xy": ("z", "x", "y")}.items():
        inb = ((coords[zc] >= 0) & (coords[zc] < Z) & (coords[yc] >= 0) & (coords[yc] < Y) & (coords[xc] >= 0) & (coords[xc] < X))
        checks[name] = dict(n_out_of_bounds=int((~inb).sum()), all_in_bounds=bool(inb.all()))
    print(f"volume shape (Z,Y,X)={shape}; coords ranges z[{coords.z.min()},{coords.z.max()}] "
          f"y[{coords.y.min()},{coords.y.max()}] x[{coords.x.min()},{coords.x.max()}]")
    for k, v in checks.items():
        print(f"  {k}: out-of-bounds rows = {v['n_out_of_bounds']}")
    print("  note: with a 512x512 square field, bounds alone cannot tell x from y -> see alignment test")
    return dict(shape=list(shape), checks=checks)


def inspect_alignment(f: h5py.File, coords: pd.DataFrame, n_frames: int, chan_keys: list[str], seed: int) -> dict:
    banner("6. DO GT POINTS LAND ON BRIGHT SPOTS? (alignment / axis-order test)")
    rng = np.random.default_rng(seed)
    ts = sorted({0, n_frames // 4, n_frames // 2, (3 * n_frames) // 4, n_frames - 1})
    results = {}
    for ck in chan_keys:
        c = int(CHAN_RE.match(ck).group(1))
        ratios, ratios_sw, frac_local_max = [], [], []
        for t in ts:
            v = read_frame(f, t, c).astype(np.float64)
            pts = coords[coords.t == t]
            z, y, x = pts.z.to_numpy(), pts.y.to_numpy(), pts.x.to_numpy()
            # random points drawn inside the same bounding box as the GT points, so
            # we compare against background from the same region of the worm
            rz = rng.integers(z.min(), z.max() + 1, len(z)); ry = rng.integers(y.min(), y.max() + 1, len(z)); rx = rng.integers(x.min(), x.max() + 1, len(z))
            gt = local_mean(v, z, y, x, 2); rnd = local_mean(v, rz, ry, rx, 2); sw = local_mean(v, z, x, y, 2)
            bg = np.median(v)
            ratios.append((np.median(gt) - bg) / max(np.median(rnd) - bg, 1e-9))
            ratios_sw.append((np.median(sw) - bg) / max(np.median(rnd) - bg, 1e-9))
            # is the GT point within 3 px of the brightest pixel of its 11x11 in-plane neighbourhood?
            hit = 0
            for zi, yi, xi in zip(z, y, x):
                win = v[zi, max(0, yi - 5):yi + 6, max(0, xi - 5):xi + 6]
                py, px = np.unravel_index(np.argmax(win), win.shape)
                py += max(0, yi - 5); px += max(0, xi - 5)
                hit += (abs(py - yi) <= 3) and (abs(px - xi) <= 3)
            frac_local_max.append(hit / len(z))
        results[ck] = dict(frames=ts, gt_over_random_contrast=ratios, swapped_xy_over_random_contrast=ratios_sw,
                           frac_gt_within_3px_of_local_max=frac_local_max)
        print(f"{ck}: background-subtracted brightness at GT points / at random points: "
              f"{np.round(ratios, 1).tolist()}  (swapped x<->y: {np.round(ratios_sw, 1).tolist()})")
        print(f"{ck}: fraction of GT points within 3 px of the local in-plane maximum: {np.round(frac_local_max, 2).tolist()}")
    aligned = all(np.median(r["gt_over_random_contrast"]) > 3 for r in results.values())
    order_ok = all(np.median(r["gt_over_random_contrast"]) > np.median(r["swapped_xy_over_random_contrast"]) for r in results.values())
    print(f"-> GT points are much brighter than random points in every channel: {aligned}")
    print(f"-> labelled (z,y,x) order beats swapped order: {order_ok}")
    return dict(per_channel=results, aligned=aligned, zyx_order_confirmed=order_ok)


# ----------------------------------------------------------------------------
# optional: stream through every frame once
# ----------------------------------------------------------------------------
def per_frame_means(f: h5py.File, n_frames: int, chan_keys: list[str], out_csv: Path, tracks: pd.DataFrame | None) -> None:
    """Stream through every frame once.  Writes per-frame volume statistics to
    out_csv and, if tracks are given, the brightness of every GT track in every
    channel to <out_csv stem>_tracks.csv (long format: t, TrackID, c0, c1, ...)."""
    banner(f"7. PER-FRAME MEANS over all {n_frames} frames (streaming, one frame at a time)")
    rows, track_rows = [], []
    for t in range(n_frames):
        row = {"t": t}
        pts = tracks[tracks.t == t] if tracks is not None else None
        per_track = None
        for ck in chan_keys:
            v = f[f"t{t}/{ck}"][()]
            row[f"{ck}_mean"] = float(v.mean())
            row[f"{ck}_p99"] = float(np.percentile(v, 99))
            if pts is not None and len(pts):
                vals = local_mean(v, pts.z.to_numpy(), pts.y.to_numpy(), pts.x.to_numpy(), 2)
                row[f"{ck}_gt_mean"] = float(vals.mean())
                if per_track is None:
                    per_track = pd.DataFrame({"t": t, "TrackID": pts.TrackID.to_numpy()})
                per_track[ck] = vals
        rows.append(row)
        if per_track is not None:
            track_rows.append(per_track)
        if t % 100 == 0:
            print(f"  t={t}: " + ", ".join(f"{k}={v:.1f}" for k, v in row.items() if k != "t"))
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")
    if track_rows:
        tp = out_csv.with_name(out_csv.stem + "_tracks.csv")
        pd.concat(track_rows, ignore_index=True).to_csv(tp, index=False)
        print(f"wrote {tp}")


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5", type=Path, required=True)
    ap.add_argument("--coords", type=Path, required=True)
    ap.add_argument("--tracks", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=None, help="write a JSON summary here")
    ap.add_argument("--n-sample-frames", type=int, default=40, help="frames sampled for channel statistics")
    ap.add_argument("--n-neurons", type=int, default=40, help="tracks sampled for per-neuron time series")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--per-frame-means", type=Path, default=None, help="CSV path; streams through ALL frames (slow)")
    args = ap.parse_args(argv)

    coords = pd.read_csv(args.coords)
    tracks = pd.read_csv(args.tracks)
    summary = dict(h5=str(args.h5), coords=str(args.coords), tracks=str(args.tracks))

    with h5py.File(args.h5, "r") as f:
        s = inspect_structure(f); summary["structure"] = s
        summary["attrs"] = inspect_attrs(f)
        n = s["n_frames"]; chans = s["channel_keys"]; shape = tuple(s["all_shapes"][0])
        summary["annotations"] = inspect_annotations(coords, tracks, n)
        summary["bounds"] = inspect_bounds(coords, shape)
        summary["channels"] = inspect_channels(f, tracks, n, chans, args.n_sample_frames, args.n_neurons, args.seed)
        summary["alignment"] = inspect_alignment(f, coords, n, chans, args.seed)
        if args.per_frame_means:
            per_frame_means(f, n, chans, args.per_frame_means, tracks)

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(jsonable(summary), indent=2))
        print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
