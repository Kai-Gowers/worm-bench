#!/usr/bin/env python3
"""Inspect one calcium-bearing NWB recording from the Kato-lab dandisets used by
WormID-Bench: SK1 (DANDI 000565) or SK2 (DANDI 000472).

The files follow the same ndx-multichannel-volume layout as the EY dataset
(DANDI 000541), so most checks are imported from scripts/inspect_ey.py:

  acquisition/CalciumImageSeries/data            (T, X, Y, Z) uint16, one GCaMP channel
  processing/CalciumActivity/CalciumSeriesSegmentation/Seg_tpoint_{t}/voxel_mask
                                                 one (x, y, z, weight) row per ROI per frame
  processing/CalciumActivity/NeuronIDs/labels    biological name per ROI (many empty)
  processing/CalciumActivity/SignalDFoF/...      stored DF/F traces (T, N)   <- differs from EY
  acquisition/NeuroPALImageRaw/data              (X, Y, Z, C) structural stack, 6-7 channels
  processing/ProcessedImage/ProcessedImage/data  filtered + histogram-matched NeuroPAL stack
  (no intervals group: no stimulus timing stored)

Extra checks specific to these files: whether NeuronIDs/ImageSegmentation is a
second copy of the per-frame segmentation, the NeuroPAL channel -> RGBW mapping,
the stored DF/F traces against our own local-mean traces, and a point-based
time-index test (does segmentation t fit image t better than t-1 / t+1?).

Only a few dozen frames are read unless --per-frame-means is given, in which
case every frame is streamed once (one frame in memory at a time).

Usage:
  python -P scripts/inspect_sk.py --nwb data/sk1/<file>.nwb \
      --out-json previews/sk1/inspection_summary.json \
      [--tracks-csv previews/sk1/calcium_tracks.csv] \
      [--per-frame-means previews/sk1/per_frame_means.csv]
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
from inspect_ey import (CALCIUM, NEUROPAL_IMG, NEUROPAL_SEG, SEG, TPOINT_RE, decode, in_bounds, inspect_alignment,  # noqa: E402
                        inspect_annotations, inspect_channels, inspect_metadata, inspect_structure, load_calcium_tracks,
                        ragged_labels, read_frame, validate_with_pynwb)

TRACE = "processing/CalciumActivity/SignalDFoF/SignalCalciumImResponseSeries"
SEG_COPY = "processing/CalciumActivity/NeuronIDs/ImageSegmentation"
PROCESSED = "processing/ProcessedImage/ProcessedImage"
NEUROPAL_VOL = "general/optophysiology/NeuroPALImVol"


# ----------------------------------------------------------------------------
# 1b. structure details that differ from EY
# ----------------------------------------------------------------------------
def inspect_structure_extra(f: h5py.File) -> dict:
    banner("1b. SK-SPECIFIC STRUCTURE (stored DF/F, segmentation copy, processed NeuroPAL stack)")
    out = {}
    if TRACE in f:
        td = f[f"{TRACE}/data"]
        out["SignalDFoF"] = dict(shape=list(td.shape), dtype=str(td.dtype), unit=decode(td.attrs.get("unit", "")),
                                 rate_hz=float(f[f"{TRACE}/starting_time"].attrs["rate"]),
                                 description=decode(f[TRACE].attrs.get("description", "")),
                                 rois_table=f[f[f"{TRACE}/rois"].attrs["table"]].name if "rois" in f[TRACE] else None,
                                 n_nan=int(np.isnan(td[()]).sum()), min=float(np.nanmin(td[()])), max=float(np.nanmax(td[()])),
                                 median=float(np.nanmedian(td[()])), frac_rows_all_zero=float((td[()] == 0).all(axis=1).mean()))
        s = out["SignalDFoF"]
        print(f"SignalDFoF/SignalCalciumImResponseSeries: shape {s['shape']} {s['dtype']} unit '{s['unit']}' rate {s['rate_hz']} Hz; "
              f"values {s['min']:.3f}..{s['max']:.3f} (median {s['median']:.3f}), {s['n_nan']} NaN, {100*s['frac_rows_all_zero']:.1f}% all-zero rows; "
              f"rois -> {s['rois_table']}; '{s['description']}'")
        print("  (EY stored raw fluorescence; here only a DF/F series is stored, so absolute brightness per neuron is not in the traces)")
    else:
        print("no SignalDFoF series"); out["SignalDFoF"] = None

    # NeuronIDs/ImageSegmentation: same objects as CalciumSeriesSegmentation (HDF5 hard links) or a copy?
    if SEG_COPY in f:
        keys = [k for k in f[SEG_COPY].keys() if TPOINT_RE.match(k)]
        same_obj = all(f[f"{SEG_COPY}/{k}"].id == f[f"{SEG}/{k}"].id for k in keys[:50]) if keys else None
        same_vals = all(np.array_equal(f[f"{SEG_COPY}/{k}/voxel_mask"][()], f[f"{SEG}/{k}/voxel_mask"][()]) for k in keys[::max(1, len(keys)//20)]) if keys else None
        out["neuronids_image_segmentation"] = dict(n_tpoints=len(keys), hard_links_to_main=same_obj, same_voxel_masks_sampled=same_vals)
        print(f"NeuronIDs/ImageSegmentation: {len(keys)} Seg_tpoint groups; hard links to CalciumSeriesSegmentation (first 50 checked): {same_obj}; "
              f"identical voxel masks on a 1-in-{max(1, len(keys)//20)} sample: {same_vals}")
    if PROCESSED in f:
        d = f[f"{PROCESSED}/data"]
        out["processed_image"] = dict(shape=list(d.shape), dtype=str(d.dtype), RGBW_channels=f[f"{PROCESSED}/RGBW_channels"][()].tolist(),
                                      description=decode(f[f"{PROCESSED}/description"][()]),
                                      channel_order=[decode(c) for c in f["general/optophysiology/ProcessedImVol/order_optical_channels/channels"][()]]
                                      if "general/optophysiology/ProcessedImVol/order_optical_channels/channels" in f else None)
        p = out["processed_image"]
        print(f"ProcessedImage: shape {p['shape']} {p['dtype']} RGBW {p['RGBW_channels']} channels {p['channel_order']}; '{p['description']}'")
    if NEUROPAL_IMG in f:
        order = [decode(c) for c in f[f"{NEUROPAL_VOL}/order_optical_channels/channels"][()]]
        rgbw = f[f"{NEUROPAL_IMG}/RGBW_channels"][()].tolist()
        out["neuropal_rgbw_names"] = dict(zip("RGBW", [order[i] for i in rgbw]))
        out["neuropal_channel_order"] = order
        print(f"NeuroPALImageRaw channel order {order}; RGBW_channels {rgbw} -> {out['neuropal_rgbw_names']}")
    d = f[f"{CALCIUM}/data"]
    out["calcium_dimension_vs_shape"] = dict(dimension=f[f"{CALCIUM}/dimension"][()].tolist(), data_shape_per_frame=list(d.shape[1:]))
    print(f"CalciumImageSeries `dimension` = {out['calcium_dimension_vs_shape']['dimension']} vs per-frame data shape {list(d.shape[1:])} "
          f"(dimension lists the axes as Y, X, Z here; the data are stored X, Y, Z)")
    print(f"calcium optical channels: {[decode(c) for c in f['general/optophysiology/CalciumImVol/order_optical_channels/channels'][()]]}")
    return out


# ----------------------------------------------------------------------------
# 5. bounds (EY's version prints EY-specific text)
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
    print(f"  note: X={X} and Y={Y} differ, so bounds alone already discriminate x from y; section 6 confirms with image content")
    return dict(shape=list(data_shape), checks=checks)


# ----------------------------------------------------------------------------
# 6b. time-index test with point annotations
# ----------------------------------------------------------------------------
def time_index_test(f: h5py.File, tracks: pd.DataFrame, n_test: int = 24) -> dict:
    banner("6b. TIME-INDEX TEST: does Seg_tpoint_t fit image t better than t-1 / t+1?")
    ts_ = tracks.sort_values(["roi", "t"]); grp = ts_.groupby("roi")
    step = grp[["x", "y"]].diff(); dt = grp["t"].diff()
    disp = np.sqrt((step[dt == 1] ** 2).sum(axis=1))
    motion = pd.DataFrame({"t": ts_.t[dt == 1], "d": disp}).groupby("t")["d"].median()
    out = dict(population_median_step_voxels=dict(median=float(motion.median()), p90=float(motion.quantile(0.9)), max=float(motion.max()),
                                                 n_frames_above_1_voxel=int((motion > 1).sum()), n_frames_above_2_voxels=int((motion > 2).sum())))
    m = out["population_median_step_voxels"]
    print(f"population median in-plane step per frame: median {m['median']:.2f}, p90 {m['p90']:.2f}, max {m['max']:.1f} voxels; "
          f"{m['n_frames_above_1_voxel']} frames > 1 voxel, {m['n_frames_above_2_voxels']} frames > 2 voxels")
    cand = motion[(motion.index > 0) & (motion.index < tracks.t.max())].sort_values(ascending=False)
    cand = cand[cand > 0].head(n_test)
    if len(cand) == 0:
        print("  no frame with population motion > 0 voxels: the test cannot discriminate (worm is immobilized)")
        out["result"] = "inconclusive: no motion"; return out
    byt = {t: g for t, g in tracks.groupby("t")}
    wins, margins, rows = 0, [], []
    for t in cand.index:
        v = read_frame(f, int(t)).astype(np.float64)
        means = {}
        for off in (-1, 0, 1):
            p = byt[int(t) + off]; ok = in_bounds(p, v.shape)
            means[off] = float(local_mean(v, p.z.to_numpy()[ok], p.y.to_numpy()[ok], p.x.to_numpy()[ok], 2).mean())
        best = max(means, key=means.get)
        wins += best == 0
        margins.append((means[0] - max(means[-1], means[1])) / means[0])
        rows.append(dict(t=int(t), motion=float(cand[t]), mean_tm1=means[-1], mean_t=means[0], mean_tp1=means[1]))
    out["frames"] = rows; out["n_tested"] = len(rows); out["n_wins_for_t"] = int(wins); out["median_margin"] = float(np.median(margins))
    print(f"  on the {len(rows)} frames with the largest population motion (median step {cand.min():.2f}..{cand.max():.2f} voxels):")
    print(f"  positions from Seg_tpoint_t give the brightest 5x5 local mean on image t in {wins}/{len(rows)} cases; "
          f"median margin over the better neighbour {100*out['median_margin']:+.1f}%")
    if cand.max() < 1:
        print("  NOTE: the motion is below one voxel, so a margin near zero is expected and the test is weak")
    return out


# ----------------------------------------------------------------------------
# 6c. do the annotations keep sitting on neurons over time? (drift / bleaching test)
# ----------------------------------------------------------------------------
def _xy_shift(a: np.ndarray, b: np.ndarray, max_shift: int = 30) -> tuple[int, int, float]:
    """Integer (dy, dx) that best aligns b to a by normalised cross-correlation of mean-subtracted images."""
    from scipy.signal import fftconvolve
    a = a - a.mean(); b = b - b.mean()
    cc = fftconvolve(a, b[::-1, ::-1], mode="same")
    cy, cx = np.array(cc.shape) // 2
    win = cc[cy - max_shift:cy + max_shift + 1, cx - max_shift:cx + max_shift + 1]
    iy, ix = np.unravel_index(np.argmax(win), win.shape)
    peak = win[iy, ix] / np.sqrt((a ** 2).sum() * (b ** 2).sum())
    return int(iy - max_shift), int(ix - max_shift), float(peak)


def drift_test(f: h5py.File, tracks: pd.DataFrame, n_frames: int, n_sample: int = 12) -> dict:
    """On sampled frames: image-based global shift relative to frame 0 (smoothed MIP cross-correlation,
    per-slice z profile), versus the shift of the stored annotations, and whether the ROI points still
    sit on local maxima of the smoothed image."""
    from scipy.ndimage import gaussian_filter
    banner("6c. DO THE ANNOTATIONS FOLLOW THE NEURONS OVER TIME? (image drift vs annotation drift, local-max hit rate)")
    ts = np.unique(np.linspace(0, n_frames - 1, n_sample).round().astype(int))
    byt = {t: g for t, g in tracks.groupby("t")}
    p0 = byt[0].set_index("roi")

    def smooth(v):
        return np.stack([gaussian_filter(v[z].astype(np.float32), 1.5) for z in range(v.shape[0])])

    v0 = smooth(read_frame(f, 0)); mip0 = v0.max(0); zprof0 = v0.mean(axis=(1, 2))
    rows = []
    for t in ts:
        v = smooth(read_frame(f, int(t))); mip = v.max(0)
        dy, dx, peak = _xy_shift(mip0, mip)
        # z: shift of the per-slice mean profile that best matches frame 0 (coarse, +-3 slices)
        zp = v.mean(axis=(1, 2)); best = max(range(-3, 4), key=lambda s: np.corrcoef(zprof0[max(0, s):len(zp) + min(0, s)], zp[max(0, -s):len(zp) + min(0, -s)])[0, 1] if len(zp) - abs(s) > 3 else -1)
        p = byt[int(t)].set_index("roi").loc[p0.index]
        ann = (p[["x", "y"]] - p0[["x", "y"]]).median()
        bg = float(np.percentile(v, 50))
        hit, contrast = 0, []
        for zi, yi, xi in zip(p.z, p.y, p.x):
            win = v[zi, max(0, yi - 5):yi + 6, max(0, xi - 5):xi + 6]
            py, px = np.unravel_index(np.argmax(win), win.shape); py += max(0, yi - 5); px += max(0, xi - 5)
            hit += (abs(py - yi) <= 3) and (abs(px - xi) <= 3)
            contrast.append(v[zi, yi, xi] - bg)
        rows.append(dict(t=int(t), image_shift_dx=dx, image_shift_dy=dy, image_shift_peak=round(peak, 3), image_zshift=int(best),
                         annotation_shift_dx=float(ann.x), annotation_shift_dy=float(ann.y),
                         frac_within_3px_of_smoothed_local_max=hit / len(p), roi_contrast_median=float(np.median(contrast)),
                         roi_contrast_p90=float(np.percentile(contrast, 90)), frac_roi_contrast_above_5=float(np.mean(np.array(contrast) > 5))))
        r = rows[-1]
        print(f"  t={r['t']:5d}: image shift vs t0 (dx,dy,dz)=({dx:+d},{dy:+d},{best:+d}) ncc {peak:.2f}; annotation median shift (dx,dy)=({ann.x:+.1f},{ann.y:+.1f}); "
              f"ROI points within 3 px of smoothed local max {100*r['frac_within_3px_of_smoothed_local_max']:.0f}%; "
              f"smoothed ROI brightness above background: median {r['roi_contrast_median']:.1f}, p90 {r['roi_contrast_p90']:.1f} counts, {100*r['frac_roi_contrast_above_5']:.0f}% > 5")
    out = dict(frames=rows)
    hits = [r["frac_within_3px_of_smoothed_local_max"] for r in rows]
    out["hit_rate_first"] = hits[0]; out["hit_rate_last"] = hits[-1]; out["hit_rate_min"] = min(hits)
    out["max_abs_image_shift_xy"] = int(max(max(abs(r["image_shift_dx"]), abs(r["image_shift_dy"])) for r in rows))
    print(f"  -> local-max hit rate first/last sampled frame {100*hits[0]:.0f}% / {100*hits[-1]:.0f}%; largest image shift vs t0 {out['max_abs_image_shift_xy']} px")
    return out


# ----------------------------------------------------------------------------
# 7. stream through every frame; compare with the stored DF/F
# ----------------------------------------------------------------------------
def per_frame_means(f: h5py.File, n_frames: int, out_csv: Path, tracks: pd.DataFrame) -> dict:
    banner(f"7. PER-FRAME STATISTICS over all {n_frames} frames (streaming, one frame at a time)")
    rows, track_rows = [], []
    gmax = 0
    t0 = time.time()
    for t in range(n_frames):
        v = read_frame(f, t)
        gmax = max(gmax, int(v.max()))
        row = {"t": t, "mean": float(v.mean()), "p50": float(np.percentile(v, 50)), "p99": float(np.percentile(v, 99)), "max": int(v.max()),
               "frac_zero": float((v == 0).mean()), "frac_at_max_code": float((v >= 4095).mean())}
        pts = tracks[tracks.t == t]
        if len(pts):
            ok = in_bounds(pts, v.shape)
            vals = np.full(len(pts), np.nan)
            vals[ok] = local_mean(v, pts.z.to_numpy()[ok], pts.y.to_numpy()[ok], pts.x.to_numpy()[ok], 2)
            row["gt_mean"] = float(np.nanmean(vals)); row["n_gt_out_of_bounds"] = int((~ok).sum())
            track_rows.append(pd.DataFrame({"t": t, "roi": pts.roi.to_numpy(), "label": pts.label.to_numpy() if "label" in pts else "", "gcamp": vals}))
        rows.append(row)
        if t % 100 == 0:
            print(f"  t={t}: " + ", ".join(f"{k}={v:.4g}" for k, v in row.items() if k != "t"))
    print(f"  streamed {n_frames} frames in {time.time() - t0:.0f} s; global max voxel value = {gmax}")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pf = pd.DataFrame(rows); pf.to_csv(out_csv, index=False); print(f"wrote {out_csv}")
    tp = out_csv.with_name(out_csv.stem + "_tracks.csv")
    mine = pd.concat(track_rows, ignore_index=True); mine.to_csv(tp, index=False); print(f"wrote {tp}")
    out = dict(global_max=gmax, mean_frac_zero=float(pf.frac_zero.mean()), mean_frac_at_4095=float(pf.frac_at_max_code.mean()),
               gt_mean_first_over_median=float(pf.gt_mean.iloc[0] / pf.gt_mean.median()), gt_mean_last_over_median=float(pf.gt_mean.iloc[-1] / pf.gt_mean.median()),
               volume_mean_first_over_median=float(pf["mean"].iloc[0] / pf["mean"].median()), volume_mean_last_over_median=float(pf["mean"].iloc[-1] / pf["mean"].median()))
    print(f"  brightness at ROI positions: first frame / median = {out['gt_mean_first_over_median']:.2f}, last / median = {out['gt_mean_last_over_median']:.2f} "
          f"(whole volume: {out['volume_mean_first_over_median']:.2f}, {out['volume_mean_last_over_median']:.2f})")

    if TRACE in f:
        stored = f[f"{TRACE}/data"][()]  # (T, N) DF/F
        sp = out_csv.with_name("stored_traces.csv")
        pd.DataFrame(stored, columns=[f"roi{i}" for i in range(stored.shape[1])]).rename_axis("t").to_csv(sp); print(f"wrote {sp}")
        wide = mine.pivot(index="t", columns="roi", values="gcamp").sort_index()
        n = min(len(wide), stored.shape[0])
        cors, f0_est = [], []
        for j in range(min(wide.shape[1], stored.shape[1])):
            a = wide.iloc[:n, j].to_numpy(); s = stored[:n, j]; ok = np.isfinite(a) & np.isfinite(s)
            cors.append(np.corrcoef(a[ok], s[ok])[0, 1] if ok.sum() > 10 else np.nan)
            # if stored = (F - F0)/F0 with F ~ k * ours, then ours = k*F0*(1 + stored): slope/intercept gives k*F0 and the fit quality
            if ok.sum() > 10:
                sl, ic = np.polyfit(s[ok], a[ok], 1); f0_est.append(ic / max(sl, 1e-9) if sl > 0 else np.nan)
        cors = np.array(cors)
        # DF/F baseline convention: what fraction of the stored trace is below zero, where does the baseline sit
        out["stored_trace_comparison"] = dict(n_frames_compared=int(n), corr_median=float(np.nanmedian(cors)), corr_p10=float(np.nanpercentile(cors, 10)),
                                              corr_min=float(np.nanmin(cors)), frac_corr_above_0_8=float(np.nanmean(cors > 0.8)),
                                              stored_value_range=[float(np.nanmin(stored)), float(np.nanmax(stored))],
                                              stored_frac_negative=float((stored < 0).mean()), stored_median_per_neuron_median=float(np.nanmedian(np.nanmedian(stored, 0))),
                                              stored_min_per_neuron_median=float(np.nanmedian(np.nanmin(stored, 0))),
                                              implied_F0_over_F_intercept_ratio_median=float(np.nanmedian(f0_est)) if f0_est else None,
                                              worst_rois=[(int(wide.columns[j]), float(cors[j])) for j in np.argsort(cors)[:5]])
        sc = out["stored_trace_comparison"]
        print(f"  stored DF/F traces vs our 5x5x1 local mean at ROI positions: per-neuron correlation median {sc['corr_median']:.2f}, "
              f"p10 {sc['corr_p10']:.2f}, min {sc['corr_min']:.2f}; {100*sc['frac_corr_above_0_8']:.0f}% of neurons > 0.8")
        print(f"  stored values {sc['stored_value_range']}, {100*sc['stored_frac_negative']:.1f}% negative; per-neuron median of medians {sc['stored_median_per_neuron_median']:.3f}, "
              f"median of per-neuron minima {sc['stored_min_per_neuron_median']:.3f}; worst (roi, r): {sc['worst_rois']}")
        print("  (a linear fit ours = a*stored + b per neuron; if stored = F/F0 - 1 then b/a = 1, i.e. the implied F0 equals the intercept:"
              f" median b/a = {sc['implied_F0_over_F_intercept_ratio_median']})")
        # most active neurons by stored trace
        rng = np.nanpercentile(stored, 95, axis=0) - np.nanpercentile(stored, 5, axis=0)
        labels = tracks.drop_duplicates("roi").set_index("roi")["label"] if "label" in tracks else {}
        top = np.argsort(-rng)[:8]
        out["most_active_stored"] = [(int(j), labels.get(j, "") if hasattr(labels, "get") else "", round(float(rng[j]), 3)) for j in top]
        print(f"  most active neurons by stored DF/F (p95-p5): {out['most_active_stored']}")
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
        summary["structure_extra"] = inspect_structure_extra(f)
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
        summary["time_index"] = time_index_test(f, tracks)
        summary["drift"] = drift_test(f, tracks, n)
        if args.per_frame_means:
            summary["per_frame"] = per_frame_means(f, n, args.per_frame_means, tracks)

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(jsonable(summary), indent=2))
        print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
