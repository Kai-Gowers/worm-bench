#!/usr/bin/env python3
"""Plots from the streaming outputs of inspect_kk.py --per-frame-means (mirrors plot_ey_traces.py).

  <out>/per_frame_means.png   population brightness over time in both channels (mean inside the
                              label masks and background), the GCaMP/RFP ratio, the number of ROIs
                              per frame, and the population motion per frame from the centroids
  <out>/neuron_traces.png     the most active neurons by the stored trace: stored trace vs our
                              mask mean in GCaMP, the RFP mask mean, and the GCaMP/RFP ratio

There is no stimulus in this dataset, so nothing is shaded.

Usage:
  python -P scripts/plot_kk_traces.py previews/kk --fps 1.67
"""
import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("dir", type=Path)
ap.add_argument("--fps", type=float, default=1.67, help="frame rate stored in the file")
ap.add_argument("--n-show", type=int, default=6)
a = ap.parse_args()

pf = pd.read_csv(a.dir / "per_frame_means.csv")
tr = pd.read_csv(a.dir / "per_frame_means_tracks.csv")
st = pd.read_csv(a.dir / "stored_traces.csv", index_col=0)
st.columns = [int(c.replace("roi", "")) for c in st.columns]
tracks = pd.read_csv(a.dir / "calcium_tracks.csv")
chans = [c[: -len("_mask_mean")] for c in pf.columns if c.endswith("_mask_mean")]
g, r = chans[0], chans[1]  # channel 0 = GCaMP, channel 1 = TagRFP-T (documented in the file)
t_s = pf.t / a.fps
gm = tr.pivot(index="t", columns="roi", values=f"{g}_mask").sort_index()
rm = tr.pivot(index="t", columns="roi", values=f"{r}_mask").sort_index()
labels = tr.drop_duplicates("roi").set_index("roi")["label"]

# population motion per frame: median in-plane centroid step of all ROIs between t-1 and t
ts_ = tracks.sort_values(["roi", "t"]); grp = ts_.groupby("roi")
step = grp[["x", "y"]].diff(); dt = grp["t"].diff()
disp = np.sqrt((step[dt == 1] ** 2).sum(axis=1))
motion = pd.DataFrame({"t": ts_.t[dt == 1], "d": disp}).groupby("t")["d"].agg(["median", "max"])

fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True, constrained_layout=True)
axes[0].plot(t_s, pf[f"{g}_mask_mean"], lw=1, color="#27ae60", label="GCaMP, mean inside label masks")
axes[0].plot(t_s, pf[f"{r}_mask_mean"], lw=1, color="#c0392b", label="TagRFP-T, mean inside label masks")
axes[0].plot(t_s, pf[f"{g}_bg_mean"], lw=1, color="#27ae60", ls=":", label="GCaMP background (label 0)")
axes[0].plot(t_s, pf[f"{r}_bg_mean"], lw=1, color="#c0392b", ls=":", label="TagRFP-T background")
axes[0].set_title("population brightness over time (8-bit gray counts)", fontsize=10)
axes[1].plot(t_s, pf[f"{g}_mask_mean"] / pf[f"{r}_mask_mean"], lw=1, color="#8e44ad", label="GCaMP / TagRFP-T mask means")
axes[1].plot(t_s[: len(st)], st.mean(1) / st.mean(1).median() * (pf[f"{g}_mask_mean"] / pf[f"{r}_mask_mean"]).median(), lw=1, color="#e67e22", alpha=0.8,
             label="stored SignalRawFluor, mean over neurons (scaled to the same median)")
axes[1].set_title("GCaMP / RFP ratio of the population vs the stored traces", fontsize=10)
axes[2].plot(t_s, pf.n_rois, lw=1, color="#2c3e50", label="ROIs present in the label volume")
axes[2].set_title("number of labelled ROIs per frame", fontsize=10)
axes[3].plot(motion.index / a.fps, motion["median"], lw=1, color="#2980b9", label="median in-plane centroid step (voxels)")
axes[3].plot(motion.index / a.fps, motion["max"], lw=0.6, color="#7f8c8d", alpha=0.6, label="max step over ROIs")
axes[3].set_yscale("symlog", linthresh=1); axes[3].set_title("population motion between consecutive frames", fontsize=10)
for ax in axes:
    ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper right")
axes[3].set_xlabel(f"time (s) at the in-file {a.fps:g} Hz")
fig.savefig(a.dir / "per_frame_means.png", dpi=110); print(f"wrote {a.dir / 'per_frame_means.png'}")

# most active neurons by the stored trace (range / median), which is spontaneous activity here
act = (st.quantile(0.95) - st.quantile(0.05)) / st.median().clip(lower=1e-6)
show = act.sort_values(ascending=False).index[: a.n_show]
n = min(len(gm), len(st))
cors = {roi: np.corrcoef(gm[roi].fillna(gm[roi].median()).iloc[:n], st[roi].iloc[:n])[0, 1] for roi in show}
fig, axes = plt.subplots(len(show), 1, figsize=(12, 1.9 * len(show)), sharex=True, constrained_layout=True)
for ax, roi in zip(np.atleast_1d(axes), show):
    ax.plot(st.index / a.fps, st[roi] / st[roi].median(), lw=0.9, color="#e67e22", label="stored trace / median")
    ax.plot(gm.index / a.fps, gm[roi] / gm[roi].median(), lw=0.9, color="#27ae60", alpha=0.8, label=f"GCaMP mask mean / median (r = {cors[roi]:.2f})")
    ax.plot(rm.index / a.fps, rm[roi] / rm[roi].median(), lw=0.7, color="#c0392b", alpha=0.7, label="TagRFP-T mask mean / median")
    ax.set_ylabel(f"roi {roi}\n{labels.get(roi, '')}", fontsize=9); ax.grid(alpha=.3); ax.legend(fontsize=7, loc="upper right", ncol=3)
np.atleast_1d(axes)[0].set_title("neurons with the largest spontaneous swings in the stored traces (labels are the biological IDs stored in the file)", fontsize=10)
np.atleast_1d(axes)[-1].set_xlabel(f"time (s) at {a.fps:g} Hz")
fig.savefig(a.dir / "neuron_traces.png", dpi=110); print(f"wrote {a.dir / 'neuron_traces.png'}")
print("most active (roi, label, (p95-p5)/median):", [(int(x), labels.get(x, ""), round(float(act[x]), 2)) for x in show])
