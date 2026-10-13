#!/usr/bin/env python3
"""Plots from the streaming outputs of inspect_sk.py --per-frame-means (mirrors plot_ey_traces.py).

  <out>/per_frame_means.png   population brightness over time (5x5 mean at ROI positions and whole
                              volume, both minus the camera offset), per-frame max, the stored DF/F
                              averaged over neurons, and the population motion from the ROI positions
  <out>/neuron_traces.png     the most active neurons by the stored DF/F: stored trace vs our local
                              mean at the stored (x, y, z) and vs the max over z of that mean

Nothing is shaded: no stimulus timing is stored in these files.

Usage:
  python -P scripts/plot_sk_traces.py previews/sk1 --fps 1.04
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
ap.add_argument("--fps", type=float, default=1.04, help="frame rate stored in the file")
ap.add_argument("--n-show", type=int, default=6)
ap.add_argument("--title", default="")
a = ap.parse_args()

pf = pd.read_csv(a.dir / "per_frame_means.csv")
tr = pd.read_csv(a.dir / "per_frame_means_tracks.csv")
st = pd.read_csv(a.dir / "stored_traces.csv", index_col=0)
st.columns = [int(c.replace("roi", "")) for c in st.columns]
tracks = pd.read_csv(a.dir / "calcium_tracks.csv")
mine = tr.pivot(index="t", columns="roi", values="gcamp").sort_index()
labels = tr.drop_duplicates("roi").set_index("roi")["label"].fillna("")
t_s = pf.t / a.fps
offset = pf.p50.median()  # camera offset: the median voxel is background

ts_ = tracks.sort_values(["roi", "t"]); grp = ts_.groupby("roi")
step = grp[["x", "y"]].diff(); dt = grp["t"].diff()
disp = np.sqrt((step[dt == 1] ** 2).sum(axis=1))
motion = pd.DataFrame({"t": ts_.t[dt == 1], "d": disp}).groupby("t")["d"].agg(["median", "max"])

fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True, constrained_layout=True)
axes[0].plot(t_s, pf.gt_mean - offset, lw=1, color="#27ae60", label=f"mean 5x5 brightness at ROI positions minus offset ({offset:.0f} counts)")
axes[0].plot(t_s, pf["mean"] - offset, lw=1, color="#7f8c8d", label="whole-volume mean minus offset")
axes[0].plot(t_s, pf.p99 - offset, lw=1, color="#2c3e50", ls=":", label="volume 99th percentile minus offset")
axes[0].set_title("population brightness over time (gray counts above the camera offset)", fontsize=10)
axes[1].plot(t_s, pf["max"], lw=0.8, color="#e67e22", label="frame max (4095 = 12-bit full scale)")
axes[1].set_yscale("log"); axes[1].set_title("per-frame maximum voxel value", fontsize=10)
axes[2].plot(t_s[: len(st)], st.mean(1), lw=1, color="#c0392b", label="stored SignalDFoF, mean over neurons")
axes[2].set_title("stored DF/F traces averaged over neurons (unit attribute says 'Percentage'; values are fractions)", fontsize=10)
axes[3].plot(motion.index / a.fps, motion["median"], lw=1, color="#2980b9", label="median in-plane ROI step (voxels)")
axes[3].plot(motion.index / a.fps, motion["max"], lw=0.6, color="#7f8c8d", alpha=0.6, label="max step over ROIs")
axes[3].set_title("annotation motion between consecutive frames", fontsize=10)
for ax in axes:
    ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper right")
axes[3].set_xlabel(f"time (s) at the in-file {a.fps:g} Hz")
if a.title:
    fig.suptitle(a.title, fontsize=11)
fig.savefig(a.dir / "per_frame_means.png", dpi=110); print(f"wrote {a.dir / 'per_frame_means.png'}")

# most active neurons by the stored DF/F (p95 - p5)
act = st.quantile(0.95) - st.quantile(0.05)
show = act.sort_values(ascending=False).index[: a.n_show]
n = min(len(mine), len(st))
cors = {roi: np.corrcoef(mine[roi].fillna(mine[roi].median()).iloc[:n], st[roi].iloc[:n])[0, 1] for roi in show}
fig, axes = plt.subplots(len(show), 1, figsize=(12, 1.9 * len(show)), sharex=True, constrained_layout=True)
for ax, roi in zip(np.atleast_1d(axes), show):
    ax.plot(st.index / a.fps, st[roi], lw=0.9, color="#c0392b", label="stored DF/F")
    m = mine[roi] - offset
    ax.plot(m.index / a.fps, m / max(m.median(), 1e-6) - 1, lw=0.9, color="#27ae60", alpha=0.8,
            label=f"our 5x5 mean at stored (x,y,z) minus offset, / median - 1 (r = {cors[roi]:.2f})")
    ax.set_ylabel(f"roi {roi}\n{labels.get(roi, '')}", fontsize=9); ax.grid(alpha=.3); ax.legend(fontsize=7, loc="upper right", ncol=2)
np.atleast_1d(axes)[0].set_title("neurons with the largest swings in the stored DF/F (labels are the biological IDs stored in the file; blank = unlabelled)", fontsize=10)
np.atleast_1d(axes)[-1].set_xlabel(f"time (s) at {a.fps:g} Hz")
fig.savefig(a.dir / "neuron_traces.png", dpi=110); print(f"wrote {a.dir / 'neuron_traces.png'}")
print("most active (roi, label, p95-p5 of stored DF/F):", [(int(x), labels.get(x, ""), round(float(act[x]), 2)) for x in show])
