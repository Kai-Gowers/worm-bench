#!/usr/bin/env python3
"""Plots from the streaming outputs of inspect_ey.py --per-frame-means.

  <out>/per_frame_means.png   population brightness over time (at ROI positions and
                              whole volume), zero-voxel fraction, with the stimulus
                              windows stored in the NWB file shaded
  <out>/neuron_traces.png     the neurons with the strongest stimulus-locked response
                              in the stored traces: stored trace vs our 5x5 local mean

Usage:
  python -P scripts/plot_ey_traces.py previews/ey --fps 4 \
      --stim nacl:60.5:70.5 --stim pentanedione:120.5:130.5 --stim butanone:180.5:190.5
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
ap.add_argument("--fps", type=float, default=4.0, help="frame rate stored in the file")
ap.add_argument("--stim", action="append", default=[], help="name:start_s:stop_s (from intervals/chemical_stimuli)")
ap.add_argument("--n-show", type=int, default=6)
a = ap.parse_args()
stims = [(s.split(":")[0], float(s.split(":")[1]), float(s.split(":")[2])) for s in a.stim]

pf = pd.read_csv(a.dir / "per_frame_means.csv")
tr = pd.read_csv(a.dir / "per_frame_means_tracks.csv")
st = pd.read_csv(a.dir / "stored_traces.csv", index_col=0)
st.columns = [int(c.replace("roi", "")) for c in st.columns]
mine = tr.pivot(index="t", columns="roi", values="gcamp").sort_index()
labels = tr.drop_duplicates("roi").set_index("roi")["label"]
t_s = pf.t / a.fps


def shade(ax):
    for name, s0, s1 in stims:
        ax.axvspan(s0, s1, color="#3498db", alpha=0.18, lw=0)
        ax.text(s0, ax.get_ylim()[1], name, fontsize=7, va="top", color="#1f618d")


fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
axes[0].plot(t_s, pf.gt_mean / pf.gt_mean.median(), lw=1, color="#27ae60", label="mean 5x5 brightness at ROI positions / median")
axes[0].plot(t_s, pf["mean"] / pf["mean"].median(), lw=1, color="#7f8c8d", label="whole-volume mean / median")
axes[0].set_yscale("log"); axes[0].set_title("population brightness over time (log scale; frame 0 is far brighter than the rest)", fontsize=10)
axes[1].plot(t_s, pf.frac_zero, lw=1, color="#8e44ad", label="fraction of voxels equal to 0")
axes[1].plot(t_s, pf["max"] / 255, lw=1, color="#e67e22", label="frame max / 255")
axes[1].set_title("zero-voxel fraction and per-frame maximum", fontsize=10)
axes[2].plot(t_s[: len(st)], st.mean(1) / st.mean(1).median(), lw=1, color="#c0392b", label="stored SignalRawFluor trace, mean over neurons / median")
axes[2].set_title("stored traces (unitless) averaged over neurons", fontsize=10)
for ax in axes:
    ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper right"); shade(ax)
axes[2].set_xlabel(f"time (s) at the in-file {a.fps:g} Hz; shaded = stimulus windows stored in intervals/chemical_stimuli")
fig.savefig(a.dir / "per_frame_means.png", dpi=110); print(f"wrote {a.dir / 'per_frame_means.png'}")

# strongest stimulus-locked responders by the stored traces (baseline floored at 10% of the typical value)
floor = 0.1 * float(np.nanmedian(st.to_numpy()))
score = pd.Series(0.0, index=st.columns)
best = {}
for name, s0, s1 in stims:
    f0, f1, pre, post = int(round(s0 * a.fps)), int(round(s1 * a.fps)), int(round(20 * a.fps)), int(round(10 * a.fps))
    base = st.iloc[max(0, f0 - pre):f0].mean(); peak = st.iloc[f0:f1 + post].max()
    r = (peak - base) / base.clip(lower=floor)
    for roi, v in r.items():
        if v > score[roi]:
            score[roi] = v; best[roi] = name
show = score.sort_values(ascending=False).index[: a.n_show]
cors = {roi: np.corrcoef(mine[roi].fillna(mine[roi].median()).iloc[: len(st)], st[roi])[0, 1] for roi in show}

fig, axes = plt.subplots(len(show), 1, figsize=(12, 1.9 * len(show)), sharex=True, constrained_layout=True)
for ax, roi in zip(np.atleast_1d(axes), show):
    ax.plot(st.index / a.fps, st[roi] / st[roi].median(), lw=0.9, color="#c0392b", label="stored trace / median")
    m = mine[roi]
    ax.plot(m.index / a.fps, m / m.median(), lw=0.9, color="#27ae60", alpha=0.8, label=f"5x5 local mean / median (r = {cors[roi]:.2f})")
    ax.set_ylabel(f"roi {roi}\n{labels.get(roi, '')}", fontsize=9); ax.grid(alpha=.3); ax.legend(fontsize=7, loc="upper right", ncol=2); shade(ax)
np.atleast_1d(axes)[0].set_title("neurons with the strongest stimulus-locked rise in the stored traces (labels are the biological IDs stored in the file)", fontsize=10)
np.atleast_1d(axes)[-1].set_xlabel(f"time (s) at {a.fps:g} Hz")
fig.savefig(a.dir / "neuron_traces.png", dpi=110); print(f"wrote {a.dir / 'neuron_traces.png'}")
print("top responders (roi, label, best stimulus, (peak-base)/base):", [(int(r), labels.get(r, ""), best.get(r), round(float(score[r]), 2)) for r in show])
