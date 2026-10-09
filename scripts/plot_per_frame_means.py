#!/usr/bin/env python3
"""Plot the per-frame intensity trace written by inspect_ascent.py --per-frame-means.

Shows mean brightness at ground-truth neuron positions for each channel over
the whole recording, normalised to its own median, so stimulus-locked calcium
responses (expected in the GCaMP channel) can be seen against the stable
reference channel.

Usage: python -P scripts/plot_per_frame_means.py previews/ascent/per_frame_means.csv previews/ascent/per_frame_means.png [--fps 3.3]
"""
import argparse
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("csv"); ap.add_argument("png")
ap.add_argument("--fps", type=float, default=3.3, help="documented frame rate used for the time axis")
ap.add_argument("--stim-period", type=float, nargs=2, default=(60.0, 30.0), help="documented buffer s / stimulus s")
a = ap.parse_args()

df = pd.read_csv(a.csv)
t_s = df.t / a.fps
fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True, constrained_layout=True)
for ax, kind, title in [(axes[0], "gt_mean", "mean brightness at GT neuron positions (5x5 window), / median"),
                        (axes[1], "mean", "whole-volume mean brightness, / median")]:
    for c, lab, col in [("c0", "c0 (documented mNeptune, reference)", "#c0392b"), ("c1", "c1 (documented GCaMP6, calcium)", "#27ae60")]:
        y = df[f"{c}_{kind}"]; ax.plot(t_s, y / y.median(), lw=1, label=lab, color=col)
    ax.set_title(title, fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=8, loc="upper right")
# documented stimulus schedule: 60 s buffer then 30 s stimulus, repeated (phase relative to frame 0 unknown)
buf, stim = a.stim_period
tt = 0
while tt < t_s.max():
    for ax in axes:
        ax.axvspan(tt + buf, min(tt + buf + stim, t_s.max()), color="#3498db", alpha=.12, lw=0)
    tt += buf + stim
axes[1].set_xlabel(f"time (s), assuming documented {a.fps} Hz; shaded = documented 30 s stimulus windows IF cycles start at frame 0 (phase unverified)")
fig.savefig(a.png, dpi=110)
print(f"wrote {a.png}")
