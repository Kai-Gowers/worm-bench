#!/usr/bin/env python3
"""Per-neuron evidence for which channel carries calcium activity.

Reads the per-track trace CSV written by `inspect_ascent.py --per-frame-means`
(columns t, TrackID, c0, c1), removes slow bleaching from each trace with a
quadratic fit, and
  1. prints, per channel, the distribution over neurons of the detrended
     temporal coefficient of variation (CV), and the fraction of neurons that
     fluctuate more in c1 than in c0;
  2. plots the neurons with the largest c1/c0 CV ratio, both channels, so a
     reader can see calcium transients in one channel and not the other.

Usage: python -P scripts/plot_neuron_traces.py previews/ascent/per_frame_means_tracks.csv previews/ascent/neuron_traces.png
"""
import argparse
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("csv"); ap.add_argument("png")
ap.add_argument("--fps", type=float, default=3.3); ap.add_argument("--n-show", type=int, default=6)
a = ap.parse_args()

df = pd.read_csv(a.csv)
chans = [c for c in df.columns if c.startswith("c")]
stats = []
detr = {}
for tid, g in df.groupby("TrackID"):
    g = g.sort_values("t")
    if len(g) < 0.9 * df.t.nunique():
        continue
    row = {"TrackID": tid}
    for c in chans:
        y = g[c].to_numpy(); t = g.t.to_numpy()
        trend = np.polyval(np.polyfit(t, y, 2), t)
        r = (y - trend) / trend.mean()           # fractional fluctuation about the bleaching trend
        row[f"{c}_cv"] = r.std()
        row[f"{c}_max_dff"] = (y / trend).max() - 1  # peak fractional excursion above trend
        detr[(tid, c)] = (t, y / trend)
    stats.append(row)
st = pd.DataFrame(stats).set_index("TrackID")
if len(chans) == 2:
    a0, a1 = chans
    st["ratio"] = st[f"{a1}_cv"] / st[f"{a0}_cv"]
    print(f"neurons analysed: {len(st)} (tracks covering >=90% of frames)")
    for c in chans:
        print(f"{c}: detrended CV median {st[f'{c}_cv'].median():.3f}, p90 {st[f'{c}_cv'].quantile(.9):.3f}, max {st[f'{c}_cv'].max():.3f}; "
              f"peak excursion above trend median {st[f'{c}_max_dff'].median():.2f}, p90 {st[f'{c}_max_dff'].quantile(.9):.2f}")
    print(f"fraction of neurons with CV({a1}) > CV({a0}): {(st['ratio'] > 1).mean():.2f}; "
          f"with CV({a1}) > 2*CV({a0}): {(st['ratio'] > 2).mean():.2f}; with CV({a0}) > 2*CV({a1}): {(st['ratio'] < .5).mean():.2f}")
    # correlation between the two channels' detrended traces, per neuron (motion artefacts show in both)
    cors = [np.corrcoef(detr[(tid, a0)][1], detr[(tid, a1)][1])[0, 1] for tid in st.index]
    print(f"per-neuron correlation of detrended {a0} vs {a1} traces: median {np.median(cors):.2f}")
    show = st.sort_values("ratio", ascending=False).index[: a.n_show]
else:
    show = st.index[: a.n_show]

fig, axes = plt.subplots(len(show), 1, figsize=(12, 1.9 * len(show)), sharex=True, constrained_layout=True)
for ax, tid in zip(np.atleast_1d(axes), show):
    for c, col in zip(chans, ["#c0392b", "#27ae60"]):
        t, y = detr[(tid, c)]
        ax.plot(t / a.fps, y, lw=0.9, color=col, label=f"{c} (CV {st.loc[tid, f'{c}_cv']:.3f})")
    ax.set_ylabel(f"track {tid}", fontsize=9); ax.grid(alpha=.3); ax.legend(fontsize=7, loc="upper right", ncol=2)
np.atleast_1d(axes)[0].set_title("brightness / bleaching trend at the GT position (5x5 window): neurons with the largest c1-vs-c0 fluctuation ratio", fontsize=10)
np.atleast_1d(axes)[-1].set_xlabel(f"time (s) assuming documented {a.fps} Hz")
fig.savefig(a.png, dpi=110)
print(f"wrote {a.png}")
