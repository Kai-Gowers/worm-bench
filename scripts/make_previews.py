#!/usr/bin/env python3
"""Generate small visual previews of the ASCENT recording.

Outputs (all written under --out, which is git-ignored):
  depth_slices_t{T}_c{C}.png        3x3 grid of the 9 z-slices of one frame
  depth_slices_t{T}_c{C}_gt.png     same with ground-truth points overlaid
  mip_t{T}.png / mip_t{T}_gt.png    max-intensity projection (c0 | c1 side by side)
  mip_movie_c0c1.mp4                 short MIP movie, c0 | c1 side by side
  mip_movie_c0c1_gt.mp4              same with GT points + track IDs overlaid

Every frame is read from the HDF5 file on demand; nothing larger than one
volume is held in memory at a time.

Usage:
  python -P scripts/make_previews.py --h5 data/ascent/<movie>.h5 \
      --tracks data/ascent/InHouse_Opterra_tracks_GT.csv --out previews/ascent \
      --slice-frames 0 550 --movie-start 150 --movie-len 150 --fps 10
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import imageio.v2 as imageio
import matplotlib
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CHANNEL_LABEL = {0: "c0 (documented: mNeptune, reference)", 1: "c1 (documented: GCaMP6, calcium)"}


def read_frame(f, t, c):
    return f[f"t{t}/c{c}"][()]


def norm_limits(vol, lo=1.0, hi=99.9):
    return np.percentile(vol, lo), np.percentile(vol, hi)


def to_uint8(img, lim):
    lo, hi = lim
    return (np.clip((img.astype(np.float32) - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)


# ----------------------------------------------------------------------------
def depth_slice_grid(vol, pts, title, path, overlay):
    Z = vol.shape[0]
    ncol = int(np.ceil(np.sqrt(Z))); nrow = int(np.ceil(Z / ncol))
    lim = norm_limits(vol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.2 * nrow), constrained_layout=True)
    for z, ax in enumerate(axes.ravel()):
        ax.set_axis_off()
        if z >= Z:
            continue
        ax.imshow(vol[z], cmap="gray", vmin=lim[0], vmax=lim[1], interpolation="nearest")
        if overlay and pts is not None:
            p = pts[pts.z == z]
            ax.scatter(p.x, p.y, s=28, facecolors="none", edgecolors="#ff4d4d", linewidths=0.8)
        ax.set_title(f"z = {z}" + (f"  ({len(pts[pts.z == z])} GT pts)" if overlay and pts is not None else ""), fontsize=9)
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def mip_pair_png(v0, v1, pts, title, path, overlay, lims):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    for ax, v, c in zip(axes, (v0, v1), (0, 1)):
        ax.imshow(v.max(0), cmap="gray", vmin=lims[c][0], vmax=lims[c][1], interpolation="nearest")
        if overlay and pts is not None:
            ax.scatter(pts.x, pts.y, s=22, facecolors="none", edgecolors="#ff4d4d", linewidths=0.8)
        ax.set_title(CHANNEL_LABEL[c] + " — max projection over z", fontsize=9)
        ax.set_axis_off()
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def compose_frame(v0, v1, lims, pts, t, seconds, overlay):
    """Side-by-side MIP of both channels as an RGB uint8 image, with optional overlay."""
    a = to_uint8(v0.max(0), lims[0]); b = to_uint8(v1.max(0), lims[1])
    rgb = np.stack([np.concatenate([a, b], 1)] * 3, -1)
    im = Image.fromarray(rgb)
    d = ImageDraw.Draw(im)
    W = a.shape[1]
    d.text((4, 4), f"t={t}  (~{seconds:.1f}s)   c0", fill=(255, 255, 0))
    d.text((W + 4, 4), "c1", fill=(255, 255, 0))
    if overlay and pts is not None:
        for row in pts.itertuples():
            for off in (0, W):
                x, y = row.x + off, row.y
                d.ellipse([x - 4, y - 4, x + 4, y + 4], outline=(255, 70, 70))
                d.text((x + 5, y - 5), str(row.TrackID), fill=(255, 170, 170))
    return np.asarray(im)


def write_movie(f, tracks, out, start, length, fps, lims, overlay, frame_s):
    stem = "mip_movie_c0c1" + ("_gt" if overlay else "")
    path = out / f"{stem}.mp4"
    try:
        writer = imageio.get_writer(path, fps=fps, codec="libx264", quality=8, macro_block_size=16)
    except Exception as e:  # pragma: no cover - fallback when ffmpeg is unavailable
        print(f"mp4 writer unavailable ({e}); writing GIF instead")
        path = out / f"{stem}.gif"
        writer = imageio.get_writer(path, duration=1000 / fps, loop=0)
    with writer:
        for t in range(start, start + length):
            v0 = read_frame(f, t, 0); v1 = read_frame(f, t, 1)
            pts = tracks[tracks.t == t] if overlay else None
            writer.append_data(compose_frame(v0, v1, lims, pts, t, t * frame_s, overlay))
    print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB, {length} frames @ {fps} fps)")


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5", type=Path, required=True)
    ap.add_argument("--tracks", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--slice-frames", type=int, nargs="+", default=[0, 550])
    ap.add_argument("--movie-start", type=int, default=150)
    ap.add_argument("--movie-len", type=int, default=150)
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--frame-interval-s", type=float, default=1 / 3.3, help="documented recording period (s/frame) for time labels")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    tracks = pd.read_csv(args.tracks)

    with h5py.File(args.h5, "r") as f:
        # fixed display limits from the first slice frame so brightness changes stay visible
        lims = {c: norm_limits(read_frame(f, args.slice_frames[0], c)) for c in (0, 1)}
        for t in args.slice_frames:
            pts = tracks[tracks.t == t]
            vols = {c: read_frame(f, t, c) for c in (0, 1)}
            for c, v in vols.items():
                for overlay in (False, True):
                    suffix = "_gt" if overlay else ""
                    depth_slice_grid(v, pts, f"ASCENT t={t}, {CHANNEL_LABEL[c]}", args.out / f"depth_slices_t{t}_c{c}{suffix}.png", overlay)
            for overlay in (False, True):
                mip_pair_png(vols[0], vols[1], pts, f"ASCENT t={t} max projections" + (" with GT points" if overlay else ""),
                             args.out / f"mip_t{t}{'_gt' if overlay else ''}.png", overlay, lims)
            print(f"wrote depth slices + MIPs for t={t}")
        for overlay in (False, True):
            write_movie(f, tracks, args.out, args.movie_start, args.movie_len, args.fps, lims, overlay, args.frame_interval_s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
