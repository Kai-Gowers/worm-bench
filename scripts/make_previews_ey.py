#!/usr/bin/env python3
"""Visual previews of one EY (DANDI 000541) NWB recording, mirroring make_previews.py.

Outputs (under --out, git-ignored):
  depth_slices_t{T}_gcamp.png (+_gt)   5x5 grid of the 21 z-slices of one calcium frame
  mip_t{T}.png (+_gt)                  calcium max projection | NeuroPAL RGB max projection
  neuropal_rgb_mip.png (+_labels)      the structural NeuroPAL stack with its neuron IDs
  mip_movie_gcamp.mp4                  short calcium MIP movie, raw | with ROI points + labels
Frames are read one at a time from the NWB file.  Display limits are fixed
from one frame (--lim-frame) so brightness changes over time are preserved.

Usage:
  python -P scripts/make_previews_ey.py --nwb data/ey/<file>.nwb --out previews/ey \
      --slice-frames 0 480 --movie-start 200 --movie-len 150 --fps 10
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import imageio.v2 as imageio
import matplotlib
import numpy as np
from PIL import Image, ImageDraw

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_previews import depth_slice_grid, to_uint8  # noqa: E402
from inspect_ey import NEUROPAL_IMG, NEUROPAL_SEG, load_calcium_tracks, load_stimuli, ragged_labels, read_frame  # noqa: E402


def neuropal_rgb_mip(f: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    """RGB max projection of the NeuroPAL stack (X, Y, Z, C) -> (Y, X, 3) uint8 and the point table."""
    d = f[f"{NEUROPAL_IMG}/data"]
    rgbw = f[f"{NEUROPAL_IMG}/RGBW_channels"][()]
    chans = []
    for c in rgbw[:3]:
        mip = d[:, :, :, int(c)].max(axis=2).T.astype(np.float32)  # (Y, X)
        chans.append(to_uint8(mip, (np.percentile(mip, 1), np.percentile(mip, 99.9))))
    g = f[NEUROPAL_SEG]; vm = g["voxel_mask"][()]
    pts = np.rec.fromarrays([vm["x"], vm["y"], vm["z"], np.array(ragged_labels(g), dtype=object)], names="x,y,z,label")
    return np.stack(chans, -1), pts


def mip_pair_png(vol, rgb, pts, npts, title, path, overlay, lim):
    fig, axes = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)
    axes[0].imshow(vol.max(0), cmap="gray", vmin=lim[0], vmax=lim[1], interpolation="nearest")
    if overlay and pts is not None:
        axes[0].scatter(pts.x, pts.y, s=22, facecolors="none", edgecolors="#ff4d4d", linewidths=0.8)
    axes[0].set_title("GCaMP calcium frame, max projection over z (256 x 128 px)", fontsize=9)
    axes[1].imshow(rgb, interpolation="nearest")
    if overlay and npts is not None:
        axes[1].scatter(npts.x, npts.y, s=22, facecolors="none", edgecolors="#ffffff", linewidths=0.6)
    axes[1].set_title("NeuroPAL structural stack, RGB max projection (512 x 256 px, 2x the calcium resolution)", fontsize=9)
    for ax in axes:
        ax.set_axis_off()
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def neuropal_png(rgb, npts, path, labels):
    fig, ax = plt.subplots(figsize=(14, 7.5), constrained_layout=True)
    ax.imshow(rgb, interpolation="nearest"); ax.set_axis_off()
    ax.scatter(npts.x, npts.y, s=18, facecolors="none", edgecolors="#ffffff", linewidths=0.5)
    if labels:
        for p in npts:
            ax.text(p.x + 2, p.y - 2, p.label, color="#ffff99", fontsize=4.5)
    ax.set_title("NeuroPAL RGB max projection (R=mNeptune 2.5, G=CyOFP1, B=mTagBFP2) with stored neuron centres"
                 + (" and biological IDs" if labels else ""), fontsize=9)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def compose_frame(vol, lim, pts, t, seconds, stim_label, overlay, scale=3):
    a = to_uint8(vol.max(0), lim)
    im = Image.fromarray(a).resize((a.shape[1] * scale, a.shape[0] * scale), Image.NEAREST).convert("RGB")
    W = im.width
    canvas = Image.new("RGB", (W * 2, im.height))
    canvas.paste(im, (0, 0)); canvas.paste(im, (W, 0))
    d = ImageDraw.Draw(canvas)
    d.text((4, 4), f"t={t}  ({seconds:.2f}s, in-file 4 Hz)  {stim_label}", fill=(255, 255, 0))
    d.text((W + 4, 4), "ROI points + IDs", fill=(255, 255, 0))
    if overlay and pts is not None:
        for row in pts.itertuples():
            x, y = W + row.x * scale + scale // 2, row.y * scale + scale // 2
            col = (255, 70, 70) if row.z < vol.shape[0] else (90, 160, 255)  # blue: z outside the stack
            d.ellipse([x - 5, y - 5, x + 5, y + 5], outline=col)
            d.text((x + 6, y - 6), str(row.label), fill=(255, 200, 200))
    return np.asarray(canvas)


def write_movie(f, tracks, out, start, length, fps, lim, frame_s, stimuli):
    path = out / "mip_movie_gcamp.mp4"
    try:
        writer = imageio.get_writer(path, fps=fps, codec="libx264", quality=8, macro_block_size=16)
    except Exception as e:  # pragma: no cover
        print(f"mp4 writer unavailable ({e}); writing GIF instead")
        path = out / "mip_movie_gcamp.gif"; writer = imageio.get_writer(path, duration=1000 / fps, loop=0)
    with writer:
        for t in range(start, start + length):
            v = read_frame(f, t)
            s = t * frame_s
            stim = next((f"STIMULUS: {x['stimulus']}" for x in stimuli if x["start_s"] <= s < x["stop_s"]), "")
            writer.append_data(compose_frame(v, lim, tracks[tracks.t == t], t, s, stim, True))
    print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB, {length} frames @ {fps} fps)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nwb", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--slice-frames", type=int, nargs="+", default=[0, 480])
    ap.add_argument("--lim-frame", type=int, default=480, help="frame whose 0..99.9th percentile sets the fixed display range")
    ap.add_argument("--movie-start", type=int, default=200)
    ap.add_argument("--movie-len", type=int, default=150)
    ap.add_argument("--fps", type=int, default=10)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.nwb, "r") as f:
        tracks = load_calcium_tracks(f)
        stimuli = load_stimuli(f)
        rate = float(f["acquisition/CalciumImageSeries/starting_time"].attrs["rate"])
        vl = read_frame(f, args.lim_frame)
        lim = (0.0, float(np.percentile(vl, 99.9)))
        print(f"display range fixed to [0, {lim[1]:.0f}] from frame {args.lim_frame}")
        rgb, npts = neuropal_rgb_mip(f)
        neuropal_png(rgb, npts, args.out / "neuropal_rgb_mip.png", labels=False)
        neuropal_png(rgb, npts, args.out / "neuropal_rgb_mip_labels.png", labels=True)
        print("wrote NeuroPAL RGB previews")
        for t in args.slice_frames:
            pts = tracks[tracks.t == t]
            v = read_frame(f, t)
            for overlay in (False, True):
                suffix = "_gt" if overlay else ""
                depth_slice_grid(v, pts, f"EY sub-20190928-13 t={t} ({t / rate:.1f} s), GCaMP; {len(pts[pts.z >= v.shape[0]])} ROI(s) with z beyond the last slice not shown",
                                 args.out / f"depth_slices_t{t}_gcamp{suffix}.png", overlay)
                mip_pair_png(v, rgb, pts, npts, f"EY sub-20190928-13 t={t} max projections" + (" with ROI points" if overlay else ""),
                             args.out / f"mip_t{t}{suffix}.png", overlay, lim)
            print(f"wrote depth slices + MIPs for t={t}")
        write_movie(f, tracks, args.out, args.movie_start, args.movie_len, args.fps, lim, 1 / rate, stimuli)
    return 0


if __name__ == "__main__":
    sys.exit(main())
