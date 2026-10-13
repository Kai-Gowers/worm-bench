#!/usr/bin/env python3
"""Visual previews of one SK1/SK2 (DANDI 000565/000472) NWB recording, mirroring make_previews_ey.py.

Outputs (under --out, git-ignored):
  depth_slices_t{T}_gcamp.png (+_gt)   grid of the 12 z-slices of one calcium frame, with/without ROI points
  mip_t{T}.png (+_gt)                  calcium max projection | NeuroPAL RGB max projection
  neuropal_rgb_mip.png (+_labels)      the structural NeuroPAL stack with its neuron IDs
  mip_movie_gcamp.mp4                  short calcium MIP movie, raw | with ROI points + labels
Frames are read one at a time from the NWB file.  Display limits are fixed from
one frame (--lim-frame), from its 1st to 99.9th percentile because these cameras
record a ~90-count offset and the neurons sit only a few tens of counts above it.
Images are shown with the stored X axis horizontal (the worm lies along X).

Usage:
  python -P scripts/make_previews_sk.py --nwb data/sk1/<file>.nwb --out previews/sk1 \
      --slice-frames 0 750 --movie-start 0 --movie-len 200 --fps 10
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
from inspect_ey import NEUROPAL_IMG, NEUROPAL_SEG, decode, load_calcium_tracks, ragged_labels, read_frame  # noqa: E402
from inspect_sk import NEUROPAL_VOL  # noqa: E402


def neuropal_rgb_mip(f: h5py.File):
    """RGB max projection of the NeuroPAL stack (X, Y, Z, C) -> (Y, X, 3) uint8, the point table, channel names."""
    d = f[f"{NEUROPAL_IMG}/data"]
    rgbw = f[f"{NEUROPAL_IMG}/RGBW_channels"][()]
    order = [decode(c) for c in f[f"{NEUROPAL_VOL}/order_optical_channels/channels"][()]]
    chans, names = [], []
    for c in rgbw[:3]:
        mip = d[:, :, :, int(c)].max(axis=2).T.astype(np.float32)  # (Y, X)
        chans.append(to_uint8(mip, (np.percentile(mip, 1), np.percentile(mip, 99.9))))
        names.append(order[int(c)])
    g = f[NEUROPAL_SEG]; vm = g["voxel_mask"][()]
    pts = np.rec.fromarrays([vm["x"], vm["y"], vm["z"], np.array(ragged_labels(g), dtype=object)], names="x,y,z,label")
    return np.stack(chans, -1), pts, names


def mip_pair_png(vol, rgb, pts, npts, title, path, overlay, lim, name):
    fig, axes = plt.subplots(2, 1, figsize=(16, 8.5), constrained_layout=True)
    axes[0].imshow(vol.max(0), cmap="gray", vmin=lim[0], vmax=lim[1], interpolation="nearest")
    if overlay and pts is not None:
        axes[0].scatter(pts.x, pts.y, s=22, facecolors="none", edgecolors="#ff4d4d", linewidths=0.8)
    axes[0].set_title(f"GCaMP calcium frame, max projection over z ({vol.shape[2]} x {vol.shape[1]} px, display {lim[0]:.0f}..{lim[1]:.0f} counts)", fontsize=9)
    axes[1].imshow(rgb, interpolation="nearest")
    if overlay and npts is not None:
        axes[1].scatter(npts.x, npts.y, s=22, facecolors="none", edgecolors="#ffffff", linewidths=0.6)
    axes[1].set_title(f"NeuroPAL structural stack, RGB max projection ({rgb.shape[1]} x {rgb.shape[0]} px, same pixel size as the calcium frame)", fontsize=9)
    for ax in axes:
        ax.set_axis_off()
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def neuropal_png(rgb, npts, names, path, labels):
    fig, ax = plt.subplots(figsize=(18, 5.5), constrained_layout=True)
    ax.imshow(rgb, interpolation="nearest"); ax.set_axis_off()
    ax.scatter(npts.x, npts.y, s=18, facecolors="none", edgecolors="#ffffff", linewidths=0.5)
    if labels:
        for p in npts:
            if p.label.strip():
                ax.text(p.x + 2, p.y - 2, p.label, color="#ffff99", fontsize=5)
    ax.set_title(f"NeuroPAL RGB max projection (R={names[0]}, G={names[1]}, B={names[2]}) with stored neuron centres"
                 + (" and biological IDs (unlabelled centres have no text)" if labels else ""), fontsize=9)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def compose_frame(vol, lim, pts, t, seconds, rate, name, scale=1):
    a = to_uint8(vol.max(0), lim)
    im = Image.fromarray(a).resize((a.shape[1] * scale, a.shape[0] * scale), Image.NEAREST).convert("RGB")
    W, H = im.size
    canvas = Image.new("RGB", (W, H * 2))
    canvas.paste(im, (0, 0)); canvas.paste(im, (0, H))
    d = ImageDraw.Draw(canvas)
    d.text((4, 4), f"{name}  t={t}  ({seconds:.1f}s, in-file {rate:g} Hz)  GCaMP max projection", fill=(255, 255, 0))
    d.text((4, H + 4), "ROI points + IDs (red = unlabelled, yellow = labelled)", fill=(255, 255, 0))
    if pts is not None:
        for row in pts.itertuples():
            x, y = row.x * scale + scale // 2, row.y * scale + scale // 2 + H
            lab = str(row.label) if isinstance(row.label, str) else ""
            col = (255, 230, 90) if lab.strip() else (255, 70, 70)
            d.ellipse([x - 4, y - 4, x + 4, y + 4], outline=col)
            if lab.strip():
                d.text((x + 5, y - 6), lab, fill=(255, 230, 150))
    return np.asarray(canvas)


def write_movie(f, tracks, out, start, length, fps, lim, rate, name):
    path = out / "mip_movie_gcamp.mp4"
    try:
        writer = imageio.get_writer(path, fps=fps, codec="libx264", quality=7, macro_block_size=16)
    except Exception as e:  # pragma: no cover
        print(f"mp4 writer unavailable ({e}); writing GIF instead")
        path = out / "mip_movie_gcamp.gif"; writer = imageio.get_writer(path, duration=1000 / fps, loop=0)
    with writer:
        for t in range(start, start + length):
            v = read_frame(f, t)
            writer.append_data(compose_frame(v, lim, tracks[tracks.t == t], t, t / rate, rate, name))
    print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB, {length} frames @ {fps} fps)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nwb", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", default=None, help="short recording name for titles (default: the NWB identifier)")
    ap.add_argument("--slice-frames", type=int, nargs="+", default=[0, 750])
    ap.add_argument("--lim-frame", type=int, default=750, help="frame whose 1..99.9th percentiles set the fixed display range")
    ap.add_argument("--movie-start", type=int, default=0)
    ap.add_argument("--movie-len", type=int, default=200)
    ap.add_argument("--fps", type=int, default=10)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.nwb, "r") as f:
        name = args.name or decode(f["identifier"][()])
        tracks = load_calcium_tracks(f)
        rate = float(f["acquisition/CalciumImageSeries/starting_time"].attrs["rate"])
        vl = read_frame(f, args.lim_frame)
        lim = (float(np.percentile(vl, 1)), float(np.percentile(vl, 99.9)))
        print(f"display range fixed to [{lim[0]:.0f}, {lim[1]:.0f}] from frame {args.lim_frame}")
        rgb, npts, names = neuropal_rgb_mip(f)
        neuropal_png(rgb, npts, names, args.out / "neuropal_rgb_mip.png", labels=False)
        neuropal_png(rgb, npts, names, args.out / "neuropal_rgb_mip_labels.png", labels=True)
        print("wrote NeuroPAL RGB previews")
        for t in args.slice_frames:
            pts = tracks[tracks.t == t]
            v = read_frame(f, t)
            for overlay in (False, True):
                suffix = "_gt" if overlay else ""
                depth_slice_grid(v, pts, f"{name} t={t} ({t / rate:.1f} s), GCaMP" + ("; circles = ROI points" if overlay else ""),
                                 args.out / f"depth_slices_t{t}_gcamp{suffix}.png", overlay)
                mip_pair_png(v, rgb, pts, npts, f"{name} t={t} max projections" + (" with ROI points" if overlay else ""),
                             args.out / f"mip_t{t}{suffix}.png", overlay, lim, name)
            print(f"wrote depth slices + MIPs for t={t}")
        write_movie(f, tracks, args.out, args.movie_start, args.movie_len, args.fps, lim, rate, name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
