#!/usr/bin/env python3
"""Visual previews of one KK (DANDI 000692) NWB recording, mirroring make_previews.py / make_previews_ey.py.

Outputs (under --out, git-ignored):
  depth_slices_t{T}_{gcamp,rfp}.png (+_gt)   5x5 grid of the 25 z-slices of one frame, one channel, with/without ROI centroids
  mip_t{T}.png (+_gt)                        GCaMP | TagRFP-T max projections; _gt adds centroids (left) and label-mask outlines (right)
  neuropal_rgb_mip.png (+_labels)            the structural NeuroPAL stack with its neuron IDs
  mip_movie.mp4                              2x2 movie: GCaMP, RFP, GCaMP + centroids and IDs, RFP + mask outlines
Images are displayed with the stored X axis vertical and Y horizontal (the worm lies along Y).
Per-frame centroids come from the cache written by inspect_kk.py --tracks-csv.

Usage:
  python -P scripts/make_previews_kk.py --nwb data/kk/<file>.nwb --tracks previews/kk/calcium_tracks.csv \
      --out previews/kk --slice-frames 0 800 --movie-start 0 --movie-len 200 --fps 10
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
from inspect_ey import decode, ragged_labels  # noqa: E402
from inspect_kk import NEUROPAL_IMG, NEUROPAL_SEG, NEUROPAL_VOL, channel_names, load_calcium_tracks, read_frame, read_seg  # noqa: E402

NAME = "KK sub-231013-02"


def neuropal_rgb_mip(f: h5py.File):
    """RGB max projection of the NeuroPAL stack (X, Y, Z, C) -> (X, Y, 3) uint8 (X vertical), the point table, channel names."""
    d = f[f"{NEUROPAL_IMG}/data"]
    rgbw = f[f"{NEUROPAL_IMG}/RGBW_channels"][()]
    g = f[NEUROPAL_VOL]
    order = [decode(c) for c in g["order_optical_channels/channels"][()]] if "order_optical_channels" in g and "channels" in g["order_optical_channels"] else None
    chans, names = [], []
    for c in rgbw[:3]:
        mip = d[:, :, :, int(c)].max(axis=2).astype(np.float32)  # (X, Y)
        chans.append(to_uint8(mip, (np.percentile(mip, 1), np.percentile(mip, 99.9))))
        names.append(order[int(c)] if order else f"ch{int(c)}")
    s = f[NEUROPAL_SEG]; vm = s["voxel_mask"][()]
    pts = np.rec.fromarrays([vm["x"], vm["y"], vm["z"], np.array(ragged_labels(s), dtype=object)], names="x,y,z,label")
    return np.stack(chans, -1), pts, names


def mask_edges(seg_mip: np.ndarray) -> np.ndarray:
    """Boolean image of label boundaries in a 2-D label map."""
    e = np.zeros(seg_mip.shape, bool)
    e[1:, :] |= seg_mip[1:, :] != seg_mip[:-1, :]
    e[:, 1:] |= seg_mip[:, 1:] != seg_mip[:, :-1]
    return e & (seg_mip > 0)


def mip_pair_png(vg, vr, seg, pts, title, path, overlay, lims):
    """vg, vr: (Z, Y, X) volumes; displayed transposed so X is vertical."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 12 * 2 * vg.shape[2] / vg.shape[1] + 1.2), constrained_layout=True)
    for ax, v, lim, name in zip(axes, (vg, vr), lims, ("GCaMP (channel 0)", "TagRFP-T (channel 1, reference)")):
        ax.imshow(v.max(0).T, cmap="gray", vmin=lim[0], vmax=lim[1], interpolation="nearest")
        ax.set_title(f"{name}, max projection over z (stored X vertical {v.shape[2]} px, Y horizontal {v.shape[1]} px)", fontsize=9)
        ax.set_axis_off()
    if overlay:
        if pts is not None:
            axes[0].scatter(pts.y, pts.x, s=20, facecolors="none", edgecolors="#ff4d4d", linewidths=0.7)
        if seg is not None:
            e = mask_edges(seg.max(0).T)
            ov = np.zeros(e.shape + (4,), np.float32); ov[e] = (0.3, 1.0, 0.3, 0.9)
            axes[1].imshow(ov, interpolation="nearest")
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def neuropal_png(rgb, npts, names, path, labels):
    fig, ax = plt.subplots(figsize=(14, 14 * rgb.shape[0] / rgb.shape[1] + 1), constrained_layout=True)
    ax.imshow(rgb, interpolation="nearest"); ax.set_axis_off()
    ax.scatter(npts.y, npts.x, s=18, facecolors="none", edgecolors="#ffffff", linewidths=0.5)
    if labels:
        for p in npts:
            ax.text(p.y + 2, p.x - 2, p.label, color="#ffff99", fontsize=4.5)
    ax.set_title(f"NeuroPAL RGB max projection (R={names[0]}, G={names[1]}, B={names[2]}; stored X vertical) with stored neuron centres"
                 + (" and biological IDs" if labels else ""), fontsize=9)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def compose_frame(vg, vr, seg, lims, pts, t, seconds, scale=1.0):
    panels = []
    for v, lim in ((vg, lims[0]), (vr, lims[1])):
        a = to_uint8(v.max(0).T, lim)
        panels.append(Image.fromarray(a).resize((int(a.shape[1] * scale), int(a.shape[0] * scale)), Image.NEAREST).convert("RGB"))
    W, H = panels[0].size
    canvas = Image.new("RGB", (W * 2, H * 2))
    canvas.paste(panels[0], (0, 0)); canvas.paste(panels[1], (W, 0)); canvas.paste(panels[0], (0, H)); canvas.paste(panels[1], (W, H))
    # mask outlines on the bottom-right panel
    e = mask_edges(seg.max(0).T)
    ys, xs = np.nonzero(e)
    arr = np.asarray(canvas).copy()
    for yy, xx in zip((ys * scale).astype(int) + H, (xs * scale).astype(int) + W):
        arr[yy:yy + 1, xx:xx + 1] = (80, 255, 80)
    canvas = Image.fromarray(arr)
    d = ImageDraw.Draw(canvas)
    d.text((4, 4), f"t={t}  ({seconds:.1f}s, in-file 1.67 Hz)  GCaMP", fill=(255, 255, 0))
    d.text((W + 4, 4), "TagRFP-T (reference)", fill=(255, 255, 0))
    d.text((4, H + 4), "GCaMP + label centroids and IDs", fill=(255, 255, 0))
    d.text((W + 4, H + 4), "TagRFP-T + label-mask outlines (max over z)", fill=(255, 255, 0))
    if pts is not None:
        for row in pts.itertuples():
            x, y = row.y * scale, row.x * scale + H  # stored Y is horizontal
            d.ellipse([x - 5, y - 5, x + 5, y + 5], outline=(255, 70, 70))
            d.text((x + 6, y - 6), str(row.label), fill=(255, 200, 200))
    return np.asarray(canvas)


def write_movie(f, tracks, out, start, length, fps, lims, frame_s):
    path = out / "mip_movie.mp4"
    try:
        writer = imageio.get_writer(path, fps=fps, codec="libx264", quality=6, macro_block_size=16)
    except Exception as e:  # pragma: no cover
        print(f"mp4 writer unavailable ({e}); writing GIF instead")
        path = out / "mip_movie.gif"; writer = imageio.get_writer(path, duration=1000 / fps, loop=0)
    with writer:
        for t in range(start, start + length):
            vg, vr, seg = read_frame(f, t, 0), read_frame(f, t, 1), read_seg(f, t)
            writer.append_data(compose_frame(vg, vr, seg, lims, tracks[tracks.t == t], t, t * frame_s))
    print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB, {length} frames @ {fps} fps)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nwb", type=Path, required=True)
    ap.add_argument("--tracks", type=Path, required=True, help="calcium_tracks.csv written by inspect_kk.py")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--slice-frames", type=int, nargs="+", default=[0, 800])
    ap.add_argument("--lim-frame", type=int, default=800, help="frame whose 1..99.9th percentiles set the fixed display range per channel")
    ap.add_argument("--movie-start", type=int, default=0)
    ap.add_argument("--movie-len", type=int, default=200)
    ap.add_argument("--fps", type=int, default=10)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.nwb, "r") as f:
        chans = channel_names(f)
        tracks = load_calcium_tracks(f, args.tracks)
        tracks["x"] = tracks.x.round().astype(int); tracks["y"] = tracks.y.round().astype(int); tracks["z"] = tracks.z.round().astype(int)
        rate = float(f["acquisition/CalciumImageSeries/starting_time"].attrs["rate"])
        lims = []
        for c in range(2):
            vl = read_frame(f, args.lim_frame, c)
            lims.append((float(np.percentile(vl, 1)), float(np.percentile(vl, 99.9))))
        print(f"display ranges fixed from frame {args.lim_frame}: {dict(zip(chans, [tuple(round(x) for x in l) for l in lims]))}")
        rgb, npts, names = neuropal_rgb_mip(f)
        neuropal_png(rgb, npts, names, args.out / "neuropal_rgb_mip.png", labels=False)
        neuropal_png(rgb, npts, names, args.out / "neuropal_rgb_mip_labels.png", labels=True)
        print("wrote NeuroPAL RGB previews")
        for t in args.slice_frames:
            pts = tracks[tracks.t == t]
            vg, vr, seg = read_frame(f, t, 0), read_frame(f, t, 1), read_seg(f, t)
            for overlay in (False, True):
                suffix = "_gt" if overlay else ""
                for v, short, cname in ((vg, "gcamp", chans[0]), (vr, "rfp", chans[1])):
                    depth_slice_grid(v, pts, f"{NAME} t={t} ({t / rate:.1f} s), {cname}; circles = label centroids (rounded)" if overlay else f"{NAME} t={t} ({t / rate:.1f} s), {cname}",
                                     args.out / f"depth_slices_t{t}_{short}{suffix}.png", overlay)
                mip_pair_png(vg, vr, seg if overlay else None, pts if overlay else None,
                             f"{NAME} t={t} max projections" + (" with label centroids (top) and mask outlines (bottom)" if overlay else ""),
                             args.out / f"mip_t{t}{suffix}.png", overlay, lims)
            print(f"wrote depth slices + MIPs for t={t}")
        write_movie(f, tracks, args.out, args.movie_start, args.movie_len, args.fps, lims, 1 / rate)
    return 0


if __name__ == "__main__":
    sys.exit(main())
