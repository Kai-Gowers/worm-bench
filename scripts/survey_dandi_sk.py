#!/usr/bin/env python3
"""Survey every NWB asset of the two Kato-lab dandisets used by WormID-Bench
(SK1 = DANDI 000565, SK2 = DANDI 000472) WITHOUT downloading them.

Most assets in both dandisets hold only a NeuroPAL structural stack; a few also
hold a calcium movie (`acquisition/CalciumImageSeries`).  The survey opens each
file over HTTP with remfile + h5py and reads only headers and small datasets by
explicit path (a generic `visititems` walk is far too slow over HTTP):

  has_calcium, calcium shape/dtype/chunks/compression, frame rate, grid spacing,
  number of Seg_tpoint_* plane segmentations and ROIs, label counts,
  NeuroPAL shape / channel order / RGBW mapping, strain, session description.

About 2-15 s per asset (listing the per-time-point groups dominates).

Usage: python -P scripts/survey_dandi_sk.py [--dandiset 000565 --version 0.241009.1504]
                                           [--max N] [--out survey.json]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time

import h5py
import remfile

API = "https://api.dandiarchive.org/api"
VERSIONS = {"000565": "0.241009.1504", "000472": "0.241009.1502"}
CALCIUM = "acquisition/CalciumImageSeries"
SEG = "processing/CalciumActivity/CalciumSeriesSegmentation"
TRACE = "processing/CalciumActivity/SignalDFoF/SignalCalciumImResponseSeries"
NEUROPAL = "acquisition/NeuroPALImageRaw"
NEUROPAL_SEG = "processing/NeuroPAL/NeuroPALSegmentation/NeuroPALNeurons"


def api_json(url: str) -> dict:
    # curl rather than urllib: the host's Python lacks the CA bundle needed for dandiarchive.org
    return json.loads(subprocess.check_output(["curl", "-sS", "-f", url]))


def decode(v):
    return v.decode() if isinstance(v, bytes) else str(v)


def ragged_labels(g: h5py.Group, name: str = "ID_labels") -> list[str]:
    chars = [decode(c) for c in g[name][()]]
    idx = g[f"{name}_index"][()]
    starts = [0] + idx[:-1].tolist()
    return ["".join(chars[a:b]) for a, b in zip(starts, idx)]


def survey_asset(a: dict) -> dict:
    t0 = time.time()
    f = h5py.File(remfile.File(f"{API}/assets/{a['asset_id']}/download/"), "r")
    r = dict(path=a["path"].split("/")[0], asset_id=a["asset_id"], size_GB=round(a["size"] / 1e9, 3),
             sha256=a.get("digest", {}).get("dandi:sha2-256"))
    r["identifier"] = decode(f["identifier"][()]); r["session_start"] = decode(f["session_start_time"][()])
    r["session_description"] = decode(f["session_description"][()])
    r["strain"] = decode(f["general/subject/strain"][()]); r["growth_stage"] = decode(f["general/subject/growth_stage"][()])
    r["acquisition"] = list(f["acquisition"].keys()); r["processing"] = list(f["processing"].keys())
    r["has_intervals"] = "intervals" in f
    r["has_calcium"] = CALCIUM in f
    npal = f[f"{NEUROPAL}/data"]
    r["neuropal_shape"] = list(npal.shape); r["neuropal_dtype"] = str(npal.dtype)
    r["neuropal_grid_spacing_um"] = f["general/optophysiology/NeuroPALImVol/grid_spacing"][()].round(4).tolist()
    r["neuropal_channel_order"] = [decode(c) for c in f["general/optophysiology/NeuroPALImVol/order_optical_channels/channels"][()]]
    r["neuropal_RGBW"] = f[f"{NEUROPAL}/RGBW_channels"][()].tolist()
    if NEUROPAL_SEG in f:
        labs = ragged_labels(f[NEUROPAL_SEG])
        r["neuropal_n_rois"] = len(labs); r["neuropal_n_labelled"] = sum(1 for l in labs if l.strip())
    r["has_processed_image"] = "processing/ProcessedImage/ProcessedImage/data" in f
    if r["has_processed_image"]:
        r["processed_shape"] = list(f["processing/ProcessedImage/ProcessedImage/data"].shape)
    if r["has_calcium"]:
        ci = f[CALCIUM]
        r["calcium_shape"] = list(ci["data"].shape); r["calcium_dtype"] = str(ci["data"].dtype)
        r["calcium_chunks"] = ci["data"].chunks; r["calcium_compression"] = ci["data"].compression
        r["rate_Hz"] = float(ci["starting_time"].attrs["rate"]); r["dimension"] = ci["dimension"][()].tolist()
        r["calcium_grid_spacing_um"] = f["general/optophysiology/CalciumImVol/grid_spacing"][()].round(4).tolist()
        r["calcium_channels"] = [decode(c) for c in f["general/optophysiology/CalciumImVol/order_optical_channels/channels"][()]]
        r["calcium_description"] = decode(ci.attrs.get("description", ""))
        keys = [k for k in f[SEG].keys() if k.startswith("Seg_tpoint_")]
        tp = sorted(int(k.split("_")[-1]) for k in keys)
        r["n_seg_tpoints"] = len(tp); r["seg_tpoint_range"] = [tp[0], tp[-1]] if tp else None
        r["seg_tpoints_contiguous"] = tp == list(range(tp[0], tp[-1] + 1)) if tp else None
        s0 = f[f"{SEG}/Seg_tpoint_0"]; r["n_rois_t0"] = int(s0["id"].shape[0]); r["n_voxels_t0"] = int(s0["voxel_mask"].shape[0])
        r["calcium_activity_members"] = list(f["processing/CalciumActivity"].keys())
        if TRACE in f:
            r["trace_shape"] = list(f[f"{TRACE}/data"].shape); r["trace_unit"] = decode(f[f"{TRACE}/data"].attrs.get("unit", ""))
        lp = "processing/CalciumActivity/NeuronIDs/labels"
        if lp in f:
            labels = [decode(l) for l in f[lp][()]]
            r["n_labels"] = len(labels); r["n_nonempty_labels"] = sum(1 for l in labels if l.strip())
            r["n_unique_nonempty_labels"] = len({l for l in labels if l.strip()})
    r["secs"] = round(time.time() - t0, 1)
    return r


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dandiset", default="000565", choices=sorted(VERSIONS)); ap.add_argument("--version", default=None)
    ap.add_argument("--max", type=int, default=None, help="survey only the first N assets (by path)")
    ap.add_argument("--calcium-only", action="store_true", help="skip assets without a calcium series (decided from the file header)")
    ap.add_argument("--out", default=None, help="write the survey as JSON here")
    args = ap.parse_args(argv)
    version = args.version or VERSIONS[args.dandiset]
    assets = api_json(f"{API}/dandisets/{args.dandiset}/versions/{version}/assets/?page_size=200")["results"]
    assets = sorted(assets, key=lambda a: a["path"])[: args.max]
    print(f"DANDI {args.dandiset} v{version}: {len(assets)} assets, {sum(a['size'] for a in assets) / 1e9:.1f} GB")
    rows = []
    for a in assets:
        r = survey_asset(a); rows.append(r)
        if r["has_calcium"]:
            print(f"{r['path']}: CALCIUM T={r['calcium_shape'][0]} XYZ={r['calcium_shape'][1:]} rate={r['rate_Hz']} grid={r['calcium_grid_spacing_um']} "
                  f"seg_tpoints={r['n_seg_tpoints']} rois={r['n_rois_t0']} labelled={r.get('n_nonempty_labels')} "
                  f"neuropal={r['neuropal_shape']} strain={r['strain']} ({r['size_GB']} GB, {r['secs']} s)")
        else:
            print(f"{r['path']}: neuropal only {r['neuropal_shape']} rois={r.get('neuropal_n_rois')} labelled={r.get('neuropal_n_labelled')} "
                  f"strain={r['strain']} ({r['size_GB']} GB, {r['secs']} s)")
    if args.out:
        json.dump(rows, open(args.out, "w"), indent=1); print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
