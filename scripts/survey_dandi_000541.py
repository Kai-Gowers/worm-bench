#!/usr/bin/env python3
"""Survey every NWB asset of DANDI 000541 (EY dataset) WITHOUT downloading it.

Opens each file over HTTP with remfile + h5py and reads only HDF5 headers and
small datasets (shapes, dtypes, rate, grid spacing, number of segmentation
time points and ROIs, labels, stimulus intervals).  About 7 s per asset.

Usage: python -P scripts/survey_dandi_000541.py [--version 0.241009.1457] [--max N] [--out survey.json]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time

import h5py
import remfile

API = "https://api.dandiarchive.org/api"


def api_json(url: str) -> dict:
    # curl rather than urllib: the host's Python lacks the CA bundle needed for dandiarchive.org
    return json.loads(subprocess.check_output(["curl", "-sS", "-f", url]))


def survey_asset(a: dict) -> dict:
    t0 = time.time()
    f = h5py.File(remfile.File(f"{API}/assets/{a['asset_id']}/download/"), "r")
    r = dict(path=a["path"].split("/")[0], asset_id=a["asset_id"], size_GB=round(a["size"] / 1e9, 3))
    ci = f["acquisition/CalciumImageSeries"]
    r["calcium_shape"] = list(ci["data"].shape); r["calcium_dtype"] = str(ci["data"].dtype)
    r["calcium_chunks"] = ci["data"].chunks; r["calcium_compression"] = ci["data"].compression
    r["rate_Hz"] = float(ci["starting_time"].attrs["rate"]); r["dimension"] = ci["dimension"][()].tolist()
    r["grid_spacing_um"] = f["general/optophysiology/CalciumImVol/grid_spacing"][()].round(4).tolist()
    npal = f["acquisition/NeuroPALImageRaw/data"]
    r["neuropal_shape"] = list(npal.shape); r["neuropal_dtype"] = str(npal.dtype)
    seg = f["processing/CalciumActivity/CalciumSeriesSegmentation"]
    tp = sorted(int(k.split("_")[-1]) for k in seg.keys() if k.startswith("Seg_tpoint_"))
    r["n_seg_tpoints"] = len(tp); r["seg_tpoint_range"] = [tp[0], tp[-1]]; r["seg_tpoints_contiguous"] = tp == list(range(tp[0], tp[-1] + 1))
    s0 = seg["Seg_tpoint_0"]; r["n_rois_t0"] = int(s0["id"].shape[0]); r["n_voxels_t0"] = int(s0["voxel_mask"].shape[0])
    r["has_dnmf_seg"] = "CalciumSeriesSegmentationdNMF" in f["processing/CalciumActivity"]
    r["trace_shape"] = list(f["processing/CalciumActivity/SignalRawFluor/SignalCalciumImResponseSeries/data"].shape)
    labels = [l.decode() if isinstance(l, bytes) else str(l) for l in f["processing/CalciumActivity/NeuronIDs/labels"][()]]
    r["n_labels"] = len(labels); r["n_nonempty_labels"] = sum(1 for l in labels if l.strip()); r["n_unique_nonempty_labels"] = len({l for l in labels if l.strip()})
    st = f["intervals/chemical_stimuli"]
    r["stim"] = [(s.decode(), float(x), float(y)) for s, x, y in zip(st["stimulus"][()], st["start_time"][()], st["stop_time"][()])]
    r["neuropal_n_rois"] = int(f["processing/NeuroPAL/NeuroPALSegmentation/NeuroPALNeurons/id"].shape[0])
    r["strain"] = f["general/subject/strain"][()].decode(); r["growth_stage"] = f["general/subject/growth_stage"][()].decode()
    r["session_start"] = f["session_start_time"][()].decode(); r["sha256"] = None
    r["secs"] = round(time.time() - t0, 1)
    return r


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dandiset", default="000541"); ap.add_argument("--version", default="0.241009.1457")
    ap.add_argument("--max", type=int, default=None, help="survey only the first N assets (by path)")
    ap.add_argument("--out", default=None, help="write the survey as JSON here")
    args = ap.parse_args(argv)
    assets = api_json(f"{API}/dandisets/{args.dandiset}/versions/{args.version}/assets/?page_size=200")["results"]
    assets = sorted(assets, key=lambda a: a["path"])[: args.max]
    print(f"{len(assets)} assets, {sum(a['size'] for a in assets) / 1e9:.1f} GB")
    rows = []
    for a in assets:
        r = survey_asset(a); rows.append(r)
        print(f"{r['path']}: T={r['calcium_shape'][0]} XYZ={r['calcium_shape'][1:]} rate={r['rate_Hz']} seg_tpoints={r['n_seg_tpoints']} "
              f"rois={r['n_rois_t0']} labelled={r['n_nonempty_labels']} stim={[s[0] for s in r['stim']]} ({r['secs']} s)")
    if args.out:
        json.dump(rows, open(args.out, "w"), indent=1); print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
