#!/usr/bin/env python3
"""Survey every NWB asset of DANDI 000692 (KK dataset, Kimura lab) WITHOUT downloading it.

Opens each ~9 GB file over HTTP with remfile + h5py and reads only HDF5 headers
and small datasets.  The KK files differ from EY (000541): two calcium channels,
gzip-chunked image data, and a volumetric label-image segmentation series
(processing/CalciumActivity/CalciumSeriesSegmentation/data) instead of per-frame
PlaneSegmentation point tables.  Everything is read by explicit path; a generic
visititems walk over these files takes minutes over HTTP.

Usage: python -P scripts/survey_dandi_000692.py [--version 0.240402.2118] [--max N] [--out survey.json]
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


def s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def ds_info(d: h5py.Dataset) -> dict:
    return dict(shape=list(d.shape), dtype=str(d.dtype), chunks=list(d.chunks) if d.chunks else None, compression=d.compression)


def survey_asset(a: dict) -> dict:
    t0 = time.time()
    f = h5py.File(remfile.File(f"{API}/assets/{a['asset_id']}/download/"), "r")
    r = dict(path=a["path"].split("/")[0], asset_id=a["asset_id"], size_GB=round(a["size"] / 1e9, 3))
    ci = f["acquisition/CalciumImageSeries"]
    r["calcium"] = ds_info(ci["data"])
    r["rate_Hz"] = float(ci["starting_time"].attrs["rate"]); r["starting_time"] = float(ci["starting_time"][()])
    r["dimension"] = ci["dimension"][()].tolist()
    iv = f["general/optophysiology/CalciumImVol"]
    r["calcium_channels"] = [k for k in iv.keys() if isinstance(iv[k], h5py.Group) and "emission_lambda" in iv[k]]
    r["calcium_channel_order"] = iv["order_optical_channels"].attrs.get("channels", None) if "order_optical_channels" in iv else None
    r["grid_spacing_um"] = iv["grid_spacing"][()].round(4).tolist()
    r["calcium_description"] = s(iv["description"][()])[:200]
    npal = f["acquisition/NeuroPALImageRaw"]
    r["neuropal"] = ds_info(npal["data"]); r["neuropal_RGBW"] = npal["RGBW_channels"][()].tolist()
    r["neuropal_grid_spacing_um"] = f["general/optophysiology/NeuroPALImVol/grid_spacing"][()].round(4).tolist()
    ca = f["processing/CalciumActivity"]
    r["calcium_activity_keys"] = list(ca.keys())
    seg = ca["CalciumSeriesSegmentation"]
    r["seg"] = ds_info(seg["data"]) if "data" in seg else None
    r["seg_rate_Hz"] = float(seg["starting_time"].attrs["rate"]) if "starting_time" in seg else None
    r["seg_description"] = s(seg.attrs.get("description", ""))[:200]
    if "SegmentationVol1" in ca:
        sv = ca["SegmentationVol1"]; r["segvol1_keys"] = list(sv.keys())
        p = sv[list(sv.keys())[0]]
        r["segvol1"] = {k: ds_info(p[k]) for k in p.keys() if isinstance(p[k], h5py.Dataset)}
    nid = ca["NeuronIDs"]
    labels = [s(l) for l in nid["labels"][()]]
    r["n_labels"] = len(labels); r["n_nonempty_labels"] = sum(1 for l in labels if l.strip()); r["n_unique_nonempty_labels"] = len({l for l in labels if l.strip()})
    r["neuronids_keys"] = list(nid.keys()); r["neuronids_description"] = s(nid["description"][()])[:200] if "description" in nid else None
    tr = ca["SignalRawFluor/SignalCalciumImResponseSeries"]
    r["trace"] = ds_info(tr["data"]); r["trace_rate_Hz"] = float(tr["starting_time"].attrs["rate"]) if "starting_time" in tr else None
    nps = f["processing/NeuroPAL/NeuroPALSegmentation/NeuroPALNeurons"]
    r["neuropal_n_rois"] = int(nps["id"].shape[0]); r["neuropal_seg_keys"] = list(nps.keys())
    sub = f["general/subject"]
    r["subject"] = {k: (float(sub[k][()]) if sub[k].dtype.kind == "f" else s(sub[k][()])) for k in sub.keys()}
    r["session_description"] = s(f["session_description"][()])[:300]; r["session_start"] = s(f["session_start_time"][()])
    r["lab"] = s(f["general/lab"][()]); r["institution"] = s(f["general/institution"][()])
    r["related_publications"] = [s(x) for x in f["general/related_publications"][()]]
    r["devices"] = list(f["general/devices"].keys())
    r["secs"] = round(time.time() - t0, 1)
    return r


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dandiset", default="000692"); ap.add_argument("--version", default="0.240402.2118")
    ap.add_argument("--max", type=int, default=None, help="survey only the first N assets (by path)")
    ap.add_argument("--out", default=None, help="write the survey as JSON here")
    args = ap.parse_args(argv)
    assets = api_json(f"{API}/dandisets/{args.dandiset}/versions/{args.version}/assets/?page_size=200")["results"]
    assets = sorted(assets, key=lambda a: a["path"])[: args.max]
    print(f"{len(assets)} assets, {sum(a['size'] for a in assets) / 1e9:.1f} GB")
    rows = []
    for a in assets:
        r = survey_asset(a); rows.append(r)
        print(f"{r['path']}: calcium={r['calcium']['shape']} {r['calcium']['dtype']} rate={r['rate_Hz']} seg={r['seg']['shape'] if r['seg'] else None} "
              f"{r['seg']['dtype'] if r['seg'] else ''} traces={r['trace']['shape']} labels={r['n_nonempty_labels']}/{r['n_labels']} "
              f"npal_rois={r['neuropal_n_rois']} grid={r['grid_spacing_um']} ({r['secs']} s)", flush=True)
        if args.out:
            json.dump(rows, open(args.out, "w"), indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
