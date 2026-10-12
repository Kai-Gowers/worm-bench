# worm-bench

A benchmark-in-progress for tracking individual neurons through 3D
calcium-imaging recordings of *C. elegans*. See `CLAUDE.md` for the project
goals. The project is currently in the data discovery and inspection phase.

## Setup

Python 3.11+ with the packages in `requirements.txt`:

```bash
pip install -r requirements.txt
```

`imageio-ffmpeg` bundles an ffmpeg binary used to write MP4 previews.
`pynwb` is used only for an optional read check of NWB files, and `remfile`
only by the remote survey script.

### Data location

Raw recordings and generated previews are **never committed**. On the
processing host they live on scratch, and the repo contains git-ignored
symlinks to them:

```
data     -> /scratch/gowers/worm-bench/data
previews -> /scratch/gowers/worm-bench/previews
```

Recreate the symlinks on a fresh clone with:

```bash
ln -sfn /scratch/gowers/worm-bench/data data
ln -sfn /scratch/gowers/worm-bench/previews previews
```

Scripts never run from inside the data directory, and Python is invoked with
`-P` so that the current directory is not on the import path.

## ASCENT dataset (Zenodo 17561700)

Findings are written up in [`docs/ascent_inspection.md`](docs/ascent_inspection.md).

### 1. Download

```bash
bash scripts/download_ascent.sh            # -> data/ascent/
```

Downloads the 6.2 GB HDF5 recording and the two ground-truth CSVs and checks
their MD5 sums against the values published by Zenodo. The HDF5 file is fetched
with 16 parallel range requests (`scripts/parallel_fetch.sh`) because Zenodo
serves roughly 0.3 MB/s per connection. The script is resumable.

### 2. Inspect

```bash
python -P scripts/inspect_ascent.py \
    --h5 data/ascent/InHouse_Opterra_CH_RG_Day1_OP50Stim_1100frames.h5 \
    --coords data/ascent/InHouse_Opterra_coords_GT.csv \
    --tracks data/ascent/InHouse_Opterra_tracks_GT.csv \
    --out-json previews/ascent/inspection_summary.json
```

Reads only a few dozen frames from the HDF5 file and prints: file structure,
dimensions, dtype, frame count, in-file metadata, which channel is calcium
versus reference, annotation coverage, bounds checks, and an alignment test
between ground-truth points and bright spots in the images.

Add `--per-frame-means previews/ascent/per_frame_means.csv` to stream through
all 1100 frames once (one at a time, about 2 minutes) and record per-frame
intensities plus the brightness of every ground-truth track in both channels
(`per_frame_means_tracks.csv`). Two small scripts plot those outputs:

```bash
python -P scripts/plot_per_frame_means.py previews/ascent/per_frame_means.csv previews/ascent/per_frame_means.png
python -P scripts/plot_neuron_traces.py previews/ascent/per_frame_means_tracks.csv previews/ascent/neuron_traces.png
```

The second prints per-neuron statistics that show which channel carries
calcium activity and plots the most responsive neurons.

### 3. Previews

```bash
python -P scripts/make_previews.py \
    --h5 data/ascent/InHouse_Opterra_CH_RG_Day1_OP50Stim_1100frames.h5 \
    --tracks data/ascent/InHouse_Opterra_tracks_GT.csv \
    --out previews/ascent
```

Writes depth-slice grids, max-projection stills, and a short max-projection
movie, each with and without ground-truth points overlaid.

## EY dataset (DANDI 000541, one recording)

Findings are written up in [`docs/ey_inspection.md`](docs/ey_inspection.md),
with a comparison to ASCENT at the end.

### 0. Survey the archive without downloading

```bash
python -P scripts/survey_dandi_000541.py --out previews/ey/survey.json
```

Opens each of the 21 NWB assets over HTTP (`remfile`) and reads only headers
and small datasets: shapes, rate, number of annotated time points and ROIs,
labels, stimulus timing. About 7 s per asset, no pixel data.

### 1. Download

```bash
bash scripts/download_ey.sh sub-20190928-13       # -> data/ey/
```

Looks the asset up through the DANDI API, downloads the 1.46 GB NWB file
with a resumable `curl` (DANDI serves ~50 MB/s on one connection) and checks
the SHA-256 published in the asset metadata, which is saved beside the file.

### 2. Inspect

```bash
python -P scripts/inspect_ey.py \
    --nwb data/ey/sub-20190928-13_ses-20190928_ophys.nwb \
    --out-json previews/ey/inspection_summary.json \
    --tracks-csv previews/ey/calcium_tracks.csv \
    --per-frame-means previews/ey/per_frame_means.csv
```

Same sections as the ASCENT script (structure, in-file metadata, channel
behaviour, annotations, bounds, alignment, streaming pass) adapted to the
NWB layout: per-frame `PlaneSegmentation` tables become one long track
table, biological labels are attached, the NeuroPAL stack and stored traces
are compared with the calcium images, and responses are checked against the
stimulus intervals stored in the file. The streaming pass over all 962
frames takes about 10 s. `--skip-pynwb` skips the optional pynwb read check.

```bash
python -P scripts/plot_ey_traces.py previews/ey --fps 4 \
    --stim nacl:60.5:70.5 --stim pentanedione:120.5:130.5 --stim butanone:180.5:190.5
```

### 3. Previews

```bash
python -P scripts/make_previews_ey.py \
    --nwb data/ey/sub-20190928-13_ses-20190928_ophys.nwb --out previews/ey
```

Depth-slice grids, calcium and NeuroPAL max projections, the NeuroPAL stack
with its neuron names, and a short max-projection movie with ROI points.

## KK dataset (DANDI 000692, one recording)

Findings are written up in [`docs/kk_inspection.md`](docs/kk_inspection.md),
with a comparison to ASCENT and EY at the end.

### 0. Survey the archive without downloading

```bash
python -P scripts/survey_dandi_000692.py --out previews/kk/survey.json
```

Opens each of the 9 NWB assets (~9 GB each) over HTTP and reads only headers
and small datasets, about 1.5 s per asset. The KK files differ from EY: two
calcium channels (GCaMP and TagRFP-T), gzip-chunked images, and a
label-volume segmentation series instead of per-frame point tables, so the
survey reads everything by explicit path.

### 1. Download

```bash
bash scripts/download_kk.sh sub-231013-02       # -> data/kk/ (9.4 GB, ~3 min)
```

### 2. Inspect

```bash
python -P scripts/inspect_kk.py \
    --nwb data/kk/sub-231013-02_ses-20231013T111100_ophys.nwb \
    --out-json previews/kk/inspection_summary.json \
    --tracks-csv previews/kk/calcium_tracks.csv \
    --per-frame-means previews/kk/per_frame_means.csv
```

Same sections as the ASCENT and EY scripts. The annotations are label
volumes (voxel value i = ROI i), so the script first derives per-frame
centroids, voxel counts and bounding boxes for every label (streams the
1,646 label volumes once, ~2 min, cached in `calcium_tracks.csv`). The
two-channel calcium-vs-reference test from ASCENT applies again; the
label-volume series also allows a time-index test (which label volume best
fits image t on the frames with the most motion) that EY could not support.
The streaming pass over all frames (images and label volumes) takes about
20 min and writes mask-mean traces for both channels, which are compared
with the stored traces.

```bash
python -P scripts/plot_kk_traces.py previews/kk --fps 1.67
```

### 3. Previews

```bash
python -P scripts/make_previews_kk.py \
    --nwb data/kk/sub-231013-02_ses-20231013T111100_ophys.nwb \
    --tracks previews/kk/calcium_tracks.csv --out previews/kk \
    --movie-start 1050 --movie-len 200
```

Depth-slice grids for both channels, GCaMP and TagRFP-T max projections
with label centroids and mask outlines, the NeuroPAL stack with its neuron
names, and a 2x2 max-projection movie over a segment where the worm moves.

## Layout

```
CLAUDE.md                  project goals and conventions
docs/ascent_inspection.md  findings for the ASCENT recording
scripts/download_ascent.sh checksum-verified download from Zenodo
scripts/parallel_fetch.sh  parallel range downloader used for large files
scripts/inspect_ascent.py  frame-by-frame inspection of the HDF5 + CSVs
scripts/make_previews.py   PNG / MP4 previews
scripts/plot_per_frame_means.py   population brightness over time
scripts/plot_neuron_traces.py     per-neuron channel comparison
docs/ey_inspection.md      findings for the EY (DANDI 000541) recording
scripts/survey_dandi_000541.py    remote header-only survey of all 21 EY assets
scripts/download_ey.sh     SHA-256-verified download of one EY asset
scripts/inspect_ey.py      NWB inspection (imports helpers from inspect_ascent.py)
scripts/make_previews_ey.py       PNG / MP4 previews for the EY recording
scripts/plot_ey_traces.py         population and per-neuron traces vs stored traces
docs/kk_inspection.md      findings for the KK (DANDI 000692) recording
docs/author_questions.md   open questions to send to the EY and KK authors
scripts/survey_dandi_000692.py    remote header-only survey of all 9 KK assets
scripts/download_kk.sh     SHA-256-verified download of one KK asset
scripts/inspect_kk.py      structure, label-volume annotations, channels, alignment, time-index test, streaming stats
scripts/make_previews_kk.py       PNG / MP4 previews for the KK recording
scripts/plot_kk_traces.py         population and per-neuron traces vs stored traces
requirements.txt
```
