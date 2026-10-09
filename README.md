# worm-bench

A benchmark-in-progress for tracking individual neurons through 3D
calcium-imaging recordings of *C. elegans*. See `CLAUDE.md` for the project
goals. The project is currently in the data discovery and inspection phase.

## Setup

Python 3.11+ with the packages in `requirements.txt`:

```bash
pip install -r requirements.txt
```

The only non-trivial dependency is `imageio-ffmpeg`, which bundles an ffmpeg
binary used to write MP4 previews.

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
requirements.txt
```
