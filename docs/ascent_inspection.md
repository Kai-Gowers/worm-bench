# ASCENT in-house recording: inspection notes

Dataset: *A device-restricted C. elegans whole-brain recording dataset for
ASCENT: annotation-free neuron tracking benchmark* (Pritchard, Han, Lu;
Zenodo record 17561700, version 1.0, published 2025-11-09, license CC-BY-4.0).
Companion preprint: Han et al., *ASCENT: Annotation-free Self-supervised
Contrastive Embeddings for 3D Neuron Tracking*, bioRxiv 10.1101/2025.07.23.666425.
Code: https://github.com/lu-lab/ascent

Every statement below is tagged:

- **Verified**: measured directly from the downloaded files by `scripts/inspect_ascent.py`.
- **Documented**: stated in the Zenodo record description or the ASCENT GitHub README, not checkable from the files.
- **Estimate**: our interpretation of the evidence.
- **Unknown**: not stated anywhere we could find.

## 1. Background for newcomers

**What was imaged.** A single adult *C. elegans* worm, held still inside a
microfluidic chip (a tiny plastic channel) while a liquid stimulus is switched
on and off around its nose. The microscope looks at the worm's head, where
most of its ~300 neurons are packed into a cluster called the head ganglion.

**Two colours, two jobs.** The worm expresses two fluorescent proteins in all
its neurons:

- **GCaMP6** (green) gets brighter when calcium enters a neuron, which happens
  when the neuron fires. Its brightness over time is the *activity signal*.
- **mNeptune** (red) is always on at a roughly constant level. It is a
  *reference* that marks where each neuron is even when it is silent, which is
  what a tracker needs.

**Volumes, slices, frames.** The microscope takes a stack of 9 images at
different depths (z-slices) every time step. One such stack is a *volume* or
*frame*; the recording is 1100 volumes in a row, in both colours. So the data
are conceptually `time × channel × depth × height × width`.

**Detections vs tracks vs identities.** A *detection* is "there is a neuron at
(z, y, x) in frame t". A *track* links detections of the same neuron across
frames with a shared track ID. A *biological identity* (e.g. "AVA") is a name
from the worm's wiring diagram; this dataset provides tracks but **no**
biological identities.

## 2. Files

| File | Size | MD5 (Zenodo, verified after download) |
|---|---|---|
| `InHouse_Opterra_CH_RG_Day1_OP50Stim_1100frames.h5` | 6,236,570,944 B (5.8 GiB) | `9581892c969bb0e7f3d84eba84988133` |
| `InHouse_Opterra_coords_GT.csv` | 2,115,526 B | `37bd1234033d360dfd888b46daa369b1` |
| `InHouse_Opterra_tracks_GT.csv` | 2,548,928 B | `5a7c414b95365bb13cc2674df03d5277` |

Two further files, `*_coords_PRED.csv` and `*_tracks_PRED.csv`, hold ASCENT's
own predictions and were not downloaded; they are method output, not ground
truth.

Files live on scratch (`/scratch/gowers/worm-bench/data/ascent/`, reached via
the git-ignored `data/` symlink). Zenodo serves about 0.3 MB/s per connection
and intermittently returns HTTP 504 pages, so `scripts/download_ascent.sh`
fetches the movie with 16 parallel range requests and verifies the checksum.

## 3. Image structure

**Verified** from the HDF5 file (`inspect_ascent.py`, section 1):

| Property | Value |
|---|---|
| Layout | 1100 top-level groups `t0` … `t1099`, contiguous, nothing else at root |
| Per frame | two datasets `c0` and `c1` |
| Dataset shape | `(9, 512, 512)` for every frame and channel |
| dtype | `uint16` |
| Storage | gzip-compressed, chunks `(2, 64, 128)` |
| Uncompressed size per channel-volume | 4.7 MB (so ~10.4 GB for the whole movie; 6.2 GB on disk) |
| Pixel values | baseline ≈ 320–330 counts, neuron peaks typically 2,000–10,000; a few pixels reach 24,000 |
| HDF5 attributes | **none** at root, group, or dataset level |

**Axis order (Estimate, strongly supported).** The documentation says each
dataset is `Z × Y × X`. The shape `(9, 512, 512)` makes the first axis the
depth axis (9 slices), and the alignment test in section 7 shows that
interpreting the second axis as `y` and the third as `x` places ground-truth
points on bright neurons, whereas swapping them does not. The HDF5 file itself
carries no axis labels.

Because every frame is a separate dataset, reading one frame costs one
decompressed 4.7 MB read; the inspection script never holds more than a
couple of volumes in memory (peak RSS 232 MB for the full inspection,
including a streaming pass over all 1100 frames in about 2 minutes).

## 4. Channels

**Documented:** `c0` = mNeptune (reference), `c1` = GCaMP6 (calcium).
Channel names are **not** stored in the file.

**Verified evidence** (`inspect_ascent.py` sections 3 and 7, plus
`plot_neuron_traces.py`):

- Neurons are bright in both channels at the same places: the per-neuron mean
  brightness correlates at 0.94 between channels, as expected for two
  pan-neuronal labels.
- Both channels bleach (fade) over the recording. Measured at ground-truth
  neuron positions, `c0` falls to 63 % and `c1` to 75 % of its initial
  brightness by the end (see `previews/ascent/per_frame_means.png`).
- After removing that slow trend, the **population-average** brightness of the
  two channels is nearly identical and shares the same sharp dips at the same
  moments. Shared dips in both colours are the signature of brief **motion**
  (a neuron shifting away from its annotated position), not calcium activity.
- Where the channels differ is at the level of **individual neurons**:
  across all 96 tracks, after removing each neuron's own bleaching trend,
  the fluctuation (coefficient of variation) in `c1` is larger than in `c0`
  for 60 % of neurons; for 25 % of neurons it is more than **twice** as large,
  while **no** neuron fluctuates twice as much in `c0` as in `c1`. The 90th
  percentile of per-neuron CV is 0.25 in `c1` versus 0.10 in `c0`, and the
  largest peak excursion above trend is 91 % in `c1` versus 40 % in `c0`
  (90th percentiles). The traces of the most channel-asymmetric neurons show
  slow, multi-second brightening events in `c1` with a flat `c0`, which is what
  calcium transients look like; the median per-neuron correlation between the
  two channels' detrended traces is 0.65, reflecting shared motion and noise.

- The responding neurons rise sharply at the **same three moments**: frames
  ≈200, ≈499 and ≈797 (medians over the six most responsive tracks). The
  spacing of ≈298 frames matches the documented 90 s stimulus cycle (60 s
  buffer + 30 s food odour), which is exactly what stimulus-evoked calcium
  responses should look like. The reference channel shows nothing at those
  moments.

**Conclusion (Verified):** `c1` carries calcium activity and `c0` is the
stable reference, as documented. The population-level difference is subtle
because the worm is immobilized and only roughly a quarter of the tracked
neurons respond strongly to the stimulus; `previews/ascent/neuron_traces.png`
shows the clearest examples.

**Caveat (Estimate):** the GCaMP channel also shows a dense, granular bright
region in the upper-right of the field that is absent from the mNeptune
channel and has no annotations (see `mip_t0.png`). This is most likely gut
autofluorescence, which is common in the green channel, rather than neurons.

## 5. Voxel size and timing

| Quantity | Value | Status |
|---|---|---|
| xy pixel size | 0.243 µm | **Documented** (Zenodo + ASCENT README); not in file |
| z step | 1.5 µm | **Documented**; not in file |
| Field of view | 512 × 0.243 = 124 µm square | Estimate from documented values |
| Depth span | 8 steps × 1.5 µm = 12 µm between first and last slice | Estimate from documented values |
| Volume rate | 3.3 Hz (one 2-channel, 9-slice volume every ~0.30 s) | **Documented**; not in file. **Cross-checked**: stimulus responses recur every ≈298 frames, and 298 frames / 90 s = 3.3 Hz |
| Duration | 1100 / 3.3 ≈ 333 s ≈ 5.6 min | Estimate from documented values |
| Per-frame timestamps | none | **Verified absent** |
| Stimulus schedule | 60 s buffer then 30 s *E. coli* OP50 supernatant, repeated | **Documented**; phase **Verified** from calcium responses: onsets at frames ≈200, 499, 797 (≈60, 151, 241 s), i.e. the cycle starts at frame 0 |
| Objective / camera | 40× 0.75 NA air, Photometrics Evolve 512 Delta EM-CCD, Bruker Opterra II swept-field confocal | **Documented** |
| Strain | ZM9624 (pan-neuronal GCaMP6 + mNeptune) | **Documented** |

The anisotropy is large: one z step equals about 6 xy pixels. The 12 µm depth
span is thinner than a worm's head ganglion (roughly 20–30 µm across), which
matters for section 6.

## 6. Ground-truth annotations

**Verified** from the two CSVs (`inspect_ascent.py` section 4):

| Property | coords_GT | tracks_GT |
|---|---|---|
| Columns | `object_id,t,z,y,x` | `TrackID,ObjectID,t,z,y,x` |
| Rows | 105,599 | 105,599 |
| All values integer | yes (int64) | yes (int64) |
| Time points | 1100 (t = 0 … 1099, none missing) | same |
| Points per frame | 95–96 (median 96) | same |
| Tracks | – | 96 distinct `TrackID`s (values 0–118, with gaps) |
| Track length | – | 95 tracks span all 1100 frames; track 88 has 1099 |
| Duplicate (track, t) pairs | – | 0 |

The two files describe exactly the same points: every `ObjectID` in the
tracks file matches an `object_id` in the coords file with identical
`t, z, y, x`. The `object_id` values are sparse (maximum 124,239 for 105,599
rows), which suggests they are leftover IDs from a larger automatic detection
set (ZephIR output, per the documentation) that was pruned during manual
proofreading. Coverage is therefore **complete in time** (every frame
annotated) but **partial in space**: 96 tracked neurons, while the images show
visibly more cells (the head contains roughly 180–190 neurons). How the 96
were chosen is **Unknown**.

**Provenance (Documented):** tracks were produced by ZephIR and then manually
proofread and corrected. They are "corrected algorithmic" rather than purely
manual annotations.

**Motion (Verified from the tracks):** the median in-plane displacement of a
neuron between consecutive frames is 1 voxel (≈0.24 µm); 49 % of steps are
exactly zero; the 99th percentile is 7.6 voxels and the maximum 40 voxels
(≈10 µm); depth never changes by more than 2 slices in one step. This is a
nearly stationary worm with occasional twitches, i.e. an easy tracking case
compared with freely moving animals.

**Depth distribution (Verified, interpretation Estimate):** annotations pile
up on the top and bottom slices:

| z | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| rows (thousands) | 25.0 | 7.8 | 7.9 | 7.7 | 3.0 | 10.2 | 11.6 | 10.5 | 22.0 |

Almost half of all points sit on z = 0 or z = 8. The likely explanation is that
the 9-slice stack does not cover the full thickness of the head, so neurons
whose true centre lies above or below the stack were annotated on the nearest
edge slice. For benchmarking, z coordinates on the edge slices should be
treated as clipped, not exact.

## 7. Are the annotations in bounds and aligned with the images?

**Verified** (`inspect_ascent.py` sections 5 and 6):

- Every point satisfies 0 ≤ z ≤ 8, 0 ≤ y, x ≤ 511. Actual ranges: z 0–8,
  y 103–366, x 9–508. Zero out-of-range rows.
- Because the field is square, bounds alone cannot distinguish `y` from `x`,
  so the script compares image brightness at the annotated points with
  brightness at random points in the same region, on 5 frames spread over the
  recording. Background-subtracted brightness at GT points is **24–83× higher**
  than at random points in both channels; with `y` and `x` swapped the ratio
  drops to about 1 (no better than random).
- 70–94 % of GT points lie within 3 pixels of the brightest pixel of their
  11×11 neighbourhood in the same slice.
- Visually (`previews/ascent/*_gt.png`, `mip_movie_c0c1_gt.mp4`), red circles
  sit on neuron bodies throughout.

Conclusion: coordinates are voxel indices in the same `(z, y, x)` order as the
HDF5 datasets, with no offset or scaling, and track IDs follow the same cells
over time.

## 8. Previews

All under `previews/ascent/` (git-ignored, on scratch). Generated by
`scripts/make_previews.py` in 40 s with 294 MB peak memory.

| File | Content |
|---|---|
| `depth_slices_t{0,550}_c{0,1}.png` (+ `_gt`) | 3×3 grid of the 9 z-slices of one frame, one channel, with/without GT circles (~0.8–0.9 MB each) |
| `mip_t{0,550}.png` (+ `_gt`) | max projection over z, `c0` and `c1` side by side (~0.6 MB) |
| `mip_movie_c0c1.mp4` | frames 150–299 (≈45 s of recording) as max projections, `c0` left, `c1` right, 10 fps (27 MB) |
| `mip_movie_c0c1_gt.mp4` | same with GT circles and track IDs (28 MB) |
| `per_frame_means.csv/.png` | per-frame mean brightness, whole volume and at GT positions, both channels |
| `per_frame_means_tracks.csv` | brightness of every track in every frame and channel (long format) |
| `neuron_traces.png` | detrended traces of the neurons that differ most between channels |
| `inspection_summary.json` | machine-readable output of the inspection script |

Display intensities are stretched between the 1st and 99.9th percentile of the
first frame and kept fixed, so brightness changes over time are real.

## 9. Inventory entry (draft)

| Field | Value |
|---|---|
| Source | Zenodo 17561700 v1.0, Lu lab (Georgia Tech), ASCENT preprint |
| Original ID | `InHouse_Opterra_CH_RG_Day1_OP50Stim_1100frames` |
| Format | HDF5, `t{n}/c{m}` groups, uint16, gzip |
| Dimensions | T=1100, C=2, Z=9, Y=512, X=512 (Verified) |
| Voxel size | 0.243 × 0.243 × 1.5 µm (Documented) |
| Rate / duration | 3.3 Hz / ≈333 s (Documented / Estimate) |
| Channels | c0 mNeptune reference, c1 GCaMP6 calcium (Documented, consistent with data) |
| Condition | immobilized in microfluidic device, chemical stimulus pulses (Documented) |
| Movement | very low; median step 1 voxel, occasional twitches up to 40 voxels (Verified) |
| Preprocessing | Unknown; file described as "raw" (Documented) |
| Annotations | 96 tracked neurons × 1100 frames = 105,599 points, integer voxel coords (Verified) |
| Annotation provenance | ZephIR tracks, manually proofread (Documented) |
| Biological identities | none |
| License | CC-BY-4.0 (Documented) |
| Verification status | files downloaded, checksums match, structure and alignment verified |

## 10. Open questions and uncertainties

1. **Exact timing.** 3.3 Hz is documented and consistent with the measured
   stimulus period, but no timestamps exist in the file; whether the interval
   is perfectly constant, and whether the two channels and nine slices are
   acquired simultaneously or sequentially, is Unknown.
2. **Stimulus phase** is now Verified from the responses (onsets at frames
   ≈200, 499, 797), but the valve timing itself is not recorded in the files,
   so the onset frames carry an uncertainty of a few frames (response latency
   and flow delay).
3. **Depth coverage.** The pile-up of annotations on the edge slices suggests
   the stack does not contain the whole ganglion. Whether edge-slice z values
   are clipped or exact is Unknown.
4. **Neuron selection.** Why 96 neurons (and not all visible cells) were
   tracked is Unknown; the GT should be treated as a reliable subset, not an
   exhaustive segmentation.
5. **Track 88** is missing one frame; whether that is a proofreading decision
   or an omission is Unknown.
6. **Raw vs processed.** The record calls the file "raw volumetric imaging
   data", but we could not confirm the absence of registration, cropping or
   background subtraction. The baseline of ~320 counts looks like an
   unprocessed camera offset.
7. **Hot pixels.** A few pixels reach ~24,000 counts in `c0` while neuron peaks
   are ~8,000; likely camera hot pixels, not checked further.
8. **Preprint full text** could not be fetched (bioRxiv rate-limits automated
   access), so documented values come from the Zenodo description and the
   GitHub README only.
