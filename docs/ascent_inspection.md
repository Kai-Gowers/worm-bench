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
- **Inferred**: a conclusion drawn from measurements plus an assumption
  (e.g. that the documented frame rate is right); stated separately from the
  measurement it rests on.
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
  moments. Dips that appear in both colours at once cannot be calcium
  activity, which would affect only `c1`; a plausible explanation is brief
  motion shifting a neuron relative to its annotated position, but this was
  not checked against the track displacements (Estimate).
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
  two channels' detrended traces is 0.65, so a large shared component (motion,
  focus drift, or noise; cause not determined) remains in both channels.

- **Measured:** the six most channel-asymmetric neurons rise sharply at
  frames ≈200, ≈499 and ≈797 (medians of the three largest positive jumps per
  trace). The spacing is ≈298 frames. The reference channel shows no
  corresponding change at those frames.
- **Inferred:** at the documented 3.3 Hz, 298 frames is 90 s, the length of
  the documented 60 s buffer + 30 s stimulus cycle. The simplest reading is
  that these are stimulus-evoked responses with stimulus onsets near those
  frames. No stimulus or valve timing is recorded in the files, so this is an
  inference from the responses, not a verified fact.

**Conclusion.** Measured: `c1` is the more temporally variable channel per
neuron, and the asymmetry is one-sided (no neuron is twice as variable in
`c0`). Inferred: `c1` carries calcium activity and `c0` is the stable
reference, consistent with the documentation. The population-level difference
is small; roughly a quarter of tracked neurons show large `c1`-only
fluctuations, and `previews/ascent/neuron_traces.png` shows the clearest
examples.

**Unannotated bright region (observation, cause not determined):** the `c1`
channel shows a dense, granular bright region in the upper-right of the field
that is absent from `c0` and has no annotations (see `mip_t0.png`). Possible
explanations include gut autofluorescence, which commonly appears in the green
channel, or non-neuronal tissue; no test was done to distinguish them.

## 5. Voxel size and timing

| Quantity | Value | Status |
|---|---|---|
| xy pixel size | 0.243 µm | **Documented** (Zenodo + ASCENT README); not in file |
| z step | 1.5 µm | **Documented**; not in file |
| Field of view | 512 × 0.243 = 124 µm square | Estimate from documented values |
| Depth span | 8 steps × 1.5 µm = 12 µm between first and last slice | Estimate from documented values |
| Volume rate | 3.3 Hz (one 2-channel, 9-slice volume every ~0.30 s) | **Documented**; not in file. Consistency check (Inferred): the measured ≈298-frame spacing of calcium responses equals 90 s at 3.3 Hz, matching the documented cycle length, if the responses are stimulus-locked |
| Duration | 1100 / 3.3 ≈ 333 s ≈ 5.6 min | Estimate from documented values |
| Per-frame timestamps | none | **Verified absent** |
| Stimulus schedule | 60 s buffer then 30 s *E. coli* OP50 supernatant, repeated | **Documented**. Phase relative to frame 0 is **Inferred**, not verified: calcium responses begin at frames ≈200, 499, 797 (≈60, 151, 241 s at 3.3 Hz), which would fit a cycle starting at frame 0 |
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
annotated). Spatial completeness is **unresolved**: 96 neurons are tracked,
and a casual look at the previews shows bright cells without circles, but no
systematic visual audit was done and the documentation does not say whether
the 96 are all identifiable neurons, a chosen subset, or all cells that could
be followed reliably. A worm's head contains roughly 180–190 neurons, but
how many lie within this 9-slice stack is not known.

**Provenance (Documented):** tracks were produced by ZephIR and then manually
proofread and corrected. They are "corrected algorithmic" rather than purely
manual annotations.

**Motion (Verified from the tracks):** the median in-plane displacement of a
neuron between consecutive frames is 1 voxel (≈0.24 µm); 49 % of steps are
exactly zero; the 99th percentile is 7.6 voxels and the maximum 40 voxels
(≈10 µm); depth never changes by more than 2 slices in one step. These are
displacements of the annotated points; they indicate a nearly stationary
preparation with occasional larger jumps whose cause (worm twitch, stage or
focus drift, or annotation correction) was not determined.

**Depth distribution (counts Verified, interpretation an open hypothesis):**
annotations pile up on the top and bottom slices:

| z | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| rows (thousands) | 25.0 | 7.8 | 7.9 | 7.7 | 3.0 | 10.2 | 11.6 | 10.5 | 22.0 |

Almost half of all points sit on z = 0 or z = 8. One hypothesis is that the
9-slice stack does not cover the full thickness of the head, so neurons whose
true centre lies above or below the stack were annotated on the nearest edge
slice (clipped z). Alternatives, such as the stack genuinely ending at dense
layers of cells, or the annotation tool snapping to edge slices, have not been
ruled out. This stays an open question (section 10); no scoring rule should be
derived from it yet.

## 7. Are the annotations in bounds and aligned with the images?

Two different questions are tested here. **Spatial alignment** asks whether
each annotated (z, y, x) lands on a neuron in the frame it belongs to.
**Track-identity correctness** asks whether a given `TrackID` follows the
same neuron from frame to frame. Only the first was tested quantitatively.

**Spatial alignment (Verified**, `inspect_ascent.py` sections 5 and 6):

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
  sit on neuron bodies in the frames inspected.

Conclusion for spatial alignment: coordinates are voxel indices in the same
`(z, y, x)` order as the HDF5 datasets, with no offset or scaling.

**Track-identity correctness (not verified):** no test was run to check that
a `TrackID` stays on the same cell across frames. Indirect evidence is
consistent with it: per-step displacements are small (section 6), and the
per-track brightness traces are smooth (`neuron_traces.png`), which would be
unlikely if IDs swapped frequently between neighbouring cells. A systematic
check (e.g. watching the overlay movie track by track, or testing for ID
swaps at the large-displacement steps) remains to be done.

## 8. Previews

All under `previews/ascent/` (git-ignored, on scratch). Generated by
`scripts/make_previews.py` in 40 s with 294 MB peak memory.

| File | Content |
|---|---|
| `depth_slices_t{0,550}_c{0,1}.png` (+ `_gt`) | 3×3 grid of the 9 z-slices of one frame, one channel, with/without GT circles (~0.8–0.9 MB each) |
| `mip_t{0,550}.png` (+ `_gt`) | max projection over z, `c0` and `c1` side by side (~0.6 MB) |
| `mip_movie_c0c1.mp4` | frames 150–299 (≈45 s of recording at the documented 3.3 Hz) as max projections, `c0` left, `c1` right, played at 10 fps, i.e. **≈3× faster than real time** (15 s playback) (27 MB) |
| `mip_movie_c0c1_gt.mp4` | same with GT circles and track IDs, same ≈3× speed-up (28 MB) |
| `per_frame_means.csv/.png` | per-frame mean brightness, whole volume and at GT positions, both channels |
| `per_frame_means_tracks.csv` | brightness of every track in every frame and channel (long format) |
| `neuron_traces.png` | detrended traces of the neurons that differ most between channels |
| `inspection_summary.json` | machine-readable output of the inspection script |

Display intensities are stretched between the 1st and 99.9th percentile of the
first frame and kept fixed, so brightness changes over time reflect the data
rather than per-frame rescaling. The time labels burned into the movie assume
the documented 3.3 Hz.

## 9. Inventory entry (draft)

| Field | Value |
|---|---|
| Source | Zenodo 17561700 v1.0, Lu lab (Georgia Tech), ASCENT preprint |
| Original ID | `InHouse_Opterra_CH_RG_Day1_OP50Stim_1100frames` |
| Format | HDF5, `t{n}/c{m}` groups, uint16, gzip |
| Dimensions | T=1100, C=2, Z=9, Y=512, X=512 (Verified) |
| Voxel size | 0.243 × 0.243 × 1.5 µm (Documented) |
| Rate / duration | 3.3 Hz / ≈333 s (Documented / Estimate) |
| Channels | c0 mNeptune reference, c1 GCaMP6 calcium (Documented; per-neuron variability is consistent with this) |
| Condition | immobilized in microfluidic device, chemical stimulus pulses (Documented) |
| Movement | annotated points move a median of 1 voxel per step, with occasional jumps up to 40 voxels (Verified); cause of the jumps not determined |
| Preprocessing | Unknown; file described as "raw" (Documented) |
| Annotations | 96 tracked neurons × 1100 frames = 105,599 points, integer voxel coords (Verified); completeness relative to visible cells unresolved |
| Annotation provenance | ZephIR tracks, manually proofread (Documented) |
| Biological identities | none |
| License | CC-BY-4.0 (Documented) |
| Verification status | files downloaded, checksums match, structure and spatial alignment verified; track-identity correctness and annotation completeness not yet audited |

## 10. Open questions and uncertainties

1. **Exact timing.** 3.3 Hz is documented; no timestamps exist in the file.
   The ≈298-frame spacing of calcium responses is consistent with 3.3 Hz if
   the responses are stimulus-locked, but that is circular as a check of the
   rate. Whether the interval is constant, and whether channels and slices are
   acquired simultaneously or sequentially, is Unknown.
2. **Stimulus phase** is Inferred from the responses (onsets at frames ≈200,
   499, 797), not verified: valve timing is not recorded, response latency
   and flow delay are unknown, and the responses could in principle have
   another periodic cause.
3. **Depth coverage.** Almost half the annotations sit on the edge slices.
   Clipped z coordinates is one hypothesis; it has not been tested, and no
   evaluation rule should depend on it until it is.
4. **Annotation completeness.** Whether the 96 tracks cover all neurons that
   are identifiable in this stack, or a subset, is unresolved. Settling it
   needs either documentation from the authors or a systematic visual audit
   of unannotated bright cells.
4b. **Track-identity correctness** has not been tested directly (section 7).
5. **Track 88** is missing one frame; whether that is a proofreading decision
   or an omission is Unknown.
6. **Raw vs processed.** The record calls the file "raw volumetric imaging
   data", but we could not confirm the absence of registration, cropping or
   background subtraction. The baseline of ~320 counts looks like an
   unprocessed camera offset.
7. **Very bright isolated pixels.** A few pixels reach ~24,000 counts in `c0`
   while neuron peaks are ~8,000. Camera hot pixels are one possible
   explanation; a small very bright structure or a saturated cell are others.
   Not checked further.
8. **Preprint full text** could not be fetched (bioRxiv rate-limits automated
   access), so documented values come from the Zenodo description and the
   GitHub README only.
9. **Scope of the record.** This Zenodo record contains exactly one recording
   (one movie plus its GT and prediction CSVs). The ASCENT paper also uses the
   NeRVE freely-moving dataset and an unnamed light-sheet dataset, which are
   not included here and would have to be obtained from their own sources.
