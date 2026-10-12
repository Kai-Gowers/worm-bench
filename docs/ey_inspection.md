# EY recording sub-20190928-13 (DANDI 000541): inspection notes

Dataset: *NeuroPAL Microfluidic Chip Images and GCaMP activity* (Yemini,
Nejatbakhsh, Venkatachalam, Lin, Sprague, Samuel, Hobert; DANDI 000541,
version 0.241009.1457, published 2024-10-09, license CC-BY-4.0). The files
accompany two papers: Yemini et al. 2021, *NeuroPAL: A Multicolor Atlas for
Whole-Brain Neuronal Identification in C. elegans* (Cell), and
Nejatbakhsh et al. 2020, *Extracting neural signals from semi-immobilized
animals with deformable non-negative matrix factorization* (dNMF). This is the
"EY" dataset of WormID-Bench / WormND (21 recordings used for the tracking
task). Benchmark code: https://github.com/focolab/WormND

Companion to `docs/ascent_inspection.md`; the section numbering and the
evidence labels are the same so the two can be read side by side:

- **Verified**: measured directly from the downloaded file by `scripts/inspect_ey.py`.
- **Documented**: stated in the NWB metadata fields, the DANDI record, or the WormND repository, not checkable from the images.
- **Inferred**: a conclusion from measurements plus an assumption, stated separately from the measurement.
- **Estimate**: our interpretation of the evidence.
- **Unknown**: not stated anywhere we could find.

A difference in kind from ASCENT: the NWB container stores voxel size, frame
rate, channel names, stimulus timing, strain, microscope and biological
neuron identities *inside the file*, so many values that were "Documented
(not in file)" for ASCENT are "Documented (in file)" here. They are still
metadata written by the authors, not measurements we could repeat.

## 0. Why this asset

All 21 assets were surveyed through the DANDI API and by opening each NWB
file remotely (HDF5 headers and small datasets only, no pixel data;
`scripts/survey_dandi_000541.py`). They share one layout: one single-channel
GCaMP series of 937–963 volumes of 256 × 128 × 21 voxels at 4 Hz, a 5-channel
NeuroPAL stack, per-frame neuron positions for all but the last frame,
166–188 ROIs, and three 10 s chemical stimuli. Sizes are 1.42–1.46 GB each,
30.5 GB in total. Differences between assets are the frame count (962 is the
mode), the ROI count (median 178), how many ROIs carry a biological name
(95–100 %), and the stimulus order (three permutations).

`sub-20190928-13` was chosen because it has the modal frame count (962), a
near-median ROI count (175), every ROI named (several other files have 1–8
unnamed ROIs), and it is part of the WormND tracking split (fold 2 of their
5-fold assignment in `neuron_tracking/*/dataset_split.csv`). Nothing about it
is known to be atypical.

## 1. Background for newcomers

**What was imaged.** One young-adult worm (strain OH16230, a NeuroPAL
strain that also expresses GCaMP pan-neuronally; the GCaMP variant is not
named in the file; Documented in file), held in a microfluidic chip while three chemicals were flowed past
its nose for 10 s each. The field of view covers the whole head: the pharynx
(the feeding organ, with its own small nervous system: I1, I2, M-cells, NSM,
MC), the nerve ring, and the ganglia behind it down to the first ventral-cord
motor neurons (VB1, DB1).

**One colour for activity, five for identity.** Unlike ASCENT there is no
stable reference channel in the movie: the calcium series is GCaMP only. The
reference role is played by a separate, static **NeuroPAL** stack recorded
once in five colours. NeuroPAL is a genetic "colour barcode": each neuron
class expresses a fixed combination of blue, green/orange and red
fluorescent proteins, so a trained annotator can name every neuron from its
colour and position. The authors did that, and the names are stored in the
file for both the NeuroPAL stack and the calcium ROIs.

**Volumes, slices, frames.** 21 z-slices per volume, 962 volumes at 4 Hz
(240.5 s). Data are stored `time × x × y × z`, one channel.

**Detections, tracks, identities.** The file gives, for 961 of the 962
frames, one voxel per neuron for the same 175 neurons in the same order, so
the ROI index is a track ID and each track also has a biological identity.
That is the full set of annotation types listed in CLAUDE.md; ASCENT had
tracks but no identities.

## 2. Files

| File | Size | SHA-256 (DANDI, verified after download) |
|---|---|---|
| `sub-20190928-13/sub-20190928-13_ses-20190928_ophys.nwb` | 1,455,259,796 B (1.36 GiB) | `304f7fde113fa7250051574eb2fbee17fe4230f70ba8bd9f3631509c97b9c7f6` |

DANDI asset id `c4a52d3e-26a9-439c-8965-13320ab6cc43`. The asset metadata
JSON is saved next to the file. Download by `scripts/download_ey.sh` (24 s at
about 58 MB/s on one connection; no parallel range fetch needed, unlike
Zenodo). The file lives on scratch under the git-ignored `data/ey/`.

## 3. Image structure

**Verified** from the HDF5 file (`inspect_ey.py`, section 1):

| Property | Calcium series | NeuroPAL stack |
|---|---|---|
| Path | `acquisition/CalciumImageSeries/data` | `acquisition/NeuroPALImageRaw/data` |
| Shape | `(962, 256, 128, 21)` | `(512, 256, 21, 5)` |
| Axes (Documented in file, `dimension` = `[256, 128, 21]`) | T, X, Y, Z | X, Y, Z, C |
| dtype | `uint16` | `uint64` |
| Storage | contiguous, no chunking, no compression | same |
| Size | 1.32 GB (1.38 MB per volume) | 110 MB |
| Value range | **0–255** in every frame (global max 255) | 0–4095 (12-bit) |
| Fraction of zero voxels | 0.875 at t=0, 0.95–0.97 for most of the movie | not measured |
| `format` attribute | `"raw"`, description "Raw GCaMP series images" | "NeuroPAL structural image" |

Because the dataset is contiguous and uncompressed, one volume is one 1.38 MB
read; the streaming pass over all 962 frames takes 8 s from disk cache with
235 MB peak RSS. (ASCENT: gzip chunks, 4.7 MB per channel-volume, 2 minutes.)

**Pixel values are not raw camera counts (Verified measurement, Inferred
interpretation).** Every voxel in the calcium series is between 0 and 255
although the dtype is 16-bit, 87–97 % of voxels are exactly 0, the median
voxel is 0, and the NWB `unit` attribute says "Voxel gray counts". A camera
produces a non-zero offset and >8-bit values (ASCENT's baseline was ~320
counts). The simplest reading is that the series was background-subtracted
(or thresholded) and rescaled to 8 bits before being written to NWB, despite
the `format: raw` attribute. The NeuroPAL stack, by contrast, has a 12-bit
range and non-zero background. We did not find a statement of the processing
in the metadata (Unknown). The WormND tracking code rescales the data again
to 8 bits per frame before use.

**Spatial relation between the two stacks (Verified):** matching the 175
NeuroPAL neuron positions to the calcium ROI positions at t=0 by name gives
per-axis linear fits NeuroPAL = 1.99·calcium − 1.1 (x), 1.99·calcium − 0.5
(y), 1.05·calcium − 1.6 (z), with residual RMS of 1.9, 1.8 and 0.7 voxels.
So the calcium movie is the NeuroPAL field at half the xy resolution (2 × 2
binning or a 2× smaller camera ROI), with the same z spacing. The two stacks
have the same documented voxel size (0.267 µm) in the file, which cannot be
literally true for both; see section 10.

## 4. Channels

**Documented (in file):** the calcium series has one optical channel,
`GFP-GCaMP`, excitation 488 nm, emission 525 nm (Semrock FF02-525/40). The
NeuroPAL stack has five: mTagBFP2 (405/445), CyOFP1 (488/610), mNeptune 2.5
(561/731), TagRFP-T (561/610) and GFP-GCaMP (488/525); `RGBW_channels =
[2, 1, 0, 3]` maps red/green/blue/white display to mNeptune, CyOFP1,
mTagBFP2 and TagRFP-T.

**What ASCENT's channel test cannot do here.** ASCENT's section 4 decided
which of two movie channels is calcium by comparing their temporal
variability. With one channel there is nothing to rank against, so that test
is inapplicable. Two substitutes were run instead:

- **Per-neuron variability (Verified, same statistic as ASCENT section 3):**
  on 40 sampled frames and 40 tracks, the detrended coefficient of variation
  of the 5 × 5 local mean has median 0.73 and 90th percentile 1.40. ASCENT's
  GCaMP channel gave a median of about 0.1 on the same statistic. The
  difference is mostly the zero background here: a 5 × 5 window around a dim
  neuron averages a few counts, so small changes are large fractions.
- **Stimulus-locked responses (Verified timing, Inferred biology):** the file
  stores three stimulus windows (NaCl 60.5–70.5 s, 2,3-pentanedione
  120.5–130.5 s, butanone 180.5–190.5 s; `intervals/chemical_stimuli`).
  Traces extracted at the ROI positions rise sharply at those onsets, and the
  stored traces agree (`previews/ey/neuron_traces.png`). The strongest
  pentanedione responder in both our traces and the stored ones is the ROI
  labelled AWAL. AWA is one of the main odour-sensing neuron pairs, so a sharp
  odour response in that ROI is biologically plausible and lends some
  credibility to the identity labels; we did not check the papers for which
  neuron each stimulus was chosen to target. The timing check is independent
  of biology: responses begin within a frame of the stored onsets. NaCl onset
  is followed by rises in ROIs labelled ASIL, ASJL, AFDL, I2L and VB1; the
  biology of those responses was not assessed.

**Whole-population brightness over time (Verified,
`previews/ey/per_frame_means.png`):** brightness at the ROI positions falls
from 5× its median at frame 0 to 2× by 1 s and to about 0.5× by the end, a
fast initial decay followed by a slow one. Frames 0–2 contain a handful of
voxels at 255 (the only frames that do). The stored traces show the same
shape. The cause of the fast initial decay (light-onset response of GCaMP,
photobleaching, or the rescaling) was not determined (Estimate: typical of
blue-light onset in GCaMP recordings).

**Agreement with the stored traces (Verified):** the per-neuron Pearson
correlation between our 5 × 5 × 1 local-mean trace and the stored
`SignalRawFluor` trace has median 0.72 (34 % of neurons above 0.8, 10th
percentile 0.25). The stored traces are unitless (0.03–2.5) and are described
only as "Raw fluorescence activity"; how they were extracted (window, z
extent, normalisation) is Unknown. The agreement is good enough to say the
stored traces come from these images at these positions, and poor enough
for dim neurons that the two should not be treated as interchangeable.

## 5. Voxel size and timing

| Quantity | Value | Status |
|---|---|---|
| xy pixel size | 0.2667 µm (`grid_spacing`, both imaging volumes) | **Documented (in file)**; inconsistent between the two stacks, see section 10 |
| z step | 1.5 µm | **Documented (in file)** |
| Field of view (calcium) | 256 × 0.267 = 68 µm × 34 µm if the value applies to the calcium stack; 137 × 68 µm if it applies to the NeuroPAL stack and the calcium stack is 2× binned | Estimate; the head of an adult is about 50 µm wide and 150 µm long, which fits the second reading |
| Depth span | 20 × 1.5 = 30 µm | Estimate from documented values |
| Volume rate | 4.0 Hz (`starting_time.rate`) | **Documented (in file)** |
| Duration | 962 / 4 = 240.5 s | Estimate from documented values |
| Per-frame timestamps | none; `starting_time` = 0 and a constant rate | **Verified absent** (same as ASCENT) |
| Stimulus schedule | NaCl 160 mM, 2-butanone 1e-4, 2,3-pentanedione 1e-4, 10 s each at 60.5, 120.5, 180.5 s, order randomised per worm | **Documented (in file)**; responses at the stated onsets **Verified** |
| Microscope | Nikon Ti-e spinning-disk confocal, 60× 1.2 NA water objective | **Documented (in file)** |
| Strain, stage, temperature | OH16230, young adult, 20 °C, sex "O" (DANDI: "Other"; hermaphrodite is an Estimate) | **Documented (in file)** |
| Session date | 2019-09-28; file written 2024-05-02 by the DANDI conversion | **Documented (in file)** |

Anisotropy: one z step is 5.6 xy pixels at 0.267 µm, or 2.8 if the calcium
pixel is really 0.533 µm.

## 6. Ground-truth annotations

**Verified** from the per-frame `PlaneSegmentation` tables
(`inspect_ey.py`, section 4):

| Property | Value |
|---|---|
| Layout | 961 groups `Seg_tpoint_0 … Seg_tpoint_960`, contiguous; each a table with `id`, `voxel_mask` (x, y, z, weight) and `voxel_mask_index` |
| Rows | 168,175 = 961 × 175 |
| Voxels per ROI | exactly one in every frame (so these are centre points, not masks); weight always 1.0 |
| Coordinates | unsigned integers; x 1–254, y 16–122, z 0–24 |
| Time coverage | frames 0–960 annotated, frame 961 (the last) not |
| ROIs per frame | 175 in every frame, identical `id` list 0–174 in every frame |
| Tracks | ROI index = track; all 175 tracks span all 961 annotated frames |
| Duplicate (roi, t) | 0 |
| Two ROIs on the same voxel in the same frame | 147 rows |
| Biological labels | `NeuronIDs/labels`: 175 names, all non-empty, all unique, 76 left/right pairs; same set as the NeuroPAL stack's labels (different order) |

**Out-of-stack coordinates (Verified):** 2,019 rows (1.2 %) have z = 21–24,
beyond the last slice (z = 20). They belong to four ROIs: FLPR in 925 of 961
frames, ADFL in 614, RIVR in 426, SMDVL in 54. The x and y ranges are always
inside the image. So for FLPR the annotated position is outside the imaged
volume for almost the whole recording. This is the mirror image of ASCENT's
edge-slice pile-up: there, out-of-stack neurons seem to have been clamped to
the first or last slice; here, the tracker was allowed to extrapolate past
the stack. WormND's loader silently clamps these to the last slice. Any
evaluation on this file must decide how to treat these rows.

**Provenance (Documented, Inferred):** the dNMF paper's method tracks neurons
by fitting a deformable model to the whole volume; the file holds two
position sets, `CalciumSeriesSegmentation` and `CalciumSeriesSegmentationdNMF`,
plus two trace sets with matching names. The two position sets agree exactly
in only 10 % of rows; mean absolute difference 1.1, 1.3 and 0.75 voxels in
x, y, z, maximum 32, 60 and 13. We do not know which set is the proofread
ground truth and which is the raw dNMF output, or whether either was
manually checked (Unknown; the DANDI record does not say). WormND uses the
first set as ground truth. The NeuroPAL identities were assigned by the
authors on the NeuroPAL stack (Documented in the NeuroPAL paper).

**Motion (Verified from the first set):** the median in-plane displacement
between consecutive frames is 0 voxels; 53 % of steps are exactly zero; the
99th percentile is 1.4 voxels and the maximum 16 (≈4 µm at 0.267 µm). Depth
changes by at most 2 slices per step and in 19 % of steps. From first to
last annotated frame, tracks move a median of 1.4 voxels and at most 61.3.
Like ASCENT, a nearly stationary preparation.

**Depth distribution (Verified):** annotations occupy z = 5–19 almost
entirely (z 0–4 together hold 871 rows, z 20 holds 214), matching the
depth-slice previews where slices 0–4 and 20 are nearly empty. There is no
edge-slice pile-up.

**Completeness (unresolved):** 175 neurons are tracked. A worm's head holds
roughly 180–190 neurons (general figure, not checked for this field); the previews
show bright cells without circles in the pharynx and in the posterior part
of the field (`mip_t480_gt.png`). Whether the unannotated cells are neurons
that could not be identified, non-neuronal cells, or deliberate omissions is
Unknown, as for ASCENT.

## 7. Are the annotations in bounds and aligned with the images?

**Axis identification (Verified, sections 5 and 6):** reading the voxel-mask
`x` as index 1 of the data (size 256), `y` as index 2 (size 128) and `z` as
index 3 (size 21) places 98.8 % of rows in bounds (the rest are the
out-of-stack z rows above); swapping x and y leaves 76,431 rows out of
bounds. Unlike ASCENT's square field, the bounds alone already settle the
axis order here. With image content: on five frames spread over the
recording, the mean 5 × 5 brightness at ROI points is 19–38× that at random
points in the same bounding box, and 91–99 % of ROI windows contain a
nonzero voxel versus 9–18 % of random windows; with x and y swapped the
ratio drops to 0.3–0.7. ASCENT's background-subtracted median ratio was not
usable here because the median voxel is 0, so mean ratios are reported
instead.

- 53–60 % of ROI points lie within 3 pixels of the brightest pixel of their
  11 × 11 in-plane neighbourhood (ASCENT: 70–94 %). The lower figure is
  consistent with half-resolution images, single-voxel positions, and dim
  neurons sitting next to bright ones; it was not investigated further.
- Visually (`previews/ey/*_gt.png`, `mip_movie_gcamp.mp4`), circles sit on
  cell bodies.

**Time alignment (not verifiable here):** the file has 961 annotation sets
for 962 frames and 961 trace rows starting at time 0, so the natural reading
is `Seg_tpoint_t` ↔ frame `t` with the last frame unannotated (Inferred). A
direct test (which of annotation sets t−1, t, t+1 best fits frame t, using a
3 × 3 window) found no preference at all: 372 / 266 / 321 frames and mean
differences under 0.5 % of brightness. Because the median motion is zero,
a ±1-frame offset is invisible in the images and also irrelevant for
tracking at this motion level; it would matter for a dataset with faster
motion, and the convention should be confirmed with the authors before
reuse.

**Track-identity correctness (not verified):** as for ASCENT, no test was run
to confirm that an ROI index stays on the same cell. Indirect evidence: zero
median motion, smooth stored traces, and biologically plausible
stimulus responses in the ROI labelled AWAL.

## 8. Previews

All under `previews/ey/` (git-ignored, on scratch). Generated by
`scripts/make_previews_ey.py` in 18 s with 297 MB peak memory, and
`scripts/plot_ey_traces.py`.

| File | Content |
|---|---|
| `depth_slices_t{0,480}_gcamp.png` (+`_gt`) | 5 × 5 grid of the 21 z-slices of one calcium frame, with/without ROI circles; the title states how many ROIs have z beyond the stack and are therefore not drawn |
| `mip_t{0,480}.png` (+`_gt`) | calcium max projection over z (top) and NeuroPAL RGB max projection (bottom, 2× resolution), with/without stored positions |
| `neuropal_rgb_mip.png`, `neuropal_rgb_mip_labels.png` | the NeuroPAL stack with its 175 neuron centres and, in the second, their names |
| `mip_movie_gcamp.mp4` | frames 200–349 (50–87.5 s at the in-file 4 Hz, covering the NaCl pulse) as 3× upscaled max projections, raw on the left, with ROI circles and names on the right, 10 fps, i.e. **2.5× faster than real time**; circles turn blue when the ROI's z is outside the stack; the stimulus name is burned in during the pulse (6.6 MB) |
| `per_frame_means.csv/.png` | per-frame mean, 99th percentile, max, zero fraction, mean at ROI positions, and the stored traces averaged over neurons, with the stimulus windows shaded |
| `per_frame_means_tracks.csv` | our 5 × 5 local-mean brightness of every track in every frame (long format, NaN where the ROI is outside the stack) |
| `stored_traces.csv` | the file's `SignalRawFluor` traces, 961 × 175 |
| `calcium_tracks.csv` | the per-frame positions as one long table: t, roi, x, y, z, weight, label |
| `neuron_traces.png` | the six neurons with the strongest stimulus-locked rise in the stored traces, stored trace versus ours |
| `inspection_summary.json` | machine-readable output of the inspection script |

Display intensities are fixed at 0 to the 99.9th percentile of frame 480 for
all frames (frame 0 is much brighter than the rest and would set too high a
ceiling), so brightness changes over time reflect the data.

## 9. Inventory entry (draft)

| Field | Value |
|---|---|
| Source | DANDI 000541 v0.241009.1457, Hobert lab (Columbia) / Yemini; NeuroPAL (Cell 2021) and dNMF (2020) papers; WormID-Bench "EY" |
| Original ID | `sub-20190928-13`, identifier `20190928_13`; asset `c4a52d3e-26a9-439c-8965-13320ab6cc43` |
| Format | NWB 2.6.0 (HDF5) with the `ndx-multichannel-volume` 0.1.12 extension; opens with pynwb 4.2 in 11 s |
| Dimensions | calcium T=962, C=1, X=256, Y=128, Z=21, uint16 with 0–255 values; NeuroPAL X=512, Y=256, Z=21, C=5, uint64 with 0–4095 values (Verified) |
| Voxel size | 0.267 × 0.267 × 1.5 µm stated for both stacks (Documented in file); the calcium stack is 2× coarser in xy than the NeuroPAL stack (Verified), so one of the two statements is wrong |
| Rate / duration | 4 Hz / 240.5 s (Documented in file / Estimate) |
| Channels | calcium: GCaMP only; NeuroPAL: mTagBFP2, CyOFP1, mNeptune 2.5, TagRFP-T, GFP (Documented in file) |
| Condition | immobilised in a microfluidic chip, three 10 s chemical pulses with timing in the file (Documented in file; responses Verified) |
| Movement | median per-step displacement 0 voxels, p99 1.4, max 16 (Verified) |
| Preprocessing | calcium series background-zeroed and 8-bit-scaled despite `format: raw` (Verified values, Inferred cause); further processing Unknown |
| Annotations | 175 tracked neurons × 961 frames = 168,175 single-voxel points, plus a second dNMF-labelled set that differs in 90 % of rows (Verified); 1.2 % of rows have z beyond the stack |
| Annotation provenance | dNMF tracking (Documented); extent of manual correction Unknown |
| Biological identities | all 175 ROIs named, 76 L/R pairs, from NeuroPAL (Verified presence, Documented method) |
| Stored traces | `SignalRawFluor` and `SignalRawFluordNMF`, 961 × 175, unitless (Verified) |
| License | CC-BY-4.0 (Documented) |
| Verification status | file downloaded, SHA-256 matches, structure, axis order and spatial alignment verified; time-index convention, track-identity correctness, annotation completeness and the ground-truth status of the two position sets not verified |
| Other assets | 20 more recordings with the same layout, 1.42–1.46 GB each, 937–963 frames, 166–188 ROIs (Verified remotely, headers only) |

## 10. Open questions and uncertainties

1. **Pixel values.** The series is 0–255 with a zero median inside a uint16
   array labelled raw. What was subtracted, whether anything was clipped
   (frames 0–2 touch 255), and whether the scaling is per frame or global is
   Unknown. For a benchmark, this file provides *processed* images; the raw
   camera frames are not in the archive.
2. **Voxel size.** Both stacks claim 0.267 µm pixels, but the calcium stack is
   exactly half the NeuroPAL resolution (slopes 1.99). Either the calcium
   pixel is 0.533 µm (2 × 2 binning) or the NeuroPAL pixel is 0.133 µm. The
   NeuroPAL paper's 0.267 µm figure would make the calcium pixel 0.533 µm
   (Estimate); not verified.
3. **Which position set is ground truth.** Two sets (`…Segmentation` and
   `…SegmentationdNMF`) differ in 90 % of rows by about one voxel and
   occasionally by tens of voxels. WormND uses the first. Whether it was
   proofread by hand is Unknown.
4. **Out-of-stack z values.** Four tracks (FLPR, ADFL, RIVR, SMDVL) spend
   54–925 frames at z = 21–24. They cannot be seen in the images there. A
   benchmark must either drop those rows or state that clamping is the rule.
5. **Time-index convention.** `Seg_tpoint_t` ↔ frame `t` is assumed; a
   ±1-frame offset is undetectable here (section 7).
6. **Annotation completeness.** Visible unannotated cells exist, as in
   ASCENT; no audit was done.
7. **Track-identity correctness** was not tested directly.
8. **Initial brightness decay.** A 5× drop in the first second; cause not
   determined. Frame 0 is atypical and should not be used alone to set
   display or normalisation parameters.
9. **Stored-trace extraction method** is Unknown; our local means agree with
   a median correlation of 0.72.
10. **Timestamps.** None per frame; the 4 Hz rate is a single number in the
    file. Whether slices are acquired sequentially within a volume (likely
    for a spinning disk with a piezo) and the within-volume timing are
    Unknown.
11. **Full text of the two papers** was not consulted for this note; the
    Documented items come from the NWB fields, the DANDI record and the
    WormND repository.

## 11. Comparison with ASCENT

| | ASCENT (Zenodo 17561700) | EY sub-20190928-13 (DANDI 000541) |
|---|---|---|
| Container | plain HDF5, `t{n}/c{m}` datasets, no metadata | NWB 2.6 with extension; voxel size, rate, channels, stimuli, strain, identities inside |
| Volume | 2 ch × 9 × 512 × 512, 0.243 × 0.243 × 1.5 µm, 3.3 Hz, 1100 frames, 333 s | 1 ch × 21 × 128 × 256, 0.267 (or 0.533) × … × 1.5 µm, 4 Hz, 962 frames, 240 s |
| Depth coverage | 12 µm, thinner than the head; annotations pile up on edge slices | 30 µm, covers the head; annotations in the middle slices, four tracks extrapolated beyond the stack |
| Reference channel | yes, mNeptune in every frame | none in the movie; a static 5-colour NeuroPAL stack instead |
| Pixel values | uint16 camera counts, baseline ~320, gzip | uint16 holding 0–255, zero background, uncompressed |
| Annotations | 96 tracks × 1100 frames, ZephIR + manual proofreading, no identities | 175 tracks × 961 frames, dNMF, proofreading status unknown, all with biological names |
| Motion | median 1 voxel per step, max 40 | median 0, max 16 |
| Stimulus | documented schedule, timing inferred from responses | timing stored in the file, responses verified at the stated onsets |
| Spatial alignment | verified (contrast 24–83×) | verified (contrast 19–38×, nonzero-window test) |
| Per-frame access | one gzip dataset per frame-channel | one contiguous slab per frame; fastest of the two |
| Open issues shared | timestamps, completeness, track-identity correctness | same, plus pixel-value provenance and voxel-size inconsistency |

For a tracking benchmark both are low-motion, immobilised recordings. The
EY file adds biological identities and a larger, fully-covered head, but its
images are already processed and single-channel; ASCENT keeps raw counts and
a reference channel but covers only a slab of the head.
