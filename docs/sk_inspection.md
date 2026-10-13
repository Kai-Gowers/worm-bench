# SK1 and SK2 recordings (DANDI 000565 and 000472): inspection notes

Two recordings from the Kato lab (FOCO lab, UCSF), one from each of the two
dandisets that WormID-Bench calls SK1 and SK2:

| | SK1 | SK2 |
|---|---|---|
| Dandiset | 000565 "C. elegans whole-brain neuroPAL and immobilized calcium imaging", v0.241009.1504 | 000472 "NeuroPAL volumetric images", v0.241009.1502 |
| Recording | `sub-20230322-20-16-50` | `sub-2022-04-26-w01-NP1` |
| Asset | `36715c41-aa5c-4a05-affd-66a9a6089f7f`, 4,293,448,011 B | `8525ab3c-4a5d-4893-a2c8-4ac7ccbd87bd`, 8,412,993,256 B |
| SHA-256 (published, verified after download) | `f442901c12bb3430cbbe09969e1cb0b6775e1bd33b68a0d668389ab0a78772d1` | `66e4243260b2de7d545711dad1d51e57d2c81f00c3ab85436387fc28895d99ab` |
| Local copy | `data/sk1/` | `data/sk2/` |
| Previews | `previews/sk1/` | `previews/sk2/` |

Both dandisets are CC-BY-4.0, published 2024-10-09, contributor "Kato,
Saul", no related publication listed. Both files were written on 2024-02-26
by the same NWB export (same layout, same extension, same device string),
so they are treated here as one source with two configurations.

Every statement is tagged **Verified** (measured in the file by
`scripts/inspect_sk.py`, logs and `inspection_summary.json` under
`previews/sk*/`), **Documented** (stated in the file, the DANDI metadata, or
an external source named in the text), **Inferred** (follows from verified
facts), **Estimate**, or **Unknown**.

## 0. Why these assets

**Most assets in both dandisets have no calcium movie.** The remote survey
(`scripts/survey_dandi_sk.py`, headers only) found (Verified):

| | assets | with `CalciumImageSeries` | NeuroPAL-only |
|---|---|---|---|
| 000565 (SK1) | 21, 46.3 GB | 4 (4.29, 4.29, 13.79, 13.79 GB) | 17 (0.25-0.92 GB each; strains FC121, FC128, FC111) |
| 000472 (SK2) | 10, 18.4 GB | 2 (8.34, 8.41 GB) | 8 (0.20-0.22 GB; all OH16230) |

This matches the WormID-Bench counts "SK1: 4, SK2: 2" exactly, and all six
appear in WormND's `dataset_split.csv` with `use_for_calcium_task = 1`.
The 000472 description talks only about NeuroPAL images; the two calcium
movies in it are not mentioned anywhere in the DANDI metadata.

The four SK1 calcium assets fall into two configurations (Verified):

| asset | calcium (T, X, Y, Z) | pixel, z step | ROIs / named | size | WormND fold |
|---|---|---|---|---|---|
| sub-20230322-18-57-04 | 1500, 1000, 240, 12 | 0.3208 µm, 2.5 µm | 118 / 25 | 4.29 GB | 0 |
| **sub-20230322-20-16-50** | 1500, 1000, 240, 12 | 0.3208 µm, 2.5 µm | 117 / 31 | 4.29 GB | 2 |
| sub-20230506-15-01-45 | 1500, 2000, 480, 12 | 0.1604 µm, 3.0 µm | 93 / 23 | 13.79 GB | 2 |
| sub-20230506-15-33-51 | 1500, 2000, 480, 12 | 0.1604 µm, 3.0 µm | 103 / 29 | 13.79 GB | 4 |

The two SK2 calcium assets share one configuration (3000 × 1000 × 240 × 12,
0.3208 µm, 1.5 µm z step, 3.3 Hz): `w00` has 121 ROIs (21 named, 6-channel
NeuroPAL stack, fold 3) and **`w01`** 140 ROIs (32 named, 9-channel
NeuroPAL stack, fold 2).

Chosen: the SK1 asset with the most named neurons in the common 0.32 µm
configuration, and the SK2 asset with more ROIs and names. Both sit in
WormND fold 2, like the EY and KK picks. The high-resolution SK1
configuration (2000 × 480, 13.8 GB) was not downloaded; it is listed above
so it is not forgotten.

## 1. Background for newcomers

**What was imaged.** One young-adult hermaphrodite per file (sex "O" in
the file; hermaphrodite is an Estimate), strain OH16230, immobilized in a
microfluidic chip. Documented externally (CGC strain record for OH16230):
the genotype is `otIs670 V; otIs672`, where otIs670 is the NeuroPAL
transgene and otIs672 drives nuclear-localised GCaMP6s in all neurons
(`rab-3::NLS::GCaMP6s` plus `arrd-4::NLS::GCaMP6s`). So the calcium signal
comes from nuclei, which look like compact blobs about 2 µm across, and the
same strain was used for the EY recording (DANDI 000541), whose GCaMP
variant the EY file leaves unstated.

**Immobilized.** The worm is held in a chip, so unlike the KK worm (which
crawled in place) the head barely moves. The annotations confirm this
(section 6): the median frame-to-frame step is zero. What does move is slow:
a drift of a few pixels and about one z slice over the whole movie, plus
photobleaching, which dims the GCaMP signal over minutes.

**One channel.** Only the GCaMP channel was recorded during the movie.
There is no stable reference channel (unlike ASCENT and KK), so a tracker
has to follow nuclei whose brightness changes with activity, and quiet
neurons can fade into the background.

**Camera offset.** These cameras add a constant of roughly 93-98 gray
counts to every pixel. Background therefore sits near 95, not near zero,
and the number that matters for a neuron is "counts above 95". SK1's
neurons are only 5-15 counts above the offset for most of the movie;
SK2's are 30-75 counts above it.

**NeuroPAL.** Each file also holds a separate multicolour structural
stack of the same head, taken before or after the movie, and a list of
neuron centres in that stack with biological names for a minority of them.
The names were then transferred to the calcium ROIs (how is not stated;
section 6).

**Stored DF/F.** Instead of raw traces the files store ΔF/F0, i.e.
(fluorescence − baseline) / baseline per neuron. A value of 1 means "twice
the baseline". How the baseline was chosen is not stated.

## 2. Files

One NWB 2.6.0 file per recording, HDF5 with the `ndx-multichannel-volume`
extension embedded (Verified). Both open with h5py and with pynwb 4.2
(13 s and 25 s, mostly spent building the 1500 / 3000 PlaneSegmentation
objects). Nothing else was downloaded.

| | SK1 | SK2 |
|---|---|---|
| Calcium raw size | 8.64 GB (5.76 MB/frame) | 17.28 GB (5.76 MB/frame) |
| File size | 4.29 GB | 8.41 GB |
| gzip ratio | 2.0 | 2.1 |
| NeuroPAL raw stack | 161 MB (7 channels) | 194 MB (9 channels) |
| Streaming pass over all frames | 205 s (0.14 s/frame) | 408 s (0.14 s/frame) |

Chunks are (2, 125, 30, 3) in (T, X, Y, Z), so one frame needs 8 × 8 × 4 =
256 chunk reads and always decompresses its neighbour frame too; reading
is nevertheless fast because the chunks are small.

## 3. Image structure

**Calcium series** `acquisition/CalciumImageSeries/data` (Verified):

| | SK1 | SK2 |
|---|---|---|
| shape (T, X, Y, Z) | 1500, 1000, 240, 12 | 3000, 1000, 240, 12 |
| dtype | uint16 | uint16 |
| rate | 1.04 Hz → 1442 s (24.0 min) | 3.3 Hz → 909 s (15.2 min) |
| `scan_line_rate` attr | 2995 | 9504 |
| grid spacing (file) | 0.3208, 0.3208, 2.5 µm | 0.3208, 0.3208, 1.5 µm |
| field of view | 321 × 77 × 27.5 µm | 321 × 77 × 16.5 µm |
| `format` / unit | "raw" / "Voxel gray counts" | same |
| description | "Raw GCaMP series images"; comment "single channel GFP-GCaMP representing GCaMP signal" | same |

The `dimension` dataset reads [240, 1000, 12] while the data are stored
(1000, 240, 12): the `dimension` field lists Y, X, Z. Coordinates in the
segmentation follow the data order (x in 0..999, y in 0..239; section 7),
so `dimension` is just inconsistent and nothing else needs swapping.

**Pixel values** (Verified, all frames streamed):

| | SK1 | SK2 |
|---|---|---|
| minimum | 24-28 | 26-33 |
| median (background) | 93-94 | 98-99 |
| 99th percentile | 121-123 | 136-148 |
| per-frame maximum | 445-1523 | 704-4095 |
| global maximum | 1523 | 4095 |
| voxels equal to 0 | none | none |
| voxels at 4095 | none | first frames only (≈6 × 10⁻⁶ of voxels at t=0) |

So the data are genuine camera counts with a ~95-count offset, 12-bit
range (SK2 clips at 4095 at the start; SK1 never gets close), and no
background subtraction or rescaling. This is the first of the four
sources whose images were not processed to 8 bits before archiving (EY and
KK both hold 0-255 values). Consecutive SK1 frames correlate at only 0.19
pixel-for-pixel: the images are shot-noise dominated and the neurons are
faint (section 4).

**NeuroPAL stack** `acquisition/NeuroPALImageRaw/data` (X, Y, Z, C) uint16,
pixel 0.3208 µm and z step 0.75 µm in both (Documented):

| | SK1 | SK2 |
|---|---|---|
| shape | 1000, 240, 48, 7 | 1000, 240, 45, 9 |
| channel order (file) | mTagBFP2, CyOFP1, CyOFP1-high filter, GFP-GCaMP, mNeptune 2.5, Tag RFP-T, mNeptune 2.5-far red | mNeptune 2.5, Tag RFP-T, CyOFP1, GFP-GCaMP, mTagBFP2, mNeptune 2.5-far red (6 names for 9 planes) |
| `RGBW_channels` | [4, 1, 0, 6] → R = mNeptune 2.5, G = CyOFP1, B = mTagBFP2, W = mNeptune far-red | [0, 2, 4, 1] → R = mNeptune 2.5, G = CyOFP1, B = mTagBFP2, W = Tag RFP-T |

The RGB mapping follows the NeuroPAL convention in both (unlike KK). The
SK2 stack has 9 channel planes but only 6 channel names; which three planes
are unnamed is Unknown. A second, processed NeuroPAL stack
(`processing/ProcessedImage`, "median filtering followed by color histogram
matching to reference NeuroPAL images") is also stored: (1000, 240, 45, 4)
in SK2 but (240, 5, 48, 1000) in SK1, an axis order and channel count
(5 planes, 4 names) that match nothing else in the file (Unknown; not used
here).

**Not present**: no `intervals` group, no stimulus table, no second
calcium channel, no raw-fluorescence traces, no label volumes.

## 4. Channels

One calcium channel: GFP-GCaMP, excitation 488 nm, emission filter Chroma
ET 525/50 (Documented). Microscope: Leica DMi8 with a Yokogawa CSU-W1 SoRA
spinning disk, 40× water-immersion 1.1 NA (Documented).

ASCENT's calcium-versus-reference test has nothing to compare against.
Per-neuron temporal variability on 40 sampled frames and 40 full-length
tracks (Verified): detrended CV median 0.075 (p90 0.28) in SK1 and 0.145
(p90 0.49) in SK2; median neuron brightness 100 counts (SK1, i.e. 7 above
background) and 151 counts (SK2, 53 above background). The whole-volume
mean is flat (temporal CV 0.003-0.005).

**Activity is clearly present** (Verified from the traces,
`previews/sk*/neuron_traces.png`): both recordings show large, slow
calcium events (several-fold rises lasting 20-60 s) in a handful of
neurons, in SK1 timed irregularly at roughly 90, 270, 480, 650, 800 and
1100 s and in SK2 at roughly 20, 110, 150, 200, 280, 430, 480, 740 s. In
SK1 the per-frame maximum voxel jumps from ~500 to ~1000 counts during
these events, and the maximum pixel lies inside a neuron each time (5-6 px
from unlabelled ROI 109, or 3 px from ROI 81), so these are neuron
responses, not stimulus light leaking into the camera. The named neurons
among the most active are AVAL, AVDL and SMDDL in SK1 and AVAR, SMDVR in
SK2; AVA is the backward-locomotion command neuron, so this pattern looks
like spontaneous reversal-like activity in an immobilized worm (Estimate).

**Stimulus: contradictory metadata** (Documented). Each file's
`session_description` is "NeuroPAL and calcium imaging of immobilized worm
with optogenetic stimulus", and SK1's `experiment_description` is
"optogenetic stimulation of single neurons", while the 000565 dandiset
description says the calcium data are "unstimulated" and no stimulus timing
is stored in either file. Whether a stimulus was delivered during these
movies is Unknown (author question 1).

## 5. Voxel size and timing

Documented in the file: 0.3208 µm pixels, 2.5 µm (SK1) / 1.5 µm (SK2) z
steps, 1.04 Hz / 3.3 Hz. Checks (Verified):

- x and y: matching named neurons between the NeuroPAL stack and calcium
  frame 0 gives slopes 0.99 / 1.00 (x) and 0.83 / 0.98 (y) with r ≥ 0.956,
  so the two stacks share the stated pixel size. SK1's lower y slope comes
  from two of its 13 matched neurons (FLPL, AVDL) sitting 12-20 px away
  from their NeuroPAL counterpart; with 13 points this is not conclusive.
- z: the same matching gives 3.81 NeuroPAL slices per calcium slice in SK2
  (r 0.98, 23 neurons) where the stated spacings predict 1.5 / 0.75 = 2.0,
  and −2.78 in SK1 (r −0.59, 13 neurons) where they predict +3.33. Taken at
  face value the SK2 calcium z step is about 2.9 µm rather than 1.5 µm, and
  the SK1 calcium stack runs in the opposite z direction to its NeuroPAL
  stack. Both are Inferred from small samples and are author question 4.
- Timing: `starting_time` is 0 with a fixed rate; no per-frame timestamps.
  Nothing in the file verifies the rate. SK2's `session_start_time` is
  midnight (time of day not recorded).

## 6. Ground-truth annotations

**Format** (Verified). `processing/CalciumActivity/CalciumSeriesSegmentation`
holds one `PlaneSegmentation` per frame, `Seg_tpoint_0` … `Seg_tpoint_{T-1}`,
contiguous, one per movie frame (every frame is annotated, unlike EY).
Each holds one `voxel_mask` row (x, y, z, weight = 1) per ROI and the same
`id` list 0..N−1 in every frame, so the ROI index is the track id, exactly
as in EY. `NeuronIDs/ImageSegmentation` contains the same `Seg_tpoint_*`
groups as HDF5 hard links to the same objects (checked on 50 groups and on
sampled voxel masks), so it is not a second annotation set. Centres are
single voxels: there are no masks and no shape information.

| | SK1 | SK2 |
|---|---|---|
| ROIs (tracks) | 117, all present in all 1500 frames | 140, all present in all 3000 frames |
| rows | 175,500 | 420,000 |
| x, y, z ranges | 315-885, 56-210, 2-11 | 33-919, 17-203, 1-10 |
| rows out of bounds | 0 | 0 |
| duplicate (roi, t) / coincident points | 0 / 0 | 0 / 0 |
| named ROIs | 31 (27 %), all unique, 2 L/R pairs | 32 (23 %), all unique, 4 L/R pairs, 9 names end in "?" |

SK1 names: ADAL ADEL ADLL AINL AIYR ALA AQR ASEL ASGL ASKL AVAL AVDL AVEL
FLPL IL1L M1 OLLL OLQDL OLQR RIAL RMDDR RMDL RMDVL RMED RMEL RMER RMEV SMDDL
SMDVL URXL URXR. SK2 names: ADER AINR? ALA AQR? AVAL AVAR AVEL? AVHR AVJL?
CEPDL FLPR? I2L? I3 M1 M4 M5? MI NSML NSMR OLQDL OLQDR RID RIVR RMDR RMDVL
RMED RMEL RMEV? RMGR SMBDL? SMDVL SMDVR. The "?" suffix is an explicit
uncertainty marker in SK2; SK1 has none (whether its names are all
confident or the marker was simply not used is Unknown, question 5).

**Motion in the annotations** (Verified):

| | SK1 | SK2 |
|---|---|---|
| in-plane step, median / p99 / max | 0 / 1.0 / 1.4 voxels | 0 / 0 / 1.4 voxels |
| steps with zero in-plane motion | 99 % | 99 % |
| z changes | never | never |
| first → last frame displacement, median / max | 3.2 / 12.2 voxels (1.0 / 3.9 µm) | 2.2 / 10.2 voxels |

So the "tracks" are near-static points that creep by single-voxel steps in
x and y and never move in z. This is consistent with an immobilized worm,
but it also means the per-frame positions carry almost no tracking
information beyond the frame-0 detection: a tracker that simply copies
frame 0 would be within one or two voxels of the stored answer for most of
the movie.

**Do the points follow the neurons?** (Verified, section 6c of the script.)
On 12 frames spread over each movie, the global image shift relative to
frame 0 was measured by cross-correlating smoothed max projections, and
compared with the median shift of the annotations:

| | image shift by the end (dx, dy, dz) | annotation shift by the end | ROI points within 3 px of a smoothed local maximum, first → last |
|---|---|---|---|
| SK1 | (−4, −3, −2) px/slices | (−3, +1) px | 62 % → 8 % |
| SK2 | (+6, 0, +1) | (0, +1) | 65 % → 25 % (minimum 9 %) |

In SK2 the image drifts 6 px in x and the annotations do not follow, so by
the second half of the movie most stored points sit 4-6 px (1.3-1.9 µm,
about one nuclear radius) off the nucleus they name, although the nuclei
themselves stay 30-75 counts above background. In SK1 the annotations
follow the drift to within 1-2 px; the collapse of the local-maximum hit
rate there is explained by the signal, not the positions: the smoothed
brightness above background at the stored points falls from a median of
15 counts (p90 99) at t=0 to 4-6 counts (p90 25-30) after t≈250 s and
stays there, and the stored z slice remains the brightest slice at every
tested frame (a ±1 slice shift always scores lower). After the first
~4 minutes most SK1 nuclei are within a few counts of the noise floor at
the stored positions.

**Time index.** The point-based time-index test (does `Seg_tpoint_t` fit
image t better than t±1?) cannot discriminate: there is no frame with
population motion above 0 voxels. The naming convention (frame 0 =
`Seg_tpoint_0`) is therefore assumed, as for EY, and the trace comparison
below peaks at lag 0 (weak support).

**NeuroPAL segmentation** (Verified): 123 (SK1) / 180 (SK2) centres in the
NeuroPAL stack, 44 / 40 of them named. Only 13 / 23 names are shared with
the calcium ROIs; 18 / 9 calcium names have no NeuroPAL counterpart. The
calcium names were therefore not simply copied from the stored NeuroPAL
segmentation (Inferred; question 6).

**Provenance.** Nothing in either file or dandiset says how the calcium
positions were produced (detection on one frame and propagation, tracking
software, manual checking). Unknown (question 2).

## 7. Are the annotations in bounds and aligned with the images?

Bounds (Verified): every row lies inside the (X, Y, Z) = (1000, 240, 12)
stack as labelled, and none do with x and y swapped (175,500 / 390,000
rows out of bounds), so the axis assignment is settled by bounds alone.

Alignment with image content, five frames per movie (Verified):

| | SK1 | SK2 |
|---|---|---|
| mean 5×5 brightness at ROI points / at random points (raw counts, offset included) | 1.3, 1.1, 1.1, 1.1, 1.1 | 2.9, 1.9, 1.8, 1.6, 1.4 |
| same with x↔y swapped | undefined (all swapped points out of bounds) | 0.96-1.0 |
| ROI points within 3 px of the raw local maximum | 0.63, 0.21, 0.22, 0.24, 0.24 | 0.63, 0.33, 0.32, 0.36, 0.37 |
| 5×5 mean minus background at the stored points, frame 0 → last (median / p90) | 15 / 102 → 4 / 27 counts | 73 / 473 → 30 / 115 counts (smoothed) |

The raw ratios are dominated by the camera offset and so understate the
contrast; the script's generic "much brighter than random" verdict (ratio
> 3) is False for both, but the background-subtracted numbers and the
previews (`mip_t0_gt.png`) show the frame-0 points sitting on nuclei in
both recordings. Frame 0 is the only frame where the annotations are
unambiguously well aligned in SK1; see section 6 for how that degrades.

**Stored DF/F versus the images** (Verified, every frame; 5×5×1 mean at
the stored position, offset not subtracted):

| | SK1 | SK2 |
|---|---|---|
| per-neuron correlation, median / p10 / min | 0.06 / −0.20 / −0.54 | 0.76 / 0.18 / −0.69 |
| neurons with r > 0.8 | 2 % | 41 % |
| SK1 after 10-frame smoothing and removing the slow trend | 0.26 (p90 0.60) | |
| best lag | 0 frames | |
| stored value range, median of per-neuron medians | −0.15 to 6.6, 0.31 | −0.36 to 11.6, 0.79 |

In SK2 the stored traces are a recognisable function of the images at the
stored positions (the most active neurons match event for event, r 0.6-0.9;
the implied baseline from a per-neuron linear fit is 1.3 × the intercept,
not 1.0, so the F0 convention is not the plain "F/F0 − 1" of the stored
positions). In SK1 they are not: only the largest events (AVAL, r 0.52)
line up, the stored traces do not show the bright first 100 s that the
images show at AVDL and SMDDL, and the per-neuron agreement is slightly
anti-correlated with brightness. The SK1 traces were evidently computed
from something other than the stored points on the stored images (larger
or different ROIs, a bleach-corrected or registered image, or a different
baseline); question 3.

## 8. Previews

Under `previews/sk1/` and `previews/sk2/` (not in git), 64-65 MB each:

- `mip_t{0,750|1500}.png` and `_gt`: GCaMP max projection over z with the
  NeuroPAL RGB projection below, with/without ROI points. The display
  range is fixed from one frame at the 1st-99.9th percentile, so neurons
  look saturated white and the ~95-count background mid-grey.
- `depth_slices_t*_gcamp[_gt].png`: the 12 z slices.
- `neuropal_rgb_mip[_labels].png`: NeuroPAL stack with its centres and names.
- `mip_movie_gcamp.mp4`: 200 frames (0-192 s SK1, 0-61 s SK2) at 10 fps,
  raw above, with ROI points and names below (yellow = named, red = unnamed).
- `per_frame_means.png`, `neuron_traces.png` (from `plot_sk_traces.py`),
  `calcium_tracks.csv` (long table of all positions), `per_frame_means*.csv`,
  `stored_traces.csv`, `inspection_summary.json`.

Visual checks (Verified by eye): in both, frame-0 points sit on nuclei;
the head lies along x with the anterior to the left (SK1) or right (SK2),
and SK2 includes the pharynx and a stretch of ventral cord along the top
edge whose cells are also annotated. In SK2 at t=1500 many cells visible
in the mid-body region carry no annotation, i.e. the 140 ROIs do not cover
every visible nucleus. SK1 nuclei are faint against noise except during
events.

**Bleaching and the first minute** (Verified from `per_frame_means.png`):
SK1's population brightness at the ROI points rises to a peak at ~50 s,
falls by two thirds between 80 and 250 s and is then flat; its first-frame
value is 1.19 × the median and the last 0.97 ×. SK2 jumps up during the
first 10 s, then decays steadily to a third of its peak by 15 minutes
(first / median 1.73, last / median 0.85); the frame maximum falls from 4095
to ~700 over the movie. A tracker evaluated on these movies faces a
non-stationary brightness.

## 9. Inventory entry (draft)

| Field | SK1 sub-20230322-20-16-50 | SK2 sub-2022-04-26-w01-NP1 | Status |
|---|---|---|---|
| Source | DANDI 000565 v0.241009.1504, Kato lab (FOCO), UCSF; experimenter Raymond Dunn | DANDI 000472 v0.241009.1502, same lab; experimenters Jackson Borchardt, Greg Bubnis | Documented |
| Publication | none listed ("related_publications" empty) | same | Documented |
| License | CC-BY-4.0 | CC-BY-4.0 | Documented |
| Files | one NWB, 4.29 GB, SHA-256 above | one NWB, 8.41 GB, SHA-256 above | Verified |
| Images | (T, X, Y, Z) = (1500, 1000, 240, 12) uint16, one GCaMP channel, gzip | (3000, 1000, 240, 12), same | Verified |
| Pixel values | camera counts, offset ~93, max 1523, no processing evident | counts, offset ~98, 12-bit, clipped at 4095 in the first frames | Verified |
| Voxel size | 0.3208 × 0.3208 × 2.5 µm (file); z step questioned | 0.3208 × 0.3208 × 1.5 µm (file); neuron matching suggests ~2.9 µm | Documented / Inferred |
| Rate, duration | 1.04 Hz, 1442 s | 3.3 Hz, 909 s | Documented |
| Condition | immobilized in microfluidic chip, 20 °C, young adult, OH16230 (NeuroPAL + pan-neuronal nuclear GCaMP6s) | same | Documented (strain genotype from CGC) |
| Stimulus | metadata contradictory; none stored | same | Unknown |
| Motion | annotations: median step 0, drift ≤ 12 voxels; images drift ≤ 5 px, ≤ 2 slices | annotations ≤ 10 voxels; images drift 6 px, 1 slice | Verified |
| Annotations | 117 single-voxel centres in all 1500 frames; z constant; track id = ROI index | 140 centres in all 3000 frames; z constant | Verified |
| Identities | 31 named (27 %) | 32 named (23 %), 9 marked "?" | Verified |
| Annotation provenance | not stated; near-static points; alignment degrades after ~250 s (signal) | not stated; alignment degrades after ~250 s (6 px image drift not followed) | Unknown / Verified |
| Traces | DF/F (1500 × 117), agrees poorly with the images | DF/F (3000 × 140), agrees well | Verified |
| NeuroPAL | 7-channel stack, 123 centres, 44 named | 9-channel stack (6 names), 180 centres, 40 named | Verified |
| WormND | fold 2, calcium task | fold 2, calcium task | Documented |
| Verification status | downloaded, checksum OK, all frames streamed | same | Verified |

## 10. Open questions and uncertainties

Phrased for sending in `docs/author_questions.md`. In brief:

1. Was an optogenetic stimulus delivered, and when (metadata contradict
   each other)?
2. How were the per-frame positions produced, and were they checked?
3. How were the stored DF/F traces computed, given the poor agreement in
   SK1 and the non-standard baseline in SK2?
4. Is the calcium z step correct, and are the NeuroPAL and calcium stacks
   acquired in the same z direction?
5. What does the "?" suffix mean, and why is it absent in SK1?
6. Why do most calcium names lack a NeuroPAL counterpart in the same file?
7. Camera bit depth, exposure and the volume rate.

Other open points (not for the authors): which three of the nine SK2
NeuroPAL planes are unnamed; the odd shape of SK1's processed NeuroPAL
stack; whether the 2000 × 480 SK1 configuration behaves like this one.

## 11. Comparison with ASCENT, EY and KK

| | ASCENT | EY | KK | SK1 | SK2 |
|---|---|---|---|---|---|
| Worm | freely moving (tracked stage) | immobilized | semi-restrained, moves | immobilized | immobilized |
| Channels in movie | 2 (GCaMP + reference) | 1 | 2 | 1 | 1 |
| Pixel values | 16-bit counts | 8-bit, background zeroed | 8-bit, clipped | 12-bit counts, offset ~93 | 12-bit counts, offset ~98 |
| Frames / rate | 1 Hz-class, minutes | 962 @ 4 Hz | 1646 @ 1.67 Hz | 1500 @ 1.04 Hz | 3000 @ 3.3 Hz |
| Annotation form | centres | centres | label volumes | centres | centres |
| Annotated frames | subset | 961 of 962 | all | all | all |
| Annotation motion | large | ~0 | large | ~0, z fixed | ~0, z fixed |
| Named neurons | some | most | 143 of 156 | 31 of 117 | 32 of 140 |
| Stored traces | raw | raw | raw (mask mean) | DF/F | DF/F |
| Traces reproducible from images + points | yes | partly | yes (r = 1.00) | no (r ≈ 0.06) | mostly (r ≈ 0.76) |
| Stimulus | none | chemical, timed | none | contradictory | contradictory |
| Alignment at frame 0 | good | good | good | good | good |
| Alignment late in movie | good | not testable | good (24/24) | points fine, signal gone | points lag a 6 px drift |

What SK adds to the collection: unprocessed camera data from a spinning
disk, a third lab and microscope, nuclear GCaMP6s in the same strain as
EY, the longest immobilized movies so far, and two configurations of pixel
size. What it lacks: a reference channel, neuron-shaped annotations,
tracking information beyond frame 0, and any documentation of the
annotation pipeline. For a tracking benchmark the honest description is
"detections on frame 0 plus a near-static propagation", and the SK2 drift
shows that even that propagation is not reliable to the pixel after a few
minutes.
