# KK recording sub-231013-02 (DANDI 000692): inspection notes

Dataset: *Whole-brain spontaneous GCaMP activity with NeuroPAL cell ID
information of semi-restricted worms* (Suzuki, Wen, Sprague, Onami, Kimura;
Kimura lab, Nagoya City University; DANDI 000692, version 0.240402.2118,
published 2024-04-02, license CC-BY-4.0, DOI 10.48324/dandi.000692/0.240402.2118).
The related publication is given as "unpublished" in every file. This is the
"KK" dataset of WormID-Bench / WormND (9 recordings used for the tracking
task). Benchmark code: https://github.com/focolab/WormND

Companion to `docs/ascent_inspection.md` and `docs/ey_inspection.md`; the
section numbering and the evidence labels are the same so the three can be
read side by side:

- **Verified**: measured directly from the downloaded file by `scripts/inspect_kk.py`.
- **Documented**: stated in the NWB metadata fields, the DANDI record, or the WormND repository, not checkable from the images.
- **Inferred**: a conclusion from measurements plus an assumption, stated separately from the measurement.
- **Estimate**: our interpretation of the evidence.
- **Unknown**: not stated anywhere we could find.

As with EY, the NWB container stores voxel size, frame rate, channel names,
strain, microscope and biological identities *inside the file*; these are
"Documented (in file)". Unlike EY, this file holds no stimulus intervals (the
session is "without Stimulation") and its annotations are label volumes, not
point tables.

## 0. Why this asset

All 9 assets were listed through the DANDI API and opened remotely
(HDF5 headers and small datasets only, about 1.5 s per file;
`scripts/survey_dandi_000692.py`). They share one layout: a two-channel
image series (GCaMP and TagRFP-T) of 1,645–1,646 volumes of 288 × 576 × 25
voxels at 1.67 Hz, a label-volume series of the same shape holding the
tracking result, one point table for a single time point, 147–159 ROIs with
127–145 real neuron names (the rest "unknown"), a 4-channel NeuroPAL stack of
200–243 × 461 × 28 voxels with 149–163 named neurons, and stored traces.
Sizes are 8.58–9.37 GB each, 80.5 GB in total. All are strain KDK92, young
adult, 20 °C, recorded on three days (2023-09-18, 09-28, 10-13), all without
stimulation. WormND uses all 9 for its tracking task.

`sub-231013-02` was chosen because it has the modal frame count (1,646), a
near-median ROI count (156, median 155), the second-highest number of unique
real names (144 of 156; the file itself gives 143 real names plus 13
"unknown"), and sits in WormND's fold 2 like the EY recording we inspected
(`neuron_identification/fDNC_for_WormND/wormnd/dataset_split.csv`,
`use_for_calcium_task = 1`). It is the largest file (9.37 GB) but the files
differ by less than 10 %. Nothing about it is known to be atypical.

## 1. Background for newcomers

**What was imaged.** One young-adult worm (strain KDK92; what it expresses
beyond "GCaMP, TagRFP-T and NeuroPAL" is not stated in the file; Documented
in file) held "semi-restricted" in a microfluidic device (DANDI description,
citing the Chronis 2007 device). "Semi-restricted" means the worm cannot
crawl away but can still bend and shift inside the channel, and the
measurements below confirm that it does: this is the first of our three
recordings with substantial motion. No stimulus was given; the activity is
spontaneous. The field covers the whole head from the tip of the nose
(IL1, IL2, OLQ, URA neurons) through the pharynx (I1–I4, M-cells, NSM, MC)
and nerve ring to the first ventral-cord motor neurons (VA1, VB2, DB1/DD1
region) and RIS; the posterior part of the field shows the diffuse
autofluorescence of the gut.

**Two colours in every frame.** Like ASCENT and unlike EY, the movie has a
stable reference channel (TagRFP-T, red) next to the calcium channel (GCaMP,
green). The file says the segmentation and tracking were done on the red
channel, and the stored traces are the green brightness inside those red
masks. A separate 4-colour NeuroPAL stack, taken once, carries the biological
names (assigned with the NeuroPAL Auto-ID software plus manual corrections,
Documented in the DANDI record).

**Volumes, slices, frames.** 25 z-slices per volume, 1,646 volumes at
1.67 Hz (16.4 min). Data are stored `time × x × y × z × channel`.

**Masks, not points.** The annotation is a `uint16` volume per frame in which
every voxel holds the ROI number it belongs to (0 = background). That gives,
for each neuron in each frame, a full 3-D mask (median 274 voxels), from
which we derive the centroid. ROI number i is the same cell in every frame
(so it is the track ID) and is the index into the name list and the stored
trace columns. ASCENT and EY gave one point per neuron per frame; here the
points are ours, the masks are the authors'.

## 2. Files

| File | Size | SHA-256 (DANDI, verified after download) |
|---|---|---|
| `sub-231013-02/sub-231013-02_ses-20231013T111100_ophys.nwb` | 9,365,258,734 B (8.72 GiB) | `50dbdc57705cc0231e6b8a3bf8862121bde865cd478be0d0ed7196747cde66f3` |

DANDI asset id `2732d548-8c7f-43e0-aec2-fd1d8edc12a1`. The asset metadata
(saved as `.asset.json` beside the file) records sex "Other", age
P3DT40260S (3 days, 11 h), strain absent at the DANDI level but KDK92 inside
the file, and the session "NeuroPAL+ Calcium Imaging without Stimulation"
starting 2023-10-13 11:11 JST. Download with `scripts/download_kk.sh`: about
3 min at ~55 MB/s on a single connection, digest checked afterwards.

The NWB file (`nwb_version` 2.6.0) embeds its specifications: core
2.6.0-alpha, hdmf-common 1.8.0, hdmf-experimental 0.5.0 and
ndx-multichannel-volume 0.1.12. pynwb 4.2 with `load_namespaces=True` reads
it in 0.4 s (Verified), but note that `CalciumSeriesSegmentation` comes back
as a `MultiChannelVolumeSeries`, not an `ImageSegmentation`, so there is no
`plane_segmentations` API; WormND's loaders special-case this dataset for
that reason.

## 3. Image structure

**Verified** from the HDF5 file (`inspect_kk.py`, section 1):

| Property | Calcium series | Label-volume series | NeuroPAL stack |
|---|---|---|---|
| Path | `acquisition/CalciumImageSeries/data` | `processing/CalciumActivity/CalciumSeriesSegmentation/data` | `acquisition/NeuroPALImageRaw/data` |
| Shape | `(1646, 288, 576, 25, 2)` | `(1646, 288, 576, 25)` | `(243, 461, 28, 4)` |
| Axes (Documented in file, `dimension` = `[288, 576, 25]`) | T, X, Y, Z, C | T, X, Y, Z | X, Y, Z, C |
| dtype | `uint16` | `uint16` | `uint16` |
| Storage | gzip level 1, chunks `(1, 36, 72, 7, 1)` | gzip, chunks `(1, 36, 72, 7)` | gzip, chunks `(31, 116, 7, 1)` |
| Size | 27.3 GB uncompressed, 9.13 GB on disk (ratio 3.0); 16.6 MB per frame | 13.7 GB uncompressed, 0.12 GB on disk | 25 MB uncompressed, 19 MB on disk |
| Value range | **0–255** in both channels (global max 255 over all frames) | 0–156 | not measured |
| Rate | 1.67 Hz, `starting_time` 0, no timestamps | same | — |
| `format` attribute | `"raw"`, "Raw GCaMP (channel 0) and RFP (channel 1) image series" | `"raw"` | — |

Reading one channel of one frame touches 256 gzip chunks and takes 0.04 s;
a full frame (both channels) 0.13 s; a label volume 0.02 s (mostly zeros).
Deriving centroids from all 1,646 label volumes takes 75 s; the streaming
pass over every frame and label volume (both channels, mask statistics)
takes 19 min with 460 MB peak memory. For comparison, EY's uncompressed
file streams in 8 s; ASCENT's gzip file in 2 min.

**Pixel values are 8-bit inside a 16-bit container (Verified measurement,
Inferred interpretation).** Every voxel is between 0 and 255. The images are
dim: whole-volume mean 3.0–3.3 counts in GCaMP and 4.4–4.8 in RFP, median
3 and 4, 99th percentile 10–14 and 15–18, background inside the worm about
3 (GCaMP) and 4.5 (RFP), neuron masks about 10–20 (GCaMP) and 16–23 (RFP).
The background is non-zero and noise-like (15–18 % of voxels are exactly
0, with values 1–8 roughly Poisson-distributed), so these look like camera
frames rescaled to 8 bits rather than thresholded data as in EY. The RFP
channel reaches 255 in every sampled frame and GCaMP reaches it early in the
recording, so the brightest voxels are clipped (Inferred). What the original
bit depth was and how the rescaling was done is Unknown (listed in
`docs/author_questions.md`).

**A six-column zero strip in the GCaMP channel (Verified).** In every
sampled frame, y columns 0–5 of channel 0 are entirely zero across all x and
z; channel 1 has no such strip and no other plane is all-zero. The simplest
reading is that the green channel was shifted by six pixels along y to
register it to the red channel and the vacated columns were zero-filled
(Inferred). It is visible as the black bar at the left of the GCaMP panels
in `previews/kk/mip_t*.png`.

**Spatial relation between the stacks (Verified):** matching the 143 named
neurons of the NeuroPAL stack to the calcium centroids at t = 0 by name
gives per-axis fits NeuroPAL = 0.97·calcium − 6.6 (x), 0.94·calcium − 86.3
(y), 1.18·calcium − 0.2 (z), with residual RMS 7.6, 6.6 and 1.7 voxels and
correlations 0.95, 0.997, 0.96 (the cross-axis correlations are all below
0.3, so no axis swap or flip). So the NeuroPAL stack is the same field at the
same pixel size, cropped (243 × 461 instead of 288 × 576, the y offset
being the crop), and the worm had moved and deformed by several microns
between the two acquisitions. This is consistent with the documented
0.32 µm pixel for both stacks, unlike the EY file.

## 4. Channels

**Documented (in file):** the calcium series has two optical channels in
order `488-525-50m` (GCaMP; excitation 488 nm, emission 525 nm, Semrock
FF03 525/50) and `561-617-73m` (TagRFP-T; 561 nm, 617 nm, Semrock FF02
617/73). The NeuroPAL stack has four: mTagBFP2 (405/445), CyOFP1
(488/617), mNeptune 2.5 (633/685) and Tag RFP-T (561/617), stored in the
order 405, 488, 633, 561 with `RGBW_channels = [0, 1, 2, 3]`. Taken
literally that puts the blue protein on the red display channel and the
far-red protein on blue, the reverse of the NeuroPAL convention; the
preview `neuropal_rgb_mip.png` follows the file and names the channels in
its title (Estimate: a metadata slip, see `docs/author_questions.md`).

**ASCENT's channel test applies here (Verified, same statistic).** On 40
sampled frames and 40 full-length tracks, the detrended coefficient of
variation of the 5 × 5 local mean at the rounded centroid is median 0.26
(90th percentile 0.60) in channel 0 and 0.155 (0.27) in channel 1, a ratio
of 1.69; 88 % of the sampled neurons vary more in channel 0. So channel 0
behaves like calcium and channel 1 like a stable reference, as the file
says. The correlation of mean neuron brightness across channels is only
0.20: which cells are bright differs between the two reporters, so the
reference channel cannot be used as a proxy for GCaMP brightness. Whole
volume means drift by 3.6 % (GCaMP) and 2.6 % (RFP) over the sampled
frames.

**Whole-population brightness over time (Verified,
`previews/kk/per_frame_means.png`):** the mean inside the masks falls from
19 to about 10 counts in GCaMP (1.7× the median at the start, 0.87× at the
end) and from 22 to 17 in RFP (1.12× to 0.85×). The RFP decline of about
25 % over 16 min is bleaching (Inferred); the larger GCaMP decline is
bleaching plus a real drop in activity, since the most active neurons are
busy in the first 350 s and quiet later (`neuron_traces.png`). The
background stays flat at 3 and 4.5 counts.

**Agreement with the stored traces (Verified):** the per-neuron Pearson
correlation between our GCaMP mask mean and the stored `SignalRawFluor`
trace is 1.000 for every neuron (minimum 0.999); a per-neuron linear fit
gives stored = 1.06–1.09 × ours with zero intercept, and the stored value is
exactly 0 in the 67 (neuron, frame) pairs where the label is absent. The
stored traces are therefore the mean GCaMP value inside these label masks,
up to a constant factor whose origin is Unknown (Estimate: computed from the
frames before the 8-bit rescaling, or from a slightly different mask). They
are not ratiometric: the GCaMP/RFP mask ratio correlates only 0.88 and the
RFP mask mean only 0.63. The `unit` attribute says "lumens", which is not
meaningful; values range 0–161.5. This is a much tighter correspondence
than EY (median 0.72) because the masks are in the file.

**Transient tracking failures (Verified presence, Documented warning):** the
DANDI description warns of "occasional transient failures in neuronal
tracking" and recommends a ±3-frame moving median on the traces. In the
stored traces 0.98 % of frame-to-frame steps change by more than 50 % and
0.28 % by more than 100 %; 70 of 156 neurons have at least one such jump.
Some of these are genuine calcium transients (the SMDVL trace rises 6× in a
few frames), so this count is an upper bound on tracking glitches.

**Spontaneous activity (Verified timing, Estimate biology):** the six
neurons with the largest swings in the stored traces are labelled RIML,
RIVL, AVAL, SMDVL, I4 and AVBR. AVA, AVB, RIM and SMD are names commonly
reported for neurons active during backward and forward movement and head
bending, which is plausible for a worm pushing against a restraint; we did
not check this further. The same events appear in our mask means (r = 1.00)
and not in the RFP channel of the same masks, so they are calcium signals,
not motion artefacts, except for AVBR where the RFP trace also fluctuates.

## 5. Voxel size and timing

| Quantity | Value | Status |
|---|---|---|
| xy pixel size | 0.32 µm (`grid_spacing`, both imaging volumes) | **Documented (in file)**; consistent between the stacks (section 3) |
| z step | 1.5 µm | **Documented (in file)** |
| Field of view (calcium) | 288 × 0.32 = 92 µm (x) by 576 × 0.32 = 184 µm (y) | Estimate from documented values; the worm lies along y |
| Depth span | 24 × 1.5 = 36 µm | Estimate from documented values |
| Volume rate | 1.67 Hz (0.6 s per volume) (`starting_time.rate`) | **Documented (in file)** |
| Duration | 1646 / 1.67 = 985.6 s = 16.4 min | Estimate from documented values |
| Per-frame timestamps | none; `starting_time` = 0 and a constant rate | **Verified absent** (same as ASCENT and EY) |
| Stimulus | none ("without Stimulation"; no `intervals` group) | **Documented (in file)** |
| Microscope | Nikon Eclipse Ti-U inverted, Yokogawa CSU-X1 SoRA spinning disk, 40× 1.3 NA objective | **Documented (in file)** |
| Strain, stage, temperature | KDK92, young adult ("YA"), 20 °C, sex "O" (DANDI: "Other"; hermaphrodite is an Estimate), born 2023-10-10 | **Documented (in file)** |
| Session date | 2023-10-13 11:11 JST; file written 2023-11-05 | **Documented (in file)** |

Anisotropy: one z step is 4.7 xy pixels. The stack is the thickest of the
three recordings (36 µm versus 30 µm for EY and 12 µm for ASCENT) and the
pixel the coarsest (0.32 µm versus 0.243 and 0.267 or 0.533).

## 6. Ground-truth annotations

**Verified** from the label-volume series (`inspect_kk.py`, section 4;
centroids cached in `previews/kk/calcium_tracks.csv`):

| Property | Value |
|---|---|
| Layout | one `uint16` volume per frame, voxel value i = ROI i, 0 = background ("Tracking results based on RFP image series"); plus `SegmentationVol1/Seg_tpoint_1`, a 156-row point table for one time point |
| Rows (label, frame) | 256,709 of 156 × 1,646 = 256,776 possible; 67 missing (0.03 %) |
| Voxels per label per frame | min 47, 10th percentile 130, median 274, 90th percentile 494, max 702; bounding box median 8 × 9 × 7 voxels (2.6 × 2.9 × 10.5 µm), never a single voxel or a single slice |
| Label values | 1–156, all present; max equals the length of `NeuronIDs/labels` |
| Time coverage | every one of the 1,646 frames has labels (151–156 per frame) |
| Tracks | ROI number = track; 151 of 156 present in every frame; the 5 others have gaps: RIS missing 42 frames (a run around t ≈ 680 where the frame count dips to 151), VA1 8, VB2 8, AVFR 6, RMGL 3 |
| Masks touching the stack edge | 5.5 % of rows touch z = 0, 3.7 % touch z = 24, 0.04 % touch the xy border |
| Centroid depth | z = 1.1–23.0, spread over the whole stack with no edge pile-up |
| Biological labels | `NeuronIDs/labels`: 156 entries, 143 real names (all unique, 56 left/right pairs) and 13 "unknown"; the NeuroPAL stack's own segmentation has 163 named neurons including all 143 |

**Provenance (Documented, Inferred):** the DANDI description says cells were
segmented and tracked with 3DeeCellTracker (Wen et al. 2021) on the RFP
channel and names were assigned with NeuroPAL Auto-ID plus manual
corrections. Whether the label volumes themselves were proofread is not
stated (Unknown). The warning about transient tracking failures suggests
they were not fully corrected (Estimate). WormND's loaders reduce these
volumes to centroids (mean label coordinate, missing frames interpolated or
copied from the previous frame), the same reduction used here.

**Motion (Verified from centroids):** this worm moves. The median in-plane
displacement between consecutive frames is 1.0 voxel (0.32 µm), the 90th
percentile 3.4, the 99th 10.5 (3.4 µm) and the maximum 75 (24 µm); only
31 % of steps are below half a voxel. Depth changes by up to 2 slices per
step at the 99th percentile and 9 at most. Taking the median over all ROIs
in each frame, 230 of the 1,645 frame pairs (14 %) have the whole
population moving more than 2 voxels, with the largest population jump
(58 voxels) at t = 1138–1139. Over the recording each track covers a region
of median 127 voxels (41 µm) in xy and 4 slices in z, although the first and
last frames differ by only 11 voxels in median: the worm shifts back and
forth along the channel rather than drifting away. For comparison ASCENT's
median step was 1 voxel with a maximum of 40 and EY's median 0 with a
maximum of 16; neither moved its whole head by tens of microns.

**Time-index convention (Verified, a test EY could not support):** on the
24 frames with the largest population motion (median step ≥ 8 voxels), the
RFP mask mean minus background computed with label volume t beats both
volume t−1 and t+1 in 24 of 24 cases, by a median margin of 36 %. Label
volume t describes image t. The separate point table
`SegmentationVol1/Seg_tpoint_1` matches label volume **0**: all 156 points
lie inside their own label (ROI i inside value i + 1) at t = 0, with a
median distance of 1.0 voxel to our centroid, whereas at t = 1 only 86 % do
and at t = 2 only 66 %. So "time point 1" is 1-based and names the first
frame (Verified for this file). The stored traces also start at frame 0.

**Completeness (unresolved):** 156 neurons are tracked, 143 named, against
the NeuroPAL stack's 163 named neurons and a head total of roughly 180–190
(general figure). The previews show bright cells without labels in the
posterior part of the field. Whether they are unidentifiable neurons or
non-neuronal cells is Unknown, as for ASCENT and EY.

## 7. Are the annotations in bounds and aligned with the images?

**Bounds (Verified):** centroids derived from a label volume with the same
(X, Y, Z) shape as the images are in bounds by construction; the real
question is whether the label volume's axes are the image's axes, which only
image content can answer.

**Axis identification and alignment (Verified, both channels, five frames
spread over the recording):** the mean 5 × 5 brightness at the rounded
centroids is 4.4–7.4× that at random points in the same bounding box in
GCaMP and 5.5–7.0× in RFP; with x and y swapped the ratio is 0.7–1.2, i.e.
chance. The mean inside each mask is 2.6–4.6× (GCaMP) and 3.2–4.0× (RFP) the
background mean. 58–68 % (GCaMP) and 71–78 % (RFP) of centroids lie within
3 pixels of the brightest pixel of their 11 × 11 in-plane neighbourhood. The
contrast ratios are lower than ASCENT's 24–83× because the 8-bit images put
neurons at 10–25 counts over a background of 3–5, and because random points
inside the dense head region often fall on other cells; the swapped-axes
control makes the conclusion unambiguous. ASCENT's background-subtracted
median ratio was not usable (the median random window equals the
background) and mean ratios are reported instead, as for EY.

- Visually (`previews/kk/*_gt.png`, `mip_movie.mp4`), circles and green
  mask outlines sit on cell bodies in both channels, including after the
  worm has moved (frames 1050–1249 in the movie).
- The masks are drawn on the red channel and the time index is verified
  (section 6), so the GCaMP values extracted through them are temporally
  aligned, subject to the six-column shift noted in section 3 having been
  applied before the masks were drawn (Inferred from the stored traces
  matching our extraction).

**Track-identity correctness (not verified directly):** no test confirms
that label i stays on the same cell through the large movements. Indirect
evidence: all 156 labels persist through the 58-voxel jump at t = 1138, the
stored traces are smooth apart from the jumps counted in section 4, and the
authors themselves warn of occasional failures. A tracking benchmark using
this file as ground truth should expect some label swaps around the
high-motion frames (Estimate).

## 8. Previews

All under `previews/kk/` (git-ignored, on scratch, 95 MB). Generated by
`scripts/make_previews_kk.py` in 65 s with 450 MB peak memory, and
`scripts/plot_kk_traces.py`. Images are shown with the stored x axis
vertical and y horizontal so the worm lies left to right (nose at the left).

| File | Content |
|---|---|
| `depth_slices_t{0,800}_{gcamp,rfp}.png` (+`_gt`) | 5 × 5 grid of the 25 z-slices of one frame and channel, with/without centroid circles |
| `mip_t{0,800}.png` (+`_gt`) | GCaMP (top) and TagRFP-T (bottom) max projections; `_gt` adds centroids on the GCaMP panel and label-mask outlines on the RFP panel |
| `neuropal_rgb_mip.png`, `neuropal_rgb_mip_labels.png` | the NeuroPAL stack with its 163 neuron centres and, in the second, their names; colours follow the file's `RGBW_channels` (see section 4) |
| `mip_movie.mp4` | frames 1050–1249 (629–748 s, covering the largest movements) as a 2 × 2 panel: GCaMP, RFP, GCaMP with centroids and names, RFP with mask outlines; 10 fps, i.e. **6× faster than real time** (39 MB; the noisy 8-bit images compress poorly) |
| `per_frame_means.csv/.png` | per-frame mean, 99th percentile, max, mask mean and background mean for both channels, ROI count, and the population motion per frame |
| `per_frame_means_tracks.csv` | our mask mean and 5 × 5 local mean of every ROI in every frame, both channels (long format, NaN where the label is absent) |
| `stored_traces.csv` | the file's `SignalRawFluor` traces, 1646 × 156 |
| `calcium_tracks.csv` | derived per-frame table: t, roi, label value, centroid x y z, voxel count, bounding box, name (24 MB) |
| `neuron_traces.png` | the six neurons with the largest swings in the stored traces: stored trace, our GCaMP mask mean, our RFP mask mean |
| `inspection_summary.json` | machine-readable output of the inspection script |

Display intensities are fixed per channel at the 1st to 99.9th percentile of
frame 800 for all frames.

## 9. Inventory entry (draft)

| Field | Value |
|---|---|
| Source | DANDI 000692 v0.240402.2118, Kimura lab (Nagoya City University) with Wen (RIKEN) and Sprague (UCSF); publication "unpublished"; WormID-Bench "KK" |
| Original ID | `sub-231013-02`, identifier `20231013-11-11-0`, subject `231013_02`; asset `2732d548-8c7f-43e0-aec2-fd1d8edc12a1` |
| Format | NWB 2.6.0 (HDF5) with `ndx-multichannel-volume` 0.1.12; opens with pynwb 4.2 in 0.4 s; segmentation is a `MultiChannelVolumeSeries`, not an `ImageSegmentation` |
| Dimensions | calcium T=1646, C=2, X=288, Y=576, Z=25, uint16 holding 0–255; label volumes same T, X, Y, Z, uint16 0–156; NeuroPAL X=243, Y=461, Z=28, C=4, uint16 (Verified) |
| Voxel size | 0.32 × 0.32 × 1.5 µm, same for both stacks (Documented in file; consistent with the measured overlap) |
| Rate / duration | 1.67 Hz / 985.6 s (Documented in file / Estimate) |
| Channels | GCaMP (488/525) and TagRFP-T (561/617) in every frame; NeuroPAL mTagBFP2, CyOFP1, mNeptune 2.5, TagRFP-T (Documented in file) |
| Condition | semi-restricted in a microfluidic device, no stimulus, spontaneous activity (Documented) |
| Movement | median per-step displacement 1.0 voxel, p99 10.5, max 75; whole-population jumps > 2 voxels in 14 % of frames; tracks cover a median 41 µm region (Verified) |
| Preprocessing | 8-bit values in a uint16 container, brightest voxels clipped at 255, GCaMP channel shifted 6 columns with zero fill, despite `format: raw` (Verified values, Inferred cause); original bit depth and registration Unknown |
| Annotations | 156 tracked neurons as full 3-D masks in all 1,646 frames (256,709 masks, median 274 voxels, 67 missing) plus one point table for frame 0 (Verified) |
| Annotation provenance | 3DeeCellTracker on the RFP channel (Documented); manual correction Unknown; authors warn of transient tracking failures |
| Biological identities | 143 of 156 ROIs named, 13 "unknown", 56 L/R pairs, from NeuroPAL Auto-ID plus manual correction (Verified presence, Documented method) |
| Stored traces | `SignalRawFluor` 1646 × 156 = 1.06–1.09 × mean GCaMP inside the masks, 0 where the mask is missing (Verified) |
| License | CC-BY-4.0 (Documented) |
| Verification status | file downloaded, SHA-256 matches, structure, axis order, spatial alignment and time-index convention verified; track-identity correctness, annotation completeness, pixel-value provenance and manual-correction status not verified |
| Other assets | 8 more recordings with the same layout, 8.58–9.00 GB, 1,645–1,646 frames, 147–159 ROIs, 127–145 real names (Verified remotely, headers only) |

## 10. Open questions and uncertainties

1. **Pixel values.** 8-bit values in a 16-bit container, clipped at 255 in
   the brightest cells, with a six-column zero strip in GCaMP. The original
   bit depth, the rescaling rule and the channel registration are Unknown.
   For a benchmark the file provides *rescaled* images; a request for the
   original frames is in `docs/author_questions.md`.
2. **Manual correction of the tracking.** The DANDI text documents the
   tracker and warns of failures; whether any frame was proofread is
   Unknown. Until answered, these are algorithm-generated annotations, to be
   kept apart from the proofread ASCENT tracks in the inventory.
3. **Track-identity correctness** was not tested directly, and this is the
   recording where it matters most because the worm moves by tens of microns.
4. **Annotation completeness.** 156 of perhaps 180–190 head neurons; the
   NeuroPAL stack names 163. No audit was done.
5. **Stored-trace scale factor.** Stored traces are 6–9 % above our mask
   means with r = 1.000; the cause is Unknown and harmless for tracking.
6. **NeuroPAL colour mapping.** `RGBW_channels = [0, 1, 2, 3]` over the
   stored order 405/488/633/561 nm reverses the usual colours; a slip or
   intended, Unknown. It does not affect the identities, which are stored as
   names.
7. **Timestamps.** None per frame; the 1.67 Hz rate is a single number.
   Within-volume slice timing is Unknown, and matters more here than for
   ASCENT or EY because the worm can move during the 0.6 s a volume takes.
8. **Strain details and publication.** KDK92 is not described in the file
   and the publication is "unpublished"; the GCaMP variant and promoters
   are Unknown.
9. **Which frame is "time point 1".** Verified as frame 0 for this file;
   assumed, not checked, for the other 8.

## 11. Comparison with ASCENT and EY

| | ASCENT (Zenodo 17561700) | EY sub-20190928-13 (DANDI 000541) | KK sub-231013-02 (DANDI 000692) |
|---|---|---|---|
| Container | plain HDF5, no metadata | NWB 2.6 + extension; metadata in file | NWB 2.6 + extension; metadata in file; no stimulus table |
| Volume | 2 ch × 9 × 512 × 512, 0.243 µm, 3.3 Hz, 1100 frames, 333 s | 1 ch × 21 × 128 × 256, 0.267 (or 0.533) µm, 4 Hz, 962 frames, 240 s | 2 ch × 25 × 576 × 288, 0.32 µm, 1.67 Hz, 1646 frames, 986 s |
| Depth coverage | 12 µm slab, edge pile-up | 30 µm, four tracks extrapolated beyond the stack | 36 µm, masks spread over the stack, 9 % touch an end slice |
| Reference channel | mNeptune in every frame | none in the movie; static NeuroPAL stack | TagRFP-T in every frame, plus a static NeuroPAL stack |
| Pixel values | uint16 camera counts, baseline ~320, gzip | uint16 holding 0–255, zero background, uncompressed | uint16 holding 0–255, background 3–5, clipped at 255, gzip (3.0×) |
| Annotations | 96 tracks × 1100 frames, points, ZephIR + manual proofreading, no identities | 175 tracks × 961 frames, single-voxel points, dNMF, proofreading unknown, all named | 156 tracks × 1646 frames, full 3-D masks (median 274 voxels), 3DeeCellTracker, proofreading unknown, 143 named |
| Motion (per-step median / max, voxels) | 1 / 40 | 0 / 16 | 1 / 75, with whole-head jumps of up to 58 voxels; tracks cover a median 41 µm |
| Stimulus | documented schedule, timing inferred | timing in file, responses verified | none; spontaneous activity in locomotion-associated neurons |
| Time-index convention | inferred | inferred (unverifiable, no motion) | **verified** (24 / 24 frames) |
| Spatial alignment | verified (contrast 24–83×) | verified (19–38×) | verified (4–7× mean ratio, swapped axes at chance) |
| Stored traces | none | unitless, r = 0.72 with our extraction | mask means, r = 1.000 with our extraction |
| Per-frame access | one gzip dataset per frame-channel (4.7 MB) | one contiguous 1.4 MB slab | 256 gzip chunks per frame-channel, 0.04 s |
| Open issues shared | timestamps, completeness, track-identity | same, plus pixel provenance, voxel-size inconsistency | same, plus pixel provenance, manual-correction status |

For a tracking benchmark KK adds what ASCENT and EY lack: a worm that
really moves, full 3-D masks instead of points, a reference channel *and*
biological names, and 16 min of recording. Against that, its images are
8-bit and clipped, its annotations are algorithm output with a documented
failure rate and no stated proofreading, and each file is 9 GB. The three
recordings now span the motion range from immobilised (EY) through
slightly drifting (ASCENT) to semi-restricted (KK), which is the axis the
benchmark most needs to cover.
