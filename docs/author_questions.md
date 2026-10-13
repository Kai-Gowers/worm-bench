# Questions for dataset authors

Open points that could not be settled from the files and that matter for
reusing the recordings in a tracking benchmark. Each is phrased so it can be
sent as is. The evidence behind each question is in the linked inspection
note.

## EY (DANDI 000541, Yemini / Hobert lab) — see `ey_inspection.md`

Recording inspected: `sub-20190928-13` (asset `c4a52d3e-26a9-439c-8965-13320ab6cc43`).

1. **Time-index convention.** `CalciumSeriesSegmentation` holds
   `Seg_tpoint_0 … Seg_tpoint_960` for a 962-frame movie. Does
   `Seg_tpoint_t` describe frame `t` of `CalciumImageSeries` (so the last
   frame is unannotated), or frame `t+1`? The worm does not move enough for
   us to tell from the images.
2. **Which position set is the proofread one.** The file holds
   `CalciumSeriesSegmentation` and `CalciumSeriesSegmentationdNMF` (and two
   matching trace sets). They differ in about 90 % of rows, usually by about
   one voxel, occasionally by tens of voxels. Was either set checked or
   corrected by hand, and which should be treated as ground truth?
3. **Calcium voxel size.** Both imaging volumes state 0.2667 µm pixels, but
   the calcium stack (256 × 128) is exactly half the NeuroPAL stack
   (512 × 256) in x and y when matched neuron by neuron. Is the calcium
   pixel 0.533 µm (2 × 2 binning), or is the NeuroPAL pixel 0.133 µm?
4. **Pixel-value processing.** The calcium series is stored as uint16 but
   holds only values 0–255, with 87–97 % of voxels exactly zero, although the
   `format` attribute says `raw`. What processing (background subtraction,
   thresholding, rescaling) was applied before the NWB export, and are the
   original camera frames available?
5. **Out-of-stack z values.** Four ROIs (FLPR, ADFL, RIVR, SMDVL) have
   z = 21–24 for many frames, beyond the 21-slice stack. Should these rows be
   treated as "outside the imaged volume" or clamped to the last slice?

## KK (DANDI 000692, Kimura lab) — see `kk_inspection.md`

Recording inspected: `sub-231013-02` (asset `2732d548-8c7f-43e0-aec2-fd1d8edc12a1`).

1. **Pixel-value processing.** Both channels are stored as uint16 but hold
   only 0–255 (8-bit) values with a background of about 3–5 counts, and the
   GCaMP channel has six all-zero columns at the start of the y axis. Were the
   images converted to 8 bits and were the channels registered to each other
   before export? Are the original 16-bit frames available?
2. **Manual correction.** The DANDI description says segmentation and
   tracking were done with 3DeeCellTracker and warns of occasional transient
   tracking failures. Were the label volumes proofread or corrected by hand
   at any time point, and if so which?
3. **Reference frame.** `SegmentationVol1/Seg_tpoint_1` matches label volume
   index 0 (every stored point lies inside its label there). Is "time point
   1" the first frame of the movie (1-based naming), as this suggests?
4. **NeuroPAL colour mapping.** `RGBW_channels = [0, 1, 2, 3]` with the
   channel order 405/445, 488/617, 633/685, 561/617 puts the 405 nm (mTagBFP2)
   channel on red and the 633 nm (mNeptune) channel on blue, the reverse of
   the NeuroPAL convention. Is the stored mapping intended?
5. **Strain and publication.** The strain is given as KDK92 and the related
   publication as "unpublished". Is there a reference (or a strain
   description) that states the GCaMP variant, promoters and the NeuroPAL
   background?

## SK1 and SK2 (DANDI 000565 and 000472, Kato lab) — see `sk_inspection.md`

Recordings inspected: `sub-20230322-20-16-50` (000565, asset
`36715c41-aa5c-4a05-affd-66a9a6089f7f`) and `sub-2022-04-26-w01-NP1` (000472,
asset `8525ab3c-4a5d-4893-a2c8-4ac7ccbd87bd`).

1. **Stimulus.** `session_description` says "immobilized worm with optogenetic
   stimulus" and the SK1 `experiment_description` says "optogenetic stimulation
   of single neurons", but the dandiset description says the calcium data are
   unstimulated and no stimulus timing is stored. Was any stimulus delivered
   during these movies, and if so when and to which neuron?
2. **Annotation provenance.** The per-frame ROI positions move by at most one
   voxel between frames and never change z. Were the positions detected once
   (on which frame or on the NeuroPAL stack?) and then propagated, or tracked
   frame by frame? Were any positions checked by hand?
3. **Stored DF/F.** The stored `SignalDFoF` traces agree well with the images
   in SK2 but poorly in SK1 (median per-neuron correlation 0.06 against a
   5×5 mean at the stored positions). How were the traces computed (ROI
   shape, background and F0 definition, bleach correction, any smoothing),
   and from which positions?
4. **Calcium z spacing.** Matching named neurons between the NeuroPAL stack
   and calcium frame 0 gives a z scale of 3.8 NeuroPAL slices per calcium
   slice in SK2, where the stored spacings (0.75 and 1.5 µm) predict 2.0. In
   SK1 the sign of the z relation is reversed. Is the calcium z step really
   1.5 µm (SK2) and 2.5 µm (SK1), and are the two stacks acquired in the same
   z direction?
5. **Uncertain identities.** SK2 labels end in "?" for 9 of the 32 named
   neurons (e.g. `AINR?`). Should these be treated as unlabelled or as best
   guesses, and does the same convention apply to SK1 (where no "?" occurs)?
6. **Labels in calcium vs NeuroPAL.** In SK1, 18 of the 31 named calcium ROIs
   have no counterpart among the 44 named NeuroPAL centres, and in SK2 9 of
   32. Were the calcium IDs assigned independently of the NeuroPAL
   segmentation stored in the file?
7. **Camera settings.** Both movies carry a ~95-count camera offset and SK2
   saturates at 4095 in the first frames. Could you confirm the camera bit
   depth, exposure time per slice and the volume rate (1.04 Hz SK1,
   3.3 Hz SK2)?
