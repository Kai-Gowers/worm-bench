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
