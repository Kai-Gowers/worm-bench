# C. elegans Calcium-Imaging Tracking Benchmark

## Project overview

We are developing a comprehensive benchmark for tracking individual
neurons through 3D calcium-imaging recordings of C. elegans.

WormID-Bench is our starting reference. Our goal is to incorporate its
tracking datasets and expand the collection across laboratories,
microscopes, experimental conditions, and worm movement patterns.

The project is currently in the data discovery and inspection phase.
Benchmark design, evaluation, and model development will follow once
we understand the available recordings and annotations.

## Research goals

We want to create a resource that makes it easier to develop and
evaluate neuron-tracking methods across diverse imaging conditions.

Potential contributions include:
- Broader coverage than existing tracking benchmarks.
- Access to valuable recordings that are not readily downloadable.
- Reliable, documented tracking annotations.
- Consistent data loading and annotation formats.
- Reproducible evaluation and meaningful generalization tests.
- Baseline methods and analysis of where existing approaches fail.

Simply combining existing datasets does not automatically establish
research novelty. The eventual contribution depends on the collection's
coverage, annotation quality, evaluation design, and scientific utility.

## Data of interest

Our primary interest is actual 3D calcium image time series with
individual neurons visible.

Extracted activity traces, static structural images, and supplementary
demonstration videos can provide useful context, but do not replace the
underlying image recordings needed for a tracking benchmark.

A recording typically contains repeated 3D volumes, conceptually:
time × channel × depth × height × width.

Available information may include:
- Calcium fluorescence images, whose brightness changes with activity.
- Stable reference fluorescence for locating neurons.
- NeuroPAL structural images for biological neuron identification.
- Neuron center coordinates or segmentation masks.
- Track IDs connecting the same neuron across time.
- Biological neuron identities.
- Extracted activity traces.

Track IDs and biological identities are different: a neuron can be
tracked consistently without knowing its biological name.

Image quality and annotation quality are also separate considerations.
A useful movie may lack annotations, while a challenging movie may have
excellent ground truth.

## Collection approach

We will inspect representative recordings before collecting entire
datasets or choosing a common storage format.

The initial examples are:
- ASCENT, which provides a clearly documented recording and tracking
  annotations.
- A recording from WormID-Bench, to understand its data organization
  and establish compatibility with the existing benchmark.

As the collection grows, we will distinguish:
- Verified recordings from unverified data leads.
- Original images from filtered or geometrically transformed images.
- Total time points from annotated time points.
- Manual annotations from algorithm-generated or corrected results.
- Independent recordings from copies reused across publications.

The collection should preserve original scientific information and
record the provenance of any transformations.

## Storage plan

GitHub hosts project code, documentation, and small metadata files.
Hugging Face is the intended platform for sharing the dataset.

Data processing takes place on a remote host. The workflow should
support individual recordings or bounded batches without requiring
the entire collection to fit on disk or in memory.

Small image and video previews will support visual inspection.
Full-resolution scientific image data remains the basis of the benchmark.

The standardized image format has not been selected yet. That decision
will depend on the source formats, preservation of pixel values and
metadata, compression, and efficient access to selected frames.

Redistribution depends on each source's license and sharing conditions.

## Recording inventory

The inventory will describe each recording's:
- Source dataset, laboratory, publication, version, and original ID.
- Files, sizes, and checksums where available.
- Image format, dimensions, axis order, and pixel dtype.
- Frame count, channels, voxel sizes, timing, and duration.
- Experimental conditions and degree of worm movement.
- Preprocessing history.
- Annotation types, neuron counts, and temporal coverage.
- Annotation provenance and biological identity coverage.
- License, attribution, and verification status.
- Remaining uncertainties.

Missing information remains explicitly unknown. Estimates are
distinguished from measurements verified in the files.

## Starting references

WormID-Bench reports 36 calcium recordings used for tracking:
- EY: 21 — https://dandiarchive.org/dandiset/000541
- KK: 9 — https://dandiarchive.org/dandiset/000692
- SK1: 4 — https://dandiarchive.org/dandiset/000565
- SK2: 2 — https://dandiarchive.org/dandiset/000472

These are the paper's benchmark counts, not a verified inventory of
every file currently available in those repositories.

Benchmark code:
https://github.com/focolab/WormND

Additional data leads:
- ASCENT: https://zenodo.org/records/17561700
- ZephIR: https://zenodo.org/records/10818810
- Leifer signal propagation atlas:
  https://dandiarchive.org/dandiset/001075
- Kato / Dunn 2025:
  https://dandiarchive.org/dandiset/001623
- Flavell subset associated with WormID:
  https://dandiarchive.org/dandiset/000776
- Targettrack: https://github.com/rahi-lab/targettrack
- See Elegans: https://github.com/enricoliit/SeeElegans
- 3DeeCellTracker: https://github.com/WenChentao/3DeeCellTracker

These sources require inspection for actual calcium-image availability,
annotation coverage, recording lengths, and overlap.

## Collaboration context

I have experience with Python, machine learning, and 3D cell tracking,
but I am new to C. elegans imaging and neuroscience.

Plain-language explanations of unfamiliar biology and imaging concepts
are helpful. Findings should clearly separate verified facts, estimates,
and unresolved questions.

Specific implementation tasks and the order of work will be provided
through individual prompts.
## Version control

Commit regularly: after each completed step (a new script, a documented
finding, a configuration change), not only at the end of a task. Use short,
descriptive commit messages.

Only code, documentation, configuration, and small metadata files are
tracked. Never commit downloaded recordings, annotation CSVs, generated
previews (PNG, MP4, GIF), or anything under `data/` and `previews/`; these
live on scratch storage and are listed in `.gitignore`. Check `git status`
before committing to confirm no data files are staged.
