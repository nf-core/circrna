# nf-core/circrna: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.0.1dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/circrna, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Multiple methods for BSJ detection (CIRIquant, find_circ, segemehl, mapsplice, circexplorer2, circrna_finder, DCC)
- Multiple methods for circRNA quantification (psirc-quant, CIRIquant, sum and max aggregations)
- DB and GTF-based circRNA annotation
- MiRNA target prediction (TargetScan, miRanda) and correlation analysis
- Basic statistical analyses (CircTest, CIRIquant differential expression)
- nf-tests for every local module and subworkflow (real and `-stub` tests), based on a small test dataset hosted in `nf-core/test-datasets`
- Stub blocks for all local modules

### `Changed`

- Local subworkflows now live in `subworkflows/local/<name>/main.nf` and no longer read `params` directly; only `NFCORE_CIRCRNA` in `main.nf` does
- The pipeline test no longer runs the miRNA prediction and finishes in about half the time
- CI runs module, subworkflow and pipeline tests as separate shard groups

### `Fixed`

- `exons_only = false` was treated as `true` in `ANNOTATION` (`bed2gtf`)
- `BSJ_DETECTION` failed when no STAR-based detection tool was selected
- `FLI_DETECTION` overwrote the genome FASTA input, which stalled CIRI-full and gave circtools the wrong FASTA
- `CIRIVIS_UNIFY` did not emit the sample meta map, which broke the CIRI-full subworkflow
- `MIRNA_BINDINGSITES` failed when only one of miranda or TargetScan was selected
- `jccirc` is now accepted by the `fli_tools` parameter schema

### `Dependencies`

### `Deprecated`
