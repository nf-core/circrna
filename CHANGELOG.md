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
- [#242](https://github.com/nf-core/circrna/pull/242) - nf-tests for every local module and subworkflow (real and `-stub` tests), based on a small test dataset hosted in `nf-core/test-datasets` (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - Stub blocks for all local modules (@nictru)

### `Changed`

- [#242](https://github.com/nf-core/circrna/pull/242) - Local subworkflows now live in `subworkflows/local/<name>/main.nf` and no longer read `params` directly; only `NFCORE_CIRCRNA` in `main.nf` does (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - The pipeline test no longer runs the miRNA prediction and finishes in about half the time (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - CI runs module, subworkflow and pipeline tests as separate shard groups (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - Local modules no longer read non-canonical `task.ext.*` keys (`suffix`, `max_shift`, `consider_strand`, `min_tools`, `metacols`, `has_header`). The values are process inputs now or fixed in the module, and unused `ext` settings are removed from `conf/modules.config` (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - Updated the nf-core modules `gawk`, `gnu/sort`, `bedtools/intersect`, `bedtools/sort`, `custom/gtffilter` and `seqkit/fx2tab`, which take the output suffix as an input (@nictru)

### `Fixed`

- [#242](https://github.com/nf-core/circrna/pull/242) - `exons_only = false` was treated as `true` in `ANNOTATION` (`bed2gtf`) (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `BSJ_DETECTION` failed when no STAR-based detection tool was selected (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `FLI_DETECTION` overwrote the genome FASTA input, which stalled CIRI-full and gave circtools the wrong FASTA (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `CIRIVIS_UNIFY` did not emit the sample meta map, which broke the CIRI-full subworkflow (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `MIRNA_BINDINGSITES` failed when only one of miranda or TargetScan was selected (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `jccirc` is now accepted by the `fli_tools` parameter schema (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `PYGTFTK_TABULATE` failed under docker because gtftk could not write to the home directory (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - An empty tool list (for example `--fli_tools ''`) is now parsed as an empty list. Previously it counted as one tool with an empty name, so single-end samples were rejected although no FLI tool was selected, and the check for an empty `--quantification_tools` never triggered (@nictru)

### `Dependencies`

### `Deprecated`
