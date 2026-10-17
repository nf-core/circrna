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
- [#243](https://github.com/nf-core/circrna/pull/243) - The per-tool BSJ count matrices are now published as `<tool>.joined.tsv` (tab-separated) instead of `<tool>.csv`. Counting is done by `JOIN_SAMPLES`, which now also works with a single sample (@nictru)
- [#243](https://github.com/nf-core/circrna/pull/243) - Renamed the config selectors `CIRCTOOLS:MAIN` to `CIRCTOOLS:DETECT` and `ANNOTATE_COMBINED:BED2GFF|BED2GTF` to `ANNOTATE_COMBINED:BED2(GFF|GTF)`, so that the options are applied to the intended processes only (@nictru)

### `Fixed`

- [#242](https://github.com/nf-core/circrna/pull/242) - `exons_only = false` was treated as `true` in `ANNOTATION` (`bed2gtf`) (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `BSJ_DETECTION` failed when no STAR-based detection tool was selected (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `FLI_DETECTION` overwrote the genome FASTA input, which stalled CIRI-full and gave circtools the wrong FASTA (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `CIRIVIS_UNIFY` did not emit the sample meta map, which broke the CIRI-full subworkflow (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `MIRNA_BINDINGSITES` failed when only one of miranda or TargetScan was selected (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `jccirc` is now accepted by the `fli_tools` parameter schema (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - `PYGTFTK_TABULATE` failed under docker because gtftk could not write to the home directory (@nictru)
- [#242](https://github.com/nf-core/circrna/pull/242) - An empty tool list (for example `--fli_tools ''`) is now parsed as an empty list. Previously it counted as one tool with an empty name, so single-end samples were rejected although no FLI tool was selected, and the check for an empty `--quantification_tools` never triggered (@nictru)
- [#243](https://github.com/nf-core/circrna/pull/243) - `CIRCTOOLS_DETECT` failed for single-end samples (#230, #210) (@nictru)
- [#243](https://github.com/nf-core/circrna/pull/243) - `CIRI` BED output used a 1-based start, so CIRI BSJs were dropped during annotation (#239) (@nictru)
- [#243](https://github.com/nf-core/circrna/pull/243) - The pipeline now fails with an error if none of the BSJ detection tools finds a BSJ (`FAIL_ON_EMPTY`) (#238) (@nictru)
- [#243](https://github.com/nf-core/circrna/pull/243) - Per-tool count matrices failed when only a single sample was present (#215) (@nictru)
- [#243](https://github.com/nf-core/circrna/pull/243) - `JCCIRC` never ran in the pipeline; it now receives its inputs correctly, its CPU request is clamped to at least 1, and `tools_with_reads` is validated (#231) (@nictru)

### `Dependencies`

- [#243](https://github.com/nf-core/circrna/pull/243) - Updated the nf-core module `csvtk/join` (its stub referenced an undefined `args` variable). The output channel is now called `out_file` (@nictru)

### `Deprecated`
