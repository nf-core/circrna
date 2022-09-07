# ![nf-core/circrna](docs/images/nf-core-circrna_logo_light.png#gh-light-mode-only) ![nf-core/circrna](docs/images/nf-core-circrna_logo_dark.png#gh-dark-mode-only)

**A reproducible workflow to characterize circular RNAs in RNA-Seq datasets**

[![GitHub Actions CI Status](w to characterize circular RNAs in RNA-Seq datasets**.)](https://github.com/nf-core/circrna/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/circrna/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/circrna/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?logo=Amazon%20AWS)](https://nf-co.re/circrna/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/circrna)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23circrna-4A154B?logo=slack)](https://nfcore.slack.com/channels/circrna)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/circrna** is a best-practice analysis pipeline for the quantification, miRNA target prediction and differential expression analysis of circular RNAs in paired-end RNA sequencing data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<p markdown="1" align="center">
<img src="docs/images/workflow.png" alt="workflow" width="100%" height="100%">
</p>

## Pipeline Summary

1. Input type conversion [`SamToFastq`](https://gatk.broadinstitute.org/hc/en-us/articles/360036485372-SamToFastq-Picard-)
2. Raw read quality control [`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
3. Adapter trimming + read filtering [`BBDUK`](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
4. circRNA quantification
    1. [`STAR`](https://github.com/alexdobin/STAR) -> [`CIRCexplorer2`](https://circexplorer2.readthedocs.io/en/latest/)
    2. [`STAR`](https://github.com/alexdobin/STAR) -> [`circRNA finder`](https://github.com/orzechoj/circRNA_finder)
    3. [`STAR`](https://github.com/alexdobin/STAR) -> [`DCC`](https://github.com/dieterich-lab/DCC)
    4. [`HISAT2`](http://daehwankimlab.github.io/hisat2/) -> [`CIRI2`](https://sourceforge.net/projects/ciri/files/CIRI2/) -> [`BWA`](http://bio-bwa.sourceforge.net/) -> [`CIRIquant`](https://github.com/Kevinzjy/CIRIquant)
    5. [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) -> [`find circ`](https://github.com/marvin-jens/find_circ)
    6. [`Bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) -> [`MapSplice`](http://www.netlab.uky.edu/p/bioinfo/MapSplice2)
    7. [`Segemehl`](https://www.bioinf.uni-leipzig.de/Software/segemehl/) -> [`Segemehl`](https://www.bioinf.uni-leipzig.de/Software/segemehl/)
5. circRNA filtering
    1. Filter candidate circRNAs by number of reads spanning back-splice junction
6. circRNA annotation
    1. Annotate candidates as circRNA, ciRNA or EI-circRNA
    2. Calculate mature spliced length
    3. Export mature spliced length as FASTA file
    4. Annotate parent gene, underlying transcripts.
    5. Export information as customised BED12 file
7. circRNA count matrix
    1. Combine results of quantification tools to produce counts matrix for downstream statistical analysis
    2. Require circRNAs in matrix to be called by at least *n* quantification tools (consensus filtering)
8. miRNA target prediction
    1. [`miRanda`](http://cbio.mskcc.org/miRNA2003/miranda.html)
    2. [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi)
    3. Filter results, miRNAs must be called by both tools
9. Differential expression analysis [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
10. MultiQC report [`MultiQC`](http://multiqc.info/)

Ouputs given by each step in the pipeline can be viewed in the [output documentation](https://nf-co.re/circrna/dev/output)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/circrna -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```console
   nextflow run nf-core/circrna -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --module 'circrna_discovery, mirna_prediction, differential_expression' --tool 'circexplorer2' --input 'samples.csv' --input_type 'fastq' --phenotype 'phenotype.csv'
   ```

## Documentation

The nf-core/circrna pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/circrna/usage) and [output](https://nf-co.re/circrna/output).

## Credits

`nf-core/circrna` was originally written by Barry Digby ([@BarryDigby](https://github.com/BarryDigby)) from the [National University of Ireland, Galway](http://www.nuigalway.ie/index-internal.html) as a member of Dr. Pilib Ó Broins lab with the financial support of Science Foundation Ireland (Grant number 18/CRT/6214). The authors would like to thank members of the nf-core team for their assistance reviewing the workflow.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#circrna` channel](https://nfcore.slack.com/channels/circrna) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
