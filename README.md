# ![nf-core/circrna](docs/images/nf-core-circrna_logo.png)

**workflow for the quantification, differential expression analysis and miRNA target prediction analysis of circRNAs in RNA-Seq data**.

[![GitHub Actions CI Status](https://github.com/nf-core/circrna/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/circrna/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/circrna/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/circrna/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/circrna.svg)](https://hub.docker.com/r/nfcore/circrna)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23circrna-4A154B?logo=slack)](https://nfcore.slack.com/channels/circrna)

## Introduction

**nf-core/circrna** is a bioinformatics pipeline used for the quantification, miRNA target prediction and differential expression analysis of circular RNAs in RNA sequencing data. Currently, the pipeline only supports the identification of circular RNAs in Human RNA-Seq data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Workflow

The diagram below gives an overview of the modules in `nf-core/circrna`:

<p markdown="1" align="center">
<img src="docs/images/workflow.png" alt="workflow" width="500">
</p>

## Pipeline Summary

1. Download reference genome files ([`Gencode`](https://www.gencodegenes.org/))
2. Download miRNA database files ([`miRbase`](http://www.mirbase.org/ftp.shtml), [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi))
3. Adapter trimming ([`BBDUK`](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/))
4. Read QC (([`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)))
5. Generate genome indices
6. circRNA quantification
    1. [`STAR`](https://github.com/alexdobin/STAR) -> [`CIRCexplorer2`](https://circexplorer2.readthedocs.io/en/latest/)
    2. [`STAR`](https://github.com/alexdobin/STAR) -> [`circRNA finder`](https://github.com/orzechoj/circRNA_finder)
    3. [`STAR`](https://github.com/alexdobin/STAR) -> [`DCC`](https://github.com/dieterich-lab/DCC)
    4. [`HISAT2`](http://daehwankimlab.github.io/hisat2/) -> [`CIRI2`](https://sourceforge.net/projects/ciri/files/CIRI2/) -> [`BWA`](http://bio-bwa.sourceforge.net/) -> [`CIRIquant`](https://github.com/Kevinzjy/CIRIquant)
    5. [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) -> [`find circ`](https://github.com/marvin-jens/find_circ)
    6. [`Bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) -> [`MapSplice`](http://www.netlab.uky.edu/p/bioinfo/MapSplice2)
7. miRNA target prediction ([`miRanda`](http://cbio.mskcc.org/miRNA2003/miranda.html), [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi))
8. DESeq2 differential expression analysis ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))

Ouputs given by each step in the pipeline can be viewed at the [output documentation](https://nf-co.re/circrna/dev/output)

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/circrna -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/circrna -profile <docker/singularity/podman/conda/institute> --input '*_R{1,2}.fastq.gz' --genome GRCh37
    ```

See [usage docs](https://nf-co.re/circrna/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/circrna pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/circrna/usage) and [output](https://nf-co.re/circrna/output).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

`nf-core/circrna` was originally written by Barry Digby [(@BarryDigby)](https://github.com/BarryDigby) from the [National University of Ireland, Galway](http://www.nuigalway.ie/index-internal.html) as a member of Dr. Pilib Ó Broins lab with the financial support of Science Foundation Ireland (Grant number 18/CRT/6214).


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#circrna` channel](https://nfcore.slack.com/channels/circrna) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/circrna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
