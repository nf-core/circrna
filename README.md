# ![nf-core/circrna](docs/images/nf-core-circrna_logo.png)

**workflow for the quantification, differential expression analysis and miRNA target prediction analysis of circRNAs in RNA-Seq data**.

[![GitHub Actions CI Status](https://github.com/nf-core/circrna/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/circrna/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/circrna/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/circrna/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/circrna.svg)](https://hub.docker.com/r/nfcore/circrna)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23circrna-4A154B?logo=slack)](https://nfcore.slack.com/channels/circrna)

## Introduction

**nf-core/circrna** is a bioinformatics pipeline used for the quantification, miRNA target prediction and differential expression analysis of circRNAs present in RNA sequencing data (currently supporting total RNA-Seq paired end sequencing data, mapped to *H. sapiens* Gencode reference genomes GRCh37, GRCh38 v34).

The pipleline has been developed in a modular fashion, permitting the user to select miRNA target prediction, differential expression analysis (or both) in addition to circRNA quantification to facilitate hypotheses surrounding circRNAs involvement in the competing endogenous RNA network.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline Summary

By default, `nf-core/circrna` utilises all 3 analysis modules: `circrna_discovery, mirna_prediction, differential_expression`.

The creation of reference genome indices and aligners used depends on the circRNA quantification tools selected, given by `--tool`.

### 1. Download Reference Genome Files

- Download Gencode GRCh37/GRCh38 *H. sapiens* reference genome, GTF files.
- Create customised annotation text file ([`gtfToGenePred`](https://anaconda.org/bioconda/ucsc-gtftogenepred))

- Download [`miRbase`](http://www.mirbase.org/ftp.shtml)  mature miRNA sequences.

- Download [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi) miRNA sequences and family conservation information.

### 2. Create Reference Index Files

- [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/) reference genome index.

- [`bwa`](https://sourceforge.net/projects/bio-bwa/files/) reference genome indices (`ciriquant` in `--tool`).

- [`HISAT2`](http://daehwankimlab.github.io/hisat2/download/) reference genome indices (`ciriquant` in `--tool`) .

- [`STAR`](https://github.com/alexdobin/STAR/releases) reference genome indices (`circexplorer2`, `circrna_finder`, `dcc` in `--tool`).

- [`Bowtie`](https://sourceforge.net/projects/bowtie-bio/) reference genome indices (`mapsplice` in `--tool`).

- [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) reference genome indices (`find_circ` in `--tool`).

### 3. Miscellaneous circRNA Tool Requirements

- Split reference genome FASTA file per chromosome (`mapsplice`, `find_circ` in `--tool` ).

- Create `CIRIquant` input `.yml` file (`ciriquant` in `--tool`).

### 4. Stage Input Data

- Accept `BAM` (converted to FASTQ pairs using [`picard`](https://sourceforge.net/projects/picard/)) or FASTQ  input data given as a path or `.csv` input file.

### 5. Quality Control

- Perform [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) , [`MultiQC`](http://multiqc.info/) on raw input data.

- Optional adapter removal + read trimming performed by [`BBDUK`](https://sourceforge.net/projects/bbmap/).

- Perform [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) , [`MultiQC`](http://multiqc.info/) on trimmed reads.

### 6. circRNA Quantification

- Align RNA-Seq reads and perform circRNA quantification ([`ciriquant`](https://sourceforge.net/projects/ciri/files/CIRIquant/), [`circrna_finder`](https://github.com/orzechoj/circRNA_finder), [`circexplorer2`](https://circexplorer2.readthedocs.io/en/latest/tutorial/setup/#installation), [`find_circ`](https://github.com/marvin-jens/find_circ), [`dcc`](https://github.com/dieterich-lab/DCC), [`mapsplice`](https://anaconda.org/bioconda/mapsplice)).
- Raw quantification tool output.
- Filtered quantification tool output as BED6 file.
- Create circRNA count matrix.

### 7. circRNA Filtering

- Filter to remove:
  - circRNAs with low evidence of reads aligned to back splice junction (read count < 2).
  - circRNAs called by only one tool (when 2+ quantification tools selected).

### 8. circRNA Annotation

- Remove unwanted biotypes from Gencode GTF file.

- Calculate circRNA mature spliced length ([`BEDTools`](https://sourceforge.net/projects/bedtools/)).
- Annotate circRNA exon-intron boundaries, mark as `circRNA`, `ciRNA`, `EIciRNA`.

- Identify circRNA parent gene ([`BEDTools`](https://sourceforge.net/projects/bedtools/)).

- Create circRNA FASTA files ([`BEDTools`](https://sourceforge.net/projects/bedtools/)).
- Create circRNA BED12 files ([`BEDTools`](https://sourceforge.net/projects/bedtools/)).
- Create master annotation file reporting circRNA ID, circRNA type, mature length, parent gene, strand.

### 9. miRNA Target Prediction

- Identify miRNA response elements in mature circRNA sequence using [`miRanda`](https://anaconda.org/bioconda/miranda) and [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi).
- Create filtered miRNA targets file for each circRNA.
- Create circos plot of circRNA exons / filtered MRE sites.
- Raw [`miRanda`](https://anaconda.org/bioconda/miranda) output.
- Raw [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi) output.

### 10. miRNA Target Filtering

- Filter to remove:
  - 6mers from `TargetScan` output.
  - miRNAs with MFE <= -20.00 Kcal/Mol.
  - Duplicate miRNA IDs targeting same circRNA MRE site (keep miRNA ID with highest score).

### 11. Differential Expression Analysis

- Perform RNA-Seq quantification using [`StringTie`](https://ccb.jhu.edu/software/stringtie/).

- Perform circRNA differential expression analysis using RNA-Seq library size factors to normalise circRNA count matrix, QC plots, PCA and clustering  ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)).
- Create circRNA expression plots, circRNA-parent gene expression plots ([`R`](https://www.r-project.org/)).
- Create differentially expressed circRNA master file reporting circRNA ID, circRNA Type, mature length, parent gene, strand, Log2FC, pvalue, padj, parent gene description.

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

nf-core/circrna was originally written by [Barry Digby](https://github.com/BarryDigby).

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
