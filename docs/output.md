# nf-core/circrna: Output

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/circrna/output](https://nf-co.re/circrna/output)

## Introduction

This documentation describes the output of `nf-core/circrna` for the test dataset which runs all 3 modules in the workflow: `circRNA discovery` , `miRNA prediction` and `differential expression` analysis of circular RNAs in RNA-Seq data.  

The processes listed below will fall into one of 4 output directories produced by `nf-core/circrna`:

```console
|-- results/
       |-- quality_control
       |-- circrna_discovery
       |-- mirna_prediction
       |-- differential_expression
```

## Pipeline Overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Download files](#download-files)
  * [Reference files](#reference-files) - Download reference files
  * [miRNA databases](#mirna-databases) - Download mature miRNA sequences
* [Quality Control](#quality-control)
  * [BAM to Fastq](#bam-to-fastq) - Convert BAM to fastq
  * [BBDUK](#bbduk) - Adapter trimming, quality and length filtering
  * [DESeq2](#deseq2) - Quality control plots from `DESeq2` analysis  
  * [FastQC](#fastqc) - Raw read QC
  * [MultiQC](#multiqc) - Consolidated fastqc reports
* [Genome index files](#genome-index-files)
* [circRNA quantification](#circrna-quantification)
  * [Miscellaneous requirements](#miscellaneous-requirements) - Generate tool specific requirements
  * [CIRCexplorer2](#circexplorer2) - Annotation of circular RNAs from STAR 2 pass mode Chimeric out files
  * [circRNA finder](#circrna-finder) - Identify circular RNAs from STAR 2 pass mode SAM files
  * [CIRIquant](#ciriquant) - De novo identification of circular RNAs using circular pseudo reference
  * [DCC](#dcc) - Identify circular RNAs from STAR 2 pass mode utilising joint *and* separate read pairs
  * [find circ](#find-circ) - Identify circular RNAs in unmapped anchors from Bowtie2
  * [MapSplice](#mapsplice) - Identify circular RNAs in unmapped anchors from Bowtie
  * [STAR](#star) - Characterise back-splice junctions in all samples using 2 pass mode
* [circRNA annotation](#circrna-annotation)
  * [Annotated circRNAs](#annotated-circrnas) - Basic circRNA information
  * [BED12 files](#bed12-files) - Individual circRNA coordinates in BED12 format
  * [count matrix](#count-matrix) - circRNA counts matrix
  * [fasta](#fasta) - mature spliced circRNA fasta sequence
* [miRNA target prediction](#mirna-target-prediction)
  * [miranda](#miranda) - Raw output from miRanda
  * [targetscan](#targetscan) - Raw output from TargetScan
  * [miRNA targets](#mirna-targets) - Filtered outputs from miRanda, TargetScan for each circRNA
  * [circos plots](#circos-plots) - circos plot of circRNA - miRNA filtered predictions
* [Differential expression analysis](#differential-expression-analysis)
  * [circRNA](#circrna) - Output directory for circRNA DESeq2 analysis
  * [Boxplots](#boxplots) - Boxplots of differentially expressed circRNAs
  * [RNA-Seq](#rna-seq) - Output directory for RNA-Seq DESeq2 analysis

## Download Files

### Reference Files

<details markdown="1">
<summary>Output files</summary>

* `circrna_discovery/reference/`
  * `*.fa`: Gencode reference FASTA file.
  * `*.gtf`: Gencode reference GTF file.
  * `*.txt`: Customised reference text annotation file.

</details>

`nf-core/circrna` has been designed exclusively with [gencode](https://www.gencodegenes.org/) reference files due to their ubiquitous compatibility with circRNA quantification tools. For this reason, ENSEMBL and UCSC reference files are not recommended. The user can specify which genome version to use (`GRCh37/GRCh38`) via the `--genome_version` parameter, see [parameter documentation](https://nf-co.re/circrna/dev/parameters#reference-genome-options) for details.

### miRNA Databases

<details markdown="1">
<summary>Output files</summary>

* `mirna_prediction/assets/`
  * `hsa_mature.fa`: mature *H. sapiens* miRNA sequences in FASTA format for `miRanda` compatibility.
  * `hsa_miR.txt`: mature *H. sapiens* miRNA sequences in tab delimited format with species ID information for `TargetScan` compatibility.  

</details>

Mature miRNA sequences are downloaded from [miRbase](http://www.mirbase.org/ftp.shtml) and [TargetScan](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi) as inputs for  `miRanda` and `targetscan.pl` miRNA prediction tools, respectively.

## Quality Control

### Bam to Fastq

<details markdown="1">
<summary>Output files</summary>

* `quality_control/preprocessing/bamtofastq/`
  * `*_R{1,2}.fq.gz`: Paired end fastq files, generated using `VALIDATION_STRINGENCY=LENIENT`.

</details>

`nf-core/circrna` can accept input BAM files generated from paired end sequencing reads (e.g `TCGA`) by invoking [picard](https://broadinstitute.github.io/picard/) `SamToFastq`, converting BAM files to paired end fastq files.

### BBDUK

<details markdown="1">
<summary>Output files</summary>

* `quality_control/preprocessing/BBDUK/`
  * `*_r{1,2}.trim.fq.gz`: Processed paired end fastq files.

</details>

[BBDUK](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) (DUK - "Decontamination Using Kmers") is capable of performing adapter trimming, quality trimming/filtering and read length filtering (refer to BBDUK [parameter documentation](https://nf-co.re/circrna/dev/parameters#read-trimming--adapter-removal)) for the quality control of sequencing reads. `nf-core/circrna` will automatically output gzipped fastq files from `BBDUK` to minimise data usage.

### DESeq2

<details markdown="1">
<summary>Output files</summary>

* `quality_control/DESeq2_QC`
  * `circRNA/`
    * `DESeq2_condition_PCA.pdf`: PCA plot of PC1 vs. PC2 displaying the highest amount of variation within the response variable `condition`.
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/circRNA/DESeq2_condition_PCA.png" alt="circRNA PCA" width="500">
    </p>

    * `DESeq2_dispersion.pdf`: Plot of re-fitted genes + gene outliers after shrinkage estimation performed by gene-wide maximum likelihood estimates (red curve) & maximum a posteriori estimates of dispersion.
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/circRNA/DESeq2_dispersion.png" alt="circRNA dispersion" width="500">
    </p>

    * `DESeq2_sample_dendogram.pdf`: Dendogram displaying sample distances using [pvclust](https://cran.r-project.org/web/packages/pvclust/index.html).
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/circRNA/DESeq2_sample_dendogram.png" alt="circRNA dendo" width="500">
    </p>

    * `DESeq2_sample_heatmap.pdf`: Heatmap displaying Manhattan distance between samples.
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/circRNA/DESeq2_sample_heatmap.png" alt="circRNA samplehm" width="500">
    </p>

  * `RNA-Seq/`
    * `DESeq2_condition_PCA.pdf`: PCA plot of PC1 vs. PC2 displaying the highest amount of variation within the response variable `condition`.
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/RNA-Seq/DESeq2_condition_PCA.png" alt="circRNA PCA" width="500">
    </p>

    * `DESeq2_dispersion.pdf`: Plot of re-fitted genes + gene outliers after shrinkage estimation performed by gene-wide maximum likelihood estimates (red curve) & maximum a posteriori estimates of dispersion.
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/RNA-Seq/DESeq2_dispersion.png" alt="circRNA dispersion" width="500">
    </p>

    * `DESeq2_sample_dendogram.pdf`: Dendogram displaying sample distances using [pvclust](https://cran.r-project.org/web/packages/pvclust/index.html).
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/RNA-Seq/DESeq2_sample_dendogram.png" alt="circRNA dendo" width="500">
    </p>

    * `DESeq2_sample_heatmap.pdf`: Heatmap displaying Manhattan distance between samples.
    <p markdown="1" align="center">
    <img src="images/output/DESeq2_QC/RNA-Seq/DESeq2_sample_heatmap.png" alt="circRNA samplehm" width="500">
    </p>

</details>

`nf-core/circrna` outputs quality control plots of normalised *log2* expression data from `DESeq2` to assess heterogeneity in the experiment samples. These plots can be useful to assess sample-sample similarity and to identify potential batch effects within the experiment. Plots are generated for both circRNAs and RNA-Seq data when the differential expression analysis module has been selected by the user (see `--module` [documentation](https://nf-co.re/circrna/dev/parameters#pipeline-options)).

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `quality_control/fastqc/`
  * `raw/`
    * `*{html,zip}`: Output files from `fastqc` for unprocessed RNA-Seq data.
  * `trimmed/`
    * `*{html,zip}`: Output files from `fastqc` for RNA-Seq reads processed by `BBDUK`.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads and the per base sequence content (%T/A/G/C). Information about adapter contamination and other over-represented sequences are also displayed.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `quality_control/multiqc/`
  * `Raw_Reads_MultiQC.html`: Summary reports of unprocessed RNA-Seq reads.
  * `Trimmed_Reads_MultiQC.html`: Summary reports of processed RNA-Seq reads.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. `nf-core` outputs HTML reports for sequencing read quality control.

## Genome Index Files

<details markdown="1">
<summary>Output files</summary>

* `circrna_discovery/index/`
  * `bowtie/`: Directory containing `Bowtie` indices.
  * `bowtie2/`: Directory containing `Bowtie2` indices.
  * `bwa/`: Directory containing `BWA` indices.
  * `hisat2/`: Directory containing `HISAT2` indices.
  * `samtools`: Directory containing `SAMtools` index file.
  * `STAR`: Directory containing `STAR` indices.

</details>

`nf-core/circrna` will automatically generate genome index files depending on the circRNA quantification tools selected for the `circrna_discovery` module. In subsequent workflow runs, users should specify the paths to the genome index files generated to save compute resources - please check `nf-core/circrna` [parameter documentation](https://nf-co.re/circrna/dev/parameters#reference-genome-options) for guidance on how to do this.

*N.B:* Index files must be regenerated if the user re-runs the analysis using different genome version.

## circRNA Quantification

### Miscellaneous Requirements

<details markdown="1">
<summary>Output files</summary>

* `circrna_discovery/reference/chromosomes/`
  * `*.fa`: Individual FASTA files per chromosome.
* `circrna_discovery/ciriquant/`
  * `travis.yml`: Example of `.yml` file below which is automatically generated for the user.

  ```bash
  name: ciriquant
  tools:
   bwa:  /opt/conda/envs/nf-core-circrna-1.0dev/bin/bwa
   hisat2:  /opt/conda/envs/nf-core-circrna-1.0dev/bin/hisat2
   stringtie:  /opt/conda/envs/nf-core-circrna-1.0dev/bin/stringtie
   samtools: /opt/conda/envs/nf-core-circrna-1.0dev/bin/samtools
  reference:
   fasta: /data/bdigby/results/circrna_discovery/reference/GRCh37.fa
   gtf: /data/bdigby/results/circrna_discovery/reference/GRCh37.gtf
   bwa_index: /data/bdigby/results/circrna_discovery/index/bwa/GRCh37
   hisat_index: /data/bdigby/results/circrna_discovery/index/hisat2/GRCh37
  ```

</details>

`CIRIquant` requires a `.yaml` file specifying the containerised paths of `BWA`, `SAMtools`, `HISAT2`, `StringTie` executables in addition to the absolute paths for the reference FASTA, GTF files and `BWA`, `HISAT2` genome index files.

`MapSplice` requires the reference FASTA file to be split into individual FASTA files per chromosome.

### CIRCexplorer2

<details markdown="1">
<summary>Output files</summary>

* `circrna_discovery/tool_outputs/circexplorer2/${sample_id}/`
  * `*.STAR.junction.bed`: Intermediate file generated by `CIRCexplorer2 parse` module, identifying STAR fusion junctions for downstream annotation.
  * `*.txt`: Output files generated by `CIRCexplorer2 annotate` module, based on BED 12 format containing circRNA genomic location information, exon cassette composition and an additional 6 columns specifying circRNA annotations.  Full descriptions of the 18 columns can be found in the `CIRCexplorer2` [documentation](https://circexplorer2.readthedocs.io/en/latest/modules/annotate/#output).

* `circrna_discovery/filtered_outputs/circexplorer2/`
  * `*_circexplorer2.bed`: Parsed `CIRCexplorer2` outputs in minimal BED file format. Low confidence circRNAs (BSJ reads < 2) have been removed.

</details>

[CIRCexplorer2](https://circexplorer2.readthedocs.io/en/latest/) uses `*.Chimeric.out.junction` files generated from `STAR` 2 pass mode to extract back-splice junction sites using the `CIRCexplorer2 parse` module. Following this, `CIRCexplorer2 annotate` performs re-alignment of reads to the back-splice junction sites to determine the precise positions of downstream donor and upstream acceptor splice sites. Back-splice junction sites are subsequently updated and annotated using the customised annotation text file.

### circRNA finder

<details markdown="1">
<summary>Output files</summary>

* `circrna_discovery/tool_outputs/circrna_finder/${sample_id}/`
  * `*.Chimeric.out.sorted.{bam,bam.bai}`: (Sorted and indexed) bam file with all chimeric reads identified by STAR. The circRNA junction spanning reads are a subset of these.
  * `*.filteredJunctions.bed`: A bed file with **all** circular junctions found by the pipeline. The score column indicates the number reads spanning each junction.
  * `*.s_filteredJunctions.bed`: A bed file with those junctions in `*.filteredJunctions.bed` that are flanked by GT-AG splice sites. The score column indicates the number reads spanning each junction.
  * `*.s_filteredJunctions_fw.bed`:  A bed file with the same circular junctions as in file (b), but here the score column gives the average number of forward spliced reads at both splice sites around each circular junction.

* `circrna_discovery/filtered_outputs/circrna_finder/`
  * `*_circrna_finder.bed`: Parsed `circrna_finder` outputs in minimal BED file format. Low confidence circRNAs (BSJ reads < 2) have been removed.

</details>

[circRNA finder](https://github.com/orzechoj/circRNA_finder) uses `*.Chimeric.out.sam`, `*.Chimeric.out.junction` & `*.SJ.out.tab` files to identify circular RNAs in RNA-Seq data.

### CIRIquant

<details markdown="1">
<summary>Output files</summary>

* `circrna_discovery/tool_outputs/ciriquant/${sample_id}/`
  * `*.log`: A `CIRIerror.log` file which should be empty, and a `${sample_id}.log` file which contains the output log of `CIRIquant`.
  * `*.bed`: `CIRI2` output file in BED 6 format.
  * `*.gtf`: Output file from `CIRIquant` in GTF format. Full description of the columns available in the `CIRIquant` [documentation](https://ciriquant-cookbook.readthedocs.io/en/latest/quantification.html#output-format).
* `circrna_discovery/tool_outputs/ciriquant/${sample_id}/align/`
  * `*.sorted.{bam, bam.bai}`: (Sorted and indexed) bam file from `HISAT2` alignment of RNA-Seq reads.
* `circrna_discovery/tool_outputs/ciriquant/${sample_id}/circ/`
  * `*.ciri`: `CIRI2` output file.
  * `*_denovo.sorted.{bam, bam.bai}`: (Sorted and indexed) bam file from `BWA` alignment of candidate circular reads to the pseudo reference.
  * `*_index.*.ht2`: `BWA` index files of the pseudo reference.
  * `*_index.fa`: Reference FASTA file of candidate circular reads.
* `circrna_discovery/filtered_outputs/ciriquant/`
  * `*_ciriquant.bed`: Parsed `CIRIquant` outputs in minimal BED file format. Low confidence circRNAs (BSJ reads < 2) have been removed.

</details>

[CIRIquant](https://github.com/Kevinzjy/CIRIquant) operates by aligning RNA-Seq reads using `HISAT2` and [CIRI2](https://sourceforge.net/projects/ciri/files/CIRI2/) to identify putative circRNAs. Next, a pseudo reference index is generated using `bwa index` by concatenating the two full-length sequences of the putative back-splice junction regions. Candidate circular reads are re-aligned against this pseudo reference using `bwa mem`, and back-splice junction reads are determined if they can be linearly and completely aligned to the putative back-splice junction regions.

### DCC

<details markdown="1">
<summary>Output files</summary>

* `/circrna_discovery/tool_outputs/dcc/${sample_id}/`
  * `*CircCoordinates`: Circular RNA annotations in BED format. Full description of the columns are available in the `DCC` [documentation](https://github.com/dieterich-lab/DCC#output-files-generated-by-dcc).
  * `*CircRNACount`: A table containing read counts for circRNAs detected.
  * `mate1/`: Output directory of STAR 2nd pass alignment for R1.
  * `mate2/`: Output directory of STAR 2nd pass alignment for R2.
* `circrna_discovery/filtered_outputs/dcc/`
  * `*_dcc.bed`: Parsed `DCC` outputs in minimal BED file format. Low confidence circRNAs (BSJ reads < 2) have been removed.

</details>

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - Read quality control
* [MultiQC](#multiqc) - Aggregate report describing results from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files:**

* `fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
* `fastqc/zips/`
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.
