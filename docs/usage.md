# nf-core/circrna: Usage

It is recommended that first time users run `nf-core/circrna` with the minimal test dataset either locally or on a HPC, referring to the [output documentation](https://nf-co.re/circrna/dev/output) before running a full analysis.

```console
nextflow run nf-core/circrna -profile test
```

Run the test dataset on a HPC:

```console
nextflow run nf-core/circrna -profile test,<docker/singularity/podman/institute>
```

## Running the pipeline

A typical command for running the pipeline is as follows:

```console
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --input_type 'fastq'
```

By default, `nf-core/circrna` runs the circRNA discovery analysis module using `CIRCexplorer2`. The above command will perform circRNA quantification using these tools on ENSEMBL GRCh37 reference annotation files as defined in the iGenomes config.

### Updating the pipeline

To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/circrna
```

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since.

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/circrna releases page](https://github.com/nf-core/circrna/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Input specifications

Input data can be passed to `nf-core/circrna` in two possible ways using the `--input` parameter.

### `--input "<path>"`

The simplest way to pass input data to `nf-core/circrna` is by providing the path to the input data with a suitable wildcard glob pattern:

#### fastq

```console
--input "/data/*_r{1,2}.fastq.gz"
```

##### bam

```console
--input "/data/*.bam"
```

> Beware that providing a path to input data will result in samples being named according to the common tuple key based on the glob pattern supplied. Take this into consideration when designing your phenotype file for differential expression analysis.

### `--input samples.csv`

Alternatively, the user may wish to provide a CSV file containing the absolute paths to input fastq/bam files.

The headers of the CSV file must be: `Sample_ID,Read1,Read2,Bam`.

> This approach is recommended for most real life situations, where in-house sequencing facilities file naming convention requires the user to manually match file names to metadata. The below input files use `TCGA` identifiers as proof of concept.

Valid examples for fastq/bam input data in a CSV file is given below:

| Sample_ID        | Read1                                                                   | Read2                                                                   | Bam |
| ---------------- | :---------------------------------------------------------------------- | ----------------------------------------------------------------------- | --- |
| TCGA-EJ-7783-11A | /data/f4c1b2b1-ba1f-4355-a1ac-3e952cf351a5_gdc_realn_rehead_R1.fastq.gz | /data/f4c1b2b1-ba1f-4355-a1ac-3e952cf351a5_gdc_realn_rehead_R2.fastq.gz | NA  |
| TCGA-G9-6365-11A | /data/8a36555b-9e27-40ee-a8df-4b15d6580a02_gdc_realn_rehead_R1.fastq.gz | /data/8a36555b-9e27-40ee-a8df-4b15d6580a02_gdc_realn_rehead_R2.fastq.gz | NA  |
| TCGA-EJ-7782-11A | /data/8b3d4a3d-2bfa-48f8-b31f-901f49a5bf6b_gdc_realn_rehead_R1.fastq.gz | /data/8b3d4a3d-2bfa-48f8-b31f-901f49a5bf6b_gdc_realn_rehead_R2.fastq.gz | NA  |
| TCGA-CH-5772-01A | /data/b6546f66-3c13-4390-9643-d1fb3d660a2f_gdc_realn_rehead_R1.fastq.gz | /data/b6546f66-3c13-4390-9643-d1fb3d660a2f_gdc_realn_rehead_R2.fastq.gz | NA  |
| TCGA-EJ-5518-01A | /data/afbbc370-5970-43d3-b9f8-f40f8e649bb6_gdc_realn_rehead_R1.fastq.gz | /data/afbbc370-5970-43d3-b9f8-f40f8e649bb6_gdc_realn_rehead_R2.fastq.gz | NA  |
| TCGA-KK-A8I4-01A | /data/81254692-ee1e-4985-bd0a-4929eed4c620_gdc_realn_rehead_R1.fastq.gz | /data/81254692-ee1e-4985-bd0a-4929eed4c620_gdc_realn_rehead_R2.fastq.gz | NA  |

---

| Sample_ID        | Read1 | Read2 | Bam                                                             |
| :--------------- | ----- | ----- | :-------------------------------------------------------------- |
| TCGA-EJ-7783-11A | NA    | NA    | /data/f4c1b2b1-ba1f-4355-a1ac-3e952cf351a5_gdc_realn_rehead.bam |
| TCGA-G9-6365-11A | NA    | NA    | /data/8a36555b-9e27-40ee-a8df-4b15d6580a02_gdc_realn_rehead.bam |
| TCGA-EJ-7782-11A | NA    | NA    | /data/8b3d4a3d-2bfa-48f8-b31f-901f49a5bf6b_gdc_realn_rehead.bam |
| TCGA-CH-5772-01A | NA    | NA    | /data/b6546f66-3c13-4390-9643-d1fb3d660a2f_gdc_realn_rehead.bam |
| TCGA-EJ-5518-01A | NA    | NA    | /data/afbbc370-5970-43d3-b9f8-f40f8e649bb6_gdc_realn_rehead.bam |
| TCGA-KK-A8I4-01A | NA    | NA    | /data/81254692-ee1e-4985-bd0a-4929eed4c620_gdc_realn_rehead.bam |

> Do not leave any cell empty in the CSV file.

### `--phenotype`

When running the differential expression analysis module, an input `phenotype.csv` file is required to specify levels for `DESeq2`. At a minimum, the user must supply one column of levels for `DESeq2` which **must be called condition**. This should be the primary contrast of interest in your experiment (e.g case vs. control). If additional columns are supplied to the phenotype file, they will be controlled for in the linear mixed model. A brief proof of concept is given below in R notation:

```R
colnames(phenotype)
 [1] 'Sample_ID' 'condition'

print(dds$design)
 [1] ' ~ condition'
```

```R
colnames(phenotype)
 [1] 'Sample_ID' 'condition' 'replicates' 'location'

print(dds$design)
 [1] ' ~ location + replicates + condition'
```

It is recommended to use an input CSV file in conjunction with your phenotype file as the `Sample_ID` column **must match** the first column of the `phenotype.csv` file.

A valid example of a `phenotype.csv` file (matching the TCGA example input CSV file above) is given below:

| Sample_ID        | condition |
| ---------------- | --------- |
| TCGA-EJ-7783-11A | control   |
| TCGA-G9-6365-11A | control   |
| TCGA-EJ-7782-11A | control   |
| TCGA-CH-5772-01A | tumor     |
| TCGA-EJ-5518-01A | tumor     |
| TCGA-KK-A8I4-01A | tumor     |

## Analysis modules

`nf-core/circrna` provides 3 analysis modules to the user:

1. circRNA quantification & annotation.
2. miRNA target prediction.
3. Differential circRNA expression analysis.

### circRNA discovery

The core module of `nf-core/circrna`, the user can utilise the most popular circRNA quantification tools to fully characterise the circRNA profile in samples. Currently, supported tools include `CIRCexplorer2`, `circRNA finder`, `CIRIquant`, `DCC`, `find circ` , `MapSplice` & `Segemehl` however, the authors of `nf-core/circrna` welcome contributions from authors of novel quantification tools to keep the workflow current.

By default, `nf-core/circrna` runs the circRNA discovery analysis module.

```console
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --input_type 'fastq' \
    --module 'circrna_discovery'
```

To view the outputs of the module, please see the output [documentation](https://nf-co.re/circrna/dev/output#circrna-quantification).

> Please note that this module must be included for every run of the workflow

#### Tool selection

The user may use one, all or any combination of circRNA quantification tools listed above in the analysis. To select which tools to use for the analysis, specify the `--tool` parameter in the configuration profile or pass it via the command line when running the workflow:

```console
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --input_type 'fastq' \
    --module 'circrna_discovery' \
    --tool 'ciriquant, dcc, find_circ'
```

> When providing multiple tools, separate each entry with a comma.

#### circRNA filtering

`nf-core/circrna` offers robust filtering of each called circRNA to reduce the number of spurious calls within the dataset.

##### BSJ reads

The user can specify the minimum number of reads spanning the back-splice junction site required for a circRNA to be considered for further analysis. circRNAs with counts below this value will be filtered to remove from the results.

To apply this filtering method, specify the `--bsj_reads` parameter in the configuration profile or pass it via the command line when running the workflow:

```console
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --input_type 'fastq' \
    --phenotype 'phenotype.csv' \
    --module 'circrna_discovery' \
    --tool 'ciriquant, dcc, find_circ' \
    --bsj_reads 2
```

Disable the filter by setting the value to 0.

##### Multiple tool filtering

When more than one tool has been provided using the `--tool` parameter, the user can specify the minimum number of tools circRNAs must be called by using `--tool_filter`. Setting this parameter to 0 or 1 will result in the union being output, i.e no filtering is applied. Setting this parameter to 2 will output circRNAs that have been called by at least 2 quantification tools and so on.

> The integer provided to the parameter must be less than or equal to the number of quantification tools provided to `--tool`.

To apply this filtering method, specify the `--tool_filter` parameter in the configuration profile or pass it via the command line when running the workflow:

```console
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --input_type 'fastq' \
    --module 'circrna_discovery' \
    --tool 'ciriquant, dcc, find_circ' \
    --bsj_reads 2 \
    --tool_filter 2
```

> This filtering method is reflected in the circRNA count matrix. Per tool circRNA annotations are subject to back-splice read filtering only.

### miRNA prediction

The second module of `nf-core/circrna`, `mirna_prediction` analyses the mature spliced sequences of circRNAs to test for the presence of miRNA response elements using both `miRanda` and `TargetScan`. Results from both tools are consolidated and filtering methods are applied to produce robust miRNA target predictions of circRNAs in the dataset.

To invoke the module, specify the `--module` parameter via the configuration profile or pass it via the command line when running the workflow:

```console
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --input_type 'fastq' \
    --module 'circrna_discovery, mirna_prediction'
```

To view the outputs of the module, please see the output [documentation](https://nf-co.re/circrna/dev/output#mirna-prediction).

### Differential circRNA analysis

The third and final module of `nf-core/circrna` performs differential expression analysis of circRNAs, returning `DESeq2` result outputs, plots and diagnostic plots for the user. In order to run this module, it is essential that your `phenotype.csv` file is in the correct format - please refer to the input [specifications](https://nf-co.re/circrna/dev/usage#differential-expression-analysis).

To invoke the module, specify the `--module` parameter via the configuration profile or pass it via the command line when running the workflow:

```console
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --input_type 'fastq' \
    --phenotype 'phenotype.csv' \
    --module 'circrna_discovery, differential_expression'
```

To view the outputs of the module, please see the output [documentation](https://nf-co.re/circrna/dev/output#differential-expression-analysis).

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
  - Pulls software from Docker Hub: [`nfcore/circrna`](https://hub.docker.com/r/nfcore/circrna/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  - Pulls software from Docker Hub: [`nfcore/circrna`](https://hub.docker.com/r/nfcore/circrna/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
  - Pulls software from Docker Hub: [`nfcore/circrna`](https://hub.docker.com/r/nfcore/circrna/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
  - Pulls software from Docker Hub: [`nfcore/circrna`](https://hub.docker.com/r/nfcore/circrna/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
  - Pulls software from Docker Hub: [`nfcore/circrna`](https://hub.docker.com/r/nfcore/circrna/)
- `conda`
  - Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  - Pulls most software from [Bioconda](https://bioconda.github.io/)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

To find the exact name of a process you wish to modify the compute resources, check the live-status of a nextflow run displayed on your terminal or check the nextflow error for a line like so: `Error executing process > 'bwa'`. In this case the name to specify in the custom config file is `bwa`.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
