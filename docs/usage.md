# nf-core/circrna: Usage

It is recommended that first time users run `nf-core/circrna` with the minimal test dataset either locally or on a HPC, referring to the [output documentation](https://nf-co.re/circrna/dev/output) before running a full analysis.

```bash
nextflow run nf-core/circrna -profile test,<docker/singularity/podman/institute>
```

# Running the pipeline

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/circrna --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/circrna -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/circrna
```

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since.

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/circrna releases page](https://github.com/nf-core/circrna/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

# Input specifications

Input data can be passed to `nf-core/circrna` using a CSV file containing the absolute paths to input fastq files.

The headers of the CSV file must be: `sample,fastq_1,fastq_2`.

Valid examples for fastq input data in a CSV file is given below:

| sample           | fastq_1                                                                 | fastq_2                                                                 |
| ---------------- | :---------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| TCGA-EJ-7783-11A | /data/f4c1b2b1-ba1f-4355-a1ac-3e952cf351a5_gdc_realn_rehead_R1.fastq.gz | /data/f4c1b2b1-ba1f-4355-a1ac-3e952cf351a5_gdc_realn_rehead_R2.fastq.gz |
| TCGA-G9-6365-11A | /data/8a36555b-9e27-40ee-a8df-4b15d6580a02_gdc_realn_rehead_R1.fastq.gz | /data/8a36555b-9e27-40ee-a8df-4b15d6580a02_gdc_realn_rehead_R2.fastq.gz |
| TCGA-EJ-7782-11A | /data/8b3d4a3d-2bfa-48f8-b31f-901f49a5bf6b_gdc_realn_rehead_R1.fastq.gz | /data/8b3d4a3d-2bfa-48f8-b31f-901f49a5bf6b_gdc_realn_rehead_R2.fastq.gz |
| TCGA-CH-5772-01A | /data/b6546f66-3c13-4390-9643-d1fb3d660a2f_gdc_realn_rehead_R1.fastq.gz | /data/b6546f66-3c13-4390-9643-d1fb3d660a2f_gdc_realn_rehead_R2.fastq.gz |
| TCGA-EJ-5518-01A | /data/afbbc370-5970-43d3-b9f8-f40f8e649bb6_gdc_realn_rehead_R1.fastq.gz | /data/afbbc370-5970-43d3-b9f8-f40f8e649bb6_gdc_realn_rehead_R2.fastq.gz |
| TCGA-KK-A8I4-01A | /data/81254692-ee1e-4985-bd0a-4929eed4c620_gdc_realn_rehead_R1.fastq.gz | /data/81254692-ee1e-4985-bd0a-4929eed4c620_gdc_realn_rehead_R2.fastq.gz |

> Do not leave any cell empty in the CSV file.

## Phenotype file

When running the differential expression analysis module via the `--module differential_expression` parameter, an input `phenotype.csv` file is required to specify levels for `DESeq2`. At a minimum, the user must supply one column of levels for `DESeq2` which **must be called condition**. This should be the primary contrast of interest in your experiment (e.g case vs. control). If additional columns are supplied to the phenotype file, they will be controlled for in the linear mixed model. A brief proof of concept is given below in R notation:

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

It is recommended to construct your input CSV file in conjunction with your phenotype file as the first column denoting sample names **must match** the first column of the `phenotype.csv` file.

A valid example of a `phenotype.csv` file (matching the TCGA example input CSV file above) is given:

| Sample_ID        | condition |
| ---------------- | --------- |
| TCGA-EJ-7783-11A | control   |
| TCGA-G9-6365-11A | control   |
| TCGA-EJ-7782-11A | control   |
| TCGA-CH-5772-01A | tumor     |
| TCGA-EJ-5518-01A | tumor     |
| TCGA-KK-A8I4-01A | tumor     |

# Analysis modules

`nf-core/circrna` provides 3 analysis modules to the user:

1. circRNA quantification & annotation.
2. miRNA target prediction.
3. Differential circRNA expression analysis.

## circRNA discovery

The core module of `nf-core/circrna`, a user can utilise up to seven circRNA quantification tools to fully characterise the circRNA profile in samples. Currently, supported tools include `CIRCexplorer2`, `circRNA finder`, `CIRIquant`, `DCC`, `find circ` , `MapSplice` & `Segemehl` however, the authors of `nf-core/circrna` welcome contributions from authors of novel quantification tools to keep the workflow current.

By default, `nf-core/circrna` runs the circRNA discovery analysis module.

```bash
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --module 'circrna_discovery'
```

To view the outputs of the module, please see the output [documentation](https://nf-co.re/circrna/dev/output#circrna-quantification).

> Please note that this module must be included for every run of the workflow

### Tool selection

The user may use one, all or any combination of circRNA quantification tools listed above in the analysis. To select which tools to use for the analysis, specify the `--tool` parameter in the configuration profile or pass it via the command line when running the workflow:

```bash
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --module 'circrna_discovery' \
    --tool 'ciriquant,dcc,find_circ'
```

> When providing multiple tools, separate each entry with a comma.

### circRNA filtering

`nf-core/circrna` offers robust filtering of each called circRNA to reduce the number of spurious calls within the dataset.

#### BSJ reads

The user can specify the minimum number of reads spanning the back-splice junction site required for a circRNA to be considered for further analysis. circRNAs with counts below this value will be filtered to remove from the results.

To apply this filtering method, specify the `--bsj_reads` parameter in the configuration profile or pass it via the command line when running the workflow:

```bash
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --phenotype 'phenotype.csv' \
    --module 'circrna_discovery' \
    --tool 'ciriquant, dcc, find_circ' \
    --bsj_reads 2
```

Disable the filter by setting the value to 0.

#### Multiple tool filtering

When more than one tool has been provided using the `--tool` parameter, the user can specify the minimum number of tools circRNAs must be called by using `--tool_filter`. Setting this parameter to 0 or 1 will result in the union being output, i.e no filtering is applied. Setting this parameter to 2 will output circRNAs that have been called by at least 2 quantification tools and so on.

> The integer provided to the parameter must be less than or equal to the number of quantification tools provided to `--tool`.

To apply this filtering method, specify the `--tool_filter` parameter in the configuration profile or pass it via the command line when running the workflow:

```bash
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --module 'circrna_discovery' \
    --tool 'ciriquant, dcc, find_circ' \
    --bsj_reads 2 \
    --tool_filter 2
```

> This filtering method is reflected in the circRNA count matrix. Per tool circRNA annotations are subject to back-splice read filtering only.

#### Handling duplicate circRNAs

In the event a circRNA has been called by more than one quantification tool, the user can specify which aggregate function to apply to the duplicated circRNA. The accepted values are 'mean' and 'max', which are passed to the workflow using the `--duplicates_fun` parameter.

## miRNA prediction

The second module of `nf-core/circrna`, `mirna_prediction` analyses the mature spliced sequences of circRNAs to test for the presence of miRNA response elements using both `miRanda` and `TargetScan`. Results from both tools are consolidated and filtering methods are applied to produce robust miRNA target predictions of circRNAs in the dataset.

To invoke the module, specify the `--module` parameter via the configuration profile or pass it via the command line when running the workflow:

```bash
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --module 'circrna_discovery, mirna_prediction'
```

To view the outputs of the module, please see the output [documentation](https://nf-co.re/circrna/dev/output#mirna-prediction).

## Differential circRNA analysis

The third and final module of `nf-core/circrna` performs differential expression analysis of circRNAs, returning `DESeq2` result outputs, plots and diagnostic plots for the user. In order to run this module, it is essential that your `phenotype.csv` file is in the correct format - please refer to the input [specifications](https://nf-co.re/circrna/dev/usage#differential-expression-analysis).

To invoke the module, specify the `--module` parameter via the configuration profile or pass it via the command line when running the workflow:

```bash
nextflow run nf-core/circrna \
    -profile <docker/singularity/podman/institute> \
    --genome 'GRCh37' \
    --input 'samples.csv' \
    --phenotype 'phenotype.csv' \
    --module 'circrna_discovery, differential_expression'
```

To view the outputs of the module, please see the output [documentation](https://nf-co.re/circrna/dev/output#differential-expression-analysis).

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
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
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

Command error:
.command.sh: line 9: 30 Killed STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
/home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

````

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
````

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
