
#### `nf-core lint` overall result: Passed :white_check_mark:

Posted for pipeline commit 0f88762

```diff
+| ✅ 151 tests passed       |+
!| ❗  51 tests had warnings |!
```

<details>

### :heavy_exclamation_mark: Test warnings:

* [Test #4](https://nf-co.re/errors#4) - Config `process.container` looks wrong. Should be `nfcore/circrna:dev` but is `barryd237/circrna:dev`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions AWS full test should test full datasets: `./.github/workflows/awsfulltest.yml`
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-apeglm=1.8.0`, `1.12.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-biomart=2.42.0`, `2.46.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-deseq2=1.26.0`, `1.30.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-enhancedvolcano=1.4.0`, `1.8.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-ihw=1.14.0`, `1.18.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-org.hs.eg.db=3.10.0`, `3.12.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-pcatools=1.2.0`, `2.2.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-rhdf5=2.30.0`, `2.34.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bioconductor-tximport=1.14.0`, `1.18.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bowtie=1.2.3`, `1.3.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `bowtie2=2.3.5.1`, `2.4.2` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `ca-certificates=2019.11.28`, `2020.12.5` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `certifi=2019.11.28`, `2020.12.5` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `curl=7.68.0`, `7.71.1` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `hisat2=2.2.0`, `2.2.1` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `htslib=1.10.2`, `1.11` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `perl=5.26.2`, `5.32.0.1` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `picard=2.23.8`, `2.24.1` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `pip=20.0.2`, `21.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `python=2.7.15`, `3.9.1` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `conda-forge::r-circlize=0.4.11`, `0.4.12` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `conda-forge::r-dplyr=1.0.2`, `1.0.3` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `conda-forge::r-ggplot2=3.3.2`, `3.3.3` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `conda-forge::r-gplots=3.1.0`, `3.1.1` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `samtools=1.10`, `1.11` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `star=2.7.6a`, `2.7.7a` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `stringtie=2.1.1`, `2.1.4` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `viennarna=2.4.15`, `2.4.17` available
* [Test #8](https://nf-co.re/errors#8) - Conda dep outdated: `wheel=0.34.2`, `0.36.2` available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 1.2.1, 1.4.0 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 2.6.9, 2.7.2 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 1.16.4, 1.19.5 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 0.15.2, 0.16.0.1 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 5.1.1, 5.4.1 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 0.20.3, 0.24.1 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 1.2.2, 1.6.0 available
* [Test #10](https://nf-co.re/errors#10) - TODO string in `README.md`: _Update the example "typical command" below used to run the pipeline_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `README.md`: _Add a brief overview of what the pipeline does and how it works_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `README.md`: _Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi._
* [Test #10](https://nf-co.re/errors#10) - TODO string in `base.config`: _Check the defaults for all processes_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `base.config`: _Customise requirements for specific processes._
* [Test #10](https://nf-co.re/errors#10) - TODO string in `test_full.config`: _Specify the paths to your full test data ( on nf-core/test-datasets or directly in repositories, e.g. SRA)_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `test_full.config`: _Give any required params for the test so that command line flags are not needed_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `awstest.yml`: _You can customise CI pipeline run tests as required_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `awsfulltest.yml`: _You can customise AWS full pipeline tests as required_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `ci.yml`: _You can customise CI pipeline run tests as required_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `usage.md`: _Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website._
* [Test #10](https://nf-co.re/errors#10) - TODO string in `output.md`: _Write this documentation describing your workflow's output_
* [Test #10](https://nf-co.re/errors#10) - TODO string in `scrape_software_versions.py`: _Add additional regexes for new tools in process get_software_versions_

### :white_check_mark: Tests passed:

* [Test #1](https://nf-co.re/errors#1) - File found: `nextflow.config`
* [Test #1](https://nf-co.re/errors#1) - File found: `nextflow_schema.json`
* [Test #1](https://nf-co.re/errors#1) - File found: `LICENSE` or `LICENSE.md` or `LICENCE` or `LICENCE.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `README.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `CHANGELOG.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `docs/README.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `docs/output.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `docs/usage.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/branch.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/ci.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/linting.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `main.nf`
* [Test #1](https://nf-co.re/errors#1) - File found: `environment.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `Dockerfile`
* [Test #1](https://nf-co.re/errors#1) - File found: `conf/base.config`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/awstest.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/awsfulltest.yml`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `Singularity`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `parameters.settings.json`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `bin/markdown_to_html.r`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `.github/workflows/push_dockerhub.yml`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `.travis.yml`
* [Test #3](https://nf-co.re/errors#3) - Licence check passed
* [Test #2](https://nf-co.re/errors#2) - Dockerfile check passed
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.name`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.nextflowVersion`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.description`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.version`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.homePage`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `timeline.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `trace.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `report.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `dag.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.cpus`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.memory`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.time`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `params.outdir`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `params.input`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.mainScript`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `timeline.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `trace.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `report.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `dag.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.container`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.version`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.nf_required_version`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.container`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.singleEnd`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.igenomesIgnore`
* [Test #4](https://nf-co.re/errors#4) - Config `timeline.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `report.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `trace.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `dag.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `manifest.name` began with `nf-core/`
* [Test #4](https://nf-co.re/errors#4) - Config variable `manifest.homePage` began with https://github.com/nf-core/
* [Test #4](https://nf-co.re/errors#4) - Config `dag.file` ended with `.svg`
* [Test #4](https://nf-co.re/errors#4) - Config variable `manifest.nextflowVersion` started with >= or !>=
* [Test #4](https://nf-co.re/errors#4) - Config `manifest.version` ends in `dev`: `'1.0dev'`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions 'branch' workflow is triggered for PRs to master: `./.github/workflows/branch.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions 'branch' workflow looks good: `./.github/workflows/branch.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions CI is triggered on expected events: `./.github/workflows/ci.yml`
* [Test #5](https://nf-co.re/errors#5) - CI is building the correct docker image: `docker build --no-cache . -t barryd237/circrna:dev`
* [Test #5](https://nf-co.re/errors#5) - CI is pulling the correct docker image: docker pull barryd237/circrna:dev
* [Test #5](https://nf-co.re/errors#5) - CI is tagging docker image correctly: docker tag barryd237/circrna:dev barryd237/circrna:dev
* [Test #5](https://nf-co.re/errors#5) - Continuous integration checks minimum NF version: `./.github/workflows/ci.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions linting workflow is triggered on PR and push: `./.github/workflows/linting.yml`
* [Test #5](https://nf-co.re/errors#5) - Continuous integration runs Markdown lint Tests: `./.github/workflows/linting.yml`
* [Test #5](https://nf-co.re/errors#5) - Continuous integration runs nf-core lint Tests: `./.github/workflows/linting.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions AWS test is triggered on workflow_dispatch: `./.github/workflows/awstest.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions AWS full test is triggered only on published release and workflow_dispatch: `./.github/workflows/awsfulltest.yml`
* [Test #6](https://nf-co.re/errors#6) - README Nextflow minimum version badge matched config. Badge: `20.04.0`, Config: `20.04.0`
* [Test #6](https://nf-co.re/errors#6) - README had a bioconda badge
* [Test #8](https://nf-co.re/errors#8) - Conda environment name was correct (nf-core-circrna-1.0dev)
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bbtools=37.62`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `bbtools=37.62`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-apeglm=1.8.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-biomart=2.42.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-deseq2=1.26.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-enhancedvolcano=1.4.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-ihw=1.14.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-org.hs.eg.db=3.10.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-pcatools=1.2.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-rhdf5=2.30.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bioconductor-tximport=1.14.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bowtie=1.2.3`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bowtie2=2.3.5.1`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bwa=0.7.17`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `bwa=0.7.17`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `bzip2=1.0.8`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `bzip2=1.0.8`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `ca-certificates=2019.11.28`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `certifi=2019.11.28`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `circexplorer2=2.3.8`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `circexplorer2=2.3.8`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `curl=7.68.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `hisat2=2.2.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `htslib=1.10.2`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `_libgcc_mutex=0.1`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `_libgcc_mutex=0.1`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `mapsplice=2.2.1`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `mapsplice=2.2.1`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `miranda=3.3a`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `miranda=3.3a`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `_openmp_mutex=4.5`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `_openmp_mutex=4.5`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `perl=5.26.2`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `picard=2.23.8`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `pip=20.0.2`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `python=2.7.15`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `conda-forge::r-argparser=0.6`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `conda-forge::r-argparser=0.6`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `conda-forge::r-circlize=0.4.11`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `conda-forge::r-dplyr=1.0.2`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `conda-forge::r-ggplot2=3.3.2`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `conda-forge::r-gplots=3.1.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `conda-forge::r-pheatmap=1.0.12`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `conda-forge::r-pheatmap=1.0.12`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `r-pvclust=2.2_0`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `r-pvclust=2.2_0`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `r-rcolorbrewer=1.1_2`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `r-rcolorbrewer=1.1_2`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `samtools=1.10`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `star=2.7.6a`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `stringtie=2.1.1`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `tophat=2.1.1`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `tophat=2.1.1`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `ucsc-genepredtobed=377`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `ucsc-genepredtobed=377`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `ucsc-genepredtogtf=377`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `ucsc-genepredtogtf=377`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `ucsc-gtftogenepred=377`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `ucsc-gtftogenepred=377`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `viennarna=2.4.15`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `wheel=0.34.2`
* [Test #8](https://nf-co.re/errors#8) - Conda dep had pinned version number: `zlib=1.2.11`
* [Test #8](https://nf-co.re/errors#8) - Conda package is the latest available: `zlib=1.2.11`
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: argparse==1.2.1
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: ciriquant==1.1.1
* [Test #8](https://nf-co.re/errors#8) - PyPi package is latest available: 1.1.1
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: numexpr==2.6.9
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: numpy==1.16.4
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: pysam==0.15.2
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: pyyaml==5.1.1
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: scikit-learn==0.20.3
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: scipy==1.2.2
* [Test #9](https://nf-co.re/errors#9) - Found all expected strings in Dockerfile file
* [Test #12](https://nf-co.re/errors#12) - Name adheres to nf-core convention
* [Test #13](https://nf-co.re/errors#13) - Did not find any cookiecutter template strings (106 files)
* [Test #14](https://nf-co.re/errors#14) - Schema lint passed
* [Test #14](https://nf-co.re/errors#14) - Schema title + description lint passed
* [Test #15](https://nf-co.re/errors#15) - Schema matched params returned from nextflow config

### Run details:

* nf-core/tools version 1.12.1
* Run at `2021-01-26 13:57:53`

</details>
