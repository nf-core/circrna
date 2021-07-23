# ![nf-core/circrna](docs/images/nf-core-circrna_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/circrna/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/circrna/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/circrna/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/circrna/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/circrna/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23circrna-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/circrna)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/circrna** is a bioinformatics best-practice analysis pipeline for circRNA quantification, differential expression analysis and miRNA target prediction of RNA-Seq data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/circrna/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Workflow

The diagram below gives an overview of the modules in `nf-core/circrna`:

<p markdown="1" align="center">
<img src="docs/images/workflow.png" alt="workflow" width="100%">
</p>

## Pipeline Summary

1. Download reference genome files ([`Gencode`](https://www.gencodegenes.org/))
2. Download miRNA database files ([`miRbase`](http://www.mirbase.org/ftp.shtml), [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi))
3. Adapter trimming ([`BBDUK`](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/))
4. Read QC ([`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
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
>>>>>>> dev

## Quick Start

1. Install [`Nextflow`](https://nf-co.re/usage/installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/circrna -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```console
    nextflow run nf-core/circrna -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --module 'circrna_discovery, mirna_prediction, differential_expression' --tool 'circexplorer2' --input 'samples.csv' --input_type 'fastq' --phenotype 'phenotype.csv'
    ```

## Documentation

The nf-core/circrna pipeline comes with documentation about the pipeline [usage](https://nf-co.re/circrna/usage), [parameters](https://nf-co.re/circrna/parameters) and [output](https://nf-co.re/circrna/output).

## Credits

`nf-core/circrna` was originally written by Barry Digby ([@BarryDigby](https://github.com/BarryDigby)) from the [National University of Ireland, Galway](http://www.nuigalway.ie/index-internal.html) as a member of Dr. Pilib Ó Broins lab with the financial support of Science Foundation Ireland (Grant number 18/CRT/6214).

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#circrna` channel](https://nfcore.slack.com/channels/circrna) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/circrna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, references of tools and data used in this pipeline are as follows:

* **BBDUK** Bushnell, B. (Unpublished). Download: [https://sourceforge.net/projects/bbmap/](https://sourceforge.net/projects/bbmap/)
* **bedtools** Quinlan, A.R. & Hall, I.M., (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics , 26(6), pp.841–842. Available at: [http://dx.doi.org/10.1093/bioinformatics/btq033](http://dx.doi.org/10.1093/bioinformatics/btq033). Download: [https://github.com/arq5x/bedtools2/releases](https://github.com/arq5x/bedtools2/)
* **Bowite** Langmead, B., Trapnell, C., Pop, M. et al., (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10, R25. Availabe at: [https://doi.org/10.1186/gb-2009-10-3-r25](https://doi.org/10.1186/gb-2009-10-3-r25). Download: [https://sourceforge.net/projects/bowtie-bio/](https://sourceforge.net/projects/bowtie-bio/)
* **Bowtie2**  Langmead, B. & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), p. 357–359. Available at: [10.1038/nmeth.1923](https:/dx.doi.org/10.1038/nmeth.1923). Download: [http://bowtie-bio.sourceforge.net/bowtie2/index.shtml](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* **bwa** Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics , 25(14), 1754–1760. Available at: [https://doi.org/10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324). Download: [http://bio-bwa.sourceforge.net/bwa.shtml](http://bio-bwa.sourceforge.net/bwa.shtml).
* **CIRCexplorer2** Zhang XO, Dong R, Zhang Y, Zhang JL, Luo Z, Zhang J, Chen LL, Yang L. (2016). Diverse alternative back-splicing and alternative splicing landscape of circular RNAs. Genome Res. 2016 Sep;26(9):1277-87. Available at: [https://doi.org/10.1101/gr.202895.115](https://doi.org/10.1101/gr.202895.115). Download: [https://circexplorer2.readthedocs.io/en/latest/tutorial/setup/](https://circexplorer2.readthedocs.io/en/latest/tutorial/setup/)
* **circRNA finder** Westholm, J.O., Lai, E.C., et al. (2016). Genome-wide Analysis of Drosophila Circular RNAs Reveals Their Structural and Sequence Properties and Age-Dependent Neural Accumulation Westholm et al. Cell Reports. Available at: [https://doi.org/10.1016/j.celrep.2014.10.062](https://doi.org/10.1016/j.celrep.2014.10.062). Download: [https://github.com/orzechoj/circRNA_finder](https://github.com/orzechoj/circRNA_finder)
* **CIRIquant** Zhang, J., Chen, S., Yang, J. et al. (2020). Accurate quantification of circular RNAs identifies extensive circular isoform switching events. Nat Commun 11, 90. Available at: [https://doi.org/10.1038/s41467-019-13840-9](https://doi.org/10.1038/s41467-019-13840-9). Download: [https://github.com/bioinfo-biols/CIRIquant](https://github.com/bioinfo-biols/CIRIquant)
* **DCC** Jun Cheng, Franziska Metge, Christoph Dieterich, (2016). Specific identification and quantification of circular RNAs from sequencing data, Bioinformatics, 32(7), 1094–1096. Available at: [https://doi.org/10.1093/bioinformatics/btv656](https://doi.org/10.1093/bioinformatics/btv656). Download: [https://github.com/dieterich-lab/DCC](https://github.com/dieterich-lab/DCC)
* **find circ** Memczak, S., Jens, M., Elefsinioti, A., Torti, F., Krueger, J., Rybak, A., Maier, L., Mackowiak, S. D., Gregersen, L. H., Munschauer, M., Loewer, A., Ziebold, U., Landthaler, M., Kocks, C., le Noble, F., & Rajewsky, N. (2013). Circular RNAs are a large class of animal RNAs with regulatory potency. Nature, 495(7441), 333–338. Available at: [https://doi.org/10.1038/nature11928](https://doi.org/10.1038/nature11928). Download: [https://github.com/marvin-jens/find_circ](https://github.com/marvin-jens/find_circ)
* **HISAT2** Kim, D., Paggi, J.M., Park, C. et al. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). Available at: [https://doi.org/10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4). Download: [http://daehwankimlab.github.io/hisat2/download/](http://daehwankimlab.github.io/hisat2/download/)
* **MapSplice2** Wang, K., Liu J., et al. (2010) MapSplice: Accurate mapping of RNA-seq reads for splice junction discovery, Nucleic Acids Research, 38(18), 178. Avaialable at: [https://doi.org/10.1093/nar/gkq622](https://doi.org/10.1093/nar/gkq622). Download: [http://www.netlab.uky.edu/p/bioinfo/MapSplice2Download](http://www.netlab.uky.edu/p/bioinfo/MapSplice2Download)
* **miRanda** Enright, A.J., John, B., Gaul, U. et al. (2003). MicroRNA targets in Drosophila. Genome Biol 5, R1. Available at: [https://doi.org/10.1186/gb-2003-5-1-r1](https://doi.org/10.1186/gb-2003-5-1-r1). Download: [http://cbio.mskcc.org/miRNA2003/miranda.html](http://cbio.mskcc.org/miRNA2003/miranda.html).
* **Picard** Broad Institute (Unpublished). Download: [http://broadinstitute.github.io/picard/](http://broadinstitute.github.io/picard/)
* **R**: R Core Team (2020). R: A language and environment for statistical computing. R Foundation for
  Statistical Computing, Vienna, Austria. Download: [https://www.R-project.org/](https://www.R-project.org/).
  * **biomaRt** Durinck S, Spellman PT, Birney E, Huber W. (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Nat Protoc. 4(8):1184-91. Available at: [https://doi.org/10.1038/nprot.2009.97](https://doi.org/10.1038/nprot.2009.97). Download: [https://bioconductor.org/packages/release/bioc/html/biomaRt.html](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
  * **circlize** Zuguang Gu, Lei Gu, Roland Eils, Matthias Schlesner, Benedikt Brors (2014). circlize implements and enhances circular visualization in R , Bioinformatics, 30,(19) 2811–2812. Available at: [https://doi.org/10.1093/bioinformatics/btu393](https://doi.org/10.1093/bioinformatics/btu393). Download: [https://cran.r-project.org/web/packages/circlize/index.html](https://cran.r-project.org/web/packages/circlize/index.html)
  * **DESeq2** Love, M.I., Huber, W. & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550. Available at: [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8). Download: [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  * **EnhancedVolcano** Blighe K, Rana S, Lewis M (2020). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. Download: [https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html)
  * **ggplot2** Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, Download: [https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org).
  * **ggpubr** Kassambara A. (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. Download: [https://rpkgs.datanovia.com/ggpubr/](https://rpkgs.datanovia.com/ggpubr/)
  * **ihw** Ignatiadis, N., Klaus, B., Zaugg, J. et al. (2016). Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. Nat Methods 13, 577–580. Available at: [https://doi.org/10.1038/nmeth.3885](https://doi.org/10.1038/nmeth.3885). Download: [https://bioconductor.org/packages/release/bioc/html/IHW.html](https://bioconductor.org/packages/release/bioc/html/IHW.html)
  * **PCAtools** Blighe K, Lun A (2020). PCAtools: PCAtools: Everything Principal Components Analysis. Download: [https://bioconductor.org/packages/release/bioc/html/PCAtools.html](https://bioconductor.org/packages/release/bioc/html/PCAtools.html)
  * **pheatmap** Kolde, R. (2019) Pretty Heatmaps. Download: [https://cran.r-project.org/package=pheatmap](https://cran.r-project.org/package=pheatmap)
  * **pvclust** Suzuki R., Shimodaira H., (2006). Pvclust: an R package for assessing the uncertainty in hierarchical clustering, Bioinformatics, 22(12), 1540–1542. Available at: [https://doi.org/10.1093/bioinformatics/btl117](https://doi.org/10.1093/bioinformatics/btl117). Download: [https://cran.r-project.org/web/packages/pvclust/index.html](https://cran.r-project.org/web/packages/pvclust/index.html)
  * **tximport** Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4. Avaialable at: [https://doi.org/10.12688/f1000research.7563.1](https://doi.org/10.12688/f1000research.7563.1). Download: [http://bioconductor.org/packages/release/bioc/html/tximport.html](http://bioconductor.org/packages/release/bioc/html/tximport.html)
* **SAMtools** Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics , 25(16), 2078–2079. [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352). Download: [http://www.htslib.org/](http://www.htslib.org/)
* **STAR** Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29(1), 15–21. Available at: [https://doi.org/10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635). Download: [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)
* **StringTie** Pertea, M., Pertea, G., Antonescu, C. et al. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol 33, 290–295. Available at: [https://doi.org/10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122). Download: [https://ccb.jhu.edu/software/stringtie/](https://ccb.jhu.edu/software/stringtie/)
* **TargetScan** Agarwal V, Bell GW, Nam JW, Bartel DP. (2015). Predicting effective microRNA target sites in mammalian mRNAs. Elife, 4:e05005. Available at: [https://doi.org/10.7554/elife.05005](https://doi.org/10.7554/elife.05005). Download: [http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi)
* **ViennaRNA** Lorenz, R., Bernhart, S.H., Höner zu Siederdissen, C. et al. (2011). ViennaRNA Package 2.0. Algorithms Mol Biol 6, 26. Available at: [https://doi.org/10.1186/1748-7188-6-26](https://doi.org/10.1186/1748-7188-6-26). Download: [https://www.tbi.univie.ac.at/RNA/#download](https://www.tbi.univie.ac.at/RNA/#download)

## Data References

This repository generated test data using:

* **CIRI_simulator.pl** Gao, Y., Wang, J. & Zhao, F. (2015). CIRI: an efficient and unbiased algorithm for de novo circular RNA identification. Genome Biol 16, 4. Available at: [https://doi.org/10.1186/s13059-014-0571-3](https://doi.org/10.1186/s13059-014-0571-3). Download: [https://sourceforge.net/projects/ciri/](https://sourceforge.net/projects/ciri/)
