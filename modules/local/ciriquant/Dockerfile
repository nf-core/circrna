FROM mambaorg/micromamba@sha256:a0a4da2d5315f0d5b9d98196f2e9ecf50152c750bed4023bb127baade7e8faca

USER root

RUN micromamba install --yes --name base --channel conda-forge --channel bioconda \
    python=3.14.7 \
    bwa=0.7.19 \
    hisat2=2.2.3 \
    samtools=1.24 \
    stringtie=3.0.3 \
    pysam=0.24.1 \
    numpy=2.5.4 \
    pyyaml=6.0.3 \
    scipy=1.18.1 \
    scikit-learn=1.9.1 \
    bioconductor-edger=4.10.6 \
    r-optparse=1.8.2 \
    r-statmod=1.5.2 \
    openjdk=25.0.2 \
    perl \
    packaging \
    git=2.56.0 \
    pip=26.2.1 \
    && micromamba clean --all --yes

RUN micromamba run --name base pip install --no-deps \
    git+https://github.com/nictru/CIRIquant.git@e4916ca7b3370cef54d76ca162be64792d8c1b16

# edgeR 4.10 removed the old decideTestsDGE alias; decideTests dispatches the
# same DGELRT method. Keep this adaptation explicit alongside the fork pin.
RUN micromamba run --name base python -c \
    'from pathlib import Path; import sysconfig; p = Path(sysconfig.get_path("purelib")) / "libs/CIRI_DE.R"; source = p.read_text(); assert source.count("decideTestsDGE(") == 2; p.write_text(source.replace("decideTestsDGE(", "decideTests("))'

# Nextflow and Apptainer bypass the micromamba activation entrypoint.
ENV PATH="/opt/conda/bin:${PATH}"

USER $MAMBA_USER
