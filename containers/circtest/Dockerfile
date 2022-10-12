FROM nfcore/base:1.13.3
LABEL authors="Barry Digby" \
    description="python 3 container for nf-core/circrna"

WORKDIR ./
COPY environment.yml ./
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/circtest/bin:$PATH

RUN R -e "install.packages('remotes', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "remotes::install_github('dieterich-lab/CircTest')"
