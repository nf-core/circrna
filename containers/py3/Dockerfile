FROM nfcore/base:1.13.3
LABEL authors="Barry Digby" \
    description="python 3 container for nf-core/circrna"

WORKDIR ./
COPY environment.yml ./
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/py3/bin:$PATH

#DCC
WORKDIR /usr/src/app
RUN wget --no-check-certificate https://github.com/dieterich-lab/DCC/archive/v0.5.0.tar.gz
RUN tar -xvf v0.5.0.tar.gz
WORKDIR /usr/src/app/DCC-0.5.0
# remove --user or else scripts installed to /root/
RUN python setup.py install
