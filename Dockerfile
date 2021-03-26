FROM nfcore/base:1.13.3
LABEL authors="Barry Digby" \
      description="Docker image containing all software requirements for the nf-core/circrna pipeline"

# install main packages:
RUN apt-get update; apt-get clean all;

RUN apt-get install --yes build-essential \
                          gcc-multilib \
                          apt-utils \
			                    curl \
                          perl \
                  			  zip \
                  			  unzip \
                          expat \
                          libexpat-dev

RUN apt-get install --yes cpanminus

RUN apt-get install --yes libxml-libxml-perl \
                  			  libxml-dom-xpath-perl \
                  			  libxml-libxml-simple-perl \
                  			  libxml-dom-perl

RUN cpanm CPAN::Meta Statistics::Lite Bio::TreeIO

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-circrna-1.0.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-circrna-1.0.0 > nf-core-circrna-1.0.0.yml

# comment CIRIquant line that attempts to run os.chmod()
RUN sed -i '126s/^/#/' /opt/conda/envs/nf-core-circrna-1.0.0/lib/python2.7/site-packages/CIRIquant/main.py

#find_circ
WORKDIR /usr/src/app
RUN wget --no-check-certificate http://www.circbase.org/download/find_circ.tar.gz
RUN tar -xvf find_circ.tar.gz
RUN cp *.py /opt/conda/envs/nf-core-circrna-1.0.0/bin

#circRNA_finder
RUN wget --no-check-certificate https://github.com/orzechoj/circRNA_finder/archive/v1.2.tar.gz
RUN tar -xvf v1.2.tar.gz
WORKDIR /usr/src/app/circRNA_finder-1.2
RUN cp *.pl /opt/conda/envs/nf-core-circrna-1.0.0/bin
RUN chmod 777 filterCirc.awk && cp filterCirc.awk /opt/conda/envs/nf-core-circrna-1.0.0/bin

## TargetScan Executables
RUN curl --output ./targetscan_70.zip http://www.targetscan.org/vert_72/vert_72_data_download/targetscan_70.zip
RUN unzip targetscan_70.zip
RUN curl --output ./targetscan_70_BL_PCT.zip http://www.targetscan.org/vert_72/vert_72_data_download/targetscan_70_BL_PCT.zip
RUN unzip targetscan_70_BL_PCT.zip
RUN curl --output ./TargetScan7_context_scores.zip http://www.targetscan.org/vert_72/vert_72_data_download/TargetScan7_context_scores.zip
RUN unzip TargetScan7_context_scores.zip
RUN chmod 777 targetscan_70.pl && mv targetscan_70.pl /opt/conda/envs/nf-core-circrna-1.0.0/bin
RUN chmod 777 TargetScan7_BL_PCT/targetscan_70_BL_bins.pl && mv TargetScan7_BL_PCT/targetscan_70_BL_bins.pl /opt/conda/envs/nf-core-circrna-1.0.0/bin
RUN chmod 777 TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl && mv TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl /opt/conda/envs/nf-core-circrna-1.0.0/bin
RUN chmod 777 TargetScan7_context_scores/targetscan_70_context_scores.pl && mv TargetScan7_context_scores/targetscan_70_context_scores.pl /opt/conda/envs/nf-core-circrna-1.0.0/bin
RUN chmod 777 TargetScan7_context_scores/targetscan_count_8mers.pl && mv TargetScan7_context_scores/targetscan_count_8mers.pl /opt/conda/envs/nf-core-circrna-1.0.0/bin

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
