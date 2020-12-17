########################################################################
# Title: Dockerfile for ngscnv - verion 1.0
# Author: Tomas W Fitzgerald
# Date: 21/07/2020
# Email: tomas@ebi.ac.uk
########################################################################

########################################################################
# Stage 1: Install hts-nim-tools with miniconda
########################################################################
FROM continuumio/miniconda3 AS builder

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -y hts-nim-tools
RUN /opt/conda/bin/conda clean --yes --force-pkgs-dirs

########################################################################
# Stage 2
########################################################################
FROM rstudio/r-base:3.6-bionic

COPY --from=builder /opt/conda /opt/conda

########################################################################
# Install samtools, bcftools, python3 and dependancies
########################################################################
ENV APP_NAME=samtools
ENV VERSION=1.9
ENV APPS=/software/applications
ENV DEST=$APPS/$APP_NAME/
ENV PATH=/opt/conda/bin/:/resources/:$APPS/$APP_NAME/$VERSION/bin:$APPS/bcftools/$VERSION/bin:$PATH
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.8
RUN apt-get update && apt-get install -y ncurses-dev && \
											apt-get install -y bzip2 && \
											apt-get install -y xz-utils && \ 
											apt-get install -y zlib1g-dev

RUN curl -L https://github.com/samtools/htslib/archive/$VERSION.tar.gz | tar xz && \
    curl -L https://github.com/samtools/samtools/archive/$VERSION.tar.gz | tar xz && \
    curl -L https://github.com/samtools/bcftools/archive/$VERSION.tar.gz | tar xz && \
    mv htslib-$VERSION htslib && \ 
    cd bcftools-$VERSION && \
    make -j HTSDIR=../htslib && \
    make prefix=$APPS/bcftools/$VERSION install && \
    cd .. && \
    cd samtools-$VERSION && \
    make -j HTSDIR=../htslib && \
    make prefix=$APPS/$APP_NAME/$VERSION install && \
    cd ../ && \
    rm -rf bcftools-$VERSION

########################################################################
# Install ngscnv - complie and linking
# Note: this package is not yet public so its a local copy and complie
########################################################################

ENV CNV_NAME=$APPS/ngscnv
ENV CNV_APP=$CNV_NAME/ngscnv

COPY ./ngscnv /$CNV_APP
RUN gcc -O3 -Wall -c -I/samtools-1.9 -I/htslib/htslib -rdynamic -o /$CNV_NAME/cnvdata.o -L/htslib/htslib -L/samtools-1.9 /$CNV_APP/cnvdata.c -lhts -lbam -lpthread -lz -lm 
RUN g++ -I /$CNV_APP -c /$CNV_APP/gc.cpp -o /$CNV_NAME/gc.o 
RUN g++ -I /$CNV_APP -c /$CNV_APP/cwavef.cpp -o /$CNV_NAME/cwavef.o
RUN g++ -I /$CNV_APP -c /$CNV_APP/adms.cpp -o /$CNV_NAME/adms.o 
RUN g++ -I /$CNV_APP -c /$CNV_APP/rd.cpp -o /$CNV_NAME/rd.o 
RUN g++ -I /$CNV_APP -c /$CNV_APP/fh.cpp -o /$CNV_NAME/fh.o 
RUN g++ -I /$CNV_APP -c /$CNV_APP/util.cpp -o /$CNV_NAME/util.o 
RUN g++ -I /$CNV_APP -c -I/samtools-1.9 -I/htslib/htslib -I/htslib/htslib /$CNV_APP/readdata.cpp -lhts -lbam -lpthread -lz -lm -o /$CNV_NAME/readdata.o 
RUN g++ -I /$CNV_APP -c -I/samtools-1.9 -I/htslib/htslib /$CNV_APP/ngscnv.cpp -lhts -lbam -lpthread -lz -lm -o /$CNV_NAME/ngscnv.o
RUN g++ -Isrc -I/samtools-1.9 -I/htslib -I/htslib/htslib -rdynamic -o /$CNV_NAME/ngs /$CNV_NAME/ngscnv.o /$CNV_NAME/gc.o /$CNV_NAME/cwavef.o /$CNV_NAME/fh.o /$CNV_NAME/util.o /$CNV_NAME/adms.o /$CNV_NAME/readdata.o /$CNV_NAME/cnvdata.o -L/htslib -L/samtools-1.9 -lhts -lbam -lpthread -lz -lm -llzma -lbz2

########################################################################
# Install helper R packages - Rbin
# Note: this could be pulled directly from github CRAN etc
# prefer to install locally for the moment
########################################################################

COPY src/kernlab_0.9-29.tar.gz /
RUN R CMD INSTALL kernlab_0.9-29.tar.gz

COPY src/MASS_7.3-51.6.tar.gz /
RUN R CMD INSTALL MASS_7.3-51.6.tar.gz

COPY src/segmented_1.2-0.tar.gz /
RUN R CMD INSTALL segmented_1.2-0.tar.gz

COPY src/survival_3.2-3.tar.gz /
RUN R CMD INSTALL survival_3.2-3.tar.gz

COPY src/mixtools_1.2.0.tar.gz /
RUN R CMD INSTALL mixtools_1.2.0.tar.gz

COPY src/ViteRbi_1.0.tar.gz /
RUN R CMD INSTALL ViteRbi_1.0.tar.gz

COPY src/Rbin_1.1.tar.gz /
RUN R CMD INSTALL Rbin_1.1.tar.gz

########################################################################
# Remove package sources
########################################################################

RUN rm -rf /$CNV_APP /Rbin_1.1.tar.gz /ViteRbi_1.0.tar.gz /samtools-$VERSION /htslib /segmented_1.2-0.tar.gz /kernlab_0.9-29.tar.gz /MASS_7.3-51.6.tar.gz /survival_3.2-3.tar.gz /mixtools_1.2.0.tar.gz 

########################################################################
# Setup the mount locations and copy some resources
# Note: resources are very dataset specific
#	- TODO: auto resource generation functionality 
########################################################################

RUN mkdir -p /resources
RUN mkdir -p /input
RUN mkdir -p /output
RUN mkdir -p /ref
ENV REF_PATH=/ref/%2s/%2s/%s

COPY src/run /resources
COPY src/run.R /resources
COPY src/cnest.py /resources
RUN chmod +x /resources/cnest.py

ENTRYPOINT ["python3.8", "/resources/cnest.py"]

########################################################################
# Run Notes:
# This docker image can be excuted with two mounting points on your server - input_location and output_location
#	The wrapper script expects sorted and indexed BAM files as the initial input (for step1)
# There are currently a total of seven steps that should be run independantly and each step should fully complete before excuting the next step
#
# Usage Examples:
#
# Step1 execution
# docker run -v <YOUR INPUT DIR>:/input_location -v  <YOUR OUTPUT DIR>:/output_location  -it cnv:1.0 ./resources/run step1 cnv-project
#
# Step2 execution
# docker run -v <YOUR INPUT DIR>:/input_location -v  <YOUR OUTPUT DIR>:/output_location  -it cnv:1.0 ./resources/run step2 cnv-project <BAMFILE> <SAMPLE_ID>
#
# Step3 execution 
# docker run -v <YOUR INPUT DIR>:/input_location -v  <YOUR OUTPUT DIR>:/output_location  -it cnv:1.0 ./resources/run step3 cnv-project
#
# Step4
# docker run -v <YOUR INPUT DIR>:/input_location -v  <YOUR OUTPUT DIR>:/output_location  -it cnv:1.0 ./resources/run step4 cnv-project <SAMPLE_INDEX>
#
# Step5 - 7 run exactly the as step4 i.e. using the step tag (e.g. step7) and <SAMPLE_INDEX> variable
# e.g. for step7:
# docker run -v <YOUR INPUT DIR>:/input_location -v  <YOUR OUTPUT DIR>:/output_location  -it cnv:1.0 ./resources/run step7 cnv-project <SAMPLE_INDEX>
#
########################################################################