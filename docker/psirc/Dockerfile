FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    git \
    apt-utils \
    autoconf \
    libcurl4-openssl-dev \
    pkg-config \
    libssl-dev \
    make \
    cmake \
    libhdf5-serial-dev \
    libbz2-dev \
    liblzma-dev \
    g++ \
    && rm -rf /var/lib/apt/lists/*

COPY install.sh /tmp/install.sh

WORKDIR /tmp

RUN bash install.sh