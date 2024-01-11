#!/usr/bin/bash

git clone https://github.com/Christina-hshi/psirc.git
cd psirc/psirc-quant

# install htslib
cd ext/htslib/
autoheader
autoconf
./configure
make
make install
cd ../..

# Install FAlite.pm
mkdir -p /usr/local/lib/site_perl
cp ../scripts/FAlite.pm /usr/local/lib/site_perl

# make release
mkdir release && cd release
cmake ..
make psirc-quant
make install
