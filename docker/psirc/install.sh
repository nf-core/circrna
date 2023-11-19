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

# make release
mkdir release && cd release
cmake ..
make psirc-quant
make install
