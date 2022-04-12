#!/bin/bash

git clone https://github.com/hahnlab/CAFE.git

cd CAFE

git checkout env_test

autoconf
./configure
make prep
make test
test/runtests

