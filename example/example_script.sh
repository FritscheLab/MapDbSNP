#!/bin/bash

Rscript ./script/positionsFromDBSNP.r \
--input=./example/example_input.txt \
--ID=ID \
--build=hg19 \
--dbsnp-version=155 \
--outdir=./example \
--prefix=example \
--data-dir=./data \
--cpus=16
