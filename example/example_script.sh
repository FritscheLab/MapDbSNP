#!/bin/bash

BB_FILE="./data/dbSnp155_hg38.bb" # Set to your BigBed file to enable fast lookup
BB_ARG=()
if [ -f "$BB_FILE" ]; then
  BB_ARG=( "--bb-file=${BB_FILE}" )
fi

CMD=(
  Rscript ./script/positionsFromDBSNP.r
  --input=./example/example_input.txt
  --ID=ID
  --build=hg19
  --build=hg38
  --dbsnp-version=155
  --outdir=./example
  --prefix=example
  --data-dir=./data
  --cpus=16
)

CMD+=( "${BB_ARG[@]}" )

"${CMD[@]}"
