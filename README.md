# MapDbSNP

Tools to add genomic positions to files that contain dbSNP IDs. The pipeline downloads dbSNP from the UCSC genome browser, filters to the required columns, splits the reference for faster lookups, and then maps IDs in parallel.

## Requirements

- R packages: `data.table`, `optparse`, `parallel`, `here`
- Command line tools: `split`, `gzip` (optional: `pigz` for faster decompression)

Install the R dependencies with:

```bash
Rscript -e 'install.packages(c("data.table","optparse","parallel","here"))'
```

## Preparing reference data (recommended)

Download and preprocess dbSNP once, then reuse across runs:

```bash
Rscript ./script/prepare_reference_data.R \
  --build=both \
  --dbsnp-version=155 \
  --data-dir=./data \
  --cpus=8
```

This fetches dbSNP 155 for hg19 and hg38 (~90–100 GB total after splitting) and the RsMerge archive, storing everything under `./data`. Use `--build=hg19` or `--build=hg38` to limit downloads, and `--split-lines` to adjust chunk size.

## Usage

```bash
Rscript ./script/positionsFromDBSNP.r [options]
```

Key options:

- `--input` path to file with dbSNP IDs (e.g., summary statistics)
- `--ID` column name containing dbSNP IDs (default: `ID`)
- `--build` genome build: `hg19` or `hg38`
- `--dbsnp-version` dbSNP release to use (`151` or `155`, default: `155`)
- `--data-dir` directory for reference data (default: `./data`)
- `--outdir` output directory
- `--prefix` prefix for output file name (defaults to input filename)
- `--cpus` CPUs to use for parallel lookups
- `--skip` skip this many lines in the input file
- `--prepare-only` download reference data and exit

## Example

```bash
Rscript ./script/positionsFromDBSNP.r \
  --input=./example/example_input.txt \
  --ID=ID \
  --build=hg19 \
  --dbsnp-version=155 \
  --outdir=./example \
  --prefix=example \
  --data-dir=./data \
  --cpus=16
```
