# MapDbSNP

Tools to add genomic positions to files that contain dbSNP IDs. The pipeline downloads dbSNP from the UCSC Genome Browser (UCSC/NCBI dbSNP mirrors), filters to the required columns, splits the reference for faster lookups, and then maps IDs in parallel. It supports dbSNP releases 151, 153, and 155 (default: 155).

## Requirements

- R packages: `data.table`, `optparse`, `parallel`, `here`
- Command line tools: `split`, `gzip` (optional: `pigz` for faster decompression), `aria2c` (optional, multi-connection downloads)
- Optional but recommended: UCSC `bigBedNamedItems` utility + dbSNP BigBed (`dbSnp155.bb`) — defaults target to hg38/dbSNP155

Install the R dependencies with:

```bash
Rscript -e 'install.packages(c("data.table","optparse","parallel","here"))'
```

Install `aria2c` (optional, for faster downloads):

- Conda: `conda install -c conda-forge aria2`
- Homebrew (macOS): `brew install aria2`
- Debian/Ubuntu: `sudo apt-get install -y aria2`
- RHEL/CentOS: `sudo yum install -y aria2`

## BigBed fast path (recommended)

Using the UCSC BigBed file skips the 90–100 GB text download and parallel `awk` scan.

1) Download the utility (place it in `./script/` or your `PATH`):
```bash
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedNamedItems -O ./script/bigBedNamedItems
chmod +x ./script/bigBedNamedItems
```

2) Download the dbSNP BigBed for your build (defaults to `./data/dbSnp<version>_<build>.bb`):
- hg38 (default): `http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb` (or `dbSnp153.bb` / `dbSnp151.bb` if you set `--dbsnp-version`)
- hg19: `http://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp155.bb` (or `dbSnp153.bb` / `dbSnp151.bb`)
Make sure the file you download matches the `--build` you use. The script will auto-download to `./data` if missing, unless `--no-bb` is set.
Tip: if you have bandwidth, `aria2c` can speed this up:
```bash
aria2c -x8 -s8 -o dbSnp155_hg38.bb http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb
mv dbSnp155_hg38.bb ./data/
```

3) Run with `--bb-file`:
```bash
Rscript ./script/positionsFromDBSNP.r \
  --input=./example/example_input.txt \
  --ID=ID \
  --build=hg38 \
  --dbsnp-version=155 \
  --bb-file=./data/dbSnp155_hg38.bb \
  --outdir=./example \
  --prefix=example_bb \
  --data-dir=./data
```
The script will still download `RsMergeArch.bcp` for ID updates if it is not present in `./data`.

## Preparing reference data (text pipeline)

Use this only if you are not using the BigBed fast path. It downloads and splits the full text dumps.

Download and preprocess dbSNP once, then reuse across runs:

```bash
Rscript ./script/prepare_reference_data.R \
  --build=both \
  --dbsnp-version=155 \
  --data-dir=./data \
  --cpus=8
```

This fetches dbSNP 155 for hg19 and hg38 (~90–100 GB total after splitting) and the RsMerge archive, storing everything under `./data`. Use `--build=hg19` or `--build=hg38` to limit downloads, and `--split-lines` to adjust chunk size.
Note: downloads automatically prefer `aria2c` (multi-connection) when available; otherwise they fall back to `download.file`.
Warning: the text/awk path is legacy and slow for large inputs; prefer the BigBed fast path whenever possible.

## Usage

```bash
Rscript ./script/positionsFromDBSNP.r [options]
```

Key options:

- `--input` path to file with dbSNP IDs (e.g., summary statistics)
- `--ID` column name containing dbSNP IDs (default: `ID`)
- `--build` genome build: `hg19` or `hg38`
- `--dbsnp-version` dbSNP release to use (`151`, `153`, or `155`; default: `155`)
- `--bb-file` path to dbSNP BigBed file (if set, text-based lookup is skipped; defaults to a downloaded `./data/dbSnp<version>_<build>.bb` when available)
- `--no-bb` disable the BigBed fast path and force text lookup
- `--data-dir` directory for reference data (default: `./data`)
- `--outdir` output directory
- `--prefix` prefix for output file name (defaults to input filename)
- `--cpus` CPUs to use for parallel lookups
- `--skip` skip this many lines in the input file
- `--prepare-only` download reference data and exit

## Example (text pipeline)

```bash
Rscript ./script/positionsFromDBSNP.r \
  --input=./example/example_input.txt \
  --ID=ID \
  --build=hg38 \
  --dbsnp-version=155 \
  --outdir=./example \
  --prefix=example \
  --data-dir=./data \
  --cpus=16
```

## References

- UCSC Genome Browser downloads (dbSNP tables): https://hgdownload.soe.ucsc.edu/goldenPath/
- UCSC gbdb BigBed sources for dbSNP: https://hgdownload.soe.ucsc.edu/gbdb/
- UCSC bigBedNamedItems utility: http://hgdownload.cse.ucsc.edu/admin/exe/
- NCBI dbSNP RsMerge archive: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/
