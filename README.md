# MapDbSNP

Tools to add genomic positions to files that contain dbSNP IDs. The pipeline downloads dbSNP from the UCSC Genome Browser (UCSC/NCBI dbSNP mirrors), filters to the required columns, splits the reference for faster lookups, and then maps IDs in parallel. It supports dbSNP releases 151, 153, and 155 (default: 155).

## Requirements

- R packages: `data.table`, `optparse`, `parallel`, `here`
- Command line tools: `split`, `gzip` (optional: `pigz` for faster decompression), `aria2c` (optional, multi-connection downloads)
- Optional but recommended: UCSC `bigBedNamedItems` utility + dbSNP BigBed (`dbSnp155.bb`) — defaults target to hg38/dbSNP155

## Mamba environment (recommended)

Create and activate the environment once:

```bash
mamba env create -f environment.yml
mamba activate mapdbsnp
```

If you update `environment.yml` later:

```bash
mamba env update -f environment.yml --prune
mamba activate mapdbsnp
```

The env includes `aria2`, `pigz`, R, and required R packages.

Linux/HPC users should also install `bigBedNamedItems` from Bioconda into the same env:

```bash
mamba activate mapdbsnp
mamba install -n mapdbsnp -c bioconda ucsc-bigbednameditems -y
which bigBedNamedItems
```

If you are not using mamba, install the R dependencies with:

```bash
Rscript -e 'install.packages(c("data.table","optparse","parallel","here"))'
```

Install `aria2c` (optional, for faster downloads):

- Conda/Mamba: `conda install -c conda-forge aria2` (or `mamba install -c conda-forge aria2`)
- Homebrew (macOS): `brew install aria2`
- Debian/Ubuntu: `sudo apt-get install -y aria2`
- RHEL/CentOS: `sudo yum install -y aria2`

## BigBed fast path (recommended)

Using the UCSC BigBed file skips the 90–100 GB text download and parallel `awk` scan.
BigBed mode streams the input in chunks and auto-sizes chunk length by input rows and workers (`--chunk-size=0`, default), targeting roughly one chunk per worker for lower overhead and predictable parallel fan-out. Chunk processing is parallelized, defaulting to `--cpus` workers.
By default, non-primary contigs (hap/alt/random-like chromosome names) are excluded from output.
If an rsID maps to multiple primary positions, the first mapping (chromosome order `1-22,X,Y,MT`, then lowest position) is kept in the main output, and the remaining mappings are written to `<prefix>_multiPos_dbSNP<version>_<build>.txt`.
Output coordinates are reported in 1-based genomic coordinates: `POS` is always the 1-based leftmost mapped position (including indels). `POS0` is UCSC 0-based and is only emitted in the multi-position diagnostic file.
Auxiliary outputs are build-specific: `<prefix>_noMatch_dbSNP<version>_<build>.txt` and `<prefix>_multiPos_dbSNP<version>_<build>.txt`.
During BigBed runs, the script prints stage timing and chunk progress (`done/total`, `%`, elapsed time, ETA) plus heartbeat updates while chunks are still running.

1) Ensure `bigBedNamedItems` is available:
- Linux/HPC (recommended with mamba):
  ```bash
  mamba activate mapdbsnp
  mamba install -n mapdbsnp -c bioconda ucsc-bigbednameditems -y
  which bigBedNamedItems
  ```
- macOS (Apple Silicon, manual UCSC binary):  
  `curl -L http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bigBedNamedItems -o ./script/bigBedNamedItems && chmod +x ./script/bigBedNamedItems`
- macOS (Intel, manual UCSC binary):  
  `curl -L http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigBedNamedItems -o ./script/bigBedNamedItems && chmod +x ./script/bigBedNamedItems`
- Linux (x86_64, manual UCSC binary):  
  `curl -L http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedNamedItems -o ./script/bigBedNamedItems && chmod +x ./script/bigBedNamedItems`

If you are on an HPC (or older enterprise Linux) and see errors like `GLIBC_2.29 not found` / `GLIBC_2.33 not found` / `GLIBC_2.34 not found`, use the Bioconda build in the `mapdbsnp` env. The script resolves `bigBedNamedItems` from `PATH` first.

```bash
mamba activate mapdbsnp
mamba install -n mapdbsnp -c bioconda ucsc-bigbednameditems -y
which bigBedNamedItems
file "$(which bigBedNamedItems)"
ldd "$(which bigBedNamedItems)" | head
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
If you do not have `aria2c`, use `curl`:
```bash
curl -L http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb -o ./data/dbSnp155_hg38.bb
```

3) Run with `--bb-file`:
```bash
mamba activate mapdbsnp
Rscript ./script/positionsFromDBSNP.r \
  --input=./example/example_input.txt \
  --ID=ID \
  --build=hg38 \
  --dbsnp-version=155 \
  --cpus=24 \
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
mamba activate mapdbsnp
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
mamba activate mapdbsnp
Rscript ./script/positionsFromDBSNP.r [options]
```

Key options:

- `--input` path to file with dbSNP IDs (e.g., summary statistics)
- `--ID` column name containing dbSNP IDs (default: `ID`)
- `--build` genome build: `hg19` or `hg38`
- `--dbsnp-version` dbSNP release to use (`151`, `153`, or `155`; default: `155`)
- `--bb-file` path to dbSNP BigBed file (if set, text-based lookup is skipped; relative names are also searched under `--data-dir`)
- `--chunk-size` rows per chunk for BigBed streaming mode (`0` = auto from row count and workers; default: `0`)
- `--bb-workers` optional override for BigBed chunk workers (`0` means use `--cpus`; default: `0`)
- `--include-alt-chrom` include non-primary contigs in output (default: off)
- `--no-bb` disable the BigBed fast path and force text lookup
- `--data-dir` directory for reference data (default: `./data`)
- `--outdir` output directory
- `--prefix` prefix for output file name (defaults to input filename)
- `--cpus` CPUs to use for parallel lookups
- `--skip` skip this many lines in the input file
- `--prepare-only` download reference data and exit

## Example (text pipeline)

```bash
mamba activate mapdbsnp
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

## Start-to-finish quickstarts

### macOS (Apple Silicon)
```bash
# Create environment once, then activate it in each new shell
mamba env create -f environment.yml
mamba activate mapdbsnp

# Get UCSC BigBed tool (Apple Silicon) and make executable
curl -L http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bigBedNamedItems -o ./script/bigBedNamedItems
chmod +x ./script/bigBedNamedItems

# Download dbSNP BigBed (hg38/dbSNP155)
aria2c -x8 -s8 -o dbSnp155_hg38.bb http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb
mv dbSnp155_hg38.bb ./data/

# Run example
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

### macOS (Intel)
```bash
# Create environment once, then activate it in each new shell
mamba env create -f environment.yml
mamba activate mapdbsnp

curl -L http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigBedNamedItems -o ./script/bigBedNamedItems
chmod +x ./script/bigBedNamedItems
aria2c -x8 -s8 -o dbSnp155_hg38.bb http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb
mv dbSnp155_hg38.bb ./data/
Rscript ./script/positionsFromDBSNP.r --input=./example/example_input.txt --ID=ID --build=hg38 --dbsnp-version=155 --bb-file=./data/dbSnp155_hg38.bb --outdir=./example --prefix=example_bb --data-dir=./data
```

### Linux (x86_64)
```bash
# Create environment once, then activate it in each new shell
mamba env create -f environment.yml
mamba activate mapdbsnp

# Install BigBed utility in the environment and use it from PATH
mamba install -n mapdbsnp -c bioconda ucsc-bigbednameditems -y
which bigBedNamedItems

# Download dbSNP BigBed (hg38/dbSNP155)
aria2c -x8 -s8 -o dbSnp155_hg38.bb http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb
mv dbSnp155_hg38.bb ./data/

# Run example
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

### Linux HPC quickstart (mamba, avoids common GLIBC issues)
```bash
# Create environment once, then activate it in each new shell
mamba env create -f environment.yml
mamba activate mapdbsnp

# Install BigBed utility in the environment and use it from PATH
mamba install -n mapdbsnp -c bioconda ucsc-bigbednameditems -y
which bigBedNamedItems
file "$(which bigBedNamedItems)"
ldd "$(which bigBedNamedItems)" | head

# Download dbSNP BigBed (hg38/dbSNP155)
aria2c -x8 -s8 -o dbSnp155_hg38.bb http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb
mv dbSnp155_hg38.bb ./data/

# Run example
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

## References

- UCSC Genome Browser downloads (dbSNP tables): https://hgdownload.soe.ucsc.edu/goldenPath/
- UCSC gbdb BigBed sources for dbSNP: https://hgdownload.soe.ucsc.edu/gbdb/
- UCSC bigBedNamedItems utility: http://hgdownload.cse.ucsc.edu/admin/exe/
- NCBI dbSNP RsMerge archive: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/
