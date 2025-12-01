suppressPackageStartupMessages({
  library("data.table")
  library("optparse")
  library("parallel")
  library("here")
})

# Extend timeout to 60 min for download
options(timeout = 3600)

source(here::here("script/function.runParallel.r"))
source(here::here("script/reference_data.R"))

option_list <- list(
  make_option("--input", type = "character", default = "", help = "file with summary statistics"),
  make_option("--ID", type = "character", default = "ID", help = "column name with SNP ID"),
  make_option("--build", type = "character", default = "hg19", help = "Genome Build, hg19 or hg38"),
  make_option("--dbsnp-version", type = "character", default = "155", help = "dbSNP release to use (151 or 155)"),
  make_option("--data-dir", type = "character", default = "", help = "Directory for reference data (default: ./data)"),
  make_option("--outdir", type = "character", default = "", help = "Output directory"),
  make_option("--prefix", type = "character", default = "", help = "Prefix for output file name without path"),
  make_option("--cpus", type = "integer", default = 6, help = "CPUs"),
  make_option("--skip", type = "integer", default = 0, help = "Skip lines"),
  make_option("--prepare-only", action = "store_true", default = FALSE, help = "Only download/prepare reference data and exit")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(t(data.frame(opt)))

ID <- opt$ID
input <- opt$input
outdir <- opt$outdir
prefix <- opt$prefix
cpus <- opt$cpus
build <- tolower(opt$build)
version <- as.character(opt$dbsnp_version)
data_dir <- if (nzchar(opt$data_dir)) opt$data_dir else here::here("data")

if (!build %in% SUPPORTED_BUILDS) stop(sprintf("Unsupported build '%s'", build))
if (!version %in% SUPPORTED_DBSNP_VERSIONS) stop(sprintf("Unsupported dbSNP version '%s'", version))

if (isTRUE(opt$prepare_only)) {
  ensure_reference_data(build, version, data_dir, cpus = cpus)
  ensure_rsmerge(data_dir)
  message("Reference data prepared. Exiting because --prepare-only was set.")
  quit(save = "no")
}

if (!nzchar(input)) stop("Please provide an --input file")
if (!file.exists(input)) stop("Input file not found")
if (!nzchar(outdir)) stop("Please provide an --outdir for outputs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(prefix)) prefix <- tools::file_path_sans_ext(basename(input))

setDTthreads(cpus)

reference <- ensure_reference_data(build, version, data_dir, cpus = cpus)
RsMerge <- ensure_rsmerge(data_dir)

if (opt$skip > 0) {
  stripped <- tempfile(fileext = ".txt")
  system(sprintf("tail -n +%s %s > %s", opt$skip + 1, shQuote(input), shQuote(stripped)))
  input <- stripped
}

header <- fread(input, nrow = 2)
rscolumn <- which(names(header) == ID)
if (length(rscolumn) != 1) {
  stop(sprintf("Column '%s' not found in input header", ID))
}

# part 1: update outdated rsIDs
message("Checking / replacing outdated dbSNP IDs")
awk_merge <- paste(
  "awk", sprintf("-v col=%s", rscolumn),
  "-f", shQuote(here::here("script", "RsMerge_awk.txt")),
  shQuote(RsMerge),
  shQuote(input)
)
updated1 <- fread(cmd = awk_merge, sep = "\t", header = TRUE, showProgress = FALSE)

temp1 <- tempfile(fileext = ".txt")
snpids <- as.character(updated1[[ID]])
snpids <- snpids[grep("^rs", snpids)]
writeLines(snpids, temp1)

# part 2: extract positions from dbSNP
message(sprintf("Extracting positions from dbSNP%s (%s)", version, build))
dbsnp <- reference$split_files
outfiles <- file.path(tempdir(), paste0(basename(dbsnp), "_", seq_along(dbsnp), ".txt"))
awk_extract <- paste("awk -f", shQuote(here::here("script", "Extract_SNPs_dbSNP_awk.txt")))
cmdLines <- sprintf(
  "%s %s %s > %s",
  awk_extract,
  shQuote(temp1),
  shQuote(dbsnp),
  shQuote(outfiles)
)

runParallel(cmdLines, min(cpus, 64))

message("Process and combine input / output")
outfiles <- outfiles[file.exists(outfiles) & file.info(outfiles)$size > 0]

if (length(outfiles) > 0) {
  snppos <- rbindlist(lapply(outfiles, fread, header = FALSE, col.names = c("CHROM", "POS0", "POS", ID)))
} else {
  snppos <- data.table::data.table(CHROM = character(), POS0 = integer(), POS = integer(), ID = character())
}

# zero based start positions
setnames(updated1, c("CHROM", "POS0", "POS"), c("CHROM_old", "POS0_old", "POS_old"), skip_absent = TRUE)

if (nrow(snppos) == 0) {
  warning("No matching SNPs found in reference")
}

updated2 <- merge(snppos, updated1, by = ID)
noMatch <- which(!updated1[[ID]] %in% snppos[[ID]])
if (length(noMatch) > 0) {
  fwrite(updated1[noMatch, ], file.path(outdir, sprintf("%s_noMatch_dbSNP%s.txt", prefix, version)), sep = "\t", quote = FALSE)
}

# Use zero based positions for larger indels to match VCF nomenclature
indels <- which(updated2$POS - updated2$POS0 > 1)
if (length(indels) > 0) updated2[indels, POS := POS0]
updated2[, POS0 := NULL]

suppressWarnings(updated2 <- updated2[order(as.numeric(CHROM), POS), ])
fwrite(updated2, file.path(outdir, sprintf("%s_dbSNP%s_%s.txt", prefix, version, build)), sep = "\t", quote = FALSE)
