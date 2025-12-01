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
  make_option("--build", type = "character", default = "hg38", help = "Genome Build, hg19 or hg38"),
  make_option("--dbsnp-version", type = "character", default = "155", help = "dbSNP release to use (151 or 155)"),
  make_option("--bb-file", type = "character", default = "", help = "Path to dbSNP BigBed file (dbSnp155.bb). Overrides text-based lookup when provided."),
  make_option("--no-bb", action = "store_true", default = FALSE, help = "Disable BigBed fast path and force text-based lookup"),
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

get_opt <- function(x, name, default = NULL) {
  candidates <- unique(c(name, gsub("\\.", "_", name), gsub("_", ".", name)))
  for (n in candidates) {
    if (!is.null(x[[n]])) {
      val <- x[[n]]
      if (length(val) > 0) return(val)
    }
  }
  default
}

ID <- opt$ID
input <- opt$input
outdir <- opt$outdir
prefix <- opt$prefix
cpus <- opt$cpus
build <- tolower(opt$build)
version <- as.character(get_opt(opt, "dbsnp.version", "155"))
data_dir <- get_opt(opt, "data.dir", here::here("data"))
skip <- as.integer(get_opt(opt, "skip", 0))

if (!build %in% SUPPORTED_BUILDS) stop(sprintf("Unsupported build '%s'", build))
if (!version %in% SUPPORTED_DBSNP_VERSIONS) stop(sprintf("Unsupported dbSNP version '%s'", version))

ensure_dir(data_dir)

find_bigbed_tool <- function() {
  candidate <- here::here("script", "bigBedNamedItems")
  if (file.exists(candidate) && file.access(candidate, 1) == 0) return(candidate)
  sys <- Sys.which("bigBedNamedItems")
  if (nzchar(sys)) return(sys)
  stop("bigBedNamedItems not found. Download it from UCSC (see README) or place it in ./script/.")
}

bb_file <- ""
no_bb <- isTRUE(get_opt(opt, "no.bb", FALSE))
if (!no_bb) {
  supplied_bb <- get_opt(opt, "bb.file", "")
  if (nzchar(supplied_bb)) {
    if (!file.exists(supplied_bb)) stop(sprintf("BigBed file not found: %s", supplied_bb))
    bb_file <- supplied_bb
  } else {
    bb_file <- tryCatch(ensure_bigbed(build, version, data_dir), error = function(e) {
      message("BigBed unavailable, falling back to text-based lookup: ", e$message)
      ""
    })
  }
}

prepare_only <- isTRUE(get_opt(opt, "prepare.only", FALSE))

if (prepare_only) {
  if (!nzchar(bb_file)) {
    ensure_reference_data(build, version, data_dir, cpus = cpus)
  }
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

reference <- if (!nzchar(bb_file)) ensure_reference_data(build, version, data_dir, cpus = cpus) else NULL
RsMerge <- ensure_rsmerge(data_dir)

if (skip > 0) {
  stripped <- tempfile(fileext = ".txt")
  system(sprintf("tail -n +%s %s > %s", skip + 1, shQuote(input), shQuote(stripped)))
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

work_dir <- tempfile(pattern = "mapdbsnp_", tmpdir = tempdir())
ensure_dir(work_dir)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

temp1 <- tempfile(fileext = ".txt", tmpdir = work_dir)
snpids <- as.character(updated1[[ID]])
snpids <- snpids[grep("^rs", snpids)]
writeLines(snpids, temp1)

snppos <- data.table::data.table(CHROM = character(), POS0 = integer(), POS = integer(), ID = character())

if (nzchar(bb_file)) {
  message(sprintf("Extracting positions from BigBed: %s", bb_file))
  bb_tool <- find_bigbed_tool()
  bb_out <- file.path(work_dir, "dbsnp_bb_hits.bed")
  cmd_bb <- sprintf("%s -nameFile %s %s %s",
                    shQuote(bb_tool),
                    shQuote(bb_file),
                    shQuote(temp1),
                    shQuote(bb_out))
  status <- system(cmd_bb)
  if (!identical(status, 0L)) stop("bigBedNamedItems failed. Ensure the binary is installed and executable.")
  if (file.exists(bb_out) && file.info(bb_out)$size > 0) {
    snppos <- fread(bb_out, select = 1:4, header = FALSE, col.names = c("CHROM", "POS0", "POS", ID))
    snppos[, CHROM := gsub("^chr", "", CHROM)]
  }
} else {
  # part 2: extract positions from dbSNP text files
  message(sprintf("Extracting positions from dbSNP%s (%s) via legacy text pipeline (slow). Prefer BigBed for speed.", version, build))
  dbsnp <- reference$split_files
  outfiles <- file.path(work_dir, paste0(basename(dbsnp), "_", seq_along(dbsnp), ".txt"))
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
  }
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
