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
  make_option("--cpus", type = "integer", default = 4, help = "CPUs"),
  make_option("--chunk-size", type = "integer", default = 1000000, help = "Rows per chunk for BigBed streaming mode (default: 1,000,000)"),
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
      if (length(val) > 0) {
        return(val)
      }
    }
  }
  default
}

ID <- opt$ID
input <- opt$input
outdir <- opt$outdir
prefix <- opt$prefix
cpus <- opt$cpus
chunk_size <- as.integer(get_opt(opt, "chunk.size", 1000000))
build <- tolower(opt$build)
version <- as.character(get_opt(opt, "dbsnp.version", "155"))
data_dir <- get_opt(opt, "data.dir", here::here("data"))
skip <- as.integer(get_opt(opt, "skip", 0))

if (!build %in% SUPPORTED_BUILDS) stop(sprintf("Unsupported build '%s'", build))
if (!version %in% SUPPORTED_DBSNP_VERSIONS) stop(sprintf("Unsupported dbSNP version '%s'", version))
if (is.na(chunk_size) || chunk_size < 1L) stop("--chunk-size must be >= 1")

ensure_dir(data_dir)

find_bigbed_tool <- function() {
  sys <- Sys.which("bigBedNamedItems")
  if (nzchar(sys)) {
    return(sys)
  }
  candidate <- here::here("script", "bigBedNamedItems")
  if (file.exists(candidate) && file.access(candidate, 1) == 0) {
    return(candidate)
  }
  stop("bigBedNamedItems not found in PATH or ./script/. Install it (see README).")
}

renamed_input_colnames <- function(cols) {
  out <- cols
  out[out == "CHROM"] <- "CHROM_old"
  out[out == "POS0"] <- "POS0_old"
  out[out == "POS"] <- "POS_old"
  out
}

sort_output_by_coords <- function(infile, outfile, work_dir) {
  hdr <- names(fread(infile, nrows = 0, showProgress = FALSE))
  chrom_col <- match("CHROM", hdr)
  pos_col <- match("POS", hdr)
  if (is.na(chrom_col) || is.na(pos_col)) stop("Cannot sort output: missing CHROM or POS columns")

  awk_script <- file.path(work_dir, "sort_mapdbsnp.awk")
  writeLines(c(
    "BEGIN {FS=OFS=\"\\t\"}",
    "NR==1 {next}",
    "{chr=$c; pos=$p; key=(chr ~ /^[0-9]+$/ ? chr+0 : 999999); print key, pos, $0}"
  ), awk_script)

  sorted_body <- file.path(work_dir, "sorted_body.tsv")
  sort_cmd <- sprintf(
    "awk -v c=%d -v p=%d -f %s %s | LC_ALL=C sort -t$'\\t' -k1,1n -k2,2n | cut -f3- > %s",
    chrom_col,
    pos_col,
    shQuote(awk_script),
    shQuote(infile),
    shQuote(sorted_body)
  )
  status <- system(sort_cmd)
  if (!identical(status, 0L)) stop("Sorting mapped output failed")

  writeLines(paste(hdr, collapse = "\t"), outfile)
  if (file.exists(sorted_body) && file.info(sorted_body)$size > 0) {
    status <- system(sprintf("cat %s >> %s", shQuote(sorted_body), shQuote(outfile)))
    if (!identical(status, 0L)) stop("Writing sorted mapped output failed")
  }
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

work_dir <- tempfile(pattern = "mapdbsnp_", tmpdir = tempdir())
ensure_dir(work_dir)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

output_main <- file.path(outdir, sprintf("%s_dbSNP%s_%s.txt", prefix, version, build))
output_nomatch <- file.path(outdir, sprintf("%s_noMatch_dbSNP%s.txt", prefix, version))

if (nzchar(bb_file)) {
  message(sprintf("Extracting positions from BigBed in chunks of %s rows: %s", format(chunk_size, big.mark = ","), bb_file))
  updated_file <- file.path(work_dir, "updated_input.tsv")
  status <- system(sprintf("%s > %s", awk_merge, shQuote(updated_file)))
  if (!identical(status, 0L)) stop("Failed to update rsIDs before BigBed lookup")

  input_cols <- names(fread(updated_file, nrows = 0, showProgress = FALSE))
  if (!ID %in% input_cols) stop(sprintf("Column '%s' not found after rsID update", ID))
  renamed_cols <- renamed_input_colnames(input_cols)

  split_prefix <- file.path(work_dir, "updated_chunk_")
  split_cmd <- sprintf(
    "tail -n +2 %s | split --suffix-length=5 --numeric-suffixes --lines=%s - %s",
    shQuote(updated_file),
    chunk_size,
    shQuote(split_prefix)
  )
  status <- system(split_cmd)
  if (!identical(status, 0L)) stop("Splitting updated input into chunks failed")
  chunk_files <- sort(list.files(work_dir, pattern = "^updated_chunk_[0-9]+$", full.names = TRUE))

  bb_tool <- find_bigbed_tool()
  matched_stream_file <- file.path(work_dir, "matched_unsorted.tsv")
  nomatch_stream_file <- file.path(work_dir, "nomatch.tsv")
  first_match <- TRUE
  first_nomatch <- TRUE
  total_chunks <- length(chunk_files)

  if (total_chunks == 0L) {
    out_cols <- c(ID, "CHROM", "POS", setdiff(renamed_cols, ID))
    empty <- as.data.table(setNames(replicate(length(out_cols), logical(0), simplify = FALSE), out_cols))
    fwrite(empty, output_main, sep = "\t", quote = FALSE)
  } else {
    for (i in seq_along(chunk_files)) {
      chunk <- fread(chunk_files[[i]], header = FALSE, col.names = input_cols, showProgress = FALSE)
      if (nrow(chunk) == 0L) next

      setnames(chunk, c("CHROM", "POS0", "POS"), c("CHROM_old", "POS0_old", "POS_old"), skip_absent = TRUE)
      chunk_ids <- as.character(chunk[[ID]])
      query_ids <- unique(chunk_ids[grep("^rs", chunk_ids)])
      snppos <- data.table(CHROM = character(), POS0 = integer(), POS = integer(), ID = character())
      matched <- NULL

      if (length(query_ids) > 0L) {
        ids_file <- file.path(work_dir, sprintf("ids_%05d.txt", i))
        bb_out <- file.path(work_dir, sprintf("bb_hits_%05d.bed", i))
        writeLines(query_ids, ids_file)

        cmd_bb <- sprintf(
          "%s -nameFile %s %s %s",
          shQuote(bb_tool),
          shQuote(bb_file),
          shQuote(ids_file),
          shQuote(bb_out)
        )
        status <- system(cmd_bb)
        if (!identical(status, 0L)) stop(sprintf("bigBedNamedItems failed on chunk %d", i))

        if (file.exists(bb_out) && file.info(bb_out)$size > 0) {
          snppos <- fread(bb_out, select = 1:4, header = FALSE, col.names = c("CHROM", "POS0", "POS", ID), showProgress = FALSE)
          snppos[, CHROM := gsub("^chr", "", CHROM)]
        }
      }

      if (nrow(snppos) > 0L) {
        matched <- merge(snppos, chunk, by = ID)
        indels <- which(matched$POS - matched$POS0 > 1)
        if (length(indels) > 0) matched[indels, POS := POS0]
        matched[, POS0 := NULL]
        fwrite(matched, matched_stream_file, sep = "\t", quote = FALSE, append = !first_match, col.names = first_match)
        first_match <- FALSE

        matched_ids <- unique(as.character(snppos[[ID]]))
        misses <- chunk[!chunk_ids %chin% matched_ids]
      } else {
        misses <- chunk
      }

      if (nrow(misses) > 0L) {
        fwrite(misses, nomatch_stream_file, sep = "\t", quote = FALSE, append = !first_nomatch, col.names = first_nomatch)
        first_nomatch <- FALSE
      }

      if (i %% 5L == 0L || i == total_chunks) {
        message(sprintf("Processed chunk %d/%d", i, total_chunks))
      }

      rm(chunk, chunk_ids, query_ids, snppos, matched, misses)
      gc(verbose = FALSE)
    }

    if (file.exists(matched_stream_file) && file.info(matched_stream_file)$size > 0) {
      sort_output_by_coords(matched_stream_file, output_main, work_dir)
    } else {
      out_cols <- c(ID, "CHROM", "POS", setdiff(renamed_cols, ID))
      empty <- as.data.table(setNames(replicate(length(out_cols), logical(0), simplify = FALSE), out_cols))
      fwrite(empty, output_main, sep = "\t", quote = FALSE)
      warning("No matching SNPs found in reference")
    }
  }

  if (file.exists(nomatch_stream_file) && file.info(nomatch_stream_file)$size > 0) {
    file.copy(nomatch_stream_file, output_nomatch, overwrite = TRUE)
  }
} else {
  updated1 <- fread(cmd = awk_merge, sep = "\t", header = TRUE, showProgress = FALSE)
  temp1 <- tempfile(fileext = ".txt", tmpdir = work_dir)
  snpids <- as.character(updated1[[ID]])
  snpids <- snpids[grep("^rs", snpids)]
  writeLines(snpids, temp1)

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
  snppos <- data.table(CHROM = character(), POS0 = integer(), POS = integer(), ID = character())
  if (length(outfiles) > 0) {
    snppos <- rbindlist(lapply(outfiles, fread, header = FALSE, col.names = c("CHROM", "POS0", "POS", ID)))
  }

  # zero based start positions
  setnames(updated1, c("CHROM", "POS0", "POS"), c("CHROM_old", "POS0_old", "POS_old"), skip_absent = TRUE)

  if (nrow(snppos) == 0) {
    warning("No matching SNPs found in reference")
  }

  updated2 <- merge(snppos, updated1, by = ID)
  noMatch <- which(!updated1[[ID]] %in% snppos[[ID]])
  if (length(noMatch) > 0) {
    fwrite(updated1[noMatch, ], output_nomatch, sep = "\t", quote = FALSE)
  }

  # Use zero based positions for larger indels to match VCF nomenclature
  indels <- which(updated2$POS - updated2$POS0 > 1)
  if (length(indels) > 0) updated2[indels, POS := POS0]
  updated2[, POS0 := NULL]

  suppressWarnings(updated2 <- updated2[order(as.numeric(CHROM), POS), ])
  fwrite(updated2, output_main, sep = "\t", quote = FALSE)
}
