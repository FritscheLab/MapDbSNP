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
  make_option("--chunk-size", type = "integer", default = 0, help = "Rows per chunk for BigBed streaming mode (0 = auto by input size and workers)"),
  make_option("--bb-workers", type = "integer", default = 0, help = "Parallel workers for BigBed chunk processing (0 = use --cpus)"),
  make_option("--include-alt-chrom", action = "store_true", default = FALSE, help = "Include non-primary contigs (hap/alt/random). Default excludes them."),
  make_option("--skip", type = "integer", default = 0, help = "Skip lines"),
  make_option("--prepare-only", action = "store_true", default = FALSE, help = "Only download/prepare reference data and exit")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(t(data.frame(opt)))

get_opt <- function(x, name, default = NULL) {
  candidates <- unique(c(
    name,
    gsub("\\.", "_", name),
    gsub("_", ".", name),
    gsub("[._]", "-", name),
    gsub("-", ".", name),
    gsub("-", "_", name)
  ))
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
bb_workers <- as.integer(get_opt(opt, "bb.workers", 0))
build <- tolower(opt$build)
version <- as.character(get_opt(opt, "dbsnp.version", "155"))
data_dir <- get_opt(opt, "data.dir", "")
if (!nzchar(data_dir)) data_dir <- here::here("data")
include_alt_chrom <- isTRUE(get_opt(opt, "include.alt.chrom", FALSE))
skip <- as.integer(get_opt(opt, "skip", 0))

if (!build %in% SUPPORTED_BUILDS) stop(sprintf("Unsupported build '%s'", build))
if (!version %in% SUPPORTED_DBSNP_VERSIONS) stop(sprintf("Unsupported dbSNP version '%s'", version))
if (is.na(chunk_size) || chunk_size < 0L) stop("--chunk-size must be >= 0")
if (is.na(bb_workers) || bb_workers < 0L) stop("--bb-workers must be >= 0")
if (bb_workers == 0L) bb_workers <- max(1L, as.integer(cpus))

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

format_seconds <- function(seconds) {
  sec <- max(0, as.numeric(seconds))
  hh <- floor(sec / 3600)
  mm <- floor((sec %% 3600) / 60)
  ss <- floor(sec %% 60)
  sprintf("%02d:%02d:%02d", hh, mm, ss)
}

progress_line <- function(done, total, start_time, label = "BigBed progress") {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  pct <- if (total > 0L) 100 * done / total else 100
  rate <- if (elapsed > 0) done / elapsed else 0
  remaining <- if (rate > 0) (total - done) / rate else NA_real_
  eta_txt <- if (is.na(remaining) || !is.finite(remaining)) "--:--:--" else format_seconds(remaining)
  sprintf(
    "%s: %s/%s chunks (%.1f%%) elapsed=%s ETA=%s",
    label,
    format(done, big.mark = ","),
    format(total, big.mark = ","),
    pct,
    format_seconds(elapsed),
    eta_txt
  )
}

sort_output_by_coords <- function(infile, outfile, work_dir) {
  t0 <- Sys.time()
  message("Sorting final mapped output by CHROM/POS...")
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
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("Finished sorting mapped output in %s", format_seconds(elapsed)))
}

count_file_lines <- function(path) {
  out <- system2("wc", args = c("-l", path), stdout = TRUE, stderr = TRUE)
  if (length(out) < 1L) stop(sprintf("Failed to count lines for %s", path))
  toks <- strsplit(trimws(out[1]), "\\s+")[[1]]
  n <- suppressWarnings(as.integer(toks[1]))
  if (is.na(n)) stop(sprintf("Could not parse line count for %s", path))
  n
}

is_primary_chrom <- function(chr) {
  x <- toupper(as.character(chr))
  grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", x)
}

chrom_rank <- function(chr) {
  x <- toupper(as.character(chr))
  rank <- rep.int(999999L, length(x))
  n <- suppressWarnings(as.integer(x))
  is_auto <- !is.na(n) & n >= 1L & n <= 22L
  rank[is_auto] <- n[is_auto]
  rank[x == "X"] <- 23L
  rank[x == "Y"] <- 24L
  rank[x %in% c("M", "MT")] <- 25L
  rank
}

split_multi_position_ids <- function(snppos, id_col) {
  if (nrow(snppos) == 0L) {
    return(list(primary = snppos, duplicates = snppos))
  }

  dt <- copy(snppos)
  dt[, CHROM_RANK := chrom_rank(CHROM)]
  setorderv(dt, c(id_col, "CHROM_RANK", "CHROM", "POS", "POS0"))
  is_dup <- duplicated(dt[[id_col]])
  dup_cols <- c(id_col, "CHROM", "POS0", "POS")
  duplicates <- dt[is_dup, ..dup_cols]
  primary <- dt[!is_dup][, CHROM_RANK := NULL]
  list(primary = primary, duplicates = duplicates)
}

bb_file <- ""
no_bb <- isTRUE(get_opt(opt, "no.bb", FALSE))
if (!no_bb) {
  supplied_bb <- get_opt(opt, "bb.file", "")
  if (nzchar(supplied_bb)) {
    bb_candidates <- unique(c(supplied_bb, file.path(data_dir, supplied_bb)))
    bb_existing <- bb_candidates[file.exists(bb_candidates)]
    if (length(bb_existing) == 0L) {
      stop(sprintf("BigBed file not found. Tried: %s", paste(bb_candidates, collapse = ", ")))
    }
    bb_file <- bb_existing[1]
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
output_multipos <- file.path(outdir, sprintf("%s_multiPos_dbSNP%s.txt", prefix, version))

if (nzchar(bb_file)) {
  message(sprintf("Extracting positions from BigBed: %s", bb_file))
  updated_file <- file.path(work_dir, "updated_input.tsv")
  message("Preparing rsID-updated temporary input...")
  t_update <- Sys.time()
  status <- system(sprintf("%s > %s", awk_merge, shQuote(updated_file)))
  if (!identical(status, 0L)) stop("Failed to update rsIDs before BigBed lookup")
  message(sprintf("Finished rsID update in %s", format_seconds(as.numeric(difftime(Sys.time(), t_update, units = "secs")))))

  input_cols <- names(fread(updated_file, nrows = 0, showProgress = FALSE))
  if (!ID %in% input_cols) stop(sprintf("Column '%s' not found after rsID update", ID))
  renamed_cols <- renamed_input_colnames(input_cols)

  message("Counting rows for chunk planning...")
  total_data_rows <- max(0L, count_file_lines(updated_file) - 1L)
  desired_workers <- max(1L, bb_workers)
  if (.Platform$OS.type == "windows" && desired_workers > 1L) desired_workers <- 1L
  effective_chunk_size <- as.integer(chunk_size)
  if (effective_chunk_size == 0L) {
    if (total_data_rows > 0L) {
      # Auto mode: aim for ~1 chunk per worker for predictable fan-out and lower overhead.
      target_chunks <- max(1L, desired_workers)
      effective_chunk_size <- as.integer(ceiling(total_data_rows / target_chunks))
      effective_chunk_size <- max(1L, effective_chunk_size)
    } else {
      effective_chunk_size <- 1L
    }
  }
  message(sprintf(
    "BigBed chunk settings: rows=%s workers=%d chunk_size=%s",
    format(total_data_rows, big.mark = ","),
    desired_workers,
    format(effective_chunk_size, big.mark = ",")
  ))

  split_prefix <- file.path(work_dir, "updated_chunk_")
  message("Splitting input into chunks...")
  t_split <- Sys.time()
  split_cmd <- sprintf(
    "tail -n +2 %s | split --suffix-length=5 --numeric-suffixes --lines=%s - %s",
    shQuote(updated_file),
    effective_chunk_size,
    shQuote(split_prefix)
  )
  status <- system(split_cmd)
  if (!identical(status, 0L)) stop("Splitting updated input into chunks failed")
  chunk_files <- sort(list.files(work_dir, pattern = "^updated_chunk_[0-9]+$", full.names = TRUE))
  message(sprintf(
    "Created %s chunk(s) in %s",
    format(length(chunk_files), big.mark = ","),
    format_seconds(as.numeric(difftime(Sys.time(), t_split, units = "secs")))
  ))

  bb_tool <- find_bigbed_tool()
  matched_stream_file <- file.path(work_dir, "matched_unsorted.tsv")
  total_chunks <- length(chunk_files)
  worker_count <- if (total_chunks > 0L) max(1L, min(total_chunks, desired_workers)) else 0L
  nomatch_parts <- character()
  multipos_parts <- character()
  total_matches <- 0L
  total_nomatch <- 0L
  total_filtered_alt <- 0L
  total_multi_pos <- 0L

  if (total_chunks == 0L) {
    out_cols <- c(ID, "CHROM", "POS", setdiff(renamed_cols, ID))
    empty <- as.data.table(setNames(replicate(length(out_cols), logical(0), simplify = FALSE), out_cols))
    fwrite(empty, output_main, sep = "\t", quote = FALSE)
  } else {
    process_chunk <- function(i) {
      data.table::setDTthreads(1L)
      chunk <- fread(chunk_files[[i]], header = FALSE, col.names = input_cols, showProgress = FALSE)
      if (nrow(chunk) == 0L) {
        return(list(
          matched_file = "",
          nomatch_file = "",
          multipos_file = "",
          matched_n = 0L,
          nomatch_n = 0L,
          filtered_alt_n = 0L,
          multi_pos_n = 0L
        ))
      }

      setnames(chunk, c("CHROM", "POS0", "POS"), c("CHROM_old", "POS0_old", "POS_old"), skip_absent = TRUE)
      chunk_ids <- as.character(chunk[[ID]])
      query_ids <- unique(chunk_ids[grep("^rs", chunk_ids)])
      snppos <- data.table(CHROM = character(), POS0 = integer(), POS = integer(), ID = character())
      matched <- NULL
      filtered_alt <- 0L
      multi_pos <- data.table()
      matched_file <- file.path(work_dir, sprintf("matched_%05d.tsv", i))
      nomatch_file <- file.path(work_dir, sprintf("nomatch_%05d.tsv", i))
      multipos_file <- file.path(work_dir, sprintf("multipos_%05d.tsv", i))

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
          snppos[, CHROM := gsub("^chr", "", CHROM, ignore.case = TRUE)]
          if (!include_alt_chrom) {
            keep_primary <- is_primary_chrom(snppos$CHROM)
            filtered_alt <- sum(!keep_primary)
            if (filtered_alt > 0L) snppos <- snppos[keep_primary]
          } else {
            filtered_alt <- 0L
          }
          if (nrow(snppos) > 0L) {
            split_pos <- split_multi_position_ids(snppos, ID)
            snppos <- split_pos$primary
            multi_pos <- split_pos$duplicates
            if (nrow(multi_pos) > 0L) {
              fwrite(multi_pos, multipos_file, sep = "\t", quote = FALSE, col.names = FALSE)
            }
          }
        }
      }

      if (nrow(snppos) > 0L) {
        matched <- merge(snppos, chunk, by = ID)
        indels <- which(matched$POS - matched$POS0 > 1)
        if (length(indels) > 0) matched[indels, POS := POS0]
        matched[, POS0 := NULL]
        fwrite(matched, matched_file, sep = "\t", quote = FALSE, col.names = FALSE)

        matched_ids <- unique(as.character(snppos[[ID]]))
        misses <- chunk[!chunk_ids %chin% matched_ids]
      } else {
        misses <- chunk
      }

      if (nrow(misses) > 0L) {
        fwrite(misses, nomatch_file, sep = "\t", quote = FALSE, col.names = FALSE)
      }

      matched_n <- if (!is.null(matched)) nrow(matched) else 0L
      nomatch_n <- nrow(misses)
      rm(chunk, chunk_ids, query_ids, snppos, matched, misses)
      gc(verbose = FALSE)
      list(
        matched_file = matched_file,
        nomatch_file = nomatch_file,
        multipos_file = multipos_file,
        matched_n = matched_n,
        nomatch_n = nomatch_n,
        filtered_alt_n = as.integer(filtered_alt),
        multi_pos_n = as.integer(nrow(multi_pos))
      )
    }

    message(sprintf("Running BigBed chunk processing with %d worker(s) across %d chunk(s)", worker_count, total_chunks))
    idx <- seq_along(chunk_files)
    chunk_results <- vector("list", total_chunks)
    t_chunks <- Sys.time()
    done_chunks <- 0L

    if (worker_count == 1L) {
      for (i in idx) {
        chunk_results[[i]] <- process_chunk(i)
        done_chunks <- done_chunks + 1L
        message(progress_line(done_chunks, total_chunks, t_chunks))
      }
    } else {
      heartbeat_sec <- 15
      next_submit <- 1L
      active_jobs <- list()
      last_heartbeat <- Sys.time()

      while (done_chunks < total_chunks) {
        while (length(active_jobs) < worker_count && next_submit <= total_chunks) {
          chunk_i <- idx[[next_submit]]
          job <- parallel::mcparallel(process_chunk(chunk_i), silent = TRUE)
          active_jobs[[as.character(job$pid)]] <- list(job = job, idx = chunk_i)
          next_submit <- next_submit + 1L
        }

        if (length(active_jobs) == 0L) break

        collected <- parallel::mccollect(lapply(active_jobs, `[[`, "job"), wait = FALSE)
        made_progress <- FALSE
        if (!is.null(collected) && length(collected) > 0L) {
          for (pid in names(collected)) {
            chunk_i <- active_jobs[[pid]]$idx
            res_i <- collected[[pid]]
            active_jobs[[pid]] <- NULL
            if (inherits(res_i, "try-error")) {
              stop(sprintf("Chunk %d failed: %s", chunk_i, as.character(res_i)))
            }
            chunk_results[[chunk_i]] <- res_i
            done_chunks <- done_chunks + 1L
            made_progress <- TRUE
          }
        }

        now <- Sys.time()
        if (made_progress) {
          message(sprintf(
            "%s (active=%d)",
            progress_line(done_chunks, total_chunks, t_chunks),
            length(active_jobs)
          ))
          last_heartbeat <- now
        } else if (as.numeric(difftime(now, last_heartbeat, units = "secs")) >= heartbeat_sec) {
          message(sprintf(
            "%s (active=%d, waiting)",
            progress_line(done_chunks, total_chunks, t_chunks, label = "BigBed heartbeat"),
            length(active_jobs)
          ))
          last_heartbeat <- now
        }

        if (!made_progress) Sys.sleep(1)
      }
    }

    matched_parts <- unlist(lapply(chunk_results, `[[`, "matched_file"), use.names = FALSE)
    matched_parts <- matched_parts[nzchar(matched_parts) & file.exists(matched_parts) & file.info(matched_parts)$size > 0]
    nomatch_parts <- unlist(lapply(chunk_results, `[[`, "nomatch_file"), use.names = FALSE)
    nomatch_parts <- nomatch_parts[nzchar(nomatch_parts) & file.exists(nomatch_parts) & file.info(nomatch_parts)$size > 0]
    multipos_parts <- unlist(lapply(chunk_results, `[[`, "multipos_file"), use.names = FALSE)
    multipos_parts <- multipos_parts[nzchar(multipos_parts) & file.exists(multipos_parts) & file.info(multipos_parts)$size > 0]
    total_matches <- sum(as.integer(unlist(lapply(chunk_results, `[[`, "matched_n"), use.names = FALSE)))
    total_nomatch <- sum(as.integer(unlist(lapply(chunk_results, `[[`, "nomatch_n"), use.names = FALSE)))
    total_filtered_alt <- sum(as.integer(unlist(lapply(chunk_results, `[[`, "filtered_alt_n"), use.names = FALSE)))
    total_multi_pos <- sum(as.integer(unlist(lapply(chunk_results, `[[`, "multi_pos_n"), use.names = FALSE)))

    if (length(matched_parts) > 0L) {
      writeLines(paste(c(ID, "CHROM", "POS", setdiff(renamed_cols, ID)), collapse = "\t"), matched_stream_file)
      for (part in matched_parts) {
        status <- system(sprintf("cat %s >> %s", shQuote(part), shQuote(matched_stream_file)))
        if (!identical(status, 0L)) stop("Failed to append matched chunk output")
      }
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

  if (length(nomatch_parts) > 0L) {
    writeLines(paste(renamed_cols, collapse = "\t"), output_nomatch)
    for (part in nomatch_parts) {
      status <- system(sprintf("cat %s >> %s", shQuote(part), shQuote(output_nomatch)))
      if (!identical(status, 0L)) stop("Failed to append no-match chunk output")
    }
  }
  if (length(multipos_parts) > 0L) {
    multipos_concat <- file.path(work_dir, "multipos_concat.tsv")
    for (part in multipos_parts) {
      status <- system(sprintf("cat %s >> %s", shQuote(part), shQuote(multipos_concat)))
      if (!identical(status, 0L)) stop("Failed to append multi-position chunk output")
    }
    multipos_sorted <- file.path(work_dir, "multipos_sorted.tsv")
    status <- system(sprintf(
      "LC_ALL=C sort -t$'\\t' -k1,1 -k2,2 -k3,3n -k4,4n -u %s > %s",
      shQuote(multipos_concat),
      shQuote(multipos_sorted)
    ))
    if (!identical(status, 0L)) stop("Failed to sort/deduplicate multi-position output")
    writeLines(paste(c(ID, "CHROM", "POS0", "POS"), collapse = "\t"), output_multipos)
    status <- system(sprintf("cat %s >> %s", shQuote(multipos_sorted), shQuote(output_multipos)))
    if (!identical(status, 0L)) stop("Failed to write multi-position output")
  }
  message(sprintf(
    "BigBed chunk processing complete: matched=%s, noMatch=%s, filtered_non_primary=%s, multi_position=%s",
    format(total_matches, big.mark = ","),
    format(total_nomatch, big.mark = ","),
    format(total_filtered_alt, big.mark = ","),
    format(total_multi_pos, big.mark = ",")
  ))
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
    if (!include_alt_chrom && nrow(snppos) > 0L) {
      keep_primary <- is_primary_chrom(snppos$CHROM)
      filtered_alt <- sum(!keep_primary)
      if (filtered_alt > 0L) {
        snppos <- snppos[keep_primary]
        message(sprintf("Filtered %s non-primary contig records", format(filtered_alt, big.mark = ",")))
      }
    }
    if (nrow(snppos) > 0L) {
      split_pos <- split_multi_position_ids(snppos, ID)
      snppos <- split_pos$primary
      if (nrow(split_pos$duplicates) > 0L) {
        fwrite(unique(split_pos$duplicates), output_multipos, sep = "\t", quote = FALSE)
      }
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
    fwrite(updated1[noMatch, ], output_nomatch, sep = "\t", quote = FALSE)
  }

  # Use zero based positions for larger indels to match VCF nomenclature
  indels <- which(updated2$POS - updated2$POS0 > 1)
  if (length(indels) > 0) updated2[indels, POS := POS0]
  updated2[, POS0 := NULL]

  suppressWarnings(updated2 <- updated2[order(as.numeric(CHROM), POS), ])
  fwrite(updated2, output_main, sep = "\t", quote = FALSE)
}
