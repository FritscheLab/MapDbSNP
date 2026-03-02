SUPPORTED_BUILDS <- c("hg19", "hg38")
SUPPORTED_DBSNP_VERSIONS <- c("151", "153", "155")

REFERENCE_URLS <- list(
  hg19 = list(
    "151" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz",
    "153" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp153.txt.gz",
    "155" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp155.txt.gz"
  ),
  hg38 = list(
    "151" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz",
    "153" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp153.txt.gz",
    "155" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp155.txt.gz"
  )
)

RS_MERGE_URL <- "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz"

# BigBed sources (much smaller than text dumps)
BIGBED_URLS <- list(
  hg19 = list(
    "153" = "http://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153.bb",
    "155" = "http://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp155.bb",
    "151" = "http://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp151.bb"
  ),
  hg38 = list(
    "153" = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb",
    "155" = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb",
    "151" = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp151.bb"
  )
)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

run_shell_cmd <- function(cmd, pipefail = FALSE) {
  if (pipefail) {
    bash <- Sys.which("bash")
    if (!nzchar(bash)) stop("bash is required for pipeline execution with pipefail")
    # Use `system()` here; on some platforms `system2()` does not reliably enforce `-o pipefail`.
    status <- system(sprintf("%s -o pipefail -c %s", shQuote(bash), shQuote(cmd)))
  } else {
    status <- system(cmd)
  }
  if (!identical(status, 0L)) {
    stop(sprintf("Shell command failed (status=%s): %s", status, cmd))
  }
  invisible(status)
}

download_if_missing <- function(url, destfile, connections = 8) {
  if (file.exists(destfile)) return(invisible(destfile))
  ensure_dir(dirname(destfile))
  message(sprintf("Downloading %s", url))

  aria <- Sys.which("aria2c")
  if (nzchar(aria)) {
    dest_dir <- normalizePath(dirname(destfile), winslash = "/", mustWork = FALSE)
    # Use multi-connection download when available
    args <- c(
      sprintf("-x%d", connections),
      sprintf("-s%d", connections),
      "-k1M",
      "-m3",
      "--allow-overwrite=true",
      "--auto-file-renaming=false",
      "-d", dest_dir,
      "-o", basename(destfile),
      url
    )
    status <- system2(aria, args = args, stdout = "", stderr = "")
    if (identical(status, 0L) && file.exists(destfile)) return(invisible(destfile))
    message("aria2c download failed or incomplete, falling back to download.file()")
  }

  download.file(url = url, destfile = destfile, mode = "wb", quiet = FALSE)
  invisible(destfile)
}

choose_decompressor <- function(cpus) {
  pigz <- Sys.which("pigz")
  if (nzchar(pigz)) {
    threads <- max(1L, min(as.integer(cpus), 16L))
    return(sprintf("pigz -p %s -dc", threads))
  }
  "gzip -dc"
}

split_prefix <- function(filtered_path) {
  sub("\\.txt$", "", filtered_path)
}

split_files_for <- function(filtered_path, data_dir) {
  prefix <- basename(split_prefix(filtered_path))
  list.files(
    data_dir,
    pattern = paste0("^", prefix, "[0-9]+$"),
    full.names = TRUE
  )
}

validate_split_reference <- function(filtered_path, split_files) {
  if (!file.exists(filtered_path)) {
    stop(sprintf(
      "Found split dbSNP chunks but missing filtered reference '%s'. Remove stale chunks and rebuild.",
      filtered_path
    ))
  }
  if (length(split_files) == 0L) stop("No split dbSNP files found")
  info <- file.info(split_files)
  if (any(is.na(info$size)) || any(info$size <= 0L)) {
    stop("Found zero-byte or unreadable dbSNP split files. Remove stale chunks and rebuild.")
  }
  filtered_info <- file.info(filtered_path)
  if (is.na(filtered_info$size) || filtered_info$size <= 0L) {
    stop("Filtered dbSNP reference is empty or unreadable. Remove stale files and rebuild.")
  }
}

ensure_reference_data <- function(build,
                                  version,
                                  data_dir = here::here("data"),
                                  cpus = 2,
                                  split_lines = 10000000,
                                  remove_download = TRUE) {
  build <- tolower(build)
  version <- as.character(version)

  if (!build %in% SUPPORTED_BUILDS) {
    stop(sprintf("Unsupported build '%s'. Choose from: %s", build, paste(SUPPORTED_BUILDS, collapse = ", ")))
  }
  if (!version %in% SUPPORTED_DBSNP_VERSIONS) {
    stop(sprintf("Unsupported dbSNP version '%s'. Choose from: %s", version, paste(SUPPORTED_DBSNP_VERSIONS, collapse = ", ")))
  }
  if (!nzchar(Sys.which("split"))) stop("Required tool 'split' not found in PATH")
  if (!nzchar(Sys.which("gzip")) && !nzchar(Sys.which("pigz"))) stop("Required tool 'gzip' or 'pigz' not found in PATH")

  ensure_dir(data_dir)

  filtered_path <- file.path(data_dir, sprintf("snp%s_%s_filtered.txt", version, build))
  split_prefix_path <- split_prefix(filtered_path)
  split_files <- split_files_for(filtered_path, data_dir)

  if (length(split_files) > 0) {
    validate_split_reference(filtered_path, split_files)
    return(list(filtered_path = filtered_path, split_files = sort(split_files)))
  }

  url <- REFERENCE_URLS[[build]][[version]]
  if (is.null(url)) stop(sprintf("No reference URL configured for build=%s, version=%s", build, version))

  gz_path <- file.path(data_dir, basename(url))
  download_if_missing(url, gz_path)

  if (!file.exists(filtered_path)) {
    decompressor <- choose_decompressor(cpus)
    filter_cmd <- paste(
      decompressor, shQuote(gz_path),
      "| cut -f 2-5",
      # Keep primary chromosomes (1-22, X, Y, M/MT); runtime filtering decides whether to drop non-primary contigs.
      "| grep -v -e Un_ -e hap -e random -e fix -e alt",
      "| sed 's/chr//g' >",
      shQuote(filtered_path)
    )
    message(sprintf("Filtering dbSNP %s %s to %s", version, build, filtered_path))
    run_shell_cmd(filter_cmd, pipefail = TRUE)
  }

  split_cmd <- paste(
    "split --suffix-length=3 --numeric-suffixes",
    sprintf("--lines=%s", split_lines),
    shQuote(filtered_path),
    shQuote(split_prefix_path)
  )
  message(sprintf("Splitting filtered file into chunks (prefix: %s)", split_prefix_path))
  run_shell_cmd(split_cmd)

  if (remove_download && file.exists(gz_path)) file.remove(gz_path)

  split_files <- split_files_for(filtered_path, data_dir)
  if (length(split_files) == 0) stop("No split dbSNP files were created")
  validate_split_reference(filtered_path, split_files)

  list(filtered_path = filtered_path, split_files = sort(split_files))
}

ensure_bigbed <- function(build,
                          version,
                          data_dir = here::here("data"),
                          download = TRUE) {
  build <- tolower(build)
  version <- as.character(version)
  ensure_dir(data_dir)

  # Prefer build-specific filename to avoid collisions if both builds are stored
  target <- file.path(data_dir, sprintf("dbSnp%s_%s.bb", version, build))

  # Backward-compatible fallback to legacy name without build
  legacy <- file.path(data_dir, sprintf("dbSnp%s.bb", version))
  if (file.exists(target)) return(target)
  if (file.exists(legacy)) {
    stop(sprintf(
      "Found legacy BigBed '%s' without build in the filename. Refusing to guess whether it is for %s. Rename it to '%s' or pass --bb-file explicitly.",
      legacy,
      build,
      basename(target)
    ))
  }

  if (!download) return("")
  url <- BIGBED_URLS[[build]][[version]]
  if (is.null(url)) stop(sprintf("No BigBed URL configured for build=%s version=%s", build, version))

  message(sprintf("Downloading dbSNP%s BigBed for %s", version, build))
  download_if_missing(url, target)
  target
}

ensure_rsmerge <- function(data_dir = here::here("data")) {
  ensure_dir(data_dir)
  if (!nzchar(Sys.which("gzip")) && !nzchar(Sys.which("pigz"))) stop("Required tool 'gzip' or 'pigz' not found in PATH")
  merge_file <- file.path(data_dir, "RsMergeArch.bcp")
  merge_gz <- paste0(merge_file, ".gz")

  if (!file.exists(merge_file)) {
    download_if_missing(RS_MERGE_URL, merge_gz)
    message("Decompressing RsMerge archive")
    decompressor <- if (nzchar(Sys.which("pigz"))) "pigz" else "gzip"
    status <- system(paste(decompressor, "-d -f", shQuote(merge_gz)))
    if (!identical(status, 0L)) stop("Decompressing RsMerge archive failed")
  }

  merge_file
}
