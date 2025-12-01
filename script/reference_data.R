SUPPORTED_BUILDS <- c("hg19", "hg38")
SUPPORTED_DBSNP_VERSIONS <- c("151", "155")

REFERENCE_URLS <- list(
  hg19 = list(
    "151" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz",
    "155" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp155.txt.gz"
  ),
  hg38 = list(
    "151" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz",
    "155" = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp155.txt.gz"
  )
)

RS_MERGE_URL <- "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz"

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

download_if_missing <- function(url, destfile) {
  if (file.exists(destfile)) return(invisible(destfile))
  message(sprintf("Downloading %s", url))
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
      "| grep -v -e Un_ -e hap -e random -e Y -e fix -e alt -e M\\t",
      "| sed 's/chr//g' >",
      shQuote(filtered_path)
    )
    message(sprintf("Filtering dbSNP %s %s to %s", version, build, filtered_path))
    status <- system(filter_cmd)
    if (!identical(status, 0L)) stop("Filtering dbSNP reference failed")
  }

  split_cmd <- paste(
    "split --suffix-length=3 --numeric-suffixes",
    sprintf("--lines=%s", split_lines),
    shQuote(filtered_path),
    shQuote(split_prefix_path)
  )
  message(sprintf("Splitting filtered file into chunks (prefix: %s)", split_prefix_path))
  status <- system(split_cmd)
  if (!identical(status, 0L)) stop("Splitting dbSNP reference failed")

  if (remove_download && file.exists(gz_path)) file.remove(gz_path)

  split_files <- split_files_for(filtered_path, data_dir)
  if (length(split_files) == 0) stop("No split dbSNP files were created")

  list(filtered_path = filtered_path, split_files = sort(split_files))
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
