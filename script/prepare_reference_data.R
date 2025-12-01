suppressPackageStartupMessages({
  library("optparse")
  library("here")
})

source(here::here("script/reference_data.R"))

option_list <- list(
  make_option("--build", type = "character", default = "hg19", help = "Genome build to prepare: hg19, hg38, or both"),
  make_option("--dbsnp-version", type = "character", default = "155", help = "dbSNP release to prepare (151 or 155)"),
  make_option("--data-dir", type = "character", default = "", help = "Directory for reference data (default: ./data)"),
  make_option("--cpus", type = "integer", default = 6, help = "CPUs to use during filtering/splitting"),
  make_option("--split-lines", type = "integer", default = 10000000, help = "Lines per chunk when splitting dbSNP")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

data_dir <- if (nzchar(opt$data_dir)) opt$data_dir else here::here("data")
builds <- tolower(opt$build)

if (builds == "both") {
  builds <- SUPPORTED_BUILDS
} else {
  builds <- tolower(builds)
  if (!builds %in% SUPPORTED_BUILDS) {
    stop(sprintf("Unsupported build '%s'. Choose from hg19, hg38, or both.", builds))
  }
}

message(sprintf("Preparing dbSNP%s reference for: %s", opt$dbsnp_version, paste(builds, collapse = ", ")))
for (b in builds) {
  ensure_reference_data(
    build = b,
    version = opt$dbsnp_version,
    data_dir = data_dir,
    cpus = opt$cpus,
    split_lines = opt$split_lines
  )
}

ensure_rsmerge(data_dir)
message("Finished preparing reference data.")
