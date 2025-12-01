runParallel <- function(cmdLines, resources = 0.5) {
  require("parallel")
  if (resources <= 0) stop("Misspecified resources")
  if (length(cmdLines) == 0) stop("No jobs specified")

  total_cores <- parallel::detectCores()

  if (resources > 0 && resources < 1) {
    ncores <- max(1, floor(total_cores * resources))
  } else {
    ncores <- as.integer(min(total_cores, resources))
  }

  ncores <- max(1, min(length(cmdLines), ncores))

  if (.Platform$OS.type == "windows") {
    # mclapply falls back to serial on Windows; run sequentially for predictability
    ncores <- 1
  }

  if (ncores == 1) {
    message(sprintf("Running %s job(s) on 1 core", length(cmdLines)))
    for (cmd in cmdLines) {
      status <- system(cmd)
      if (!identical(status, 0L)) stop(sprintf("Job failed: %s", cmd))
    }
    message("Done")
    return(invisible())
  }

  okfiles <- file.path(tempdir(), sprintf("Check_%05d", seq_along(cmdLines)))
  cmdLines <- paste0(cmdLines, "; touch ", shQuote(okfiles))

  message(sprintf("Running %s jobs on %s cores", length(cmdLines), ncores))
  errors <- unlist(parallel::mclapply(cmdLines, system, mc.cores = ncores, mc.preschedule = FALSE))
  failedjobs <- length(which(errors > 0 | !file.exists(okfiles)))
  file.remove(okfiles)

  if (failedjobs > 0) stop(paste(failedjobs, "jobs failed"))
  message("Done")
}
