#!/usr/bin/env Rscript

# Reproduce the manuscript analyses in dependency order.
#
# Usage:
#   Rscript scripts/run_all.R [data_file] [output_dir] [permutations] [bootstraps]
#
# Repeatability is not run here because it requires the separate raw
# remeasurement archive and measurement-order file. Its maintained script is
# packaged as repeatability_analysis.R and accepts those locations explicitly.

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "data/trait_data.csv"
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "analysis_outputs/weta_trait_analysis"
}
permutations <- if (length(args) >= 3L) args[[3L]] else "9999"
bootstraps <- if (length(args) >= 4L) args[[4L]] else "9999"

if (!file.exists(data_file)) stop("Trait data not found: ", data_file)

command <- commandArgs(trailingOnly = FALSE)
file_argument <- command[grepl("^--file=", command)]
script_dir <- if (length(file_argument) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_argument)))
} else {
  normalizePath("scripts")
}

rscript <- file.path(R.home("bin"), "Rscript")

run_stage <- function(script, stage_args) {
  path <- file.path(script_dir, script)
  if (!file.exists(path)) stop("Analysis stage not found: ", path)
  message("Running ", script)
  # Quote every path/argument because the canonical data filename contains
  # parentheses and project paths may contain spaces on another machine.
  status <- system2(rscript, shQuote(c(path, stage_args)))
  if (!identical(status, 0L)) {
    stop(script, " failed with exit status ", status)
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_stage("core_trait_analysis.R", c(data_file, output_dir))
run_stage("reviewer_reanalysis.R", c(data_file, output_dir))
run_stage(
  "female_reference_multivariate_reanalysis.R",
  c(data_file, output_dir, permutations, bootstraps)
)

manifest <- data.frame(
  item = c(
    "run_completed_utc",
    "input_file",
    "input_md5",
    "output_directory",
    "residual_randomizations",
    "bootstrap_replicates",
    "R_version"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    normalizePath(data_file),
    unname(tools::md5sum(data_file)),
    normalizePath(output_dir),
    permutations,
    bootstraps,
    R.version.string
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  manifest,
  file.path(output_dir, "analysis_run_manifest.csv"),
  row.names = FALSE
)

message("All manuscript analyses completed: ", normalizePath(output_dir))
