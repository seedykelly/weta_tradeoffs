#!/usr/bin/env Rscript

# Reproduce the manuscript analyses in dependency order.
#
# Usage:
#   Rscript scripts/run_all.R [data_file] [output_dir] [permutations] [bootstraps]
#     [rarefactions]
# Defaults: 9999 permutations, 9999 bootstraps, 2000 rarefactions.
#
# Repeatability is not run here because it requires the separate raw
# remeasurement archive and measurement-order file. Its maintained script is
# packaged as repeatability_analysis.R and accepts those locations explicitly.

# Resolve defaults from this script, including when sourced from another directory.
.script_file <- local({
  source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(source_files)) {
    tail(source_files, 1L)[[1L]]
  } else {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg) != 1L) stop("Run this file with Rscript or source().")
    path <- sub("^--file=", "", file_arg[[1L]])
    # Rscript may encode spaces in --file as "~+~".
    if (!file.exists(path)) path <- gsub("~+~", " ", path, fixed = TRUE)
    path
  }
})
.script_file <- normalizePath(.script_file, mustWork = TRUE)
source(file.path(dirname(.script_file), "workflow_helpers.R"), local = TRUE)
project_root <- dirname(dirname(.script_file))

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else file.path(project_root, "data", "trait_data.csv")
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  file.path(project_root, "analysis_outputs", "weta_trait_analysis")
}
if (length(args) > 5L) stop("Expected at most five arguments; see usage.")
permutations <- parse_replicates(
  if (length(args) >= 3L) args[[3L]] else 9999L, "Permutations", 999L
)
bootstraps <- parse_replicates(
  if (length(args) >= 4L) args[[4L]] else 9999L, "Bootstraps", 999L
)
rarefactions <- parse_replicates(
  if (length(args) >= 5L) args[[5L]] else 2000L, "Rarefactions"
)

if (!file.exists(data_file)) stop("Trait data not found: ", data_file)

script_dir <- dirname(.script_file)

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
run_stage("reviewer_reanalysis.R", c(data_file, output_dir, permutations, rarefactions))
run_stage(
  "female_reference_multivariate_reanalysis.R",
  c(data_file, output_dir, permutations, bootstraps)
)
# Runs last because it consolidates the covariance sign audit from the
# upstream correlation table and adds the instar-matched decomposition and
# the precision reporting.
run_stage(
  "instar_matched_precision_analysis.R",
  c(data_file, output_dir)
)

manifest <- data.frame(
  item = c(
    "run_completed_utc",
    "input_file",
    "input_md5",
    "output_directory",
    "residual_randomizations",
    "bootstrap_replicates",
    "pca_parallel_permutations",
    "covariance_matrix_permutations",
    "covariance_rarefactions",
    "random_skewers_replicates",
    "core_rng_seed",
    "reviewer_model_rng_seed",
    "precision_model_rng_seed",
    "R_version"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    normalizePath(data_file),
    unname(tools::md5sum(data_file)),
    normalizePath(output_dir),
    permutations,
    bootstraps,
    permutations,
    permutations,
    rarefactions,
    10000L,
    20260921L,
    2301L,
    2302L,
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
