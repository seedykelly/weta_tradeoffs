#!/usr/bin/env Rscript
# Fast regression checks. All generated data stay in a temporary directory.
# Usage: Rscript tests/test_workflow.R
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_file <- sub("^--file=", "", file_arg[[1L]])
if (!file.exists(test_file)) test_file <- gsub("~+~", " ", test_file, fixed = TRUE)
project <- dirname(dirname(normalizePath(test_file, mustWork = TRUE)))
source(file.path(project, "scripts", "workflow_helpers.R"))
scratch <- tempfile("weta-workflow-tests-")
dir.create(scratch)
old_wd <- setwd(scratch)
rscript <- file.path(R.home("bin"), "Rscript")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
stopifnot(identical(parse_replicates("999", "Test", 999L), 999L))
for (bad in c("998", "999.5", "NA", "Inf", "no", "2147483648")) {
  expect_error(parse_replicates(bad, "Test", 999L), "whole number")
}
expect_error(parse_replicates(character(), "Test"), "whole number")
expect_error(parse_replicates(c(999, 1000), "Test"), "whole number")
stopifnot(validate_unique_ids(data.frame(ID = c("a", "b")), "fixture"))
expect_error(validate_unique_ids(data.frame(ID = c("a", "a")), "fixture"), "duplicate")
expect_error(validate_unique_ids(data.frame(ID = c("a", NA)), "fixture"), "missing")
expect_error(validate_unique_ids(data.frame(ID = c("a", " ")), "fixture"), "blank")
expect_error(validate_unique_ids(data.frame(ID = c("a", " b")), "fixture"), "spaces")
expect_error(validate_unique_ids(data.frame(other = "a"), "fixture"), "no ID")

run <- function(script, args, expected = 0L, pattern = NULL) {
  log <- tempfile(tmpdir = scratch, fileext = ".log")
  status <- system2(rscript, shQuote(c(script, args)), stdout = log, stderr = log)
  output <- paste(readLines(log, warn = FALSE), collapse = "\n")
  if (status != expected || (!is.null(pattern) && !grepl(pattern, output))) {
    stop("Unexpected result for ", basename(script), ":\n", output)
  }
}
script <- function(name) file.path(project, "scripts", name)
data_file <- file.path(project, "data", "trait_data.csv")
bad_output <- file.path(scratch, "invalid counts")
run(script("run_all.R"), c(data_file, bad_output, "999.5", "999"), 1L, "whole number")
run(script("run_all.R"), c(data_file, bad_output, "999", "999", "0"), 1L, "whole number")
run(script("reviewer_reanalysis.R"), c(data_file, bad_output, "Inf"), 1L, "whole number")
run(script("female_reference_multivariate_reanalysis.R"),
    c(data_file, bad_output, "999", "999.5"), 1L, "whole number")
stopifnot(!dir.exists(bad_output))

# Rebuild every cleaned measurement file without overwriting the originals.
clean <- file.path(scratch, "clean measurements")
raw <- file.path(project, "data", "raw", "measurements")
for (trait in c("foreleg", "midleg", "hindleg", "head", "pronotum")) {
  run(script(paste0("data_cleaning_", trait, ".R")), c(raw, clean))
}
original_clean <- file.path(project, "data", "clean", "measurements")
stopifnot(file.copy(file.path(original_clean, "id_sex.csv"), clean))
joined <- file.path(scratch, "joined traits.csv")
run(script("join_cleaned_files.R"), c(clean, joined))
result <- read.csv(joined)
baseline <- read.csv(file.path(project, "data", "clean", "trait_data.csv"))
result <- result[order(result$ID), ]
baseline <- baseline[order(baseline$ID), ]
rownames(result) <- rownames(baseline) <- NULL
stopifnot(isTRUE(all.equal(result, baseline, tolerance = 1e-12)))

# A duplicate in ANY source must fail before an output is written.
for (name in c("forelegs_clean.csv", "midlegs_clean.csv", "hindlegs_clean.csv",
               "head_clean.csv", "pronotum_clean.csv", "id_sex.csv")) {
  path <- file.path(clean, name)
  original <- read.csv(path)
  write.csv(rbind(original, original[1L, ]), path, row.names = FALSE, na = "")
  target <- file.path(scratch, paste0("invalid-", name))
  run(script("join_cleaned_files.R"), c(clean, target), 1L, "duplicate IDs")
  stopifnot(!file.exists(target))
  write.csv(original, path, row.names = FALSE, na = "")
}

# Missing male head measurements must not silently assign the ninth instar.
head_path <- file.path(clean, "head_clean.csv")
head <- read.csv(head_path)
sex <- read.csv(file.path(clean, "id_sex.csv"))
male <- intersect(head$ID, sex$ID[sex$sex == "m"])[[1L]]
head$head_length[head$ID == male] <- NA_real_
write.csv(head, head_path, row.names = FALSE, na = "")
run(script("join_cleaned_files.R"), c(clean, joined))
missing_head_result <- read.csv(joined, na.strings = "")
stopifnot(is.na(missing_head_result$morph[missing_head_result$ID == male]))

# A relocated copy must use its own raw data and default output, not this project.
copy_root <- file.path(scratch, "relocated project")
copy_scripts <- file.path(copy_root, "scripts")
copy_raw <- file.path(copy_root, "data", "raw", "measurements")
dir.create(copy_scripts, recursive = TRUE)
dir.create(copy_raw, recursive = TRUE)
stopifnot(all(file.copy(script(c("data_cleaning_pronotum.R", "workflow_helpers.R")), copy_scripts)))
stopifnot(file.copy(file.path(raw, "pronotum_1.csv"), copy_raw))
copy_script <- file.path(copy_scripts, "data_cleaning_pronotum.R")
run(copy_script, character())
stopifnot(file.exists(file.path(copy_root, "data", "clean", "measurements", "pronotum_clean.csv")))
# Also exercise source() from a directory unrelated to the project.
env <- new.env(parent = globalenv())
invisible(capture.output(source(copy_script, local = env)))
stopifnot(identical(env$project_root, normalizePath(copy_root)))
setwd(old_wd)
cat("Workflow regression tests passed. Temporary outputs:", scratch, "\n")
