# Regression checks for QMD display helpers; no analyses or renders are run.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))

check_tables <- function() {
  extracted <- tempfile(fileext = ".R")
  on.exit(unlink(extracted), add = TRUE)
  knitr::purl(file.path(project_dir, "manuscript.qmd"),
              output = extracted, documentation = 0, quiet = TRUE)
  expressions <- parse(extracted)
  env <- new.env(parent = globalenv())
  wanted <- c("morph_label", "normalize_table_labels", "n_integration_traits")
  assignments <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted
  }, as.list(expressions))
  stopifnot(length(assignments) == length(wanted))
  env$n_distinct <- dplyr::n_distinct
  env$integration_pc_loadings <- readr::read_csv(
    file.path(project_dir, "analysis_outputs", "weta_trait_analysis", "tables",
              "pca_loadings.csv"), show_col_types = FALSE
  )
  for (expr in assignments) eval(expr, env)

  input <- data.frame(
    Group = c("eighth", "ninth", "tenth", "female"),
    Contrast = c("tenth - eighth", "ninth - midpoint(eighth,tenth)",
                 "tenth - female", "eighth - female"),
    Value = c("88.8%", "88.8–98.3%", "1.10–1.14", "8"),
    Quantity = c("8 traits", "10 comparisons", "9 rows", "female reference"),
    estimate = c(88.8, 1.10, 9, 10)
  )
  actual <- env$normalize_table_labels(input)
  stopifnot(
    identical(actual$Value, input$Value),
    identical(actual$Quantity, input$Quantity),
    identical(actual$estimate, input$estimate),
    identical(actual$Group,
              c("Eighth instar", "Ninth instar", "Tenth instar", "Female")),
    actual$Contrast[1] == "Tenth instar - Eighth instar",
    env$n_integration_traits == 5L
  )
  message("QMD table regression checks passed.")
}

check_tables()
