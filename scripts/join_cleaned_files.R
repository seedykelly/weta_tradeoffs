library(tidyverse)

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

# Usage: Rscript scripts/join_cleaned_files.R [clean_measurement_dir] [output_file]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Expected at most source_dir and output_file.")
source_dir <- if (length(args) >= 1L) args[[1L]] else {
  file.path(project_root, "data", "clean", "measurements")
}
output_file <- if (length(args) >= 2L) args[[2L]] else {
  file.path(dirname(source_dir), "trait_data.csv")
}

# Import cleaned measurement files
forelegs <- read_csv(
  file.path(source_dir, "forelegs_clean.csv"),
  show_col_types = FALSE
)

midlegs <- read_csv(
  file.path(source_dir, "midlegs_clean.csv"),
  show_col_types = FALSE
)

hindlegs <- read_csv(
  file.path(source_dir, "hindlegs_clean.csv"),
  show_col_types = FALSE
)

head <- read_csv(
  file.path(source_dir, "head_clean.csv"),
  show_col_types = FALSE
)

pronotum <- read_csv(
  file.path(source_dir, "pronotum_clean.csv"),
  show_col_types = FALSE
)

id_sex <- read_csv(
  file.path(source_dir, "id_sex.csv"),
  show_col_types = FALSE
)


# Validate every input BEFORE joining: duplicate keys multiply observations.
inputs <- list(pronotum = pronotum, head = head, forelegs = forelegs,
               midlegs = midlegs, hindlegs = hindlegs, id_sex = id_sex)
iwalk(inputs, validate_unique_ids)
if (!"sex" %in% names(id_sex) || anyNA(id_sex$sex) ||
    any(!id_sex$sex %in% c("f", "m"))) {
  stop("id_sex must contain sex coded as 'f' or 'm' for every ID.")
}
measurement_ids <- unique(unlist(lapply(inputs[names(inputs) != "id_sex"], `[[`, "ID")))
missing_sex <- setdiff(measurement_ids, id_sex$ID)
if (length(missing_sex)) {
  stop("Missing sex metadata for IDs: ", paste(missing_sex, collapse = ", "))
}

trait_data <- list(
  pronotum,
  head,
  forelegs,
  midlegs,
  hindlegs
) %>%
  reduce(full_join, by = "ID") %>%
  left_join(id_sex, by = "ID") %>%
  relocate(ID, sex)

trait_data %>%
  print(n=320)

glimpse(trait_data)

nrow(trait_data)
n_distinct(trait_data$ID)

list(
  pronotum = pronotum,
  head = head,
  forelegs = forelegs,
  midlegs = midlegs,
  hindlegs = hindlegs,
  id_sex = id_sex
) %>%
  imap_dfr(
    ~ tibble(
      file = .y,
      rows = nrow(.x),
      unique_IDs = n_distinct(.x$ID),
      duplicates = nrow(.x) - n_distinct(.x$ID)
    )
  )


# ============================================================
# Average left and right measurements for bilateral traits
# ============================================================
#
# rowMeans(..., na.rm = TRUE) returns the available side when
# only one side was measured, and NaN when both are missing.
# ============================================================

trait_data <- trait_data |>
  mutate(
    forefemur = rowMeans(pick(left_forefemur,  right_forefemur),  na.rm = TRUE),
    foretibia = rowMeans(pick(left_foretibia,  right_foretibia),  na.rm = TRUE),
    ear       = rowMeans(pick(left_ear,        right_ear),        na.rm = TRUE),
    midfemur  = rowMeans(pick(left_midfemur,   right_midfemur),   na.rm = TRUE),
    midtibia  = rowMeans(pick(left_midtibia,   right_midtibia),   na.rm = TRUE),
    hindfemur = rowMeans(pick(left_hindfemur,  right_hindfemur),  na.rm = TRUE),
    hindtibia = rowMeans(pick(left_hindtibia,  right_hindtibia),  na.rm = TRUE),
    eye       = rowMeans(pick(left_eye_length, right_eye_length), na.rm = TRUE)
  ) |>
  select(
    ID, sex,
    pronotum,
    head_length, head_width,
    forefemur, foretibia,
    midfemur,  midtibia,
    hindfemur, hindtibia,
    ear, eye
  )

glimpse(trait_data)


# ------------------------------------------------------------
# Add morph column
# ------------------------------------------------------------

trait_data <- trait_data |>
  mutate(
    morph = case_when(
      sex == "f"                        ~ "female",
      is.na(head_length)                 ~ NA_character_,
      sex == "m" & head_length < 18.50579 ~ "eighth",
      sex == "m" & head_length > 24.15225 ~ "tenth",
      sex == "m"                        ~ "ninth"
    )
  ) |>
  relocate(morph, .after = sex)


# ------------------------------------------------------------
# Save final combined dataset
# ------------------------------------------------------------

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(trait_data, output_file, na = "")
