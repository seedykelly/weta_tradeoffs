# ============================================================
# Import and clean pronotum measurements
# File: pronotum_1.csv
# ============================================================

library(tidyverse)


# ------------------------------------------------------------
# 1. File locations
# ------------------------------------------------------------

raw_dir <- paste0(
  "/Users/uqam/Documents/Research_Admin/research projects/",
  "Trait compensation/data/raw/measurements"
)

clean_dir <- file.path(
  dirname(dirname(raw_dir)),
  "clean",
  "measurements"
)

dir.create(
  clean_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Identify the pronotum file
# ------------------------------------------------------------

pronotum_files <- list.files(
  path = raw_dir,
  pattern = "^pronotum_[0-9]+\\.csv$",
  full.names = TRUE
) |>
  str_sort(numeric = TRUE)

if (length(pronotum_files) == 0) {
  stop(
    "No pronotum file was found in:\n",
    raw_dir
  )
}

if (length(pronotum_files) > 1) {
  stop(
    "More than one pronotum file was found:\n",
    paste(basename(pronotum_files), collapse = "\n"),
    "\n\nThe code expects exactly one pronotum file."
  )
}

pronotum_file <- pronotum_files[[1]]

message(
  "Pronotum file found:\n",
  basename(pronotum_file)
)


# ------------------------------------------------------------
# 3. Import the pronotum file
# ------------------------------------------------------------

pronotum_raw <- read_csv(
  pronotum_file,
  na = c("", "NA", "NaN", "-", "N/A"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 4. Check that required columns are present
# ------------------------------------------------------------

required_columns <- c(
  "Label",
  "Length"
)

missing_columns <- setdiff(
  required_columns,
  names(pronotum_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing from ",
    basename(pronotum_file),
    ": ",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 5. Extract ID and clean the measurement
# ------------------------------------------------------------

pronotum <- pronotum_raw |>
  select(
    Label,
    Length
  ) |>
  mutate(
    source_row = row_number(),
    
    # Example:
    # ACGN-Pronotum.jpg
    # becomes:
    # ACGN
    ID = str_extract(
      Label,
      "^[^-]+"
    ),
    
    ID = str_to_upper(
      str_trim(ID)
    ),
    
    pronotum = parse_double(
      Length,
      na = c("", "NA", "NaN", "-", "N/A")
    )
  ) |>
  select(
    ID,
    pronotum
  ) |>
  arrange(ID)


# ------------------------------------------------------------
# 6. Check that every label produced a valid ID
# ------------------------------------------------------------

invalid_pronotum_labels <- pronotum |>
  filter(
    is.na(ID) |
      ID == ""
  )

if (nrow(invalid_pronotum_labels) > 0) {
  
  print(invalid_pronotum_labels)
  
  stop(
    "Some pronotum labels could not be converted to IDs. ",
    "Examine 'invalid_pronotum_labels'."
  )
}


# ------------------------------------------------------------
# 7. Check for duplicated individuals
# ------------------------------------------------------------

duplicate_pronotum_IDs <- pronotum |>
  count(
    ID,
    name = "n_measurements"
  ) |>
  filter(n_measurements > 1)

if (nrow(duplicate_pronotum_IDs) > 0) {
  
  print(duplicate_pronotum_IDs)
  
  stop(
    "Some individuals have more than one pronotum measurement. ",
    "Examine 'duplicate_pronotum_IDs'."
  )
}


# ------------------------------------------------------------
# 8. Summarize missing measurements
# ------------------------------------------------------------

pronotum_missing_summary <- pronotum |>
  summarise(
    n_individuals = n(),
    missing_pronotum = sum(is.na(pronotum)),
    complete_individuals = sum(!is.na(pronotum))
  )

print(pronotum_missing_summary)


# List individuals with missing pronotum measurements

pronotum_missing <- pronotum |>
  filter(
    is.na(pronotum)
  )

print(pronotum_missing)


# ------------------------------------------------------------
# 9. Optional range check
# ------------------------------------------------------------
#
# Values are flagged but are not altered or removed.
# Adjust the upper threshold if necessary.
# ------------------------------------------------------------

pronotum_suspect_values <- pronotum |>
  filter(
    !is.na(pronotum) &
      (
        pronotum <= 0 |
          pronotum > 15
      )
  )

print(pronotum_suspect_values)


# ------------------------------------------------------------
# 10. View the cleaned data
# ------------------------------------------------------------

glimpse(pronotum)

print(
  pronotum,
  n = 20
)


# ------------------------------------------------------------
# 11. Save the cleaned pronotum data
# ------------------------------------------------------------

pronotum_output_file <- file.path(
  clean_dir,
  "pronotum_clean.csv"
)

write_csv(
  pronotum,
  pronotum_output_file,
  na = ""
)

message(
  "Cleaned pronotum data saved to:\n",
  pronotum_output_file
)