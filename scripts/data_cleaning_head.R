# ============================================================
# Import and clean head measurements
# Files: head_1.csv, head_2.csv, etc.
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
# 2. Find all head files
# ------------------------------------------------------------

head_files <- list.files(
  path = raw_dir,
  pattern = "^head_[0-9]+\\.csv$",
  full.names = TRUE
) |>
  str_sort(numeric = TRUE)

if (length(head_files) == 0) {
  stop(
    "No files matching 'head_[number].csv' were found in:\n",
    raw_dir
  )
}

message(
  "Found ",
  length(head_files),
  " head files:\n",
  paste(basename(head_files), collapse = "\n")
)


# ------------------------------------------------------------
# 3. Function to import one head file
# ------------------------------------------------------------

read_head_file <- function(file) {
  
  dat <- read_csv(
    file,
    na = c("", "NA", "NaN", "-", "N/A"),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "Label",
    "Length"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "The following required columns are missing from ",
      basename(file),
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  dat |>
    select(
      Label,
      Length
    ) |>
    mutate(
      source_file = basename(file),
      source_row = row_number(),
      
      # Example:
      # ACGN-Head.jpg:1698-2669
      # becomes:
      # ACGN
      ID = str_extract(
        Label,
        "^[^-]+"
      ),
      
      ID = str_to_upper(
        str_trim(ID)
      ),
      
      Length = parse_double(
        Length,
        na = c("", "NA", "NaN", "-", "N/A")
      )
    ) |>
    relocate(
      source_file,
      source_row,
      ID,
      Label
    )
}


# ------------------------------------------------------------
# 4. Import and combine the head files
# ------------------------------------------------------------

head_raw <- map_dfr(
  head_files,
  read_head_file
)

glimpse(head_raw)


# ------------------------------------------------------------
# 5. Check that every label produced a valid ID
# ------------------------------------------------------------

invalid_head_labels <- head_raw |>
  filter(
    is.na(ID) |
      ID == ""
  )

if (nrow(invalid_head_labels) > 0) {
  
  print(invalid_head_labels)
  
  stop(
    "Some labels could not be converted to IDs. ",
    "Examine 'invalid_head_labels'."
  )
}


# ------------------------------------------------------------
# 6. Check the number of rows per individual within each file
# ------------------------------------------------------------
#
# The expected order is:
#
# 1. Head length
# 2. Head width
# 3. Left eye length
# 4. Right eye length
#
# A measurement can be NA, provided that its row still exists.
# ------------------------------------------------------------

head_counts <- head_raw |>
  count(
    source_file,
    ID,
    name = "n_measurements"
  )

incorrect_head_counts <- head_counts |>
  filter(
    n_measurements != 4
  )

if (nrow(incorrect_head_counts) > 0) {
  
  print(incorrect_head_counts)
  
  problematic_head_rows <- head_raw |>
    semi_join(
      incorrect_head_counts,
      by = c("source_file", "ID")
    ) |>
    arrange(
      source_file,
      ID,
      source_row
    )
  
  print(problematic_head_rows)
  
  stop(
    "\nSome individuals do not have exactly four rows.\n",
    "Examine 'incorrect_head_counts' and ",
    "'problematic_head_rows' before continuing."
  )
}


# ------------------------------------------------------------
# 7. Check whether an ID occurs in more than one head file
# ------------------------------------------------------------

duplicate_head_IDs <- head_raw |>
  distinct(
    source_file,
    ID
  ) |>
  count(
    ID,
    name = "n_files"
  ) |>
  filter(
    n_files > 1
  )

if (nrow(duplicate_head_IDs) > 0) {
  
  print(duplicate_head_IDs)
  
  duplicate_head_rows <- head_raw |>
    semi_join(
      duplicate_head_IDs,
      by = "ID"
    ) |>
    arrange(
      ID,
      source_file,
      source_row
    )
  
  print(duplicate_head_rows)
  
  stop(
    "\nSome individuals occur in more than one head file.\n",
    "Examine 'duplicate_head_IDs' and ",
    "'duplicate_head_rows' before continuing."
  )
}


# ------------------------------------------------------------
# 8. Assign the four measurements to traits
# ------------------------------------------------------------

head_long <- head_raw |>
  group_by(
    source_file,
    ID
  ) |>
  arrange(
    source_row,
    .by_group = TRUE
  ) |>
  mutate(
    measurement_number = row_number(),
    
    trait = case_when(
      measurement_number == 1 ~ "head_length",
      measurement_number == 2 ~ "head_width",
      measurement_number == 3 ~ "left_eye_length",
      measurement_number == 4 ~ "right_eye_length"
    ),
    
    value = Length
  ) |>
  ungroup()


# ------------------------------------------------------------
# 9. Convert to one row per individual
# ------------------------------------------------------------

head_measurements <- head_long |>
  select(
    ID,
    trait,
    value
  ) |>
  pivot_wider(
    names_from = trait,
    values_from = value
  ) |>
  select(
    ID,
    head_length,
    head_width,
    left_eye_length,
    right_eye_length
  ) |>
  arrange(ID)


# View cleaned data

glimpse(head_measurements)

print(
  head_measurements,
  n = 20
)


# ------------------------------------------------------------
# 10. Summarize missing measurements
# ------------------------------------------------------------

head_missing_summary <- head_measurements |>
  summarise(
    n_individuals = n(),
    
    missing_head_length = sum(
      is.na(head_length)
    ),
    
    missing_head_width = sum(
      is.na(head_width)
    ),
    
    missing_left_eye_length = sum(
      is.na(left_eye_length)
    ),
    
    missing_right_eye_length = sum(
      is.na(right_eye_length)
    ),
    
    complete_individuals = sum(
      complete.cases(
        head_length,
        head_width,
        left_eye_length,
        right_eye_length
      )
    )
  )

print(head_missing_summary)


# Individuals with at least one missing measurement

head_missing_measurements <- head_measurements |>
  filter(
    if_any(
      -ID,
      is.na
    )
  )

print(head_missing_measurements)


# ------------------------------------------------------------
# 11. Optional range checks
# ------------------------------------------------------------
#
# These observations are only flagged. They are not changed
# or removed.
# ------------------------------------------------------------

head_suspect_values <- head_measurements |>
  filter(
    (
      !is.na(head_length) &
        (head_length <= 0 | head_length > 25)
    ) |
      (
        !is.na(head_width) &
          (head_width <= 0 | head_width > 20)
      ) |
      (
        !is.na(left_eye_length) &
          (left_eye_length <= 0 | left_eye_length > 10)
      ) |
      (
        !is.na(right_eye_length) &
          (right_eye_length <= 0 | right_eye_length > 10)
      )
  )

print(head_suspect_values)


# ------------------------------------------------------------
# 12. Save cleaned head measurements
# ------------------------------------------------------------

head_output_file <- file.path(
  clean_dir,
  "head_clean.csv"
)

write_csv(
  head_measurements,
  head_output_file,
  na = ""
)

message(
  "Cleaned head data saved to:\n",
  head_output_file
)