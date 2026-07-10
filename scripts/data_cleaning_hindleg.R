# ============================================================
# Import and clean hindleg measurements
# Files: hindleg_1.csv, hindleg_2.csv, etc.
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
# 2. Find all hindleg files
# ------------------------------------------------------------

hindleg_files <- list.files(
  path = raw_dir,
  pattern = "^hindleg_[0-9]+\\.csv$",
  full.names = TRUE
) |>
  str_sort(numeric = TRUE)

if (length(hindleg_files) == 0) {
  stop(
    "No files matching 'hindleg_[number].csv' were found in:\n",
    raw_dir
  )
}

message(
  "Found ",
  length(hindleg_files),
  " hindleg files:\n",
  paste(basename(hindleg_files), collapse = "\n")
)


# ------------------------------------------------------------
# 3. Function to import one hindleg file
# ------------------------------------------------------------

read_hindleg_file <- function(file) {
  
  dat <- read_csv(
    file,
    na = c("", "NA", "NaN", "-", "N/A"),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  required_columns <- c("Label", "Length")
  
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
    select(Label, Length) |>
    mutate(
      source_file = basename(file),
      source_row = row_number(),
      
      # Extract the individual ID from the beginning of Label.
      #
      # ACGN-3rd leg-left.jpg:2531-2079
      # becomes ACGN
      ID = str_extract(Label, "^[^-]+"),
      ID = str_to_upper(str_trim(ID)),
      
      # Extract side directly from Label
      side = case_when(
        str_detect(
          Label,
          regex(
            "3rd\\s+legs?\\s*-\\s*left",
            ignore_case = TRUE
          )
        ) ~ "left",
        
        str_detect(
          Label,
          regex(
            "3rd\\s+legs?\\s*-\\s*right",
            ignore_case = TRUE
          )
        ) ~ "right",
        
        TRUE ~ NA_character_
      ),
      
      # Convert the Fiji Length column to numeric
      Length = parse_double(
        Length,
        na = c("", "NA", "NaN", "-", "N/A")
      )
    ) |>
    relocate(
      source_file,
      source_row,
      ID,
      side,
      Label
    )
}


# ------------------------------------------------------------
# 4. Import and combine all hindleg files
# ------------------------------------------------------------

hindleg_raw <- map_dfr(
  hindleg_files,
  read_hindleg_file
)

glimpse(hindleg_raw)


# ------------------------------------------------------------
# 5. Check ID and side extraction
# ------------------------------------------------------------

invalid_hindleg_labels <- hindleg_raw |>
  filter(
    is.na(ID) |
      ID == "" |
      is.na(side)
  )

if (nrow(invalid_hindleg_labels) > 0) {
  
  print(invalid_hindleg_labels)
  
  stop(
    "Some labels could not be converted to an ID and side. ",
    "Examine 'invalid_hindleg_labels'."
  )
}


# ------------------------------------------------------------
# 6. Check whether an ID-side occurs in multiple files
# ------------------------------------------------------------

duplicate_hindleg_sides <- hindleg_raw |>
  distinct(
    source_file,
    ID,
    side
  ) |>
  count(
    ID,
    side,
    name = "n_files"
  ) |>
  filter(n_files > 1)

if (nrow(duplicate_hindleg_sides) > 0) {
  
  print(duplicate_hindleg_sides)
  
  stop(
    "Some ID-side combinations occur in more than one file. ",
    "Examine 'duplicate_hindleg_sides'."
  )
}


# ------------------------------------------------------------
# 7. Count measurements for each available leg
# ------------------------------------------------------------
#
# Normally:
#
# Row 1 = femur length
# Row 2 = tibia length
#
# An entire side can be absent.
#
# A side with only one row is retained. Its first measurement
# is interpreted as femur length and tibia is recorded as NA.
# ------------------------------------------------------------

hindleg_counts <- hindleg_raw |>
  count(
    source_file,
    ID,
    side,
    name = "n_measurements"
  )

# More than two rows cannot be interpreted safely

excess_hindleg_rows <- hindleg_counts |>
  filter(n_measurements > 2)

if (nrow(excess_hindleg_rows) > 0) {
  
  print(excess_hindleg_rows)
  
  stop(
    "Some hindlegs have more than two measurement rows. ",
    "Examine 'excess_hindleg_rows'."
  )
}


# Legs containing only one measurement

incomplete_hindlegs <- hindleg_counts |>
  filter(n_measurements == 1)

if (nrow(incomplete_hindlegs) > 0) {
  
  message(
    "\nThe following legs contain only one measurement.\n",
    "The measurement will be treated as femur length, and ",
    "tibia length will be set to NA."
  )
  
  print(incomplete_hindlegs)
}


# ------------------------------------------------------------
# 8. Assign femur and tibia measurements
# ------------------------------------------------------------

hindleg_long <- hindleg_raw |>
  group_by(
    source_file,
    ID,
    side
  ) |>
  arrange(
    source_row,
    .by_group = TRUE
  ) |>
  mutate(
    measurement_number = row_number(),
    
    trait = case_when(
      side == "left" &
        measurement_number == 1 ~ "left_hindfemur",
      
      side == "left" &
        measurement_number == 2 ~ "left_hindtibia",
      
      side == "right" &
        measurement_number == 1 ~ "right_hindfemur",
      
      side == "right" &
        measurement_number == 2 ~ "right_hindtibia"
    ),
    
    value = Length
  ) |>
  ungroup()


# ------------------------------------------------------------
# 9. Convert to one row per individual
# ------------------------------------------------------------

hindleg_traits <- c(
  "left_hindfemur",
  "left_hindtibia",
  "right_hindfemur",
  "right_hindtibia"
)

hindlegs <- hindleg_long |>
  select(
    ID,
    trait,
    value
  ) |>
  pivot_wider(
    names_from = trait,
    values_from = value
  )


# Add any trait column that is completely absent from all files

missing_trait_columns <- setdiff(
  hindleg_traits,
  names(hindlegs)
)

if (length(missing_trait_columns) > 0) {
  
  for (column_name in missing_trait_columns) {
    hindlegs[[column_name]] <- NA_real_
  }
}


# Arrange columns in the requested order

hindlegs <- hindlegs |>
  select(
    ID,
    all_of(hindleg_traits)
  ) |>
  arrange(ID)


# View the cleaned data

glimpse(hindlegs)

print(
  hindlegs,
  n = 20
)


# ------------------------------------------------------------
# 10. Identify individuals missing an entire side
# ------------------------------------------------------------

hindleg_side_presence <- hindleg_raw |>
  distinct(
    ID,
    side
  ) |>
  mutate(present = TRUE) |>
  complete(
    ID,
    side = c("left", "right"),
    fill = list(present = FALSE)
  ) |>
  pivot_wider(
    names_from = side,
    values_from = present,
    names_prefix = "has_"
  ) |>
  arrange(ID)

hindleg_missing_sides <- hindleg_side_presence |>
  filter(
    !has_left |
      !has_right
  )

print(hindleg_missing_sides)


# ------------------------------------------------------------
# 11. Summarize missing measurements
# ------------------------------------------------------------

hindleg_missing_summary <- hindlegs |>
  summarise(
    n_individuals = n(),
    
    missing_left_hindfemur  = sum(is.na(left_hindfemur)),
    missing_left_hindtibia  = sum(is.na(left_hindtibia)),
    missing_right_hindfemur = sum(is.na(right_hindfemur)),
    missing_right_hindtibia = sum(is.na(right_hindtibia)),
    
    complete_individuals = sum(
      complete.cases(
        left_hindfemur,
        left_hindtibia,
        right_hindfemur,
        right_hindtibia
      )
    )
  )

print(hindleg_missing_summary)


# Individuals with at least one missing hindleg measurement

hindleg_missing_measurements <- hindlegs |>
  filter(
    if_any(
      -ID,
      is.na
    )
  )

print(hindleg_missing_measurements)


# ------------------------------------------------------------
# 12. Optional range check
# ------------------------------------------------------------
#
# Values are flagged but are not removed or changed.
# Adjust the maximum if necessary.
# ------------------------------------------------------------

hindleg_suspect_values <- hindlegs |>
  filter(
    if_any(
      -ID,
      ~ !is.na(.x) & (.x <= 0 | .x > 40)
    )
  )

print(hindleg_suspect_values)


# ------------------------------------------------------------
# 13. Save cleaned hindleg data
# ------------------------------------------------------------

hindleg_output_file <- file.path(
  clean_dir,
  "hindlegs_clean.csv"
)

write_csv(
  hindlegs,
  hindleg_output_file,
  na = ""
)

message(
  "Cleaned hindleg data saved to:\n",
  hindleg_output_file
)