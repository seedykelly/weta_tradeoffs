# ============================================================
# Import and clean midleg measurements
# Files: second_leg_1.csv, second_leg_2.csv, etc.
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
# 2. Find all midleg files
# ------------------------------------------------------------

midleg_files <- list.files(
  path = raw_dir,
  pattern = "^second_leg_[0-9]+\\.csv$",
  full.names = TRUE
) |>
  str_sort(numeric = TRUE)

if (length(midleg_files) == 0) {
  stop(
    "No files matching 'second_leg_[number].csv' were found in:\n",
    raw_dir
  )
}

message(
  "Found ",
  length(midleg_files),
  " midleg files:\n",
  paste(basename(midleg_files), collapse = "\n")
)


# ------------------------------------------------------------
# 3. Function for importing one midleg CSV
# ------------------------------------------------------------

read_midleg_file <- function(file) {
  
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
      
      # Extract the individual identifier from the beginning
      # of the label.
      #
      # Example:
      # ACIU-2nd leg-left.jpg:1454-2285
      # becomes:
      # ACIU
      ID = str_extract(Label, "^[^-]+"),
      ID = str_to_upper(str_trim(ID)),
      
      # Extract side from the label.
      #
      # This handles:
      # "2nd leg-left"
      # "2nd leg- left"
      # "2nd legs-left"
      side = case_when(
        str_detect(
          Label,
          regex(
            "2nd\\s+legs?\\s*-\\s*left",
            ignore_case = TRUE
          )
        ) ~ "left",
        
        str_detect(
          Label,
          regex(
            "2nd\\s+legs?\\s*-\\s*right",
            ignore_case = TRUE
          )
        ) ~ "right",
        
        TRUE ~ NA_character_
      ),
      
      # Convert the measurement to numeric
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
# 4. Import and combine all midleg files
# ------------------------------------------------------------

midleg_raw <- map_dfr(
  midleg_files,
  read_midleg_file
)

glimpse(midleg_raw)


# ------------------------------------------------------------
# 5. Check that every label produced an ID and side
# ------------------------------------------------------------

invalid_midleg_labels <- midleg_raw |>
  filter(
    is.na(ID) |
      ID == "" |
      is.na(side)
  )

if (nrow(invalid_midleg_labels) > 0) {
  
  print(invalid_midleg_labels)
  
  stop(
    "Some labels could not be converted to an ID and leg side. ",
    "Examine 'invalid_midleg_labels'."
  )
}


# ------------------------------------------------------------
# 6. Check whether an ID-side combination occurs in
#    more than one file
# ------------------------------------------------------------
#
# An ID can have a left and right leg, but the same side should
# not occur in multiple files.
# ------------------------------------------------------------

duplicate_midleg_sides <- midleg_raw |>
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

if (nrow(duplicate_midleg_sides) > 0) {
  
  print(duplicate_midleg_sides)
  
  stop(
    "Some ID-side combinations occur in more than one file. ",
    "Examine 'duplicate_midleg_sides'."
  )
}


# ------------------------------------------------------------
# 7. Check the number of rows per leg
# ------------------------------------------------------------
#
# Each available leg should have two rows:
#
# 1. Femur length
# 2. Tibia length
#
# An entire left or right leg may be absent. That is allowed.
#
# A measurement value may also be NA, but the row must still
# be present.
# ------------------------------------------------------------

midleg_counts <- midleg_raw |>
  count(
    source_file,
    ID,
    side,
    name = "n_measurements"
  )

incorrect_midleg_counts <- midleg_counts |>
  filter(n_measurements != 2)

if (nrow(incorrect_midleg_counts) > 0) {
  
  print(incorrect_midleg_counts)
  
  stop(
    "Some available legs do not have exactly two measurement rows. ",
    "Examine 'incorrect_midleg_counts'."
  )
}


# ------------------------------------------------------------
# 8. Assign femur and tibia measurements
# ------------------------------------------------------------
#
# Within each ID and side:
#
# Row 1 = femur
# Row 2 = tibia
# ------------------------------------------------------------

midleg_long <- midleg_raw |>
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
        measurement_number == 1 ~ "left_midfemur",
      
      side == "left" &
        measurement_number == 2 ~ "left_midtibia",
      
      side == "right" &
        measurement_number == 1 ~ "right_midfemur",
      
      side == "right" &
        measurement_number == 2 ~ "right_midtibia"
    ),
    
    value = Length
  ) |>
  ungroup()


# ------------------------------------------------------------
# 9. Convert to one row per individual
# ------------------------------------------------------------

midleg_traits <- c(
  "left_midfemur",
  "left_midtibia",
  "right_midfemur",
  "right_midtibia"
)

midlegs <- midleg_long |>
  select(
    ID,
    trait,
    value
  ) |>
  pivot_wider(
    names_from = trait,
    values_from = value
  )


# Add any trait columns that might be completely absent
# from the entire collection of files.

missing_trait_columns <- setdiff(
  midleg_traits,
  names(midlegs)
)

if (length(missing_trait_columns) > 0) {
  
  for (column_name in missing_trait_columns) {
    midlegs[[column_name]] <- NA_real_
  }
}


# Put columns in the requested order

midlegs <- midlegs |>
  select(
    ID,
    all_of(midleg_traits)
  ) |>
  arrange(ID)


# View the cleaned data

glimpse(midlegs)

print(
  midlegs,
  n = 20
)


# ------------------------------------------------------------
# 10. Identify individuals missing an entire leg
# ------------------------------------------------------------
#
# This uses the presence of labels rather than the numerical
# measurements. Therefore, a leg with two rows containing NA
# is distinguished from a leg that was not photographed or
# measured at all.
# ------------------------------------------------------------

midleg_side_presence <- midleg_raw |>
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

midleg_missing_sides <- midleg_side_presence |>
  filter(
    !has_left |
      !has_right
  )

print(midleg_missing_sides)


# ------------------------------------------------------------
# 11. Summarize missing measurements
# ------------------------------------------------------------

midleg_missing_summary <- midlegs |>
  summarise(
    n_individuals = n(),
    
    missing_left_midfemur  = sum(is.na(left_midfemur)),
    missing_left_midtibia  = sum(is.na(left_midtibia)),
    missing_right_midfemur = sum(is.na(right_midfemur)),
    missing_right_midtibia = sum(is.na(right_midtibia)),
    
    complete_individuals = sum(
      !is.na(left_midfemur) &
        !is.na(left_midtibia) &
        !is.na(right_midfemur) &
        !is.na(right_midtibia)
    )
  )

print(midleg_missing_summary)


# Individuals with at least one missing measurement

midleg_missing_measurements <- midlegs |>
  filter(
    if_any(
      -ID,
      is.na
    )
  )

print(midleg_missing_measurements)


# ------------------------------------------------------------
# 12. Optional range check
# ------------------------------------------------------------
#
# Values are only flagged; they are not removed or changed.
# Adjust the upper limit if necessary.
# ------------------------------------------------------------

midleg_suspect_values <- midlegs |>
  filter(
    if_any(
      -ID,
      ~ !is.na(.x) & (.x <= 0 | .x > 30)
    )
  )

print(midleg_suspect_values)


# ------------------------------------------------------------
# 13. Save the cleaned midleg data
# ------------------------------------------------------------

midleg_output_file <- file.path(
  clean_dir,
  "second_leg_clean.csv"
)

write_csv(
  midlegs,
  midleg_output_file,
  na = ""
)

message(
  "Cleaned midleg data saved to:\n",
  midleg_output_file
)