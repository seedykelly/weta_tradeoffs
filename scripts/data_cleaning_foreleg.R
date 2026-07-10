# ============================================================
# Import and clean left and right foreleg measurements
# ============================================================

library(tidyverse)


# ------------------------------------------------------------
# 1. File locations
# ------------------------------------------------------------

raw_dir <- paste0(
  "/Users/uqam/Documents/Research_Admin/research projects/",
  "Trait compensation/data/raw/measurements"
)

# This creates:
# .../Trait compensation/data/clean/measurements

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
# 2. Function to import and clean one side of the foreleg
# ------------------------------------------------------------

clean_foreleg <- function(side, raw_dir, clean_dir) {
  
  side <- match.arg(
    side,
    choices = c("left", "right")
  )
  
  # Trait names for the selected side
  trait_names <- paste0(
    side,
    c("_forefemur", "_foretibia", "_ear")
  )
  
  # File-name pattern
  file_pattern <- paste0(
    "^first_leg_",
    side,
    "_[0-9]+\\.csv$"
  )
  
  
  # ----------------------------------------------------------
  # Find files
  # ----------------------------------------------------------
  
  files <- list.files(
    path = raw_dir,
    pattern = file_pattern,
    full.names = TRUE
  ) |>
    str_sort(numeric = TRUE)
  
  if (length(files) == 0) {
    stop(
      "No ",
      side,
      " foreleg files were found in:\n",
      raw_dir
    )
  }
  
  message(
    "\nFound ",
    length(files),
    " ",
    side,
    " foreleg files:\n",
    paste(basename(files), collapse = "\n")
  )
  
  
  # ----------------------------------------------------------
  # Import one CSV file
  # ----------------------------------------------------------
  
  read_foreleg_file <- function(file) {
    
    dat <- read_csv(
      file,
      na = c("", "NA", "NaN", "-", "N/A"),
      show_col_types = FALSE
    )
    
    # Check that the required Fiji columns are present
    required_columns <- c("Label", "Area", "Length")
    
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
      select(Label, Area, Length) |>
      mutate(
        source_file = basename(file),
        source_row = row_number(),
        
        # Extract everything before the first hyphen.
        #
        # Example:
        # ACGN-1st leg-left.jpg:1908-3095
        # becomes:
        # ACGN
        ID = str_extract(Label, "^[^-]+"),
        
        ID = str_to_upper(str_trim(ID)),
        
        Area = as.numeric(Area),
        Length = as.numeric(Length)
      ) |>
      relocate(
        source_file,
        source_row,
        ID,
        Label
      )
  }
  
  
  # ----------------------------------------------------------
  # Import and combine all files for this side
  # ----------------------------------------------------------
  
  raw_data <- map_dfr(
    files,
    read_foreleg_file
  )
  
  
  # ----------------------------------------------------------
  # Check for obsolete sex-prefixed labels
  # ----------------------------------------------------------
  #
  # These labels should no longer occur. Stopping here prevents
  # IDs such as FEMALE or MALE from being interpreted as the ID.
  # ----------------------------------------------------------
  
  legacy_labels <- raw_data |>
    filter(
      str_detect(
        Label,
        regex("^(FEMALE|MALE)-", ignore_case = TRUE)
      )
    ) |>
    distinct(source_file, Label)
  
  if (nrow(legacy_labels) > 0) {
    
    print(legacy_labels)
    
    stop(
      "\nLegacy FEMALE- or MALE-prefixed labels were found.\n",
      "Replace or remove the affected files before continuing."
    )
  }
  
  
  # ----------------------------------------------------------
  # Check that every row produced an ID
  # ----------------------------------------------------------
  
  invalid_labels <- raw_data |>
    filter(
      is.na(ID) |
        ID == ""
    )
  
  if (nrow(invalid_labels) > 0) {
    
    print(invalid_labels)
    
    stop(
      "Some labels could not be converted to individual IDs."
    )
  }
  
  
  # ----------------------------------------------------------
  # Check number of measurements per individual
  # ----------------------------------------------------------
  #
  # Each individual should have exactly three rows:
  #
  # 1. Femur length
  # 2. Tibia length
  # 3. Ear area
  #
  # A value can be NA, but its measurement row must still exist.
  # ----------------------------------------------------------
  
  measurement_counts <- raw_data |>
    count(
      source_file,
      ID,
      name = "n_measurements"
    )
  
  incorrect_counts <- measurement_counts |>
    filter(n_measurements != 3)
  
  if (nrow(incorrect_counts) > 0) {
    
    print(incorrect_counts)
    
    stop(
      "\nSome individuals do not have exactly three rows.\n",
      "Missing values are acceptable, but each measurement row ",
      "must still be present."
    )
  }
  
  
  # ----------------------------------------------------------
  # Assign measurement rows to traits
  # ----------------------------------------------------------
  
  long_data <- raw_data |>
    group_by(source_file, ID) |>
    arrange(source_row, .by_group = TRUE) |>
    mutate(
      measurement_number = row_number(),
      
      trait = trait_names[measurement_number],
      
      # Rows 1 and 2 are lengths.
      # Row 3 is ear area.
      value = if_else(
        measurement_number == 3,
        Area,
        Length
      )
    ) |>
    ungroup()
  
  
  # ----------------------------------------------------------
  # Check whether an ID occurs in multiple files
  # ----------------------------------------------------------
  
  duplicate_IDs <- long_data |>
    distinct(source_file, ID) |>
    count(
      ID,
      name = "n_files"
    ) |>
    filter(n_files > 1)
  
  if (nrow(duplicate_IDs) > 0) {
    
    print(duplicate_IDs)
    
    stop(
      "\nSome IDs occur in more than one ",
      side,
      " foreleg file."
    )
  }
  
  
  # ----------------------------------------------------------
  # Convert to one row per individual
  # ----------------------------------------------------------
  
  clean_data <- long_data |>
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
      all_of(trait_names)
    ) |>
    arrange(ID)
  
  
  # ----------------------------------------------------------
  # Summarize missing data
  # ----------------------------------------------------------
  
  missing_summary <- clean_data |>
    summarise(
      n_individuals = n(),
      
      across(
        all_of(trait_names),
        ~ sum(is.na(.x)),
        .names = "missing_{.col}"
      ),
      
      complete_individuals = sum(
        complete.cases(
          across(all_of(trait_names))
        )
      )
    )
  
  message(
    "\nMissing-data summary for the ",
    side,
    " foreleg:"
  )
  
  print(missing_summary)
  
  
  # ----------------------------------------------------------
  # Save the cleaned CSV
  # ----------------------------------------------------------
  
  output_file <- file.path(
    clean_dir,
    paste0(
      "first_leg_",
      side,
      "_clean.csv"
    )
  )
  
  write_csv(
    clean_data,
    output_file,
    na = ""
  )
  
  message(
    "\nCleaned ",
    side,
    " foreleg data saved to:\n",
    output_file
  )
  
  
  # Return the cleaned data to R
  clean_data
}


# ============================================================
# 3. Process the left foreleg
# ============================================================

left_foreleg <- clean_foreleg(
  side = "left",
  raw_dir = raw_dir,
  clean_dir = clean_dir
)


# Resulting columns:
#
# ID
# left_forefemur
# left_foretibia
# left_ear

glimpse(left_foreleg)

print(
  left_foreleg,
  n = 20
)


# ============================================================
# 4. Process the right foreleg
# ============================================================

right_foreleg <- clean_foreleg(
  side = "right",
  raw_dir = raw_dir,
  clean_dir = clean_dir
)


# Resulting columns:
#
# ID
# right_forefemur
# right_foretibia
# right_ear

glimpse(right_foreleg)

print(
  right_foreleg,
  n = 20
)


# ============================================================
# 5. List individuals with missing measurements
# ============================================================

left_foreleg_missing <- left_foreleg |>
  filter(
    if_any(
      -ID,
      is.na
    )
  )

right_foreleg_missing <- right_foreleg |>
  filter(
    if_any(
      -ID,
      is.na
    )
  )

left_foreleg_missing
right_foreleg_missing


# ============================================================
# 6. Combine left and right forelegs
# ============================================================
#
# full_join() retains an individual even if it was measured
# on only one side.
# ============================================================

forelegs <- left_foreleg |>
  full_join(
    right_foreleg,
    by = "ID"
  ) |>
  arrange(ID)

glimpse(forelegs)

print(
  forelegs,
  n = 20
)


# ------------------------------------------------------------
# Save the combined foreleg dataset
# ------------------------------------------------------------

write_csv(
  forelegs,
  file.path(
    clean_dir,
    "forelegs_clean.csv"
  ),
  na = ""
)