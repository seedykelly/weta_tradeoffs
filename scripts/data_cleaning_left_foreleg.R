# ============================================================
# Import and clean left foreleg measurements
# ============================================================

library(tidyverse)

# Folder containing the raw Fiji measurement files
raw_dir <- paste0(
  "/Users/uqam/Documents/Research_Admin/research projects/",
  "Trait compensation/data/raw/measurements"
)


# ------------------------------------------------------------
# 1. Find the left-foreleg CSV files
# ------------------------------------------------------------

left_foreleg_files <- list.files(
  path = raw_dir,
  pattern = "^first_leg_left_[0-9]+\\.csv$",
  full.names = TRUE
) |>
  stringr::str_sort(numeric = TRUE)

if (length(left_foreleg_files) == 0) {
  stop(
    "No files matching 'first_leg_left_[number].csv' ",
    "were found in raw_dir."
  )
}

message(
  "Found ",
  length(left_foreleg_files),
  " left-foreleg files:\n",
  paste(basename(left_foreleg_files), collapse = "\n")
)


# ------------------------------------------------------------
# 2. Function for importing one Fiji CSV file
# ------------------------------------------------------------

read_left_foreleg_file <- function(file) {
  
  read_csv(
    file,
    na = c("", "NA", "NaN", "-", "N/A"),
    show_col_types = FALSE
  ) |>
    select(Label, Area, Length) |>
    mutate(
      source_file = basename(file),
      source_row = row_number(),
      
      # Record sex where it appears in the filename
      sex = case_when(
        str_detect(Label, regex("^FEMALE-", ignore_case = TRUE)) ~ "female",
        str_detect(Label, regex("^MALE-", ignore_case = TRUE))   ~ "male",
        TRUE                                                    ~ NA_character_
      ),
      
      # Extract a unique individual identifier.
      #
      # Examples:
      # ACGN-1st leg-left.jpg      -> ACGN
      # SBJVC-1st leg-left.jpg     -> SBJVC
      # FEMALE-BF001-1st.jpg       -> FEMALE-BF001
      # MALE-BF001-1st left.jpg    -> MALE-BF001
      #
      # The sex prefix is retained when present because female
      # and male files can contain the same code, such as BF001.
      ID = case_when(
        str_detect(
          Label,
          regex("^(FEMALE|MALE)-", ignore_case = TRUE)
        ) ~ str_extract(
          Label,
          regex("^(FEMALE|MALE)-[^-]+", ignore_case = TRUE)
        ),
        
        TRUE ~ str_extract(Label, "^[^-]+")
      ),
      
      ID = str_to_upper(ID),
      
      # Ensure that measurement columns are numeric
      Area = as.numeric(Area),
      Length = as.numeric(Length)
    ) |>
    relocate(source_file, source_row, ID, sex, Label)
}


# ------------------------------------------------------------
# 3. Import and combine all left-foreleg files
# ------------------------------------------------------------

left_foreleg_raw <- map_dfr(
  left_foreleg_files,
  read_left_foreleg_file
)


# Examine the combined raw data if desired
glimpse(left_foreleg_raw)


# ------------------------------------------------------------
# 4. Check that every label produced an ID
# ------------------------------------------------------------

invalid_labels <- left_foreleg_raw |>
  filter(is.na(ID) | ID == "")

if (nrow(invalid_labels) > 0) {
  
  print(invalid_labels)
  
  stop(
    "Some labels could not be converted to individual IDs. ",
    "Examine the object 'invalid_labels'."
  )
}


# ------------------------------------------------------------
# 5. Check the number of rows per individual
# ------------------------------------------------------------
#
# Each left foreleg should have three rows, in this order:
#
# 1. Femur length
# 2. Tibia length
# 3. Ear area
#
# An NA measurement is acceptable as long as its row is present.
# This check identifies cases in which an entire row is absent or
# an individual was measured more than three times.
# ------------------------------------------------------------

left_foreleg_counts <- left_foreleg_raw |>
  count(
    source_file,
    ID,
    name = "n_measurements"
  )

incorrect_counts <- left_foreleg_counts |>
  filter(n_measurements != 3)

if (nrow(incorrect_counts) > 0) {
  
  print(incorrect_counts)
  
  stop(
    "Some individuals do not have exactly three measurement rows. ",
    "Examine the object 'incorrect_counts'."
  )
}


# ------------------------------------------------------------
# 6. Assign each row to its corresponding trait
# ------------------------------------------------------------

left_foreleg_long <- left_foreleg_raw |>
  group_by(source_file, ID) |>
  arrange(source_row, .by_group = TRUE) |>
  mutate(
    measurement_number = row_number(),
    
    trait = case_when(
      measurement_number == 1 ~ "left_forefemur",
      measurement_number == 2 ~ "left_foretibia",
      measurement_number == 3 ~ "left_ear"
    ),
    
    # Femur and tibia are stored in the Length column.
    # Ear area is stored in the Area column.
    value = case_when(
      trait == "left_ear" ~ Area,
      TRUE                ~ Length
    )
  ) |>
  ungroup()


# ------------------------------------------------------------
# 7. Check whether an ID occurs in more than one file
# ------------------------------------------------------------

duplicate_IDs <- left_foreleg_long |>
  distinct(source_file, ID) |>
  count(ID, name = "n_files") |>
  filter(n_files > 1)

if (nrow(duplicate_IDs) > 0) {
  
  print(duplicate_IDs)
  
  stop(
    "Some IDs occur in more than one left-foreleg file. ",
    "Examine the object 'duplicate_IDs'."
  )
}


# ------------------------------------------------------------
# 8. Convert to one row per individual
# ------------------------------------------------------------

left_foreleg <- left_foreleg_long |>
  select(ID, sex, trait, value) |>
  pivot_wider(
    names_from = trait,
    values_from = value
  ) |>
  select(
    ID,
    sex,
    left_forefemur,
    left_foretibia,
    left_ear
  ) |>
  arrange(ID)


# View the cleaned data
print(left_foreleg, n = 20)

glimpse(left_foreleg)


# ------------------------------------------------------------
# 9. Missing-data checks
# ------------------------------------------------------------
#
# These checks report missing values but do not remove individuals.
# ------------------------------------------------------------

left_foreleg_missing_summary <- left_foreleg |>
  summarise(
    n_individuals = n(),
    
    missing_left_forefemur = sum(is.na(left_forefemur)),
    missing_left_foretibia = sum(is.na(left_foretibia)),
    missing_left_ear       = sum(is.na(left_ear)),
    
    complete_individuals = sum(
      complete.cases(
        left_forefemur,
        left_foretibia,
        left_ear
      )
    )
  )

print(left_foreleg_missing_summary)


# Individuals with at least one missing measurement
left_foreleg_missing <- left_foreleg |>
  filter(
    if_any(
      c(
        left_forefemur,
        left_foretibia,
        left_ear
      ),
      is.na
    )
  )

print(left_foreleg_missing)


# ------------------------------------------------------------
# 10. Optional range checks
# ------------------------------------------------------------
#
# These do not delete values. They simply make it easier to find
# possible data-entry or measurement errors.
# Adjust the limits if necessary for your species.
# ------------------------------------------------------------

left_foreleg_suspect_values <- left_foreleg |>
  filter(
    (!is.na(left_forefemur) &
       (left_forefemur <= 0 | left_forefemur > 30)) |
      
      (!is.na(left_foretibia) &
         (left_foretibia <= 0 | left_foretibia > 30)) |
      
      (!is.na(left_ear) &
         (left_ear <= 0 | left_ear > 20))
  )

print(left_foreleg_suspect_values)


# ------------------------------------------------------------
# 11. Save the cleaned dataset
# ------------------------------------------------------------

clean_dir <- file.path(
  dirname(raw_dir),
  "clean"
)

dir.create(
  clean_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

left_foreleg_output_file <- file.path(
  clean_dir,
  "first_leg_left_clean.csv"
)

write_csv(
  left_foreleg,
  left_foreleg_output_file,
  na = ""
)

message(
  "Cleaned left-foreleg data saved to:\n",
  left_foreleg_output_file
)
