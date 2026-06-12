library(tidyverse)

# Folder containing the raw CSV files
raw_dir <- "/Users/uqam/Documents/Research_Admin/research projects/Trait compensation/data/raw"

# -------------------------------------------------------------------------
# 1. Import all CSV files and clean IDs
# -------------------------------------------------------------------------

weta_imported <- list.files(
  path = raw_dir,
  pattern = "\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(
    ~ read_csv(
      .x,
      col_types = cols(.default = col_character()),
      na = c("", "NA", "N/A", "-")
    ) %>%
      mutate(source_file = basename(.x))
  ) %>%
  mutate(
    side = case_when(
      str_detect(source_file, regex("left", ignore_case = TRUE))  ~ "left",
      str_detect(source_file, regex("right", ignore_case = TRUE)) ~ "right",
      TRUE ~ NA_character_
    ),
    
    Label_clean = Label %>%
      str_remove("-1st.*$") %>%
      str_replace(regex("^(MALE|FEMALE)-", ignore_case = TRUE), "\\1") %>%
      str_to_upper() %>%
      str_trim()
  ) %>%
  filter(
    !str_detect(Label_clean, regex("^(MALE|FEMALE)", ignore_case = TRUE))
  )

# -------------------------------------------------------------------------
# 2. Check that each individual has exactly three rows per side/file
# -------------------------------------------------------------------------

row_check <- weta_imported %>%
  count(source_file, side, Label_clean, name = "n_rows") %>%
  arrange(n_rows, source_file, Label_clean)

row_problems <- row_check %>%
  filter(n_rows != 3)

row_problems

# write_csv(
#   row_problems,
#   "/Users/uqam/Documents/Research_Admin/research projects/Trait compensation/data/row_problems.csv"
# )

# -------------------------------------------------------------------------
# 3. Label rows as femur, tibia, or ear
# -------------------------------------------------------------------------

weta_long <- weta_imported %>%
  group_by(source_file, side, Label_clean) %>%
  mutate(
    row_within_individual = row_number(),
    
    trait = case_when(
      row_within_individual == 1 ~ "femur",
      row_within_individual == 2 ~ "tibia",
      row_within_individual == 3 ~ "ear",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    Area_num   = parse_number(Area, na = c("", "NA", "N/A", "-")),
    Length_num = parse_number(Length, na = c("", "NA", "N/A", "-")),
    
    value = case_when(
      trait %in% c("femur", "tibia") ~ Length_num,
      trait == "ear" ~ Area_num,
      TRUE ~ NA_real_
    ),
    
    trait_side = paste(trait, side, sep = "_")
  ) %>%
  select(
    Label = Label_clean,
    source_file,
    side,
    row_within_individual,
    trait,
    trait_side,
    Area_num,
    Length_num,
    value
  )

# -------------------------------------------------------------------------
# 4. Convert to wide format
# -------------------------------------------------------------------------

weta_wide <- weta_long %>%
  select(Label, trait_side, value) %>%
  pivot_wider(
    names_from = trait_side,
    values_from = value
  ) %>%
  select(
    Label,
    femur_left, tibia_left, ear_left,
    femur_right, tibia_right, ear_right
  )

weta_wide

# -------------------------------------------------------------------------
# 5. Save cleaned wide-format dataset
# -------------------------------------------------------------------------

write_csv(
  weta_wide,
  "/Users/uqam/Documents/Research_Admin/research projects/Trait compensation/data/weta_first_leg_wide.csv"
)

weta_wide %>%
  summarise(
    n_individuals = n(),
    missing_femur_left  = sum(is.na(femur_left)),
    missing_tibia_left  = sum(is.na(tibia_left)),
    missing_ear_left    = sum(is.na(ear_left)),
    missing_femur_right = sum(is.na(femur_right)),
    missing_tibia_right = sum(is.na(tibia_right)),
    missing_ear_right   = sum(is.na(ear_right))
  )

weta_wide <- weta_wide %>%
  mutate(
    femur_abs_diff = abs(femur_left - femur_right),
    tibia_abs_diff = abs(tibia_left - tibia_right),
    ear_abs_diff   = abs(ear_left - ear_right)
  )

weta_wide %>%
  print(n=435)

weta_wide %>%
  summarise(
    mean_femur_abs_diff = mean(femur_abs_diff, na.rm = TRUE),
    sd_femur_abs_diff   = sd(femur_abs_diff, na.rm = TRUE),
    
    mean_tibia_abs_diff = mean(tibia_abs_diff, na.rm = TRUE),
    sd_tibia_abs_diff   = sd(tibia_abs_diff, na.rm = TRUE),
    
    mean_ear_abs_diff = mean(ear_abs_diff, na.rm = TRUE),
    sd_ear_abs_diff   = sd(ear_abs_diff, na.rm = TRUE)
  )

femur_diff_sorted <- weta_wide %>%
  arrange(desc(femur_abs_diff)) %>%
  select(Label, femur_left, femur_right, femur_abs_diff) %>%
  print(n=400)

tibia_diff_sorted <- weta_wide %>%
  arrange(desc(tibia_abs_diff)) %>%
  select(Label, tibia_left, tibia_right, tibia_abs_diff)

ear_diff_sorted <- weta_wide %>%
  arrange(desc(ear_abs_diff)) %>%
  select(Label, ear_left, ear_right, ear_abs_diff)

# weta_long %>%
#   filter(Label == "ZVNV") %>%
#   arrange(side, row_within_individual)
# 
# weta_wide %>%
#   filter(Label == "RCMK")








