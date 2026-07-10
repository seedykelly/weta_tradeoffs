library(tidyverse)

source_dir <- "/Users/uqam/Documents/Research_Admin/research projects/Trait compensation/data/clean/measurements"

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
      sex == "m" & head_length < 18.50579 ~ "eighth",
      sex == "m" & head_length > 24.15225 ~ "tenth",
      sex == "m"                        ~ "ninth"
    )
  ) |>
  relocate(morph, .after = sex)


# ------------------------------------------------------------
# Save final combined dataset
# ------------------------------------------------------------

write_csv(
  trait_data,
  file.path(dirname(source_dir), "trait_data.csv"),
  na = ""
)


