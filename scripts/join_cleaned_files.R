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
