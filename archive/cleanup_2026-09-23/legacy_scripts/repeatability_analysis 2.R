# ============================================================
# Measurement repeatability for weta traits
# ============================================================
#
# Three measurement rounds:
#   1 = first measurement session (from trait_data.csv / dat)
#   2 = second session (*_two.csv files)
#   3 = third session  (*_three.csv files)
#
# Bilateral traits (left/right legs, eyes, ears) are averaged
# per individual per round before repeatability is calculated,
# matching the approach in join_cleaned_files.R.
#
# Known data anomalies (flagged by messages, handled below):
#   foreleg_two:   ECJN left (2 rows, no ear), XWCJ left (1 row)
#   foreleg_three: ECJN left (2 rows), XWCJ left (1 row),
#                  PZMO right (6 rows — first 3 used)
#   midleg_two:    DOCI right (1 row), LRSS right (3 rows — first 2 used)
#   hindleg_two:   MJSX left (1 row), YZUH right (1 row)
#   hindleg_three: LAZJ right (1 row), MJSX left (1 row), YZUH right (1 row)
#   pronotum:      JUJH absent in rounds 2 and 3
# ============================================================

library(tidyverse)
library(rptR)


# ============================================================
# 1. Paths
# ============================================================

raw_dir <- paste0(
  "/Users/uqam/Documents/Research_Admin/research projects/",
  "Trait compensation/data/raw/measurements"
)

project_dir <- paste0(
  "/Users/uqam/Documents/Research_Admin/research projects/",
  "Trait compensation"
)

clean_dir <- paste0(
  "/Users/uqam/Documents/Research_Admin/research projects/",
  "Trait compensation/data/clean/measurements"
)

na_vals <- c("", "NA", "NaN", "-", "N/A")


# ============================================================
# 2. Repeatability individual list
# ============================================================

rep_metadata <- read_csv(
  file.path(project_dir, "weta_repeatability_measurement_order.csv"),
  show_col_types = FALSE
)

rep_ids <- rep_metadata$ID

message("Repeatability individuals: ", length(rep_ids))


# ============================================================
# 3. Round 1 — extract from trait_data.csv
# ============================================================

dat <- read_csv(
  file.path(dirname(clean_dir), "trait_data.csv"),
  show_col_types = FALSE
)

round1 <- dat |>
  filter(ID %in% rep_ids) |>
  select(
    ID,
    pronotum, head_length, head_width,
    forefemur, foretibia, ear,
    midfemur, midtibia,
    hindfemur, hindtibia,
    eye
  ) |>
  mutate(measurement_round = 1L)

message("Round 1 individuals found: ", nrow(round1))

missing_r1 <- setdiff(rep_ids, round1$ID)
if (length(missing_r1) > 0) {
  warning("IDs missing from round 1: ", paste(missing_r1, collapse = ", "))
}


# ============================================================
# 4. Cleaning functions for single *_two / *_three files
# ============================================================
# Each function reads one file and returns a tibble with one
# row per individual and the averaged trait value(s).
# ============================================================


# --- Pronotum ---

clean_pronotum_file <- function(path) {

  raw <- read_csv(
    path,
    na = na_vals,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) |>
    mutate(
      ID       = str_to_upper(str_trim(str_extract(Label, "^[^-]+"))),
      pronotum = parse_double(Length, na = na_vals)
    )

  dups <- raw |> count(ID) |> filter(n > 1)
  if (nrow(dups) > 0) {
    message(basename(path), ": duplicate pronotum IDs: ",
            paste(dups$ID, collapse = ", "))
  }

  raw |>
    select(ID, pronotum) |>
    arrange(ID)
}


# --- Head ---
# Row order per individual: 1 = head_length, 2 = head_width,
#                           3 = left_eye_length, 4 = right_eye_length
# eye = mean(left, right)

clean_head_file <- function(path) {

  read_csv(
    path,
    na = na_vals,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) |>
    mutate(
      source_row = row_number(),
      ID         = str_to_upper(str_trim(str_extract(Label, "^[^-]+"))),
      Length     = parse_double(Length, na = na_vals)
    ) |>
    group_by(ID) |>
    arrange(source_row, .by_group = TRUE) |>
    mutate(
      n     = row_number(),
      trait = case_when(
        n == 1 ~ "head_length",
        n == 2 ~ "head_width",
        n == 3 ~ "left_eye_length",
        n == 4 ~ "right_eye_length"
      ),
      value = Length
    ) |>
    ungroup() |>
    select(ID, trait, value) |>
    pivot_wider(names_from = trait, values_from = value) |>
    mutate(
      eye = rowMeans(
        pick(any_of(c("left_eye_length", "right_eye_length"))),
        na.rm = TRUE
      ),
      eye = if_else(is.nan(eye), NA_real_, eye)
    ) |>
    select(ID, head_length, head_width, eye) |>
    arrange(ID)
}


# --- Foreleg ---
# Combined left + right file.
# Side detected from Label: "1st leg-left" / "1st leg-right".
# Row order per side: 1 = femur length, 2 = tibia length, 3 = ear area.
# Ear uses Area column; femur and tibia use Length.
# Extra rows per ID-side beyond 3 are discarded with a message.
# Fewer than 3 rows → missing traits become NA.
# Averaged traits: forefemur, foretibia, ear.

clean_foreleg_file <- function(path) {

  raw <- read_csv(
    path,
    na = na_vals,
    show_col_types = FALSE
  ) |>
    mutate(
      source_row = row_number(),
      ID         = str_to_upper(str_trim(str_extract(Label, "^[^-]+"))),
      side       = case_when(
        str_detect(Label, regex("1st\\s+leg\\s*-\\s*left",  ignore_case = TRUE)) ~ "left",
        str_detect(Label, regex("1st\\s+leg\\s*-\\s*right", ignore_case = TRUE)) ~ "right",
        TRUE ~ NA_character_
      ),
      Area   = as.numeric(Area),
      Length = as.numeric(Length)
    )

  # Report anomalous row counts
  counts <- raw |> count(ID, side, name = "n_rows")

  extra <- counts |> filter(n_rows > 3)
  if (nrow(extra) > 0) {
    message(basename(path), ": extra foreleg rows detected (using first 3 per side):\n",
            paste(paste0("  ", extra$ID, " ", extra$side,
                         " (", extra$n_rows, " rows)"), collapse = "\n"))
  }

  short <- counts |> filter(n_rows < 3)
  if (nrow(short) > 0) {
    message(basename(path), ": incomplete foreleg data (missing trait(s) → NA):\n",
            paste(paste0("  ", short$ID, " ", short$side,
                         " (", short$n_rows, " of 3 rows)"), collapse = "\n"))
  }

  trait_order <- c("forefemur", "foretibia", "ear")

  long <- raw |>
    group_by(ID, side) |>
    arrange(source_row, .by_group = TRUE) |>
    mutate(
      n = row_number(),
      max_n = n_distinct(source_row)  # Total rows for this ID-side before filtering
    ) |>
    ungroup() |>
    group_by(ID, side) |>
    mutate(
      max_n = max(max_n)  # Propagate max to all rows in group
    ) |>
    ungroup() |>
    filter(n <= 3) |>
    group_by(ID, side) |>
    mutate(
      # Assign traits based on row position and whether we have 1, 2, 3, or 3+ rows
      # For 3+ rows: use first 3 (femur, tibia, ear)
      # For 3 rows: row 1 = femur, row 2 = tibia, row 3 = ear
      # For 2 rows: row 1 = femur, row 2 = ear
      # For 1 row: row 1 = femur
      trait_base = case_when(
        (max_n >= 3L) & n == 1L ~ "forefemur",
        (max_n >= 3L) & n == 2L ~ "foretibia",
        (max_n >= 3L) & n == 3L ~ "ear",
        max_n == 2L & n == 1L ~ "forefemur",
        max_n == 2L & n == 2L ~ "ear",
        max_n == 1L & n == 1L ~ "forefemur"
      ),
      trait = paste0(side, "_", trait_base),
      # Use Area for ear, Length for legs
      value = if_else(trait_base == "ear", Area, Length)
    ) |>
    ungroup() |>
    select(ID, trait, value) |>
    pivot_wider(names_from = trait, values_from = value)

  # Ensure all six columns exist (may be absent if no individual had that side)
  bilateral_cols <- c("left_forefemur", "right_forefemur",
                      "left_foretibia", "right_foretibia",
                      "left_ear",       "right_ear")
  for (col in setdiff(bilateral_cols, names(long))) {
    long[[col]] <- NA_real_
  }

  long |>
    mutate(
      forefemur = rowMeans(pick(left_forefemur, right_forefemur), na.rm = TRUE),
      foretibia = rowMeans(pick(left_foretibia, right_foretibia), na.rm = TRUE),
      ear       = rowMeans(pick(left_ear,       right_ear),       na.rm = TRUE),
      across(c(forefemur, foretibia, ear), ~ if_else(is.nan(.x), NA_real_, .x))
    ) |>
    select(ID, forefemur, foretibia, ear) |>
    arrange(ID)
}


# --- Midleg ---
# Combined left + right file.
# Side from Label: "2nd leg-left" / "2nd leg-right".
# Row order per side: 1 = femur, 2 = tibia.
# Sides with 1 row: femur only, tibia = NA.
# Extra rows beyond 2 discarded with a message.

clean_midleg_file <- function(path) {

  raw <- read_csv(
    path,
    na = na_vals,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) |>
    mutate(
      source_row = row_number(),
      ID         = str_to_upper(str_trim(str_extract(Label, "^[^-]+"))),
      side       = case_when(
        str_detect(Label, regex("2nd\\s+legs?\\s*-\\s*left",  ignore_case = TRUE)) ~ "left",
        str_detect(Label, regex("2nd\\s+legs?\\s*-\\s*right", ignore_case = TRUE)) ~ "right",
        TRUE ~ NA_character_
      ),
      Length = parse_double(Length, na = na_vals)
    )

  counts <- raw |> count(ID, side, name = "n_rows")

  extra <- counts |> filter(n_rows > 2)
  if (nrow(extra) > 0) {
    message(basename(path), ": extra midleg rows (using first 2 per side):\n",
            paste(paste0("  ", extra$ID, " ", extra$side,
                         " (", extra$n_rows, " rows)"), collapse = "\n"))
  }

  short <- counts |> filter(n_rows == 1)
  if (nrow(short) > 0) {
    message(basename(path), ": midleg sides with 1 row (femur only, tibia = NA):\n",
            paste(paste0("  ", short$ID, " ", short$side), collapse = "\n"))
  }

  trait_order <- c("midfemur", "midtibia")

  long <- raw |>
    group_by(ID, side) |>
    arrange(source_row, .by_group = TRUE) |>
    mutate(n = row_number()) |>
    filter(n <= 2) |>
    mutate(
      trait = paste0(side, "_", trait_order[n]),
      value = Length
    ) |>
    ungroup() |>
    select(ID, trait, value) |>
    pivot_wider(names_from = trait, values_from = value)

  bilateral_cols <- c("left_midfemur", "right_midfemur",
                      "left_midtibia", "right_midtibia")
  for (col in setdiff(bilateral_cols, names(long))) {
    long[[col]] <- NA_real_
  }

  long |>
    mutate(
      midfemur = rowMeans(pick(left_midfemur, right_midfemur), na.rm = TRUE),
      midtibia = rowMeans(pick(left_midtibia, right_midtibia), na.rm = TRUE),
      across(c(midfemur, midtibia), ~ if_else(is.nan(.x), NA_real_, .x))
    ) |>
    select(ID, midfemur, midtibia) |>
    arrange(ID)
}


# --- Hindleg ---
# Combined left + right file.
# Side from Label: "3rd leg-left" / "3rd leg-right".
# Row order per side: 1 = femur, 2 = tibia.
# Sides with 1 row: femur only, tibia = NA (same rule as midleg).
# Extra rows beyond 2 discarded with a message.

clean_hindleg_file <- function(path) {

  raw <- read_csv(
    path,
    na = na_vals,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) |>
    mutate(
      source_row = row_number(),
      ID         = str_to_upper(str_trim(str_extract(Label, "^[^-]+"))),
      side       = case_when(
        str_detect(Label, regex("3rd\\s+legs?\\s*-\\s*left",  ignore_case = TRUE)) ~ "left",
        str_detect(Label, regex("3rd\\s+legs?\\s*-\\s*right", ignore_case = TRUE)) ~ "right",
        TRUE ~ NA_character_
      ),
      Length = parse_double(Length, na = na_vals)
    )

  counts <- raw |> count(ID, side, name = "n_rows")

  extra <- counts |> filter(n_rows > 2)
  if (nrow(extra) > 0) {
    message(basename(path), ": extra hindleg rows (using first 2 per side):\n",
            paste(paste0("  ", extra$ID, " ", extra$side,
                         " (", extra$n_rows, " rows)"), collapse = "\n"))
  }

  short <- counts |> filter(n_rows == 1)
  if (nrow(short) > 0) {
    message(basename(path), ": hindleg sides with 1 row (femur only, tibia = NA):\n",
            paste(paste0("  ", short$ID, " ", short$side), collapse = "\n"))
  }

  trait_order <- c("hindfemur", "hindtibia")

  long <- raw |>
    group_by(ID, side) |>
    arrange(source_row, .by_group = TRUE) |>
    mutate(n = row_number()) |>
    filter(n <= 2) |>
    mutate(
      trait = paste0(side, "_", trait_order[n]),
      value = Length
    ) |>
    ungroup() |>
    select(ID, trait, value) |>
    pivot_wider(names_from = trait, values_from = value)

  bilateral_cols <- c("left_hindfemur", "right_hindfemur",
                      "left_hindtibia", "right_hindtibia")
  for (col in setdiff(bilateral_cols, names(long))) {
    long[[col]] <- NA_real_
  }

  long |>
    mutate(
      hindfemur = rowMeans(pick(left_hindfemur, right_hindfemur), na.rm = TRUE),
      hindtibia = rowMeans(pick(left_hindtibia, right_hindtibia), na.rm = TRUE),
      across(c(hindfemur, hindtibia), ~ if_else(is.nan(.x), NA_real_, .x))
    ) |>
    select(ID, hindfemur, hindtibia) |>
    arrange(ID)
}


# ============================================================
# 5. Build one measurement round from a suffix ("_two"/"_three")
# ============================================================

build_round <- function(suffix, round_number) {

  message("\n--- Building round ", round_number, " (suffix: ", suffix, ") ---")

  head_dat     <- clean_head_file(    file.path(raw_dir, paste0("head",     suffix, ".csv")))
  foreleg_dat  <- clean_foreleg_file( file.path(raw_dir, paste0("foreleg",  suffix, ".csv")))
  midleg_dat   <- clean_midleg_file(  file.path(raw_dir, paste0("midleg",   suffix, ".csv")))
  hindleg_dat  <- clean_hindleg_file( file.path(raw_dir, paste0("hindleg",  suffix, ".csv")))
  pronotum_dat <- clean_pronotum_file(file.path(raw_dir, paste0("pronotum", suffix, ".csv")))

  combined <- list(pronotum_dat, head_dat, foreleg_dat, midleg_dat, hindleg_dat) |>
    reduce(full_join, by = "ID") |>
    filter(ID %in% rep_ids) |>
    select(
      ID,
      pronotum, head_length, head_width,
      forefemur, foretibia, ear,
      midfemur, midtibia,
      hindfemur, hindtibia,
      eye
    ) |>
    mutate(measurement_round = round_number)

  message("  Individuals in round ", round_number, ": ", nrow(combined))

  combined
}


# ============================================================
# 6. Process rounds 2 and 3
# ============================================================

round2 <- build_round("_two",   2L)
round3 <- build_round("_three", 3L)


# ============================================================
# 7. Combine all three rounds
# ============================================================

repeat_data <- bind_rows(round1, round2, round3) |>
  mutate(measurement_round = factor(measurement_round))

# Add the pre-specified sex and morph metadata for robustness analyses while
# preserving repeat_data in its original, measurement-only form.
repeat_metadata <- rep_metadata |>
  select(ID, sex, morph) |>
  distinct()

if (anyDuplicated(repeat_metadata$ID)) {
  stop("Repeatability metadata contain duplicate IDs.", call. = FALSE)
}

repeat_data_with_meta <- repeat_data |>
  left_join(repeat_metadata, by = "ID", relationship = "many-to-one") |>
  mutate(
    sex = factor(sex),
    morph = factor(morph, levels = c("female", "eighth", "ninth", "tenth"))
  )

if (any(is.na(repeat_data_with_meta$morph))) {
  stop("Missing morph metadata for one or more repeatability individuals.",
       call. = FALSE)
}

message("\nTotal rows in repeat_data: ", nrow(repeat_data),
        " (", n_distinct(repeat_data$ID), " individuals × 3 rounds)")

# Summary of missingness by trait and round
repeat_data |>
  group_by(measurement_round) |>
  summarise(
    n = n(),
    across(
      -ID,
      ~ sum(is.na(.x)),
      .names = "NA_{.col}"
    )
  ) |>
  print()


# ============================================================
# 8. Repeatability analysis (rptR, Gaussian ICC)
# ============================================================
#
# Model: trait ~ 1 + (1 | ID)
#
# Repeatability R = σ²_ID / (σ²_ID + σ²_residual)
# where σ²_ID is between-individual variance and σ²_residual is
# within-individual (measurement error) variance.
#
# measurement_round is not included as a fixed effect because
# we are interested in total measurement repeatability across
# sessions. If systematic inter-session differences are present,
# re-run with: formula = value ~ measurement_round.
# ============================================================

trait_cols <- c(
  "pronotum", "head_length", "head_width",
  "forefemur", "foretibia", "ear",
  "midfemur", "midtibia",
  "hindfemur", "hindtibia",
  "eye"
)

set.seed(7293)

repeatability_results <- map(trait_cols, function(trait) {

  message("Running rptR for: ", trait)

  trait_data <- repeat_data |>
    select(ID, measurement_round, value = all_of(trait)) |>
    filter(!is.na(value)) |>
    mutate(ID = factor(ID))

  n_ind <- n_distinct(trait_data$ID)
  n_obs <- nrow(trait_data)
  message("  n individuals = ", n_ind, ", n observations = ", n_obs)

  rpt(
    formula  = value ~ 1 + (1 | ID),
    grname   = "ID",
    data     = trait_data,
    datatype = "Gaussian",
    nboot    = 1000,
    npermut  = 0
  )
})

names(repeatability_results) <- trait_cols


# ============================================================
# 8a. Robustness analyses
# ============================================================
#
# The primary model follows the preregistered/published-style analysis above.
# Because females and male morphs occupy different trait ranges, a pooled ICC
# can partly reflect among-morph variance. The morph-adjusted model estimates
# repeatability after removing average differences among the four groups.
# A second robustness model also removes systematic differences among rounds;
# this is a consistency ICC rather than the primary absolute-agreement ICC.
# ============================================================

run_adjusted_repeatability <- function(trait, include_round = FALSE) {

  trait_data <- repeat_data_with_meta |>
    transmute(
      ID = factor(ID),
      morph,
      measurement_round,
      value = .data[[trait]]
    ) |>
    filter(!is.na(value))

  model_formula <- if (include_round) {
    value ~ morph + measurement_round + (1 | ID)
  } else {
    value ~ morph + (1 | ID)
  }

  rpt(
    formula  = model_formula,
    grname   = "ID",
    data     = trait_data,
    datatype = "Gaussian",
    nboot    = 1000,
    npermut  = 0,
    adjusted = TRUE
  )
}

message("\nRunning morph-adjusted repeatability models")
set.seed(7293)
morph_adjusted_results <- map(
  trait_cols,
  ~ run_adjusted_repeatability(.x, include_round = FALSE)
)
names(morph_adjusted_results) <- trait_cols

message("\nRunning morph- and round-adjusted repeatability models")
set.seed(7293)
morph_round_adjusted_results <- map(
  trait_cols,
  ~ run_adjusted_repeatability(.x, include_round = TRUE)
)
names(morph_round_adjusted_results) <- trait_cols


# ============================================================
# 9. Summary table
# ============================================================

repeatability_summary <- map_dfr(
  trait_cols,
  function(trait) {

    res <- repeatability_results[[trait]]

    tibble(
      trait    = trait,
      R        = round(res$R$ID,         3),
      SE       = round(res$se$se,        3),
      CI_lower = round(res$CI_emp[1, 1], 3),
      CI_upper = round(res$CI_emp[1, 2], 3),
      LRT_P    = signif(res$P$LRT_P,     3),
      n_obs    = res$nobs,
      n_ind    = res$ngroups[["ID"]]
    )
  }
)

summarise_rpt_results <- function(results, model_name) {

  map_dfr(
    trait_cols,
    function(trait) {

      res <- results[[trait]]

      tibble(
        trait    = trait,
        model    = model_name,
        R        = round(res$R$ID,         4),
        SE       = round(res$se$se,        4),
        CI_lower = round(res$CI_emp[1, 1], 4),
        CI_upper = round(res$CI_emp[1, 2], 4),
        LRT_P    = signif(res$P$LRT_P,     4),
        n_obs    = res$nobs,
        n_ind    = res$ngroups[["ID"]]
      )
    }
  )
}

repeatability_robustness <- bind_rows(
  summarise_rpt_results(repeatability_results, "primary_unadjusted"),
  summarise_rpt_results(morph_adjusted_results, "morph_adjusted"),
  summarise_rpt_results(
    morph_round_adjusted_results,
    "morph_and_round_adjusted"
  )
)


# ============================================================
# 9a. Measurement-round diagnostics
# ============================================================
#
# Round effects are tested with a repeated-measures mixed model containing
# measurement round as a fixed effect and ID as a random intercept. Pairwise
# bias is the mean within-individual difference (later minus earlier round).
# Pairwise technical error of measurement (TEM) is
# sqrt(sum(difference^2) / (2 * n)); relative TEM is expressed as a percentage
# of the mean across the two rounds. Random TEM removes the mean pairwise bias.
# ============================================================

repeat_long <- repeat_data_with_meta |>
  pivot_longer(
    cols = all_of(trait_cols),
    names_to = "trait",
    values_to = "value"
  )

measurement_round_means <- repeat_long |>
  filter(!is.na(value)) |>
  group_by(trait, measurement_round) |>
  summarise(
    n = n(),
    mean = mean(value),
    SD = sd(value),
    .groups = "drop"
  )

round_pairs <- tribble(
  ~round_earlier, ~round_later,
  "1",            "2",
  "1",            "3",
  "2",            "3"
)

measurement_round_pairwise <- map_dfr(
  trait_cols,
  function(trait) {

    wide <- repeat_data |>
      select(ID, measurement_round, value = all_of(trait)) |>
      pivot_wider(
        names_from = measurement_round,
        values_from = value,
        names_prefix = "round_"
      )

    pmap_dfr(
      round_pairs,
      function(round_earlier, round_later) {

        earlier <- wide[[paste0("round_", round_earlier)]]
        later   <- wide[[paste0("round_", round_later)]]
        keep    <- complete.cases(earlier, later)
        earlier <- earlier[keep]
        later   <- later[keep]
        difference <- later - earlier
        n_pairs <- length(difference)

        mean_bias <- mean(difference)
        bias_SE <- sd(difference) / sqrt(n_pairs)
        critical_t <- qt(0.975, df = n_pairs - 1)
        pair_mean <- mean(c(earlier, later))
        TEM <- sqrt(sum(difference^2) / (2 * n_pairs))
        random_TEM <- sd(difference) / sqrt(2)

        tibble(
          trait = trait,
          comparison = paste0("round_", round_later, "_minus_round_", round_earlier),
          n_pairs = n_pairs,
          mean_bias = mean_bias,
          bias_SE = bias_SE,
          bias_CI_lower = mean_bias - critical_t * bias_SE,
          bias_CI_upper = mean_bias + critical_t * bias_SE,
          bias_P = if (sd(difference) > 0) {
            t.test(difference, mu = 0)$p.value
          } else {
            NA_real_
          },
          bias_percent = 100 * mean_bias / pair_mean,
          TEM = TEM,
          relative_TEM_percent = 100 * TEM / pair_mean,
          random_TEM = random_TEM,
          relative_random_TEM_percent = 100 * random_TEM / pair_mean
        )
      }
    )
  }
) |>
  group_by(trait) |>
  mutate(bias_P_Holm = p.adjust(bias_P, method = "holm")) |>
  ungroup()

round_effect_tests <- map_dfr(
  trait_cols,
  function(trait) {

    trait_data <- repeat_data_with_meta |>
      transmute(
        ID = factor(ID),
        morph,
        measurement_round,
        value = .data[[trait]]
      ) |>
      filter(!is.na(value))

    round_model <- lmerTest::lmer(
      value ~ measurement_round + (1 | ID),
      data = trait_data,
      REML = TRUE
    )

    round_anova <- anova(round_model, ddf = "Satterthwaite")
    round_row <- round_anova["measurement_round", , drop = FALSE]

    within_subject <- trait_data |>
      group_by(ID) |>
      summarise(
        within_SS = sum((value - mean(value))^2),
        within_df = n() - 1L,
        .groups = "drop"
      ) |>
      filter(within_df > 0)

    overall_TEM <- sqrt(
      sum(within_subject$within_SS) / sum(within_subject$within_df)
    )

    complete_triplets <- trait_data |>
      select(ID, measurement_round, value) |>
      pivot_wider(
        names_from = measurement_round,
        values_from = value,
        names_prefix = "round_"
      ) |>
      drop_na(round_1, round_2, round_3)

    measurement_matrix <- complete_triplets |>
      select(round_1, round_2, round_3) |>
      as.matrix()

    n_complete <- nrow(measurement_matrix)
    k_rounds <- ncol(measurement_matrix)
    grand_mean <- mean(measurement_matrix)
    subject_means <- rowMeans(measurement_matrix)
    round_means <- colMeans(measurement_matrix)

    SS_subject <- k_rounds * sum((subject_means - grand_mean)^2)
    SS_round <- n_complete * sum((round_means - grand_mean)^2)
    SS_total <- sum((measurement_matrix - grand_mean)^2)
    SS_error <- SS_total - SS_subject - SS_round

    MS_subject <- SS_subject / (n_complete - 1L)
    MS_round <- SS_round / (k_rounds - 1L)
    MS_error <- SS_error / ((n_complete - 1L) * (k_rounds - 1L))

    ICC_A1 <- (MS_subject - MS_error) /
      (
        MS_subject + (k_rounds - 1L) * MS_error +
          k_rounds * (MS_round - MS_error) / n_complete
      )

    ICC_C1 <- (MS_subject - MS_error) /
      (MS_subject + (k_rounds - 1L) * MS_error)

    tibble(
      trait = trait,
      trait_mean = mean(trait_data$value),
      round_F = round_row$`F value`,
      round_num_df = round_row$NumDF,
      round_den_df = round_row$DenDF,
      round_P = round_row$`Pr(>F)`,
      overall_TEM = overall_TEM,
      relative_TEM_percent = 100 * overall_TEM / mean(trait_data$value),
      random_TEM_after_round_adjustment = sigma(round_model),
      relative_random_TEM_percent =
        100 * sigma(round_model) / mean(trait_data$value),
      n_complete_triplets = n_complete,
      absolute_agreement_ICC_A1 = ICC_A1,
      consistency_ICC_C1 = ICC_C1
    )
  }
) |>
  mutate(
    round_P_BH = p.adjust(round_P, method = "BH"),
    round_P_Holm = p.adjust(round_P, method = "holm")
  )

round_ranges <- measurement_round_means |>
  group_by(trait) |>
  summarise(
    round_mean_min = min(mean),
    round_mean_max = max(mean),
    round_mean_range = round_mean_max - round_mean_min,
    .groups = "drop"
  )

pairwise_extremes <- measurement_round_pairwise |>
  group_by(trait) |>
  summarise(
    max_absolute_pairwise_bias = max(abs(mean_bias)),
    max_absolute_pairwise_bias_percent = max(abs(bias_percent)),
    max_pairwise_relative_TEM_percent = max(relative_TEM_percent),
    .groups = "drop"
  )

measurement_round_effects <- round_effect_tests |>
  left_join(round_ranges, by = "trait") |>
  left_join(pairwise_extremes, by = "trait") |>
  left_join(
    repeatability_robustness |>
      filter(model == "morph_and_round_adjusted") |>
      select(
        trait,
        round_adjusted_R = R,
        round_adjusted_CI_lower = CI_lower,
        round_adjusted_CI_upper = CI_upper
      ),
    by = "trait"
  ) |>
  mutate(round_mean_range_percent = 100 * round_mean_range / trait_mean)


# Test for within-session drift despite randomized measurement order. Metadata
# round 1 order corresponds to measurement round 2, and metadata round 2 order
# corresponds to measurement round 3. Slopes are expressed per 10 positions.
order_lookup <- rep_metadata |>
  select(
    ID,
    round_2_order = remeasurement_round_1_order,
    round_3_order = remeasurement_round_2_order
  ) |>
  pivot_longer(
    cols = starts_with("round_"),
    names_to = "measurement_round",
    values_to = "measurement_order"
  ) |>
  mutate(
    measurement_round = str_extract(measurement_round, "[23]"),
    measurement_round = factor(measurement_round, levels = c("2", "3"))
  )

round_1_values <- repeat_long |>
  filter(measurement_round == "1") |>
  select(ID, trait, round_1_value = value)

order_differences <- repeat_long |>
  filter(measurement_round %in% c("2", "3")) |>
  left_join(round_1_values, by = c("ID", "trait")) |>
  left_join(order_lookup, by = c("ID", "measurement_round")) |>
  filter(!is.na(value), !is.na(round_1_value), !is.na(measurement_order)) |>
  mutate(
    signed_difference = value - round_1_value,
    absolute_difference = abs(signed_difference),
    order_per_10 = measurement_order / 10
  )

measurement_order_effects <- order_differences |>
  group_by(trait, measurement_round) |>
  group_modify(
    ~ {
      signed_model <- lm(signed_difference ~ order_per_10, data = .x)
      absolute_model <- lm(absolute_difference ~ order_per_10, data = .x)
      signed_coef <- summary(signed_model)$coefficients["order_per_10", ]
      absolute_coef <- summary(absolute_model)$coefficients["order_per_10", ]

      tibble(
        n = nrow(.x),
        signed_slope_per_10_positions = signed_coef["Estimate"],
        signed_slope_SE = signed_coef["Std. Error"],
        signed_slope_P = signed_coef["Pr(>|t|)"],
        absolute_error_slope_per_10_positions = absolute_coef["Estimate"],
        absolute_error_slope_SE = absolute_coef["Std. Error"],
        absolute_error_slope_P = absolute_coef["Pr(>|t|)"]
      )
    }
  ) |>
  ungroup() |>
  mutate(
    signed_slope_P_BH = p.adjust(signed_slope_P, method = "BH"),
    absolute_error_slope_P_BH =
      p.adjust(absolute_error_slope_P, method = "BH")
  )

cat("\n============================================================\n")
cat("Measurement repeatability (ICC) — 1000-bootstrap 95% CI\n")
cat("============================================================\n")
print(repeatability_summary, n = Inf)


# ============================================================
# 10. Save results
# ============================================================

output_dir <- Sys.getenv(
  "WETA_REPEATABILITY_OUTPUT_DIR",
  unset = file.path(project_dir, "analysis_outputs")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  repeatability_summary,
  file.path(output_dir, "measurement_repeatability.csv")
)

saveRDS(
  repeatability_results,
  file.path(output_dir, "repeatability_rptR_objects.rds")
)

write_csv(
  repeatability_robustness,
  file.path(output_dir, "measurement_repeatability_robustness.csv")
)

saveRDS(
  list(
    primary_unadjusted = repeatability_results,
    morph_adjusted = morph_adjusted_results,
    morph_and_round_adjusted = morph_round_adjusted_results
  ),
  file.path(output_dir, "repeatability_robustness_rptR_objects.rds")
)

write_csv(
  measurement_round_means,
  file.path(output_dir, "measurement_round_means.csv")
)

write_csv(
  measurement_round_pairwise,
  file.path(output_dir, "measurement_round_pairwise_diagnostics.csv")
)

write_csv(
  measurement_round_effects,
  file.path(output_dir, "measurement_round_effects.csv")
)

write_csv(
  measurement_order_effects,
  file.path(output_dir, "measurement_order_effects.csv")
)

# Save the combined repeat-measurement data
write_csv(
  repeat_data,
  file.path(output_dir, "repeat_measurement_data.csv"),
  na = ""
)

message("\nOutputs saved to: ", output_dir)
