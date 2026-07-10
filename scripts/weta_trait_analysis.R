# ============================================================
# WELLINGTON TREE WETA MORPHOLOGICAL ALLOCATION
# Complete, annotated analysis workflow
# Females and eighth-, ninth-, and tenth-instar males
# ============================================================
#
# PRIMARY QUESTIONS
# 1. Do the four developmental groups differ in foreleg,
#    midleg, and hindleg length after accounting for body size?
# 2. Do group differences depend on leg identity (relative
#    locomotory allocation)?
# 3. Do allometric slopes differ among groups and leg pairs?
# 4. Is ear size integrated with foretibia length, and does that
#    relationship differ among developmental groups?
# 5. Is eye size integrated with head size, and does that
#    relationship differ among developmental groups?
#
# IMPORTANT INTERPRETIVE NOTES
# - Pronotum length is used as the measure of structural body size.
# - All continuous morphological variables are analysed on a log
#   scale. Slopes are therefore allometric exponents.
# - logP_c = 0 corresponds to a pronotum length of 7.2 mm, which
#   lies within the observed range of all four groups.
# - Back-transformed EMMs are geometric means / predicted medians,
#   not arithmetic means.
# - The joint repeated-trait model is the primary analysis of
#   relative allocation among the three leg pairs.
# ============================================================


# 0. USER SETTINGS -------------------------------------------

# Directory containing trait_data.csv
source_dir <- paste0(
  "/Users/uqam/Documents/Research_Admin/research projects/",
  "Trait compensation/data/clean"
)

data_file <- file.path(source_dir, "trait_data.csv")

# Project-level output directory
project_dir <- dirname(dirname(source_dir))
output_dir  <- file.path(project_dir, "analysis_outputs", "weta_trait_analysis")
tables_dir  <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
models_dir  <- file.path(output_dir, "models")

invisible(lapply(
  c(output_dir, tables_dir, figures_dir, models_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# Reference pronotum length for adjusted comparisons.
# This value lies within the observed body-size range of all groups.
reference_pronotum <- 7.2

# Familywise significance threshold used for planned model tests.
alpha <- 0.05

# Use treatment contrasts so that coefficients in lm() summaries
# retain female and Foreleg as the reference levels. Omnibus tests
# are conducted with nested-model comparisons or lmerTest::anova(),
# and biological interpretation is based on emmeans contrasts.
options(contrasts = c("contr.treatment", "contr.poly"))


# 1. PACKAGES ------------------------------------------------

required_packages <- c(
  "tidyverse",
  "emmeans",
  "car",
  "lme4",
  "lmerTest",
  "performance",
  "see",
  "broom"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Install the following packages before running this script:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nRun:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Keep emmeans output stable and explicit.
emmeans::emm_options(
  lmer.df = "kenward-roger",
  disable.pbkrtest = FALSE
)


# 2. HELPER FUNCTIONS ----------------------------------------

# Convert an emmeans or emtrends object to a tibble with confidence
# intervals and tests.
tidy_emm <- function(x, adjust = NULL) {
  if (is.null(adjust)) {
    out <- summary(x, infer = c(TRUE, TRUE))
  } else {
    out <- summary(x, infer = c(TRUE, TRUE), adjust = adjust)
  }
  tibble::as_tibble(out)
}

# Extract a partial F-test comparing nested lm models.
# The reduced and full models must be fitted to the same observations.
nested_f_test <- function(reduced_model, full_model, analysis_name) {
  tab <- anova(reduced_model, full_model)

  tibble(
    analysis = analysis_name,
    reduced_residual_df = tab$Res.Df[1],
    full_residual_df = tab$Res.Df[2],
    numerator_df = tab$Df[2],
    reduced_RSS = tab$RSS[1],
    full_RSS = tab$RSS[2],
    sum_of_squares = tab$`Sum of Sq`[2],
    F = tab$F[2],
    p_value = tab$`Pr(>F)`[2]
  )
}

# Back-transform EMMs from a log response.
back_transform_emm <- function(emm_object, outcome, scale_label) {
  tidy_emm(emm_object) |>
    mutate(
      outcome = outcome,
      response_scale = scale_label,
      predicted = exp(emmean),
      lower = exp(lower.CL),
      upper = exp(upper.CL)
    )
}

# Convert pairwise contrasts on a log scale to ratios.
# A positive estimate means the first group in the contrast is larger.
ratio_contrast_table <- function(
    contrast_object,
    outcome,
    adjust = "holm"
) {
  tidy_emm(contrast_object, adjust = adjust) |>
    mutate(
      outcome = outcome,
      ratio = exp(estimate),
      ratio_lower = exp(lower.CL),
      ratio_upper = exp(upper.CL),
      percent_difference = 100 * (ratio - 1)
    )
}

# Diagnostics for lm objects, retaining original row number and ID.
get_lm_diagnostics <- function(model, original_data) {
  mf <- model.frame(model)

  used_rows <- suppressWarnings(as.integer(rownames(mf)))
  if (anyNA(used_rows)) {
    used_rows <- match(rownames(mf), rownames(original_data))
  }

  broom::augment(model) |>
    mutate(
      row_id = used_rows,
      ID = original_data$ID[used_rows],
      group = original_data$group[used_rows]
    ) |>
    relocate(row_id, ID, group)
}

# Flag observations for inspection. These are screening criteria,
# not automatic deletion rules.
flag_lm_diagnostics <- function(diag_data, model) {
  n <- nobs(model)
  p <- length(coef(model))

  diag_data |>
    mutate(
      high_cooks = .cooksd > 4 / n,
      high_leverage = .hat > 2 * p / n,
      large_residual = abs(.std.resid) > 3,
      flagged = high_cooks | high_leverage | large_residual
    ) |>
    arrange(desc(.cooksd))
}

# Save a performance::check_model panel without interrupting the
# rest of the analysis if a graphics device fails.
save_check_model <- function(model, filename) {
  try({
    png(
      filename,
      width = 3300,
      height = 4200,
      res = 300
    )
    print(performance::check_model(model))
    dev.off()
  }, silent = TRUE)

  # Close an open device if an error occurred after png().
  if (dev.cur() > 1) {
    try(dev.off(), silent = TRUE)
  }
}

# Consistent plot theme.
theme_weta <- function() {
  theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}


# 3. IMPORT DATA ---------------------------------------------

if (!file.exists(data_file)) {
  stop(
    paste0(
      "Data file not found:\n", data_file,
      "\nCheck source_dir in Section 0."
    ),
    call. = FALSE
  )
}

dat_raw <- readr::read_csv(
  data_file,
  na = c("", "NA"),
  show_col_types = FALSE
)

required_columns <- c(
  "ID", "sex", "morph", "pronotum",
  "head_length", "head_width",
  "forefemur", "foretibia",
  "midfemur", "midtibia",
  "hindfemur", "hindtibia",
  "ear", "eye"
)

missing_columns <- setdiff(required_columns, names(dat_raw))
if (length(missing_columns) > 0) {
  stop(
    paste(
      "The data file is missing required columns:",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

cat("\nDATA STRUCTURE\n")
glimpse(dat_raw)


# 4. QUALITY CONTROL -----------------------------------------

# One row per individual is required.
duplicate_ids <- dat_raw |>
  count(ID, name = "n") |>
  filter(n > 1)

if (nrow(duplicate_ids) > 0) {
  print(duplicate_ids)
  stop("Duplicate IDs detected. Resolve them before analysis.", call. = FALSE)
}

expected_groups <- c("female", "eighth", "ninth", "tenth")
unexpected_groups <- setdiff(unique(na.omit(dat_raw$morph)), expected_groups)

if (length(unexpected_groups) > 0) {
  stop(
    paste(
      "Unexpected morph labels:",
      paste(unexpected_groups, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Missing-value summary.
missing_summary <- dat_raw |>
  summarise(across(everything(), ~ sum(is.na(.x)))) |>
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) |>
  arrange(desc(n_missing), variable)

# Raw ranges before transformation.
measurement_variables <- c(
  "pronotum", "head_length", "head_width",
  "forefemur", "foretibia",
  "midfemur", "midtibia",
  "hindfemur", "hindtibia",
  "ear", "eye"
)

raw_range_summary <- dat_raw |>
  summarise(
    across(
      all_of(measurement_variables),
      list(
        n = ~ sum(!is.na(.x)),
        minimum = ~ min(.x, na.rm = TRUE),
        maximum = ~ max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) |>
  pivot_longer(
    everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(n|minimum|maximum)$"
  )

# All logged measurements must be positive.
nonpositive_summary <- dat_raw |>
  summarise(
    across(
      all_of(measurement_variables),
      ~ sum(.x <= 0, na.rm = TRUE),
      .names = "{.col}"
    )
  ) |>
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_nonpositive"
  )

if (any(nonpositive_summary$n_nonpositive > 0)) {
  print(filter(nonpositive_summary, n_nonpositive > 0))
  stop(
    "Non-positive measurements detected. Log transformation is not valid.",
    call. = FALSE
  )
}

# Sex-by-morph table is retained as a QC check.
sex_morph_table <- dat_raw |>
  count(sex, morph, name = "n")

cat("\nSEX BY MORPH\n")
print(sex_morph_table)


# 5. CREATE ANALYSIS VARIABLES -------------------------------

dat <- dat_raw |>
  mutate(
    row_id = row_number(),
    group = factor(
      morph,
      levels = c("female", "eighth", "ninth", "tenth")
    ),

    # Ordered male instar number for descriptive use only.
    male_morph_ordered = case_when(
      group == "eighth" ~ 8,
      group == "ninth"  ~ 9,
      group == "tenth"  ~ 10,
      TRUE               ~ NA_real_
    ),

    # Total femur + tibia length for each leg pair.
    foreleg = forefemur + foretibia,
    midleg  = midfemur + midtibia,
    hindleg = hindfemur + hindtibia,

    # Composite linear head size: geometric mean of head length
    # and head width.
    head_size = sqrt(head_length * head_width),

    # Ear was measured as an area. Square-root transformation
    # expresses it as a linear dimension before log transformation.
    ear_linear = sqrt(ear),

    # logP_c = 0 represents reference_pronotum exactly.
    logP_c = log(pronotum) - log(reference_pronotum),

    # Logged total leg lengths.
    log_foreleg = log(foreleg),
    log_midleg  = log(midleg),
    log_hindleg = log(hindleg),

    # Logged segment lengths.
    log_forefemur = log(forefemur),
    log_foretibia = log(foretibia),
    log_midfemur  = log(midfemur),
    log_midtibia  = log(midtibia),
    log_hindfemur = log(hindfemur),
    log_hindtibia = log(hindtibia),

    # Logged sensory-supporting structures and sensory traits.
    log_head_size = log(head_size),
    log_eye = log(eye),
    log_ear_linear = log(ear_linear)
  ) |>
  mutate(
    # Centre supporting-trait covariates at their sample means.
    # This makes group EMMs refer to a typical foretibia/head size.
    log_foretibia_c = log_foretibia - mean(log_foretibia, na.rm = TRUE),
    log_head_size_c = log_head_size - mean(log_head_size, na.rm = TRUE)
  )

# Confirm that centring worked.
centering_check <- dat |>
  summarise(
    reference_pronotum = reference_pronotum,
    mean_log_foretibia_c = mean(log_foretibia_c, na.rm = TRUE),
    mean_log_head_size_c = mean(log_head_size_c, na.rm = TRUE)
  )


# 6. DESCRIPTIVE STATISTICS ----------------------------------

group_summary <- dat |>
  group_by(group) |>
  summarise(
    n = n(),

    pronotum_n = sum(!is.na(pronotum)),
    pronotum_mean = mean(pronotum, na.rm = TRUE),
    pronotum_sd = sd(pronotum, na.rm = TRUE),
    pronotum_min = min(pronotum, na.rm = TRUE),
    pronotum_max = max(pronotum, na.rm = TRUE),

    foreleg_n = sum(!is.na(foreleg)),
    foreleg_mean = mean(foreleg, na.rm = TRUE),
    foreleg_sd = sd(foreleg, na.rm = TRUE),

    midleg_n = sum(!is.na(midleg)),
    midleg_mean = mean(midleg, na.rm = TRUE),
    midleg_sd = sd(midleg, na.rm = TRUE),

    hindleg_n = sum(!is.na(hindleg)),
    hindleg_mean = mean(hindleg, na.rm = TRUE),
    hindleg_sd = sd(hindleg, na.rm = TRUE),

    ear_n = sum(!is.na(ear)),
    ear_mean = mean(ear, na.rm = TRUE),
    ear_sd = sd(ear, na.rm = TRUE),

    eye_n = sum(!is.na(eye)),
    eye_mean = mean(eye, na.rm = TRUE),
    eye_sd = sd(eye, na.rm = TRUE),

    .groups = "drop"
  )

body_size_ranges <- dat |>
  group_by(group) |>
  summarise(
    minimum = min(pronotum, na.rm = TRUE),
    maximum = max(pronotum, na.rm = TRUE),
    .groups = "drop"
  )

common_body_size <- tibble(
  common_lower = max(body_size_ranges$minimum),
  common_upper = min(body_size_ranges$maximum),
  reference_pronotum = reference_pronotum,
  reference_inside_common_range =
    reference_pronotum >= common_lower &
    reference_pronotum <= common_upper
)

cat("\nGROUP SUMMARY\n")
print(group_summary)
cat("\nCOMMON BODY-SIZE SUPPORT\n")
print(common_body_size)

if (!common_body_size$reference_inside_common_range) {
  warning(
    "The reference pronotum lies outside the range shared by all groups."
  )
}


# 7. BODY-SIZE DISTRIBUTION ----------------------------------

p_body <- ggplot(dat, aes(x = pronotum, fill = group)) +
  geom_density(alpha = 0.35, na.rm = TRUE) +
  geom_vline(
    xintercept = reference_pronotum,
    linetype = 2,
    linewidth = 0.6
  ) +
  labs(
    x = "Pronotum length (mm)",
    y = "Density",
    fill = "Group"
  ) +
  theme_weta()


# 8. LEG DATA IN LONG FORMAT ---------------------------------

legs_long <- dat |>
  select(
    row_id, ID, group, pronotum, logP_c,
    foreleg, midleg, hindleg
  ) |>
  pivot_longer(
    cols = c(foreleg, midleg, hindleg),
    names_to = "leg",
    values_to = "leg_length"
  ) |>
  mutate(
    leg = factor(
      leg,
      levels = c("foreleg", "midleg", "hindleg"),
      labels = c("Foreleg", "Midleg", "Hindleg")
    ),
    log_leg_length = log(leg_length),
    log_pronotum = log(pronotum)
  )

p_leg_allometry <- ggplot(
  legs_long,
  aes(
    x = log_pronotum,
    y = log_leg_length,
    colour = group
  )
) +
  geom_point(alpha = 0.50, na.rm = TRUE) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    na.rm = TRUE
  ) +
  facet_wrap(~ leg, scales = "free_y") +
  labs(
    x = "Log pronotum length",
    y = "Log leg length",
    colour = "Group"
  ) +
  theme_weta()


# 9. SEPARATE LEG ALLOMETRIES -------------------------------
# These models describe each leg pair separately. The joint mixed
# model in Section 10 is the primary test of relative allocation.

# Foreleg
mod_fore_common <- lm(
  log_foreleg ~ logP_c + group,
  data = dat,
  na.action = na.omit
)

mod_fore_interaction <- lm(
  log_foreleg ~ logP_c * group,
  data = dat,
  na.action = na.omit
)

# Midleg
mod_mid_common <- lm(
  log_midleg ~ logP_c + group,
  data = dat,
  na.action = na.omit
)

mod_mid_interaction <- lm(
  log_midleg ~ logP_c * group,
  data = dat,
  na.action = na.omit
)

# Hindleg
mod_hind_common <- lm(
  log_hindleg ~ logP_c + group,
  data = dat,
  na.action = na.omit
)

mod_hind_interaction <- lm(
  log_hindleg ~ logP_c * group,
  data = dat,
  na.action = na.omit
)

# Omnibus tests of common versus group-specific slopes.
leg_slope_model_tests <- bind_rows(
  nested_f_test(
    mod_fore_common,
    mod_fore_interaction,
    "Foreleg: logP_c x group"
  ),
  nested_f_test(
    mod_mid_common,
    mod_mid_interaction,
    "Midleg: logP_c x group"
  ),
  nested_f_test(
    mod_hind_common,
    mod_hind_interaction,
    "Hindleg: logP_c x group"
  )
)

cat("\nSEPARATE LEG SLOPE TESTS\n")
print(leg_slope_model_tests)

# Group-specific allometric slopes from the interaction models.
fore_slopes <- emtrends(
  mod_fore_interaction,
  ~ group,
  var = "logP_c"
)
mid_slopes <- emtrends(
  mod_mid_interaction,
  ~ group,
  var = "logP_c"
)
hind_slopes <- emtrends(
  mod_hind_interaction,
  ~ group,
  var = "logP_c"
)

leg_slopes <- bind_rows(
  tidy_emm(fore_slopes) |> mutate(leg = "Foreleg"),
  tidy_emm(mid_slopes)  |> mutate(leg = "Midleg"),
  tidy_emm(hind_slopes) |> mutate(leg = "Hindleg")
) |>
  relocate(leg)

fore_slope_pairs <- pairs(fore_slopes, adjust = "holm")
mid_slope_pairs  <- pairs(mid_slopes, adjust = "holm")
hind_slope_pairs <- pairs(hind_slopes, adjust = "holm")

leg_slope_contrasts <- bind_rows(
  tidy_emm(fore_slope_pairs, adjust = "holm") |>
    mutate(leg = "Foreleg"),
  tidy_emm(mid_slope_pairs, adjust = "holm") |>
    mutate(leg = "Midleg"),
  tidy_emm(hind_slope_pairs, adjust = "holm") |>
    mutate(leg = "Hindleg")
) |>
  relocate(leg)


# 10. PRIMARY JOINT LOCOMOTORY-ALLOCATION MODEL --------------
# Each individual contributes up to three leg measurements, so ID
# is included as a random intercept. The three-way interaction asks
# whether group-specific allometric slopes depend on leg identity.

mod_leg_allocation <- lmerTest::lmer(
  log_leg_length ~ logP_c * group * leg + (1 | ID),
  data = legs_long,
  REML = FALSE,
  na.action = na.omit
)

leg_allocation_anova <- lmerTest::anova(
  mod_leg_allocation,
  type = 3,
  ddf = "Satterthwaite"
) |>
  as.data.frame() |>
  rownames_to_column("term") |>
  as_tibble()

cat("\nJOINT LEG-ALLOCATION MODEL\n")
print(leg_allocation_anova)
print(summary(mod_leg_allocation))

# Adjusted group means at pronotum = reference_pronotum.
emm_allocation <- emmeans(
  mod_leg_allocation,
  ~ group | leg,
  at = list(logP_c = 0)
)

# All six group comparisons separately within each leg pair.
allocation_pairs <- pairs(
  emm_allocation,
  by = "leg",
  adjust = "holm"
)

# More directly useful Results contrasts: each male group versus
# females within each leg pair. Ratios > 1 mean male > female.
allocation_vs_female <- contrast(
  emm_allocation,
  method = "trt.vs.ctrl",
  ref = 1,
  by = "leg",
  adjust = "holm"
)

adjusted_leg_means <- back_transform_emm(
  emm_allocation,
  outcome = "Leg length",
  scale_label = "original length"
) |>
  mutate(
    leg = factor(leg, levels = c("Foreleg", "Midleg", "Hindleg"))
  )

leg_pairwise_contrasts <- ratio_contrast_table(
  allocation_pairs,
  outcome = "Leg length",
  adjust = "holm"
)

leg_vs_female_contrasts <- ratio_contrast_table(
  allocation_vs_female,
  outcome = "Leg length",
  adjust = "holm"
)

# Joint slope estimates and planned decompositions of the significant
# or potentially important three-way interaction.
joint_slopes <- emtrends(
  mod_leg_allocation,
  ~ group * leg,
  var = "logP_c"
)

joint_slopes_table <- tidy_emm(joint_slopes)

# Group slope comparisons within each leg: six tests per leg family.
joint_group_slope_pairs <- pairs(
  joint_slopes,
  by = "leg",
  adjust = "holm"
)

joint_group_slope_contrasts <- tidy_emm(
  joint_group_slope_pairs,
  adjust = "holm"
)

# Leg slope comparisons within each group: three tests per group.
joint_leg_slope_pairs <- pairs(
  joint_slopes,
  by = "group",
  adjust = "holm"
)

joint_leg_slope_contrasts <- tidy_emm(
  joint_leg_slope_pairs,
  adjust = "holm"
)

# Difference-in-differences: tests whether a particular group slope
# contrast differs between two leg types. Holm is applied across the
# complete set of these decompositions.
three_way_slope_contrasts <- contrast(
  joint_slopes,
  interaction = c("pairwise", "pairwise"),
  adjust = "holm"
)

three_way_slope_contrasts_table <- tidy_emm(
  three_way_slope_contrasts,
  adjust = "holm"
)

# Figure of size-adjusted total leg lengths.
p_adjusted_legs <- ggplot(
  adjusted_leg_means,
  aes(x = group, y = predicted)
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12
  ) +
  facet_wrap(~ leg, scales = "free_y") +
  labs(
    x = NULL,
    y = paste0(
      "Predicted leg length at pronotum = ",
      reference_pronotum,
      " mm"
    )
  ) +
  theme_weta()


# 11. LEG-SEGMENT ALLOMETRIES -------------------------------
# Secondary analyses decompose total-leg patterns into femur and
# tibia components. All six segment models are fitted consistently.

segments_long <- dat |>
  select(
    row_id, ID, group, pronotum, logP_c,
    forefemur, foretibia,
    midfemur, midtibia,
    hindfemur, hindtibia
  ) |>
  pivot_longer(
    cols = c(
      forefemur, foretibia,
      midfemur, midtibia,
      hindfemur, hindtibia
    ),
    names_to = "trait",
    values_to = "segment_length"
  ) |>
  mutate(
    leg = case_when(
      str_starts(trait, "fore") ~ "Foreleg",
      str_starts(trait, "mid")  ~ "Midleg",
      str_starts(trait, "hind") ~ "Hindleg"
    ),
    segment = case_when(
      str_ends(trait, "femur") ~ "Femur",
      str_ends(trait, "tibia") ~ "Tibia"
    ),
    leg = factor(leg, levels = c("Foreleg", "Midleg", "Hindleg")),
    segment = factor(segment, levels = c("Femur", "Tibia")),
    log_segment_length = log(segment_length)
  )

p_segments <- ggplot(
  segments_long,
  aes(
    x = logP_c,
    y = log_segment_length,
    colour = group
  )
) +
  geom_point(alpha = 0.45, na.rm = TRUE) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    na.rm = TRUE
  ) +
  facet_grid(segment ~ leg, scales = "free_y") +
  labs(
    x = paste0(
      "Log pronotum length centred at ",
      reference_pronotum,
      " mm"
    ),
    y = "Log segment length",
    colour = "Group"
  ) +
  theme_weta()

segment_response_map <- c(
  Forefemur = "log_forefemur",
  Foretibia = "log_foretibia",
  Midfemur = "log_midfemur",
  Midtibia = "log_midtibia",
  Hindfemur = "log_hindfemur",
  Hindtibia = "log_hindtibia"
)

segment_models_common <- purrr::map(
  segment_response_map,
  ~ lm(
    as.formula(paste(.x, "~ logP_c + group")),
    data = dat,
    na.action = na.omit
  )
)

segment_models_interaction <- purrr::map(
  segment_response_map,
  ~ lm(
    as.formula(paste(.x, "~ logP_c * group")),
    data = dat,
    na.action = na.omit
  )
)

segment_model_tests <- purrr::imap_dfr(
  segment_models_interaction,
  function(full_model, trait_name) {
    nested_f_test(
      segment_models_common[[trait_name]],
      full_model,
      paste0(trait_name, ": logP_c x group")
    )
  }
)

segment_trends <- purrr::map(
  segment_models_interaction,
  ~ emtrends(.x, ~ group, var = "logP_c")
)

segment_slopes <- purrr::imap_dfr(
  segment_trends,
  ~ tidy_emm(.x) |> mutate(segment = .y)
) |>
  relocate(segment)

segment_slope_contrasts <- purrr::imap_dfr(
  segment_trends,
  function(trend_object, trait_name) {
    tidy_emm(
      pairs(trend_object, adjust = "holm"),
      adjust = "holm"
    ) |>
      mutate(segment = trait_name)
  }
) |>
  relocate(segment)


# 12. EAR–FORETIBIA RELATIONSHIP -----------------------------
# The response is log(sqrt(ear area)), a linearized ear dimension.
# The full model asks whether ear–tibia and ear–body-size slopes
# differ among groups. Current analyses support common slopes, so
# the additive model is used for final adjusted group comparisons.

p_ear <- ggplot(
  dat,
  aes(
    x = log_foretibia,
    y = log_ear_linear,
    colour = group
  )
) +
  geom_point(alpha = 0.55, na.rm = TRUE) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    na.rm = TRUE
  ) +
  labs(
    x = "Log foretibia length",
    y = "Log linearized ear size",
    colour = "Group"
  ) +
  theme_weta()

# Full model with group-specific body-size and ear–tibia slopes.
mod_ear_full <- lm(
  log_ear_linear ~
    logP_c * group +
    log_foretibia_c * group,
  data = dat,
  na.action = na.omit
)

# Remove the group-specific ear–tibia slopes first.
mod_ear_no_tibia_interaction <- lm(
  log_ear_linear ~
    logP_c * group +
    log_foretibia_c,
  data = dat,
  na.action = na.omit
)

# Fully additive model with common body-size and ear–tibia slopes.
mod_ear_additive <- lm(
  log_ear_linear ~
    logP_c +
    log_foretibia_c +
    group,
  data = dat,
  na.action = na.omit
)

ear_model_tests <- bind_rows(
  nested_f_test(
    mod_ear_no_tibia_interaction,
    mod_ear_full,
    "Ear: group x foretibia"
  ),
  nested_f_test(
    mod_ear_additive,
    mod_ear_no_tibia_interaction,
    "Ear: group x pronotum"
  )
)

# Warn if new data no longer support the prespecified additive model.
if (ear_model_tests$p_value[1] < alpha) {
  warning(
    "The group x foretibia interaction is supported. Retain mod_ear_full."
  )
}
if (ear_model_tests$p_value[2] < alpha) {
  warning(
    paste0(
      "The group x pronotum interaction is supported. ",
      "Retain mod_ear_no_tibia_interaction."
    )
  )
}

# Type II partial tests are appropriate for the final additive model.
ear_additive_anova <- car::Anova(
  mod_ear_additive,
  type = 2
) |>
  as.data.frame() |>
  rownames_to_column("term") |>
  as_tibble()

# Common ear–tibia slope, controlling for pronotum and group.
ear_tibia_common_slope <- emtrends(
  mod_ear_additive,
  ~ 1,
  var = "log_foretibia_c"
)

# Group differences at pronotum = 7.2 mm and mean log foretibia.
ear_group_emm <- emmeans(
  mod_ear_additive,
  ~ group,
  at = list(
    logP_c = 0,
    log_foretibia_c = 0
  )
)

ear_group_pairs <- pairs(
  ear_group_emm,
  adjust = "holm"
)

ear_vs_female <- contrast(
  ear_group_emm,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "holm"
)

adjusted_ear_means <- tidy_emm(ear_group_emm) |>
  mutate(
    # Linearized ear dimension = sqrt(area).
    predicted_ear_linear = exp(emmean),
    ear_linear_lower = exp(lower.CL),
    ear_linear_upper = exp(upper.CL),

    # Return to the original area scale by squaring.
    predicted_ear_area = exp(2 * emmean),
    ear_area_lower = exp(2 * lower.CL),
    ear_area_upper = exp(2 * upper.CL)
  )

ear_group_contrasts <- ratio_contrast_table(
  ear_group_pairs,
  outcome = "Linearized ear size",
  adjust = "holm"
) |>
  mutate(
    area_ratio = ratio^2,
    area_percent_difference = 100 * (area_ratio - 1)
  )

ear_vs_female_contrasts <- ratio_contrast_table(
  ear_vs_female,
  outcome = "Linearized ear size",
  adjust = "holm"
) |>
  mutate(
    area_ratio = ratio^2,
    area_percent_difference = 100 * (area_ratio - 1)
  )


# 13. EYE–HEAD RELATIONSHIP ----------------------------------
# This mirrors the ear analysis. Eye length is related to a linear
# composite head size while pronotum controls overall body size.

p_eye <- ggplot(
  dat,
  aes(
    x = log_head_size,
    y = log_eye,
    colour = group
  )
) +
  geom_point(alpha = 0.55, na.rm = TRUE) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    na.rm = TRUE
  ) +
  labs(
    x = "Log head size",
    y = "Log eye length",
    colour = "Group"
  ) +
  theme_weta()

mod_eye_full <- lm(
  log_eye ~
    logP_c * group +
    log_head_size_c * group,
  data = dat,
  na.action = na.omit
)

mod_eye_no_head_interaction <- lm(
  log_eye ~
    logP_c * group +
    log_head_size_c,
  data = dat,
  na.action = na.omit
)

mod_eye_additive <- lm(
  log_eye ~
    logP_c +
    log_head_size_c +
    group,
  data = dat,
  na.action = na.omit
)

eye_model_tests <- bind_rows(
  nested_f_test(
    mod_eye_no_head_interaction,
    mod_eye_full,
    "Eye: group x head size"
  ),
  nested_f_test(
    mod_eye_additive,
    mod_eye_no_head_interaction,
    "Eye: group x pronotum"
  )
)

if (eye_model_tests$p_value[1] < alpha) {
  warning(
    "The group x head-size interaction is supported. Retain mod_eye_full."
  )
}
if (eye_model_tests$p_value[2] < alpha) {
  warning(
    paste0(
      "The group x pronotum interaction is supported. ",
      "Retain mod_eye_no_head_interaction."
    )
  )
}

eye_additive_anova <- car::Anova(
  mod_eye_additive,
  type = 2
) |>
  as.data.frame() |>
  rownames_to_column("term") |>
  as_tibble()

eye_head_common_slope <- emtrends(
  mod_eye_additive,
  ~ 1,
  var = "log_head_size_c"
)

eye_group_emm <- emmeans(
  mod_eye_additive,
  ~ group,
  at = list(
    logP_c = 0,
    log_head_size_c = 0
  )
)

eye_group_pairs <- pairs(
  eye_group_emm,
  adjust = "holm"
)

eye_vs_female <- contrast(
  eye_group_emm,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "holm"
)

adjusted_eye_means <- tidy_emm(eye_group_emm) |>
  mutate(
    predicted_eye_length = exp(emmean),
    eye_lower = exp(lower.CL),
    eye_upper = exp(upper.CL)
  )

eye_group_contrasts <- ratio_contrast_table(
  eye_group_pairs,
  outcome = "Eye length",
  adjust = "holm"
)

eye_vs_female_contrasts <- ratio_contrast_table(
  eye_vs_female,
  outcome = "Eye length",
  adjust = "holm"
)


# 14. MODEL DIAGNOSTICS --------------------------------------
# The joint mixed model is primary. Separate leg models are checked
# using the interaction form for forelegs and hindlegs and the
# common-slope form for midlegs, based on the current omnibus tests.

mod_fore_final <- mod_fore_interaction
mod_mid_final  <- mod_mid_common
mod_hind_final <- mod_hind_interaction

# Collinearity tables are easier to save and report than VIF plots.
collinearity_tables <- bind_rows(
  performance::check_collinearity(mod_fore_final) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(model = "Foreleg"),
  performance::check_collinearity(mod_mid_final) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(model = "Midleg"),
  performance::check_collinearity(mod_hind_final) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(model = "Hindleg"),
  performance::check_collinearity(mod_ear_additive) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(model = "Ear"),
  performance::check_collinearity(mod_eye_additive) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(model = "Eye")
) |>
  relocate(model)

# Observation-level diagnostics for the lm models.
fore_diag <- flag_lm_diagnostics(
  get_lm_diagnostics(mod_fore_final, dat),
  mod_fore_final
)
mid_diag <- flag_lm_diagnostics(
  get_lm_diagnostics(mod_mid_final, dat),
  mod_mid_final
)
hind_diag <- flag_lm_diagnostics(
  get_lm_diagnostics(mod_hind_final, dat),
  mod_hind_final
)
ear_diag <- flag_lm_diagnostics(
  get_lm_diagnostics(mod_ear_additive, dat),
  mod_ear_additive
)
eye_diag <- flag_lm_diagnostics(
  get_lm_diagnostics(mod_eye_additive, dat),
  mod_eye_additive
)

fore_diag_flagged <- filter(fore_diag, flagged)
mid_diag_flagged  <- filter(mid_diag, flagged)
hind_diag_flagged <- filter(hind_diag, flagged)
ear_diag_flagged  <- filter(ear_diag, flagged)
eye_diag_flagged  <- filter(eye_diag, flagged)

# Save graphical diagnostic panels.
save_check_model(
  mod_fore_final,
  file.path(figures_dir, "diagnostics_foreleg.png")
)
save_check_model(
  mod_mid_final,
  file.path(figures_dir, "diagnostics_midleg.png")
)
save_check_model(
  mod_hind_final,
  file.path(figures_dir, "diagnostics_hindleg.png")
)
save_check_model(
  mod_leg_allocation,
  file.path(figures_dir, "diagnostics_joint_leg_model.png")
)
save_check_model(
  mod_ear_additive,
  file.path(figures_dir, "diagnostics_ear_model.png")
)
save_check_model(
  mod_eye_additive,
  file.path(figures_dir, "diagnostics_eye_model.png")
)


# 15. SAVE RESULTS TABLES ------------------------------------

# QC and descriptive outputs.
write_csv(missing_summary, file.path(tables_dir, "qc_missing_values.csv"))
write_csv(raw_range_summary, file.path(tables_dir, "qc_raw_ranges.csv"))
write_csv(nonpositive_summary, file.path(tables_dir, "qc_nonpositive_values.csv"))
write_csv(sex_morph_table, file.path(tables_dir, "qc_sex_by_morph.csv"))
write_csv(group_summary, file.path(tables_dir, "group_descriptive_statistics.csv"))
write_csv(body_size_ranges, file.path(tables_dir, "body_size_ranges.csv"))
write_csv(common_body_size, file.path(tables_dir, "common_body_size_range.csv"))
write_csv(centering_check, file.path(tables_dir, "centering_check.csv"))

# Separate leg analyses.
write_csv(leg_slope_model_tests, file.path(tables_dir, "leg_slope_model_tests.csv"))
write_csv(leg_slopes, file.path(tables_dir, "leg_group_specific_slopes.csv"))
write_csv(leg_slope_contrasts, file.path(tables_dir, "leg_slope_pairwise_contrasts.csv"))

# Primary joint locomotory-allocation model.
write_csv(leg_allocation_anova, file.path(tables_dir, "joint_leg_model_type3_anova.csv"))
write_csv(adjusted_leg_means, file.path(tables_dir, "joint_adjusted_leg_means.csv"))
write_csv(leg_pairwise_contrasts, file.path(tables_dir, "joint_leg_all_pairwise_contrasts.csv"))
write_csv(leg_vs_female_contrasts, file.path(tables_dir, "joint_leg_males_vs_female.csv"))
write_csv(joint_slopes_table, file.path(tables_dir, "joint_leg_slopes.csv"))
write_csv(joint_group_slope_contrasts, file.path(tables_dir, "joint_group_slopes_within_leg.csv"))
write_csv(joint_leg_slope_contrasts, file.path(tables_dir, "joint_leg_slopes_within_group.csv"))
write_csv(three_way_slope_contrasts_table, file.path(tables_dir, "joint_three_way_slope_contrasts.csv"))

# Segment-level analyses.
write_csv(segment_model_tests, file.path(tables_dir, "segment_slope_model_tests.csv"))
write_csv(segment_slopes, file.path(tables_dir, "segment_group_specific_slopes.csv"))
write_csv(segment_slope_contrasts, file.path(tables_dir, "segment_slope_pairwise_contrasts.csv"))

# Ear analyses.
write_csv(ear_model_tests, file.path(tables_dir, "ear_nested_model_tests.csv"))
write_csv(ear_additive_anova, file.path(tables_dir, "ear_additive_type2_anova.csv"))
write_csv(tidy_emm(ear_tibia_common_slope), file.path(tables_dir, "ear_tibia_common_slope.csv"))
write_csv(adjusted_ear_means, file.path(tables_dir, "ear_adjusted_group_means.csv"))
write_csv(ear_group_contrasts, file.path(tables_dir, "ear_all_pairwise_group_contrasts.csv"))
write_csv(ear_vs_female_contrasts, file.path(tables_dir, "ear_males_vs_female.csv"))

# Eye analyses.
write_csv(eye_model_tests, file.path(tables_dir, "eye_nested_model_tests.csv"))
write_csv(eye_additive_anova, file.path(tables_dir, "eye_additive_type2_anova.csv"))
write_csv(tidy_emm(eye_head_common_slope), file.path(tables_dir, "eye_head_common_slope.csv"))
write_csv(adjusted_eye_means, file.path(tables_dir, "eye_adjusted_group_means.csv"))
write_csv(eye_group_contrasts, file.path(tables_dir, "eye_all_pairwise_group_contrasts.csv"))
write_csv(eye_vs_female_contrasts, file.path(tables_dir, "eye_males_vs_female.csv"))

# Diagnostics.
write_csv(collinearity_tables, file.path(tables_dir, "model_collinearity.csv"))
write_csv(fore_diag, file.path(tables_dir, "diagnostics_foreleg_all.csv"))
write_csv(mid_diag, file.path(tables_dir, "diagnostics_midleg_all.csv"))
write_csv(hind_diag, file.path(tables_dir, "diagnostics_hindleg_all.csv"))
write_csv(ear_diag, file.path(tables_dir, "diagnostics_ear_all.csv"))
write_csv(eye_diag, file.path(tables_dir, "diagnostics_eye_all.csv"))
write_csv(fore_diag_flagged, file.path(tables_dir, "diagnostics_foreleg_flagged.csv"))
write_csv(mid_diag_flagged, file.path(tables_dir, "diagnostics_midleg_flagged.csv"))
write_csv(hind_diag_flagged, file.path(tables_dir, "diagnostics_hindleg_flagged.csv"))
write_csv(ear_diag_flagged, file.path(tables_dir, "diagnostics_ear_flagged.csv"))
write_csv(eye_diag_flagged, file.path(tables_dir, "diagnostics_eye_flagged.csv"))


# 16. SAVE FIGURES -------------------------------------------

ggsave(
  file.path(figures_dir, "body_size_distributions.png"),
  p_body,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "leg_allometries.png"),
  p_leg_allometry,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "adjusted_leg_means.png"),
  p_adjusted_legs,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "segment_allometries.png"),
  p_segments,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "ear_foretibia_relationship.png"),
  p_ear,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "eye_head_relationship.png"),
  p_eye,
  width = 7,
  height = 5,
  dpi = 300
)


# 17. SAVE MODEL OBJECTS AND TEXT SUMMARY --------------------

model_objects <- list(
  fore_common = mod_fore_common,
  fore_interaction = mod_fore_interaction,
  mid_common = mod_mid_common,
  mid_interaction = mod_mid_interaction,
  hind_common = mod_hind_common,
  hind_interaction = mod_hind_interaction,
  joint_leg_allocation = mod_leg_allocation,
  segment_common = segment_models_common,
  segment_interaction = segment_models_interaction,
  ear_full = mod_ear_full,
  ear_no_tibia_interaction = mod_ear_no_tibia_interaction,
  ear_additive = mod_ear_additive,
  eye_full = mod_eye_full,
  eye_no_head_interaction = mod_eye_no_head_interaction,
  eye_additive = mod_eye_additive
)

saveRDS(
  model_objects,
  file.path(models_dir, "weta_trait_models.rds")
)

# A compact text file containing the key inferential outputs for
# drafting the Results section.
capture.output(
  {
    cat("WELLINGTON TREE WETA MORPHOLOGICAL ALLOCATION\n")
    cat("Analysis run:", format(Sys.time()), "\n\n")

    cat("REFERENCE PRONOTUM\n")
    print(reference_pronotum)
    print(common_body_size)

    cat("\nGROUP DESCRIPTIVES\n")
    print(group_summary)

    cat("\nSEPARATE LEG SLOPE MODEL TESTS\n")
    print(leg_slope_model_tests)

    cat("\nSEPARATE LEG SLOPES\n")
    print(leg_slopes)

    cat("\nSEPARATE LEG SLOPE CONTRASTS\n")
    print(leg_slope_contrasts)

    cat("\nJOINT LEG MODEL TYPE III ANOVA\n")
    print(leg_allocation_anova)

    cat("\nJOINT ADJUSTED LEG MEANS\n")
    print(adjusted_leg_means)

    cat("\nJOINT MALE-VERSUS-FEMALE LEG CONTRASTS\n")
    print(leg_vs_female_contrasts)

    cat("\nJOINT LEG SLOPES\n")
    print(joint_slopes_table)

    cat("\nJOINT GROUP SLOPE CONTRASTS WITHIN LEG\n")
    print(joint_group_slope_contrasts)

    cat("\nJOINT LEG SLOPE CONTRASTS WITHIN GROUP\n")
    print(joint_leg_slope_contrasts)

    cat("\nTHREE-WAY SLOPE CONTRASTS\n")
    print(three_way_slope_contrasts_table)

    cat("\nSEGMENT MODEL TESTS\n")
    print(segment_model_tests)

    cat("\nEAR MODEL TESTS\n")
    print(ear_model_tests)
    print(ear_additive_anova)
    print(tidy_emm(ear_tibia_common_slope))
    print(adjusted_ear_means)
    print(ear_vs_female_contrasts)

    cat("\nEYE MODEL TESTS\n")
    print(eye_model_tests)
    print(eye_additive_anova)
    print(tidy_emm(eye_head_common_slope))
    print(adjusted_eye_means)
    print(eye_vs_female_contrasts)
  },
  file = file.path(output_dir, "results_ready_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)

cat("\nANALYSIS COMPLETE\n")
cat("Outputs saved to:\n", output_dir, "\n")
