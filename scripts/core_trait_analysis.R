# ============================================================
# WELLINGTON TREE WETA MORPHOLOGICAL ALLOCATION
# Core univariate analysis workflow
# Females and eighth-, ninth-, and tenth-instar males
# Display labels: Female, Eighth instar, Ninth instar, Tenth instar
# ============================================================
#
# ANALYSIS HIERARCHY
#
# Central analyses retained in the manuscript are the joint repeated-trait
# generalized least-squares leg model, planned size-adjusted contrasts, the tenth-instar male versus female
# anterior-posterior allocation contrast, and body-size-adjusted ear and eye
# means. The one-stage raw-trait weapon-leg model is fitted in
# reviewer_reanalysis.R.
#
# Separate leg models, segment decompositions, local anatomical-coupling
# models and legacy residual-predictor models remain below for the audit trail.
# They are not treated as co-equal tests in the revised manuscript.
#
# ORIGINAL ANALYSIS QUESTIONS
# 1. Do the four developmental groups differ in foreleg,
#    midleg, and hindleg length after accounting for body size?
# 2. Do group differences depend on leg identity (relative
#    locomotory allocation), particularly anterior versus
#    posterior allocation?
# 3. Do allometric slopes differ among groups and leg pairs?
# 4. Is relative head/weapon size positively integrated with
#    leg size, as predicted by secondary sexual trait compensation?
# 5. Do groups differ in ear and eye investment at a common
#    structural body size?
# 6. Are ear and eye size locally integrated with foretibia and
#    head size, respectively, after accounting for group-specific
#    body-size allometry?
# 7. Is sensory morphology negatively associated with relative
#    weapon investment, as expected under a trade-off, or instead
#    positively integrated or developmentally decoupled?
#
# IMPORTANT INTERPRETIVE NOTES
# - Pronotum length is used as the measure of structural body size.
# - All continuous morphological variables are analysed on a log
#   scale. Slopes are therefore allometric exponents.
# - logP_c = 0 corresponds to a pronotum length of 7.2 mm, which
#   lies within the observed range shared by all four groups.
# - Back-transformed EMMs are geometric means / predicted medians,
#   not arithmetic means.
# - The joint repeated-trait GLS model is the primary analysis of
#   relative allocation among the three leg pairs.
# - Morph comparisons of ear and eye size are made at a common
#   pronotum size. Local ear-tibia and eye-head integration is
#   tested using body-adjusted within-group residual variation,
#   avoiding extrapolation to unsupported common tibia/head sizes.
# - Positive phenotypic covariance is consistent with integration
#   or compensation; negative covariance is consistent with a
#   trade-off; a slope near zero is consistent with decoupling.
# ============================================================


# 0. USER SETTINGS -------------------------------------------

# Command-line usage:
# Rscript core_trait_analysis.R [data_file] [output_dir]
args <- commandArgs(trailingOnly = TRUE)

data_file <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "data/trait_data.csv"
}

if (!file.exists(data_file)) {
  stop(
    paste0(
      "Trait data were not found at:\n",
      data_file,
      "\nPass the canonical trait CSV as the first command-line argument."
    ),
    call. = FALSE
  )
}

# Project-level output directory
project_dir <- normalizePath(getwd())
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "analysis_outputs/weta_trait_analysis"
}
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
  "nlme",
  "lme4",
  "lmerTest",
  "performance",
  "see",
  "broom",
  "flextable"
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
  lmer.df = "satterthwaite",
  disable.pbkrtest = TRUE
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

# Lower-case publication labels used whenever sex/morph group is plotted
# on the x-axis. Model and table labels remain unchanged.
group_axis_labels <- c(
  "Female" = "female",
  "Eighth instar" = "eighth instar",
  "Ninth instar" = "ninth instar",
  "Tenth instar" = "tenth instar"
)


# Planned contrasts used repeatedly throughout the analysis.
# Coefficients correspond to group order:
# Female, Eighth instar, Ninth instar, Tenth instar.
strategy_contrasts_full <- list(
  "Eighth instar - Tenth instar" = c(0, 1, 0, -1),
  "Ninth instar - midpoint(Eighth instar,Tenth instar)" =
    c(0, -0.5, 1, -0.5),
  "Tenth instar - Female" = c(-1, 0, 0, 1)
)

# Construct a contrast vector across a group-by-leg emmeans grid.
# Named weights make the code robust to the internal row order.
make_factorial_contrast <- function(
    emm_object,
    group_weights,
    leg_weights
) {
  grid <- as.data.frame(emm_object)

  if (!all(c("group", "leg") %in% names(grid))) {
    stop(
      "The emmeans grid must contain columns named 'group' and 'leg'.",
      call. = FALSE
    )
  }

  group_labels <- as.character(grid$group)
  leg_labels <- as.character(grid$leg)

  if (!all(group_labels %in% names(group_weights))) {
    stop("group_weights is missing one or more group names.", call. = FALSE)
  }
  if (!all(leg_labels %in% names(leg_weights))) {
    stop("leg_weights is missing one or more leg names.", call. = FALSE)
  }

  unname(
    group_weights[group_labels] *
      leg_weights[leg_labels]
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

# Sex-by-morph table is retained as a QC check. Use the same
# publication-ready morph labels as the manuscript outputs.
sex_morph_table <- dat_raw |>
  mutate(
    morph = dplyr::recode(
      morph,
      female = "Female",
      eighth = "Eighth instar",
      ninth = "Ninth instar",
      tenth = "Tenth instar"
    )
  ) |>
  count(sex, morph, name = "n")

cat("\nSEX BY MORPH\n")
print(sex_morph_table)


# 5. CREATE ANALYSIS VARIABLES -------------------------------

dat <- dat_raw |>
  mutate(
    row_id = row_number(),
    # Retain the raw morph codes as factor inputs, but use complete,
    # publication-ready group labels in models, tables, and figures.
    group = factor(
      morph,
      levels = c("female", "eighth", "ninth", "tenth"),
      labels = c("Female", "Eighth instar", "Ninth instar", "Tenth instar")
    ),

    # Ordered male instar number for descriptive use only.
    male_morph_ordered = case_when(
      group == "Eighth instar" ~ 8,
      group == "Ninth instar"  ~ 9,
      group == "Tenth instar"  ~ 10,
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
  group_by(group) |>
  mutate(
    # Within-group centring prevents morph differences in absolute
    # head or foretibia size from being treated as within-group
    # developmental integration.
    log_foretibia_within =
      log_foretibia - mean(log_foretibia, na.rm = TRUE),
    log_head_size_within =
      log_head_size - mean(log_head_size, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(
    # Grand-mean-centred variables are retained only for descriptive
    # checks and compatibility with earlier exploratory models.
    log_foretibia_c = log_foretibia - mean(log_foretibia, na.rm = TRUE),
    log_head_size_c = log_head_size - mean(log_head_size, na.rm = TRUE)
  )

# Confirm that centring worked.
centering_check <- dat |>
  summarise(
    reference_pronotum = reference_pronotum,
    mean_log_foretibia_c = mean(log_foretibia_c, na.rm = TRUE),
    mean_log_head_size_c = mean(log_head_size_c, na.rm = TRUE),
    maximum_absolute_group_mean_foretibia_within =
      max(abs(tapply(log_foretibia_within, group, mean, na.rm = TRUE))),
    maximum_absolute_group_mean_head_within =
      max(abs(tapply(log_head_size_within, group, mean, na.rm = TRUE)))
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



# 7A. SUPPORTING-TRAIT RANGES --------------------------------
# Pronotum overlaps among all groups, whereas head size and
# foretibia length do not. These tables document why sensory-group
# comparisons use common body size rather than a common absolute
# head or foretibia size.

supporting_trait_ranges <- dat |>
  group_by(group) |>
  summarise(
    pronotum_min = min(pronotum, na.rm = TRUE),
    pronotum_max = max(pronotum, na.rm = TRUE),
    head_size_min = min(head_size, na.rm = TRUE),
    head_size_max = max(head_size, na.rm = TRUE),
    foretibia_min = min(foretibia, na.rm = TRUE),
    foretibia_max = max(foretibia, na.rm = TRUE),
    .groups = "drop"
  )

common_support_ranges <- tibble(
  trait = c("pronotum", "head_size", "foretibia"),
  common_lower = c(
    max(supporting_trait_ranges$pronotum_min),
    max(supporting_trait_ranges$head_size_min),
    max(supporting_trait_ranges$foretibia_min)
  ),
  common_upper = c(
    min(supporting_trait_ranges$pronotum_max),
    min(supporting_trait_ranges$head_size_max),
    min(supporting_trait_ranges$foretibia_max)
  )
) |>
  mutate(
    has_common_overlap = common_lower <= common_upper
  )

cat("\nSUPPORTING-TRAIT RANGES\n")
print(supporting_trait_ranges)
print(common_support_ranges)


# 7B. BODY-ADJUSTED RELATIVE TRAIT INDICES -------------------
# Relative head size and relative foretibia size are residuals from
# group-specific allometries. They therefore quantify whether an
# individual has a larger or smaller head/foretibia than expected
# for its pronotum size and group. These variables are used for
# direct tests of compensation, integration, and trade-offs.

mod_head_allometry_common <- lm(
  log_head_size ~ logP_c + group,
  data = dat,
  na.action = na.exclude
)

mod_head_allometry_full <- lm(
  log_head_size ~ logP_c * group,
  data = dat,
  na.action = na.exclude
)

mod_foretibia_allometry_common <- lm(
  log_foretibia ~ logP_c + group,
  data = dat,
  na.action = na.exclude
)

mod_foretibia_allometry_full <- lm(
  log_foretibia ~ logP_c * group,
  data = dat,
  na.action = na.exclude
)

relative_trait_model_tests <- bind_rows(
  nested_f_test(
    mod_head_allometry_common,
    mod_head_allometry_full,
    "Head size: logP_c x group"
  ),
  nested_f_test(
    mod_foretibia_allometry_common,
    mod_foretibia_allometry_full,
    "Foretibia: logP_c x group"
  )
)

head_allometry_trends <- emtrends(
  mod_head_allometry_full,
  ~ group,
  var = "logP_c"
)

foretibia_allometry_trends <- emtrends(
  mod_foretibia_allometry_full,
  ~ group,
  var = "logP_c"
)

relative_trait_allometric_slopes <- bind_rows(
  tidy_emm(head_allometry_trends) |>
    mutate(trait = "Head size"),
  tidy_emm(foretibia_allometry_trends) |>
    mutate(trait = "Foretibia")
) |>
  relocate(trait)

dat <- dat |>
  mutate(
    relative_head_size =
      as.numeric(residuals(mod_head_allometry_full)),
    relative_foretibia_size =
      as.numeric(residuals(mod_foretibia_allometry_full))
  )

relative_size_check <- dat |>
  group_by(group) |>
  summarise(
    mean_relative_head_size =
      mean(relative_head_size, na.rm = TRUE),
    mean_relative_foretibia_size =
      mean(relative_foretibia_size, na.rm = TRUE),
    sd_relative_head_size =
      sd(relative_head_size, na.rm = TRUE),
    sd_relative_foretibia_size =
      sd(relative_foretibia_size, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nRELATIVE-TRAIT ALLOMETRY TESTS\n")
print(relative_trait_model_tests)
print(relative_size_check)


# 8. LEG DATA IN LONG FORMAT ---------------------------------

legs_long <- dat |>
  select(
    row_id, ID, group, pronotum, logP_c,
    relative_head_size,
    foreleg, midleg, hindleg
  ) |>
  pivot_longer(
    cols = c(foreleg, midleg, hindleg),
    names_to = "leg",
    values_to = "leg_length"
  ) |>
  mutate(
    ID = factor(ID),
    leg = factor(
      leg,
      levels = c("foreleg", "midleg", "hindleg"),
      labels = c("Foreleg", "Midleg", "Hindleg")
    ),
    leg_order = as.integer(leg),
    log_leg_length = log(leg_length),
    log_pronotum = log(pronotum)
  ) |>
  arrange(ID, leg_order)

if (anyDuplicated(legs_long[c("ID", "leg")])) {
  stop("Each individual must contribute at most one observation per leg pair.")
}

# Sum coding makes Type-III-style joint tests invariant to the reference level.
contrasts(legs_long$group) <- contr.sum(nlevels(legs_long$group))
contrasts(legs_long$leg) <- contr.sum(nlevels(legs_long$leg))

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
  theme_weta() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(5.5, 5.5, 18, 5.5)
  ) +
  guides(colour = guide_legend(nrow = 1))


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
# Each individual contributes up to three leg measurements. A generalized
# least-squares model estimates a fully unstructured within-individual
# covariance matrix: separate residual variances for the three leg pairs and
# freely estimated correlations among them. This avoids the compound-symmetry
# assumption imposed by a random-intercept-only model. An ID-by-leg random
# effect is not separately identifiable because there is one observation per
# individual-by-leg cell.

mod_leg_allocation <- nlme::gls(
  log_leg_length ~ logP_c * group * leg,
  data = legs_long,
  correlation = nlme::corSymm(form = ~ leg_order | ID),
  weights = nlme::varIdent(form = ~ 1 | leg),
  method = "REML",
  na.action = na.omit,
  control = nlme::glsControl(opt = "optim", msMaxIter = 1000)
)

# Two symmetric values around the centred body-size reference retain tests of
# slope terms while evaluating non-slope effects at pronotum = 7.2 mm.
zero_interval <- function(x) c(-0.01, 0.01)

leg_joint_tests <- emmeans::joint_tests(
  mod_leg_allocation,
  cov.reduce = list(logP_c = zero_interval),
  mode = "appx-satterthwaite"
) |>
  as.data.frame() |>
  tibble::as_tibble()

leg_allocation_anova <- leg_joint_tests |>
  transmute(
    term = `model term`,
    `Sum Sq` = NA_real_,
    `Mean Sq` = NA_real_,
    NumDF = df1,
    DenDF = df2,
    `F value` = F.ratio,
    `Pr(>F)` = p.value
  )

cat("\nJOINT LEG-ALLOCATION MODEL\n")
print(leg_allocation_anova)
print(summary(mod_leg_allocation))

# Document the covariance choice against homogeneous and heterogeneous
# compound-symmetry alternatives fitted to the same fixed-effects structure.
mod_leg_cs <- update(
  mod_leg_allocation,
  correlation = nlme::corCompSymm(form = ~ 1 | ID),
  weights = NULL
)
mod_leg_hetcs <- update(
  mod_leg_allocation,
  correlation = nlme::corCompSymm(form = ~ 1 | ID),
  weights = nlme::varIdent(form = ~ 1 | leg)
)

leg_covariance_model_comparison <- anova(
  mod_leg_cs,
  mod_leg_hetcs,
  mod_leg_allocation
) |>
  as.data.frame() |>
  tibble::rownames_to_column("covariance_model") |>
  tibble::as_tibble()

leg_covariance_parameters <- tibble::tibble(
  parameter = c(
    "Foreleg-midleg residual correlation",
    "Foreleg-hindleg residual correlation",
    "Midleg-hindleg residual correlation",
    "Midleg/foreleg residual SD ratio",
    "Hindleg/foreleg residual SD ratio"
  ),
  estimate = c(
    coef(mod_leg_allocation$modelStruct$corStruct, unconstrained = FALSE),
    coef(mod_leg_allocation$modelStruct$varStruct, unconstrained = FALSE)
  )
)

# Adjusted group means at pronotum = reference_pronotum.
emm_allocation <- emmeans(
  mod_leg_allocation,
  ~ group | leg,
  at = list(logP_c = 0),
  mode = "appx-satterthwaite"
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
) |>
  mutate(Holm_family = paste0("Six all-pairwise group contrasts within ", leg))

leg_vs_female_contrasts <- ratio_contrast_table(
  allocation_vs_female,
  outcome = "Leg length",
  adjust = "holm"
) |>
  mutate(Holm_family = paste0("Three male-versus-female contrasts within ", leg))

# Joint slope estimates and planned decompositions of the significant
# or potentially important three-way interaction.
joint_slopes <- emtrends(
  mod_leg_allocation,
  ~ group * leg,
  var = "logP_c",
  mode = "appx-satterthwaite"
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
) |>
  mutate(Holm_family = paste0("Six group slope contrasts within ", leg))

# Leg slope comparisons within each group: three tests per group.
joint_leg_slope_pairs <- pairs(
  joint_slopes,
  by = "group",
  adjust = "holm"
)

joint_leg_slope_contrasts <- tidy_emm(
  joint_leg_slope_pairs,
  adjust = "holm"
) |>
  mutate(Holm_family = paste0("Three leg slope contrasts within ", group))

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
) |>
  mutate(Holm_family = "All group-by-leg slope difference-in-differences")

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
  scale_x_discrete(labels = group_axis_labels) +
  theme_weta()

# Planned strategy contrasts within each leg pair.
leg_strategy_contrasts_object <- contrast(
  emm_allocation,
  method = strategy_contrasts_full,
  by = "leg",
  adjust = "holm"
)

leg_strategy_contrasts <- ratio_contrast_table(
  leg_strategy_contrasts_object,
  outcome = "Leg length",
  adjust = "holm"
) |>
  mutate(Holm_family = paste0("Three planned strategy contrasts within ", leg))

# Explicit anterior-versus-posterior allocation contrasts.
# A difference-in-differences compares group differences in the
# ratio of anterior leg size to hindleg size. Back-transformation
# yields a ratio of ratios.

emm_allocation_grid <- emmeans(
  mod_leg_allocation,
  ~ group * leg,
  at = list(logP_c = 0),
  mode = "appx-satterthwaite"
)

group_weight_sets <- list(
  "Tenth instar - Eighth instar" =
    c(
      "Female" = 0, "Eighth instar" = -1,
      "Ninth instar" = 0, "Tenth instar" = 1
    ),
  "Ninth instar - midpoint(Eighth instar,Tenth instar)" =
    c(
      "Female" = 0, "Eighth instar" = -0.5,
      "Ninth instar" = 1, "Tenth instar" = -0.5
    ),
  "Tenth instar - Female" =
    c(
      "Female" = -1, "Eighth instar" = 0,
      "Ninth instar" = 0, "Tenth instar" = 1
    )
)

leg_weight_sets <- list(
  "foreleg relative to hindleg" =
    c(Foreleg = 1, Midleg = 0, Hindleg = -1),
  "midleg relative to hindleg" =
    c(Foreleg = 0, Midleg = 1, Hindleg = -1),
  "anterior mean relative to hindleg" =
    c(Foreleg = 0.5, Midleg = 0.5, Hindleg = -1)
)

allocation_difference_methods <- list()

for (group_contrast_name in names(group_weight_sets)) {
  for (leg_contrast_name in names(leg_weight_sets)) {
    contrast_name <- paste(
      group_contrast_name,
      leg_contrast_name,
      sep = ": "
    )

    allocation_difference_methods[[contrast_name]] <-
      make_factorial_contrast(
        emm_allocation_grid,
        group_weights =
          group_weight_sets[[group_contrast_name]],
        leg_weights =
          leg_weight_sets[[leg_contrast_name]]
      )
  }
}

allocation_difference_object <- contrast(
  emm_allocation_grid,
  method = allocation_difference_methods,
  adjust = "holm"
)

allocation_difference_contrasts <- ratio_contrast_table(
  allocation_difference_object,
  outcome = "Anterior-posterior allocation",
  adjust = "holm"
) |>
  rename(
    ratio_of_ratios = ratio,
    ratio_of_ratios_lower = ratio_lower,
    ratio_of_ratios_upper = ratio_upper
  ) |>
  mutate(Holm_family = "Nine planned anterior-posterior difference-in-differences")

# Within-group anterior-to-hindleg ratios are retained for
# description. The among-group difference-in-differences above is
# the more direct test of morph-specific allocation.
allocation_within_group_methods <- list()

for (group_name in levels(dat$group)) {
  group_indicator <- setNames(
    rep(0, length(levels(dat$group))),
    levels(dat$group)
  )
  group_indicator[group_name] <- 1

  allocation_within_group_methods[[
    paste0(group_name, ": anterior mean - hindleg")
  ]] <- make_factorial_contrast(
    emm_allocation_grid,
    group_weights = group_indicator,
    leg_weights =
      c(Foreleg = 0.5, Midleg = 0.5, Hindleg = -1)
  )
}

allocation_within_group_object <- contrast(
  emm_allocation_grid,
  method = allocation_within_group_methods,
  adjust = "holm"
)

allocation_within_group_contrasts <- ratio_contrast_table(
  allocation_within_group_object,
  outcome = "Within-group anterior-hindleg allocation",
  adjust = "holm"
) |>
  mutate(Holm_family = "Four descriptive within-group anterior-hindleg contrasts")


# 10B. DIRECT HEAD-LEG COMPENSATION --------------------------
# Tomkins-style secondary sexual trait compensation predicts that
# individuals with relatively large heads/weapons should also have
# relatively large supporting legs, particularly anterior legs and
# especially in the weapon-specialized tenth-instar morph.

mod_head_leg_compensation <- nlme::gls(
  log_leg_length ~
    logP_c * group * leg +
    relative_head_size * group * leg,
  data = legs_long,
  correlation = nlme::corSymm(form = ~ leg_order | ID),
  weights = nlme::varIdent(form = ~ 1 | leg),
  method = "REML",
  na.action = na.omit,
  control = nlme::glsControl(opt = "optim", msMaxIter = 1000)
)

head_leg_compensation_anova <- emmeans::joint_tests(
  mod_head_leg_compensation,
  cov.reduce = list(
    logP_c = zero_interval,
    relative_head_size = emmeans::make.meanint(1)
  ),
  mode = "appx-satterthwaite"
) |>
  as.data.frame() |>
  tibble::as_tibble() |>
  transmute(
    term = `model term`,
    `Sum Sq` = NA_real_,
    `Mean Sq` = NA_real_,
    NumDF = df1,
    DenDF = df2,
    `F value` = F.ratio,
    `Pr(>F)` = p.value
  )

head_leg_trends <- emtrends(
  mod_head_leg_compensation,
  ~ group * leg,
  var = "relative_head_size",
  mode = "appx-satterthwaite"
)

head_leg_trends_table <- tidy_emm(head_leg_trends)

# Do head-leg integration slopes differ among groups within a leg?
head_leg_group_slope_pairs <- pairs(
  head_leg_trends,
  by = "leg",
  adjust = "holm"
)

head_leg_group_slope_contrasts <- tidy_emm(
  head_leg_group_slope_pairs,
  adjust = "holm"
) |>
  mutate(Holm_family = paste0("Six group slope contrasts within ", leg))

# Do the three leg pairs differ in their head-integration slope
# within each group?
head_leg_leg_slope_pairs <- pairs(
  head_leg_trends,
  by = "group",
  adjust = "holm"
)

head_leg_leg_slope_contrasts <- tidy_emm(
  head_leg_leg_slope_pairs,
  adjust = "holm"
) |>
  mutate(Holm_family = paste0("Three leg slope contrasts within ", group))

# Planned morph contrasts in compensation slopes within each leg.
head_leg_strategy_slope_object <- contrast(
  head_leg_trends,
  method = strategy_contrasts_full,
  by = "leg",
  adjust = "holm"
)

head_leg_strategy_slope_contrasts <- tidy_emm(
  head_leg_strategy_slope_object,
  adjust = "holm"
) |>
  mutate(Holm_family = paste0("Three planned strategy slope contrasts within ", leg))

p_head_leg_compensation <- ggplot(
  legs_long,
  aes(
    x = relative_head_size,
    y = log_leg_length,
    colour = group
  )
) +
  geom_point(alpha = 0.45, na.rm = TRUE) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    na.rm = TRUE
  ) +
  facet_wrap(~ leg, scales = "free_y") +
  labs(
    x = "Relative head size (within-group allometric residual)",
    y = "Log leg length",
    colour = "Group"
  ) +
  theme_weta() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(5.5, 5.5, 30, 5.5)
  ) +
  guides(colour = guide_legend(nrow = 1))


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


# 12. BODY-SIZE-ADJUSTED SENSORY INVESTMENT -----------------
# These are the primary morph-comparison models for sensory traits.
# All groups overlap at pronotum = 7.2 mm, so adjusted means are
# interpolations rather than extrapolations.

# 12A. Ear investment relative to body size ------------------

p_ear_body <- ggplot(
  dat,
  aes(
    x = logP_c,
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
    x = paste0(
      "Log pronotum length centred at ",
      reference_pronotum,
      " mm"
    ),
    y = "Log linearized ear size",
    colour = "Group"
  ) +
  theme_weta()

mod_ear_body_common <- lm(
  log_ear_linear ~ logP_c + group,
  data = dat,
  na.action = na.omit
)

mod_ear_body_full <- lm(
  log_ear_linear ~ logP_c * group,
  data = dat,
  na.action = na.omit
)

ear_body_model_test <- nested_f_test(
  mod_ear_body_common,
  mod_ear_body_full,
  "Ear relative to body size: logP_c x group"
)

# Use the common-slope model for adjusted group means unless the
# group-by-body-size interaction is supported.
ear_body_uses_group_specific_slopes <-
  isTRUE(ear_body_model_test$p_value < alpha)

mod_ear_body_final <- if (ear_body_uses_group_specific_slopes) {
  mod_ear_body_full
} else {
  mod_ear_body_common
}

ear_body_trends <- emtrends(
  mod_ear_body_full,
  ~ group,
  var = "logP_c"
)

ear_body_slopes <- tidy_emm(ear_body_trends)

ear_body_slope_contrasts <- tidy_emm(
  pairs(ear_body_trends, adjust = "holm"),
  adjust = "holm"
)

ear_body_emm <- emmeans(
  mod_ear_body_final,
  ~ group,
  at = list(logP_c = 0)
)

ear_body_group_pairs <- pairs(
  ear_body_emm,
  adjust = "holm"
)

ear_body_vs_female <- contrast(
  ear_body_emm,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "holm"
)

ear_body_strategy_object <- contrast(
  ear_body_emm,
  method = strategy_contrasts_full,
  adjust = "holm"
)

adjusted_ear_body_means <- tidy_emm(ear_body_emm) |>
  mutate(
    predicted_ear_linear = exp(emmean),
    ear_linear_lower = exp(lower.CL),
    ear_linear_upper = exp(upper.CL),
    predicted_ear_area = exp(2 * emmean),
    ear_area_lower = exp(2 * lower.CL),
    ear_area_upper = exp(2 * upper.CL)
  )

ear_body_group_contrasts <- ratio_contrast_table(
  ear_body_group_pairs,
  outcome = "Linearized ear size at common body size",
  adjust = "holm"
) |>
  mutate(
    area_ratio = ratio^2,
    area_percent_difference = 100 * (area_ratio - 1),
    Holm_family = "Six all-pairwise group contrasts for tympanal area"
  )

ear_body_vs_female_contrasts <- ratio_contrast_table(
  ear_body_vs_female,
  outcome = "Linearized ear size at common body size",
  adjust = "holm"
) |>
  mutate(
    area_ratio = ratio^2,
    area_percent_difference = 100 * (area_ratio - 1),
    Holm_family = "Three male-versus-female contrasts for tympanal area"
  )

ear_body_strategy_contrasts <- ratio_contrast_table(
  ear_body_strategy_object,
  outcome = "Linearized ear size at common body size",
  adjust = "holm"
) |>
  mutate(
    area_ratio = ratio^2,
    area_percent_difference = 100 * (area_ratio - 1),
    Holm_family = "Three planned strategy contrasts for tympanal area"
  )

p_adjusted_ear_body <- ggplot(
  adjusted_ear_body_means,
  aes(x = group, y = predicted_ear_area)
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = ear_area_lower, ymax = ear_area_upper),
    width = 0.12
  ) +
  labs(
    x = NULL,
    y = paste0(
      "Predicted ear area at pronotum = ",
      reference_pronotum,
      " mm"
    )
  ) +
  scale_x_discrete(labels = group_axis_labels) +
  theme_weta()


# 12B. Eye investment relative to body size ------------------

p_eye_body <- ggplot(
  dat,
  aes(
    x = logP_c,
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
    x = paste0(
      "Log pronotum length centred at ",
      reference_pronotum,
      " mm"
    ),
    y = "Log eye length",
    colour = "Group"
  ) +
  theme_weta()

mod_eye_body_common <- lm(
  log_eye ~ logP_c + group,
  data = dat,
  na.action = na.omit
)

mod_eye_body_full <- lm(
  log_eye ~ logP_c * group,
  data = dat,
  na.action = na.omit
)

eye_body_model_test <- nested_f_test(
  mod_eye_body_common,
  mod_eye_body_full,
  "Eye relative to body size: logP_c x group"
)

# Use the common-slope model for adjusted group means unless the
# group-by-body-size interaction is supported.
eye_body_uses_group_specific_slopes <-
  isTRUE(eye_body_model_test$p_value < alpha)

mod_eye_body_final <- if (eye_body_uses_group_specific_slopes) {
  mod_eye_body_full
} else {
  mod_eye_body_common
}

eye_body_trends <- emtrends(
  mod_eye_body_full,
  ~ group,
  var = "logP_c"
)

eye_body_slopes <- tidy_emm(eye_body_trends)

eye_body_slope_contrasts <- tidy_emm(
  pairs(eye_body_trends, adjust = "holm"),
  adjust = "holm"
)

eye_body_emm <- emmeans(
  mod_eye_body_final,
  ~ group,
  at = list(logP_c = 0)
)

eye_body_group_pairs <- pairs(
  eye_body_emm,
  adjust = "holm"
)

eye_body_vs_female <- contrast(
  eye_body_emm,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "holm"
)

eye_body_strategy_object <- contrast(
  eye_body_emm,
  method = strategy_contrasts_full,
  adjust = "holm"
)

adjusted_eye_body_means <- tidy_emm(eye_body_emm) |>
  mutate(
    predicted_eye_length = exp(emmean),
    eye_lower = exp(lower.CL),
    eye_upper = exp(upper.CL)
  )

eye_body_group_contrasts <- ratio_contrast_table(
  eye_body_group_pairs,
  outcome = "Eye length at common body size",
  adjust = "holm"
) |>
  mutate(Holm_family = "Six all-pairwise group contrasts for eye length")

eye_body_vs_female_contrasts <- ratio_contrast_table(
  eye_body_vs_female,
  outcome = "Eye length at common body size",
  adjust = "holm"
) |>
  mutate(Holm_family = "Three male-versus-female contrasts for eye length")

eye_body_strategy_contrasts <- ratio_contrast_table(
  eye_body_strategy_object,
  outcome = "Eye length at common body size",
  adjust = "holm"
) |>
  mutate(Holm_family = "Three planned strategy contrasts for eye length")

p_adjusted_eye_body <- ggplot(
  adjusted_eye_body_means,
  aes(x = group, y = predicted_eye_length)
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = eye_lower, ymax = eye_upper),
    width = 0.12
  ) +
  labs(
    x = NULL,
    y = paste0(
      "Predicted eye length at pronotum = ",
      reference_pronotum,
      " mm"
    )
  ) +
  scale_x_discrete(labels = group_axis_labels) +
  theme_weta()


# 13. LOCAL DEVELOPMENTAL INTEGRATION ------------------------
# Relative supporting-trait values are residuals from group-specific
# body-size allometries. Thus, slopes are estimated from observed
# within-group variation and do not require a common absolute head
# or foretibia size across groups.

# 13A. Ear-foretibia integration -----------------------------

p_ear_local <- ggplot(
  dat,
  aes(
    x = relative_foretibia_size,
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
    x = "Relative foretibia size (within-group allometric residual)",
    y = "Log linearized ear size",
    colour = "Group"
  ) +
  theme_weta()

mod_ear_local_common <- lm(
  log_ear_linear ~
    logP_c * group +
    relative_foretibia_size,
  data = dat,
  na.action = na.omit
)

mod_ear_local_full <- lm(
  log_ear_linear ~
    logP_c * group +
    relative_foretibia_size * group,
  data = dat,
  na.action = na.omit
)

ear_local_model_test <- nested_f_test(
  mod_ear_local_common,
  mod_ear_local_full,
  "Ear integration: group x relative foretibia"
)

ear_local_trends <- emtrends(
  mod_ear_local_full,
  ~ group,
  var = "relative_foretibia_size"
)

ear_local_slopes <- tidy_emm(ear_local_trends)

ear_local_slope_contrasts <- tidy_emm(
  pairs(ear_local_trends, adjust = "holm"),
  adjust = "holm"
)

ear_local_strategy_slope_object <- contrast(
  ear_local_trends,
  method = strategy_contrasts_full,
  adjust = "holm"
)

ear_local_strategy_slope_contrasts <- tidy_emm(
  ear_local_strategy_slope_object,
  adjust = "holm"
)


# 13B. Eye-head integration and weapon-eye trade-off ---------
# A positive relative-head slope indicates developmental
# integration; a negative slope supports a local head-eye
# trade-off; a slope near zero is consistent with decoupling.

p_eye_local <- ggplot(
  dat,
  aes(
    x = relative_head_size,
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
    x = "Relative head size (within-group allometric residual)",
    y = "Log eye length",
    colour = "Group"
  ) +
  theme_weta()

mod_eye_local_common <- lm(
  log_eye ~
    logP_c * group +
    relative_head_size,
  data = dat,
  na.action = na.omit
)

mod_eye_local_full <- lm(
  log_eye ~
    logP_c * group +
    relative_head_size * group,
  data = dat,
  na.action = na.omit
)

eye_local_model_test <- nested_f_test(
  mod_eye_local_common,
  mod_eye_local_full,
  "Eye integration: group x relative head size"
)

eye_local_trends <- emtrends(
  mod_eye_local_full,
  ~ group,
  var = "relative_head_size"
)

eye_local_slopes <- tidy_emm(eye_local_trends)

eye_local_slope_contrasts <- tidy_emm(
  pairs(eye_local_trends, adjust = "holm"),
  adjust = "holm"
)

eye_local_strategy_slope_object <- contrast(
  eye_local_trends,
  method = strategy_contrasts_full,
  adjust = "holm"
)

eye_local_strategy_slope_contrasts <- tidy_emm(
  eye_local_strategy_slope_object,
  adjust = "holm"
)


# 14. DIRECT WEAPON-EAR COVARIANCE ---------------------------
# This exploratory model asks whether individuals with relatively
# large heads also have relatively large or small ears, after
# controlling for group-specific body-size allometry. A negative
# slope is consistent with a weapon-sensory trade-off; a positive
# slope indicates integration; a slope near zero suggests
# developmental decoupling.

p_ear_weapon <- ggplot(
  dat,
  aes(
    x = relative_head_size,
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
    x = "Relative head size (within-group allometric residual)",
    y = "Log linearized ear size",
    colour = "Group"
  ) +
  theme_weta()

mod_ear_weapon_common <- lm(
  log_ear_linear ~
    logP_c * group +
    relative_head_size,
  data = dat,
  na.action = na.omit
)

mod_ear_weapon_full <- lm(
  log_ear_linear ~
    logP_c * group +
    relative_head_size * group,
  data = dat,
  na.action = na.omit
)

ear_weapon_model_test <- nested_f_test(
  mod_ear_weapon_common,
  mod_ear_weapon_full,
  "Ear-weapon covariance: group x relative head size"
)

ear_weapon_trends <- emtrends(
  mod_ear_weapon_full,
  ~ group,
  var = "relative_head_size"
)

ear_weapon_slopes <- tidy_emm(ear_weapon_trends)

ear_weapon_slope_contrasts <- tidy_emm(
  pairs(ear_weapon_trends, adjust = "holm"),
  adjust = "holm"
)

ear_weapon_strategy_slope_object <- contrast(
  ear_weapon_trends,
  method = strategy_contrasts_full,
  adjust = "holm"
)

ear_weapon_strategy_slope_contrasts <- tidy_emm(
  ear_weapon_strategy_slope_object,
  adjust = "holm"
)


# 15. MODEL DIAGNOSTICS --------------------------------------
# The joint allocation and head-leg compensation models are the
# primary mixed models. Observation-level diagnostics are saved for
# the principal lm models. Flagging criteria identify observations
# for inspection and do not trigger automatic deletion.

mod_fore_final <- mod_fore_interaction
mod_mid_final  <- mod_mid_common
mod_hind_final <- mod_hind_interaction

safe_collinearity <- function(model, model_name) {
  tryCatch(
    performance::check_collinearity(model) |>
      as.data.frame() |>
      as_tibble() |>
      mutate(
        model = model_name,
        error = NA_character_
      ) |>
      relocate(model),
    error = function(e) {
      tibble(
        model = model_name,
        error = conditionMessage(e)
      )
    }
  )
}

collinearity_model_list <- list(
  foreleg = mod_fore_final,
  midleg = mod_mid_final,
  hindleg = mod_hind_final,
  head_allometry = mod_head_allometry_full,
  foretibia_allometry = mod_foretibia_allometry_full,
  ear_body = mod_ear_body_final,
  eye_body = mod_eye_body_final,
  ear_local = mod_ear_local_full,
  eye_local = mod_eye_local_full,
  ear_weapon = mod_ear_weapon_full
)

collinearity_tables <- purrr::imap_dfr(
  collinearity_model_list,
  ~ safe_collinearity(.x, .y)
)

lm_diagnostic_models <- list(
  foreleg = mod_fore_final,
  midleg = mod_mid_final,
  hindleg = mod_hind_final,
  head_allometry = mod_head_allometry_full,
  foretibia_allometry = mod_foretibia_allometry_full,
  ear_body = mod_ear_body_final,
  eye_body = mod_eye_body_final,
  ear_local = mod_ear_local_full,
  eye_local = mod_eye_local_full,
  ear_weapon = mod_ear_weapon_full
)

lm_diagnostics <- purrr::map(
  lm_diagnostic_models,
  ~ flag_lm_diagnostics(
    get_lm_diagnostics(.x, dat),
    .x
  )
)

lm_diagnostics_flagged <- purrr::map(
  lm_diagnostics,
  ~ filter(.x, flagged)
)

# Save graphical diagnostic panels for lm models.
purrr::iwalk(
  lm_diagnostic_models,
  function(model_object, model_name) {
    save_check_model(
      model_object,
      file.path(
        figures_dir,
        paste0("diagnostics_", model_name, ".png")
      )
    )
  }
)

# Save graphical diagnostics for the two primary mixed models.
save_check_model(
  mod_leg_allocation,
  file.path(figures_dir, "diagnostics_joint_leg_allocation.png")
)

save_check_model(
  mod_head_leg_compensation,
  file.path(figures_dir, "diagnostics_head_leg_compensation.png")
)


# 16. SAVE RESULTS TABLES ------------------------------------

# QC, support, and descriptive outputs.
write_csv(missing_summary, file.path(tables_dir, "qc_missing_values.csv"))
write_csv(raw_range_summary, file.path(tables_dir, "qc_raw_ranges.csv"))
write_csv(nonpositive_summary, file.path(tables_dir, "qc_nonpositive_values.csv"))
write_csv(sex_morph_table, file.path(tables_dir, "qc_sex_by_morph.csv"))
write_csv(group_summary, file.path(tables_dir, "group_descriptive_statistics.csv"))
write_csv(body_size_ranges, file.path(tables_dir, "body_size_ranges.csv"))
write_csv(common_body_size, file.path(tables_dir, "common_body_size_range.csv"))
write_csv(supporting_trait_ranges, file.path(tables_dir, "supporting_trait_ranges.csv"))
write_csv(common_support_ranges, file.path(tables_dir, "common_support_ranges.csv"))
write_csv(centering_check, file.path(tables_dir, "centering_check.csv"))
write_csv(relative_trait_model_tests, file.path(tables_dir, "relative_trait_model_tests.csv"))
write_csv(relative_trait_allometric_slopes, file.path(tables_dir, "relative_trait_allometric_slopes.csv"))
write_csv(relative_size_check, file.path(tables_dir, "relative_size_check.csv"))

# Separate leg analyses.
write_csv(leg_slope_model_tests, file.path(tables_dir, "leg_slope_model_tests.csv"))
write_csv(leg_slopes, file.path(tables_dir, "leg_group_specific_slopes.csv"))
write_csv(leg_slope_contrasts, file.path(tables_dir, "leg_slope_pairwise_contrasts.csv"))

# Primary joint locomotory-allocation model.
write_csv(leg_allocation_anova, file.path(tables_dir, "joint_leg_model_type3_anova.csv"))
write_csv(leg_covariance_model_comparison, file.path(tables_dir, "joint_leg_covariance_model_comparison.csv"))
write_csv(leg_covariance_parameters, file.path(tables_dir, "joint_leg_unstructured_covariance.csv"))
write_csv(adjusted_leg_means, file.path(tables_dir, "joint_adjusted_leg_means.csv"))
write_csv(leg_pairwise_contrasts, file.path(tables_dir, "joint_leg_all_pairwise_contrasts.csv"))
write_csv(leg_vs_female_contrasts, file.path(tables_dir, "joint_leg_males_vs_female.csv"))
write_csv(leg_strategy_contrasts, file.path(tables_dir, "joint_leg_strategy_contrasts.csv"))
write_csv(allocation_difference_contrasts, file.path(tables_dir, "joint_anterior_posterior_difference_contrasts.csv"))
write_csv(allocation_within_group_contrasts, file.path(tables_dir, "joint_within_group_anterior_hindleg_contrasts.csv"))
write_csv(joint_slopes_table, file.path(tables_dir, "joint_leg_slopes.csv"))
write_csv(joint_group_slope_contrasts, file.path(tables_dir, "joint_group_slopes_within_leg.csv"))
write_csv(joint_leg_slope_contrasts, file.path(tables_dir, "joint_leg_slopes_within_group.csv"))
write_csv(three_way_slope_contrasts_table, file.path(tables_dir, "joint_three_way_slope_contrasts.csv"))

# Direct head-leg compensation.
write_csv(head_leg_compensation_anova, file.path(tables_dir, "head_leg_compensation_type3_anova.csv"))
write_csv(head_leg_trends_table, file.path(tables_dir, "head_leg_compensation_slopes.csv"))
write_csv(head_leg_group_slope_contrasts, file.path(tables_dir, "head_leg_group_slopes_within_leg.csv"))
write_csv(head_leg_leg_slope_contrasts, file.path(tables_dir, "head_leg_leg_slopes_within_group.csv"))
write_csv(head_leg_strategy_slope_contrasts, file.path(tables_dir, "head_leg_strategy_slope_contrasts.csv"))

# Segment-level analyses.
write_csv(segment_model_tests, file.path(tables_dir, "segment_slope_model_tests.csv"))
write_csv(segment_slopes, file.path(tables_dir, "segment_group_specific_slopes.csv"))
write_csv(segment_slope_contrasts, file.path(tables_dir, "segment_slope_pairwise_contrasts.csv"))

# Ear investment relative to body size.
write_csv(ear_body_model_test, file.path(tables_dir, "ear_body_slope_model_test.csv"))
write_csv(ear_body_slopes, file.path(tables_dir, "ear_body_group_specific_slopes.csv"))
write_csv(ear_body_slope_contrasts, file.path(tables_dir, "ear_body_slope_contrasts.csv"))
write_csv(adjusted_ear_body_means, file.path(tables_dir, "ear_body_adjusted_group_means.csv"))
write_csv(ear_body_group_contrasts, file.path(tables_dir, "ear_body_all_pairwise_group_contrasts.csv"))
write_csv(ear_body_vs_female_contrasts, file.path(tables_dir, "ear_body_males_vs_female.csv"))
write_csv(ear_body_strategy_contrasts, file.path(tables_dir, "ear_body_strategy_contrasts.csv"))

# Eye investment relative to body size.
write_csv(eye_body_model_test, file.path(tables_dir, "eye_body_slope_model_test.csv"))
write_csv(eye_body_slopes, file.path(tables_dir, "eye_body_group_specific_slopes.csv"))
write_csv(eye_body_slope_contrasts, file.path(tables_dir, "eye_body_slope_contrasts.csv"))
write_csv(adjusted_eye_body_means, file.path(tables_dir, "eye_body_adjusted_group_means.csv"))
write_csv(eye_body_group_contrasts, file.path(tables_dir, "eye_body_all_pairwise_group_contrasts.csv"))
write_csv(eye_body_vs_female_contrasts, file.path(tables_dir, "eye_body_males_vs_female.csv"))
write_csv(eye_body_strategy_contrasts, file.path(tables_dir, "eye_body_strategy_contrasts.csv"))

# Local sensory integration.
write_csv(ear_local_model_test, file.path(tables_dir, "ear_local_integration_model_test.csv"))
write_csv(ear_local_slopes, file.path(tables_dir, "ear_local_integration_slopes.csv"))
write_csv(ear_local_slope_contrasts, file.path(tables_dir, "ear_local_integration_slope_contrasts.csv"))
write_csv(ear_local_strategy_slope_contrasts, file.path(tables_dir, "ear_local_strategy_slope_contrasts.csv"))

write_csv(eye_local_model_test, file.path(tables_dir, "eye_local_integration_model_test.csv"))
write_csv(eye_local_slopes, file.path(tables_dir, "eye_local_integration_slopes.csv"))
write_csv(eye_local_slope_contrasts, file.path(tables_dir, "eye_local_integration_slope_contrasts.csv"))
write_csv(eye_local_strategy_slope_contrasts, file.path(tables_dir, "eye_local_strategy_slope_contrasts.csv"))

# Direct weapon-ear covariance.
write_csv(ear_weapon_model_test, file.path(tables_dir, "ear_weapon_covariance_model_test.csv"))
write_csv(ear_weapon_slopes, file.path(tables_dir, "ear_weapon_covariance_slopes.csv"))
write_csv(ear_weapon_slope_contrasts, file.path(tables_dir, "ear_weapon_covariance_slope_contrasts.csv"))
write_csv(ear_weapon_strategy_slope_contrasts, file.path(tables_dir, "ear_weapon_strategy_slope_contrasts.csv"))

# Diagnostics.
write_csv(collinearity_tables, file.path(tables_dir, "model_collinearity.csv"))

purrr::iwalk(
  lm_diagnostics,
  ~ write_csv(
    .x,
    file.path(tables_dir, paste0("diagnostics_", .y, "_all.csv"))
  )
)

purrr::iwalk(
  lm_diagnostics_flagged,
  ~ write_csv(
    .x,
    file.path(tables_dir, paste0("diagnostics_", .y, "_flagged.csv"))
  )
)


# 17. SAVE FIGURES -------------------------------------------

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
  file.path(figures_dir, "head_leg_compensation.png"),
  p_head_leg_compensation,
  width = 9,
  height = 7.5,
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
  file.path(figures_dir, "ear_body_allometry.png"),
  p_ear_body,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "ear_body_adjusted_means.png"),
  p_adjusted_ear_body,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "eye_body_allometry.png"),
  p_eye_body,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "eye_body_adjusted_means.png"),
  p_adjusted_eye_body,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "ear_local_foretibia_integration.png"),
  p_ear_local,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "eye_local_head_integration.png"),
  p_eye_local,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "ear_weapon_covariance.png"),
  p_ear_weapon,
  width = 7,
  height = 5,
  dpi = 300
)

# Confirm that every table and figure required by the QMD report
# exists before the analysis is declared complete.
report_required_outputs <- c(
  file.path(tables_dir, "group_descriptive_statistics.csv"),
  file.path(tables_dir, "common_body_size_range.csv"),
  file.path(tables_dir, "joint_leg_model_type3_anova.csv"),
  file.path(tables_dir, "joint_adjusted_leg_means.csv"),
  file.path(tables_dir, "joint_leg_males_vs_female.csv"),
  file.path(tables_dir, "joint_leg_strategy_contrasts.csv"),
  file.path(tables_dir, "joint_anterior_posterior_difference_contrasts.csv"),
  file.path(tables_dir, "joint_leg_slopes.csv"),
  file.path(tables_dir, "joint_group_slopes_within_leg.csv"),
  file.path(tables_dir, "joint_three_way_slope_contrasts.csv"),
  file.path(tables_dir, "segment_slope_model_tests.csv"),
  file.path(tables_dir, "head_leg_compensation_type3_anova.csv"),
  file.path(tables_dir, "head_leg_compensation_slopes.csv"),
  file.path(tables_dir, "head_leg_strategy_slope_contrasts.csv"),
  file.path(tables_dir, "ear_body_slope_model_test.csv"),
  file.path(tables_dir, "ear_body_adjusted_group_means.csv"),
  file.path(tables_dir, "ear_body_strategy_contrasts.csv"),
  file.path(tables_dir, "eye_body_slope_model_test.csv"),
  file.path(tables_dir, "eye_body_adjusted_group_means.csv"),
  file.path(tables_dir, "eye_body_strategy_contrasts.csv"),
  file.path(tables_dir, "ear_local_integration_model_test.csv"),
  file.path(tables_dir, "ear_local_integration_slopes.csv"),
  file.path(tables_dir, "eye_local_integration_model_test.csv"),
  file.path(tables_dir, "eye_local_integration_slopes.csv"),
  file.path(tables_dir, "ear_weapon_covariance_model_test.csv"),
  file.path(tables_dir, "ear_weapon_covariance_slopes.csv"),
  file.path(figures_dir, "adjusted_leg_means.png"),
  file.path(figures_dir, "leg_allometries.png"),
  file.path(figures_dir, "head_leg_compensation.png"),
  file.path(figures_dir, "ear_body_adjusted_means.png"),
  file.path(figures_dir, "eye_body_adjusted_means.png")
)

missing_report_outputs <- report_required_outputs[
  !file.exists(report_required_outputs)
]

if (length(missing_report_outputs) > 0) {
  stop(
    paste0(
      "The analysis completed, but the following QMD inputs are missing:\n",
      paste(missing_report_outputs, collapse = "\n")
    ),
    call. = FALSE
  )
}


# 18. SAVE MODEL OBJECTS AND TEXT SUMMARY --------------------

model_objects <- list(
  head_allometry_common = mod_head_allometry_common,
  head_allometry_full = mod_head_allometry_full,
  foretibia_allometry_common = mod_foretibia_allometry_common,
  foretibia_allometry_full = mod_foretibia_allometry_full,

  fore_common = mod_fore_common,
  fore_interaction = mod_fore_interaction,
  mid_common = mod_mid_common,
  mid_interaction = mod_mid_interaction,
  hind_common = mod_hind_common,
  hind_interaction = mod_hind_interaction,
  joint_leg_allocation = mod_leg_allocation,
  joint_leg_compound_symmetry = mod_leg_cs,
  joint_leg_heterogeneous_compound_symmetry = mod_leg_hetcs,
  head_leg_compensation = mod_head_leg_compensation,
  segment_common = segment_models_common,
  segment_interaction = segment_models_interaction,

  ear_body_common = mod_ear_body_common,
  ear_body_full = mod_ear_body_full,
  ear_body_final = mod_ear_body_final,
  eye_body_common = mod_eye_body_common,
  eye_body_full = mod_eye_body_full,
  eye_body_final = mod_eye_body_final,

  ear_local_common = mod_ear_local_common,
  ear_local_full = mod_ear_local_full,
  eye_local_common = mod_eye_local_common,
  eye_local_full = mod_eye_local_full,

  ear_weapon_common = mod_ear_weapon_common,
  ear_weapon_full = mod_ear_weapon_full
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

    cat("REFERENCE PRONOTUM AND COMMON SUPPORT\n")
    print(reference_pronotum)
    print(common_body_size)
    print(common_support_ranges)

    cat("\nGROUP DESCRIPTIVES\n")
    print(group_summary)

    cat("\nRELATIVE-TRAIT ALLOMETRIES\n")
    print(relative_trait_model_tests)
    print(relative_trait_allometric_slopes)
    print(relative_size_check)

    cat("\nSEPARATE LEG SLOPE MODEL TESTS\n")
    print(leg_slope_model_tests)

    cat("\nSEPARATE LEG SLOPES\n")
    print(leg_slopes)

    cat("\nJOINT LEG MODEL TYPE III ANOVA\n")
    print(leg_allocation_anova)

    cat("\nJOINT ADJUSTED LEG MEANS\n")
    print(adjusted_leg_means)

    cat("\nJOINT STRATEGY CONTRASTS WITHIN LEG\n")
    print(leg_strategy_contrasts)

    cat("\nANTERIOR-POSTERIOR ALLOCATION DIFFERENCES\n")
    print(allocation_difference_contrasts)

    cat("\nJOINT LEG SLOPES\n")
    print(joint_slopes_table)

    cat("\nHEAD-LEG COMPENSATION TYPE III ANOVA\n")
    print(head_leg_compensation_anova)

    cat("\nHEAD-LEG COMPENSATION SLOPES\n")
    print(head_leg_trends_table)

    cat("\nHEAD-LEG PLANNED STRATEGY SLOPE CONTRASTS\n")
    print(head_leg_strategy_slope_contrasts)

    cat("\nSEGMENT MODEL TESTS\n")
    print(segment_model_tests)

    cat("\nEAR INVESTMENT RELATIVE TO BODY SIZE\n")
    cat(
      "Adjusted means model:",
      if (ear_body_uses_group_specific_slopes) {
        "group-specific slopes"
      } else {
        "common slope"
      },
      "\n"
    )
    print(ear_body_model_test)
    print(ear_body_slopes)
    print(adjusted_ear_body_means)
    print(ear_body_strategy_contrasts)

    cat("\nEAR-FORETIBIA LOCAL INTEGRATION\n")
    print(ear_local_model_test)
    print(ear_local_slopes)
    print(ear_local_strategy_slope_contrasts)

    cat("\nEYE INVESTMENT RELATIVE TO BODY SIZE\n")
    cat(
      "Adjusted means model:",
      if (eye_body_uses_group_specific_slopes) {
        "group-specific slopes"
      } else {
        "common slope"
      },
      "\n"
    )
    print(eye_body_model_test)
    print(eye_body_slopes)
    print(adjusted_eye_body_means)
    print(eye_body_strategy_contrasts)

    cat("\nEYE-HEAD LOCAL INTEGRATION / TRADE-OFF\n")
    print(eye_local_model_test)
    print(eye_local_slopes)
    print(eye_local_strategy_slope_contrasts)

    cat("\nDIRECT WEAPON-EAR COVARIANCE\n")
    print(ear_weapon_model_test)
    print(ear_weapon_slopes)
    print(ear_weapon_strategy_slope_contrasts)
  },
  file = file.path(output_dir, "results_ready_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)

cat("\nANALYSIS COMPLETE\n")
cat("Outputs saved to:\n", output_dir, "\n")
