#!/usr/bin/env Rscript

# Instar-matched decomposition and precision reporting.
#
# This is a fourth stage of the manuscript pipeline. It runs after
# core_trait_analysis.R, reviewer_reanalysis.R and
# female_reference_multivariate_reanalysis.R, and it adds only the two
# analyses that are not already produced upstream:
#
#   1. INSTAR-MATCHED DECOMPOSITION
#      The male morphs differ in weapon size, but they also differ in the
#      number of juvenile instars completed before maturity. Within males,
#      instar number and weapon expression change together, so the
#      among-morph gradient cannot separate them. Females mature at the
#      tenth instar, as tenth-instar males do, but do not express the
#      exaggerated male head and mandibles. Tenth-instar males and females
#      are therefore matched on instar number while differing in weapon
#      expression, which holds terminal instar number constant while leaving sex-associated
#      differences and chronological development time unresolved.
#      Fitting instar number and sex as separate terms decomposes each
#      trait's divergence into a per-instar component and an
#      instar-matched sex contrast, and compares their coefficients directly.
#
#   2. PRECISION AND MINIMUM DETECTABLE EFFECTS
#      A null interaction is only interpretable alongside the smallest
#      effect the design could have detected. For every planned contrast
#      the standard error, the 95% confidence interval and the minimum
#      detectable effect at 80% power are reported.
#
# The dimensionality, PCA, parallel-analysis and covariance-matrix results
# are NOT recomputed here. Those are produced upstream by
# reviewer_reanalysis.R and female_reference_multivariate_reanalysis.R and
# are read from their tables. This stage only consolidates the covariance
# sign audit from the upstream correlation table, so that the structural
# results are reported from a single set of upstream estimates.
#
# Usage:
#   Rscript scripts/instar_matched_precision_analysis.R [data_file] [output_dir]
#
# IMPORTANT INTERPRETIVE NOTES
# - Pronotum length is the measure of structural body size and is centred
#   on 7.2 mm (logP_c = 0), which lies within the common observed support
#   of all four groups.
# - All continuous morphological variables are analysed on a log scale, so
#   slopes are allometric exponents and differences are ratios.
# - These coefficients describe a within-male per-instar gradient and an
#   instar-matched sex contrast. They do not partition causal contributions
#   of development and weaponry. Legacy weapon_component column names are
#   retained for compatibility with archived outputs only.
# - A dominant axis among group means and positive within-group covariance
#   describe different levels of variation; neither identifies a mechanism.

# Resolve defaults from this script, including when sourced from another directory.
.script_file <- local({
  source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(source_files)) {
    tail(source_files, 1L)[[1L]]
  } else {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg) != 1L) stop("Run this file with Rscript or source().")
    path <- sub("^--file=", "", file_arg[[1L]])
    # Rscript may encode spaces in --file as "~+~".
    if (!file.exists(path)) path <- gsub("~+~", " ", path, fixed = TRUE)
    path
  }
})
.script_file <- normalizePath(.script_file, mustWork = TRUE)
source(file.path(dirname(.script_file), "workflow_helpers.R"), local = TRUE)
project_root <- dirname(dirname(.script_file))

args <- commandArgs(trailingOnly = TRUE)

data_file  <- if (length(args) >= 1L) args[[1L]] else file.path(project_root, "data", "trait_data.csv")
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  file.path(project_root, "analysis_outputs", "weta_trait_analysis")
}
if (length(args) > 2L) {
  stop("This stage does not resample; supply only data_file and output_dir.")
}

if (!file.exists(data_file)) {
  stop("Trait data not found: ", data_file)
}

# The structural analyses are mandatory upstream products: the sign audit
# below is derived from the upstream correlation table, so this stage
# fails loudly rather than silently reporting nothing if they are absent.
tables_dir <- file.path(output_dir, "tables")

upstream_correlations <- file.path(tables_dir, "group_trait_correlations.csv")

if (!file.exists(upstream_correlations)) {
  stop(
    "Upstream correlation table not found: ", upstream_correlations,
    "\nRun the earlier pipeline stages first (see run_all.R)."
  )
}

reference_pronotum <- 7.2
alpha <- 0.05
# emmeans may use simulated perturbations for approximate GLS degrees of freedom.
set.seed(2302)


# 1. PACKAGES ----------------------------------------------------------------

required_packages <- c("tidyverse", "nlme", "emmeans")

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


# 2. HELPER FUNCTIONS --------------------------------------------------------

# Residual degrees of freedom for a gls fit. df.residual() is not reliable
# for these objects, so the value is taken from the fitted dimensions.
gls_df <- function(model) model$dims$N - model$dims$p

# Minimum detectable effect for a single contrast estimate: the smallest
# true effect detectable with the stated power at alpha, given the observed
# standard error and denominator degrees of freedom. The significance level
# is passed in explicitly rather than read from a same-named global, so the
# helper stays self-contained.
min_detectable_effect <- function(se, df, power = 0.80, alpha_level = 0.05) {
  crit <- qt(1 - alpha_level / 2, df) + qt(power, df)
  crit * se
}

# Back-transform a log-scale minimum detectable effect to a percentage.
mde_percent <- function(se, df) 100 * (exp(min_detectable_effect(se, df)) - 1)

# emmeans matches a user-supplied contrast vector to the grid POSITIONALLY,
# not by the names on the vector. A vector written in a different order than
# the grid is therefore applied to the wrong groups without any warning,
# which silently scrambles the contrast. group_contrast() re-orders the
# supplied coefficients into the factor level order that the grid uses.
group_contrast <- function(vals, group_levels) {
  missing_names <- setdiff(group_levels, names(vals))

  if (length(missing_names) > 0) {
    stop(
      paste("Contrast is missing coefficients for:",
            paste(missing_names, collapse = ", ")),
      call. = FALSE
    )
  }

  vals[group_levels]
}

# Verify in the fitted model that the named groups really differ as the
# decomposition assumes, so that a mislabelled factor cannot pass.
# gls_df() is defined above so that precision figures are well defined for
# generalized least-squares fits as well as for lm fits.


# 3. IMPORT AND PREPARE DATA -------------------------------------------------

dat_raw <- readr::read_csv(data_file, na = c("", "NA"), show_col_types = FALSE)

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

dat <- dat_raw |>
  mutate(
    group = factor(
      morph,
      levels = c("female", "eighth", "ninth", "tenth"),
      labels = c("Female", "Eighth instar", "Ninth instar", "Tenth instar")
    ),
    # Number of juvenile instars completed before maturity. Females mature
    # at the tenth instar, as tenth-instar males do; this is the basis of
    # the instar-matched comparison.
    instar = case_when(
      morph == "eighth" ~ 8L,
      morph == "ninth"  ~ 9L,
      morph == "tenth"  ~ 10L,
      morph == "female" ~ 10L
    ),
    sex = if_else(morph == "female", "female", "male"),
    foreleg = forefemur + foretibia,
    midleg  = midfemur + midtibia,
    hindleg = hindfemur + hindtibia,
    ear_lin = sqrt(ear),
    log_head = log(sqrt(head_length * head_width)),
    logP_c = log(pronotum) - log(reference_pronotum)
  ) |>
  mutate(
    across(
      c(foreleg, midleg, hindleg, ear_lin, eye),
      ~ log(.x),
      .names = "l_{.col}"
    )
  )

# Identifiability of the decomposition requires two things. Instar number
# must vary within males, so that the per-instar component is estimable from
# the male morphs; and sex must vary within the tenth instar, so that the
# weapon component is estimable from the instar-matched contrast. Females are
# all tenth instar by design, which is what makes them the instar-matched
# reference, so instar number is deliberately NOT required to vary within
# females. If a future revision of the data breaks either condition the
# script stops, because the two components would no longer be separately
# estimable.
male_instar_levels <- dat |> filter(sex == "male") |> distinct(instar) |> nrow()

if (male_instar_levels < 2L) {
  stop("Instar number does not vary within males; the per-instar component is not estimable.",
       call. = FALSE)
}

tenth_sexes <- dat |> filter(instar == 10L) |> distinct(sex) |> nrow()

if (tenth_sexes < 2L) {
  stop("Sex does not vary within the tenth instar; the weapon component is not estimable.",
       call. = FALSE)
}

female_instars <- dat |> filter(sex == "female") |> distinct(instar) |> pull(instar)

if (!identical(sort(female_instars), 10L)) {
  stop(
    "Females are not exclusively tenth instar, so the female group no longer ",
    "supplies an instar-matched reference.",
    call. = FALSE
  )
}

cat("\nDECOMPOSITION IDENTIFIABILITY CHECK\n")
cat("Instar number varies within males:", male_instar_levels, "levels.\n")
cat("Sex varies within the tenth instar:", tenth_sexes, "sexes.\n")
cat("Females are exclusively tenth instar, matching tenth-instar males.\n")
print(dat |> count(group, instar, sex, name = "n"))

trait_cols <- c("l_foreleg", "l_midleg", "l_hindleg", "l_ear_lin", "l_eye")


# ============================================================
# 4. INSTAR-MATCHED DECOMPOSITION
# ============================================================

# The sex coefficient is male minus female (female is the reference level),
# so it is already expressed as male relative to female, matching the sign
# convention used elsewhere in the manuscript.

decompose_trait <- function(trait_label) {
  dat_t <- dat |>
    filter(!is.na(.data[[trait_label]]), !is.na(logP_c)) |>
    mutate(instar_c = instar - 10)

  fit <- lm(
    reformulate(c("logP_c", "instar_c", "sex"), response = trait_label),
    data = dat_t
  )

  co <- summary(fit)$coefficients
  df_resid <- df.residual(fit)

  instar_est <- co["instar_c", "Estimate"]
  sex_est    <- co["sexmale", "Estimate"]
  instar_se  <- co["instar_c", "Std. Error"]
  sex_se     <- co["sexmale", "Std. Error"]

  # The two components are compared directly, by testing their difference,
  # rather than by comparing their p-values. Coerced to plain numerics so
  # that the results frame stays writable to CSV.
  V <- vcov(fit)
  diff_est <- as.numeric(instar_est - sex_est)
  diff_se  <- as.numeric(sqrt(V["instar_c", "instar_c"] +
                                V["sexmale", "sexmale"] -
                                2 * V["instar_c", "sexmale"]))
  diff_t   <- diff_est / diff_se
  diff_p   <- 2 * pt(abs(diff_t), df_resid, lower.tail = FALSE)

  tibble(
    trait = trait_label,
    n = nobs(fit),
    df_residual = df_resid,
    instar_component_log_per_instar = instar_est,
    instar_component_SE = instar_se,
    instar_component_CI_lower = instar_est - qt(0.975, df_resid) * instar_se,
    instar_component_CI_upper = instar_est + qt(0.975, df_resid) * instar_se,
    instar_component_pct_per_instar = 100 * (exp(instar_est) - 1),
    instar_component_p = co["instar_c", "Pr(>|t|)"],
    weapon_component_log = sex_est,
    weapon_component_SE = sex_se,
    weapon_component_CI_lower = sex_est - qt(0.975, df_resid) * sex_se,
    weapon_component_CI_upper = sex_est + qt(0.975, df_resid) * sex_se,
    weapon_component_ratio_male_vs_female = exp(sex_est),
    weapon_component_pct_male_vs_female = 100 * (exp(sex_est) - 1),
    weapon_component_p = co["sexmale", "Pr(>|t|)"],
    component_difference_log = diff_est,
    component_difference_SE = diff_se,
    component_difference_t = diff_t,
    component_difference_p = diff_p,
    instar_component_min_detectable_pct = mde_percent(instar_se, df_resid),
    weapon_component_min_detectable_pct = mde_percent(sex_se, df_resid)
  )
}

trait_decomposition <- map_dfr(trait_cols, decompose_trait) |>
  mutate(
    trait = str_remove(trait, "^l_"),
    instar_component_holm_p = p.adjust(instar_component_p, method = "holm"),
    weapon_component_holm_p = p.adjust(weapon_component_p, method = "holm"),
    component_difference_holm_p = p.adjust(component_difference_p, method = "holm"),
    larger_component = if_else(component_difference_log < 0,
                               "weapon", "per_instar")
  )

# Joint decomposition across the three leg pairs, with the unstructured
# within-individual covariance and separate residual variances per leg pair.
# This is the multivariate analogue of the per-trait decomposition and
# respects the non-independence of three legs measured on one individual.

leg_long <- dat |>
  select(ID, group, sex, instar, logP_c, l_foreleg, l_midleg, l_hindleg) |>
  pivot_longer(
    cols = c(l_foreleg, l_midleg, l_hindleg),
    names_to = "leg",
    values_to = "log_leg"
  ) |>
  mutate(
    leg = factor(
      leg,
      levels = c("l_foreleg", "l_midleg", "l_hindleg"),
      labels = c("Foreleg", "Midleg", "Hindleg")
    ),
    instar_c = instar - 10
  ) |>
  filter(!is.na(log_leg))

leg_joint <- gls(
  log_leg ~ (logP_c + instar_c + sex) * leg,
  data = leg_long,
  correlation = corSymm(form = ~ as.integer(leg) | ID),
  weights = varIdent(form = ~ 1 | leg),
  method = "REML",
  control = glsControl(opt = "optim", maxIter = 200, msMaxIter = 200)
)

joint_df <- gls_df(leg_joint)

joint_components <- emmeans::emtrends(leg_joint, ~ leg, var = "instar_c", infer = c(TRUE, TRUE)) |>
  as_tibble() |>
  transmute(
    leg,
    instar_component_log_per_instar = instar_c.trend,
    instar_component_SE = SE,
    instar_component_df = df,
    instar_component_pct_per_instar = 100 * (exp(instar_c.trend) - 1),
    instar_component_p = p.value
  )

joint_weapon <- emmeans::emmeans(leg_joint, ~ sex | leg) |>
  pairs() |>
  as_tibble() |>
  transmute(
    leg,
    # pairs() returns female - male, so the contrast is negated to give
    # male relative to female, matching the per-trait decomposition and the
    # sign convention used elsewhere in the manuscript.
    weapon_component_log = -estimate,
    weapon_component_SE = SE,
    weapon_component_df = df,
    weapon_component_p = p.value,
    weapon_component_pct_male_vs_female = 100 * (exp(-estimate) - 1)
  )

joint_table <- left_join(joint_components, joint_weapon, by = "leg") |>
  mutate(
    df_residual = joint_df,
    instar_min_detectable_pct = mde_percent(instar_component_SE, instar_component_df),
    weapon_min_detectable_pct = mde_percent(weapon_component_SE, weapon_component_df),
    larger_component = if_else(
      abs(instar_component_log_per_instar) > abs(weapon_component_log),
      "per_instar", "weapon"
    )
  )

cat("\nINSTAR-MATCHED DECOMPOSITION BY TRAIT\n")
print(trait_decomposition |>
        select(trait, n,
               instar_component_pct_per_instar, instar_component_holm_p,
               weapon_component_pct_male_vs_female, weapon_component_holm_p,
               component_difference_holm_p, larger_component),
      n = 10)

cat("\nJOINT DECOMPOSITION ACROSS LEG PAIRS\n")
print(joint_table)


# ============================================================
# 5. PRECISION AND MINIMUM DETECTABLE EFFECTS
# ============================================================

# Use the upstream contrast and slope tables so estimates, standard errors
# and approximate Satterthwaite degrees of freedom agree exactly with the
# main models. Precision intervals are pointwise (unadjusted), and the MDE
# describes a single two-sided test at alpha = 0.05, not a Holm family.
read_upstream <- function(filename) {
  path <- file.path(tables_dir, filename)
  if (!file.exists(path)) stop("Missing upstream precision input: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

precision_leg <- read_upstream("joint_leg_strategy_contrasts.csv") |>
  mutate(
    reverse = contrast == "Eighth instar - Tenth instar",
    estimate = if_else(reverse, -estimate, estimate),
    contrast = if_else(reverse, "Tenth instar - Eighth instar", contrast),
    contrast = str_replace(contrast, "instar,Tenth", "instar, Tenth")
  ) |>
  transmute(
    leg, contrast, estimate, SE, df,
    ratio = exp(estimate),
    percent_difference = 100 * (ratio - 1),
    ci_lower = exp(estimate - qt(0.975, df) * SE),
    ci_upper = exp(estimate + qt(0.975, df) * SE),
    min_detectable_percent = mde_percent(SE, df),
    contrast_type = "planned"
  )

allocation_labels <- c(
  "Tenth instar - Female: foreleg relative to hindleg",
  "Tenth instar - Female: midleg relative to hindleg",
  "Tenth instar - Eighth instar: foreleg relative to hindleg"
)
precision_allocation <- read_upstream(
  "joint_anterior_posterior_difference_contrasts.csv"
) |>
  filter(contrast %in% allocation_labels) |>
  transmute(
    leg = NA_character_, contrast, estimate, SE, df,
    ratio = exp(estimate),
    percent_difference = 100 * (ratio - 1),
    ci_lower = exp(estimate - qt(0.975, df) * SE),
    ci_upper = exp(estimate + qt(0.975, df) * SE),
    min_detectable_percent = mde_percent(SE, df),
    contrast_type = "difference_in_differences"
  )
stopifnot(nrow(precision_leg) == 9L, nrow(precision_allocation) == 3L)
precision_leg_table <- bind_rows(precision_leg, precision_allocation)

weapon_leg_slopes <- read_upstream("head_leg_raw_slopes.csv") |>
  transmute(
    leg, group,
    weapon_leg_slope = log_head_size.trend,
    SE, df,
    ci_lower = lower.CL,
    ci_upper = upper.CL,
    excludes_zero = (ci_lower > 0 | ci_upper < 0),
    min_detectable_slope = min_detectable_effect(SE, df)
  )

# Male weapon-ear slope precision. The male slopes are imprecise; their
# minimum detectable magnitude shows the range of covariances the data could
# actually have excluded.
ear_weapon_model <- lm(l_ear_lin ~ logP_c * group + log_head * group, data = dat)
ear_df <- df.residual(ear_weapon_model)

ear_weapon_slopes <- emmeans::emtrends(ear_weapon_model, ~ group, var = "log_head",
                            infer = c(TRUE, TRUE)) |>
  as_tibble() |>
  transmute(
    group,
    weapon_ear_slope = log_head.trend,
    SE,
    df = ear_df,
    ci_lower = log_head.trend - qt(0.975, ear_df) * SE,
    ci_upper = log_head.trend + qt(0.975, ear_df) * SE,
    excludes_zero = (ci_lower > 0 | ci_upper < 0),
    min_detectable_slope = min_detectable_effect(SE, ear_df)
  )

male_min_detectable <- ear_weapon_slopes |>
  filter(group != "Female") |>
  summarise(
    n_male_morphs = n(),
    mean_min_detectable_slope = mean(min_detectable_slope),
    median_min_detectable_slope = median(min_detectable_slope),
    max_min_detectable_slope = max(min_detectable_slope),
    min_observed_male_slope = min(weapon_ear_slope),
    max_observed_male_slope = max(weapon_ear_slope)
  )

cat("\nPLANNED LEG CONTRAST PRECISION\n")
print(precision_leg_table)
cat("\nWEAPON-LEG SLOPE PRECISION\n")
print(weapon_leg_slopes, n = 20)
cat("\nWEAPON-EAR SLOPE PRECISION\n")
print(ear_weapon_slopes)
cat("\nSMALLEST MALE WEAPON-EAR SLOPE DISTINGUISHABLE FROM ZERO\n")
print(male_min_detectable)


# ============================================================
# 6. COVARIANCE SIGN AUDIT (consolidated from upstream)
# ============================================================

# The within-group correlations are computed upstream, on group-specific
# size-adjusted residuals. They are consolidated here rather than
# recomputed, so that a single set of upstream estimates is reported.
# This is a descriptive audit of phenotypic covariance. Positive
# associations cannot by themselves demonstrate the absence of a cost,
# because variation among individuals in resource acquisition generates
# positive covariance among costly traits even when those traits compete
# for allocation (van Noordwijk & de Jong 1986).

upstream_cor <- readr::read_csv(upstream_correlations, show_col_types = FALSE)

# The multivariate analyses require a complete data matrix, so the audit
# reports the complete-case count for the five non-head traits rather than
# the full sample size.
integration_traits <- c("foreleg", "midleg", "hindleg", "ear_lin", "eye")

group_sizes <- dat |>
  filter(if_all(all_of(integration_traits), ~ !is.na(.x)), !is.na(logP_c)) |>
  count(group, name = "n_individuals")

sign_audit <- upstream_cor |>
  filter(trait_1 != trait_2) |>
  mutate(
    trait_pair = if_else(trait_1 < trait_2,
                         paste(trait_1, trait_2, sep = " - "),
                         paste(trait_2, trait_1, sep = " - "))
  ) |>
  distinct(group, trait_pair, .keep_all = TRUE) |>
  group_by(group) |>
  summarise(
    n_pairwise = n(),
    n_positive = sum(correlation > 0),
    n_negative = sum(correlation < 0),
    minimum_correlation = min(correlation),
    maximum_correlation = max(correlation),
    mean_correlation = mean(correlation),
    median_correlation = median(correlation),
    .groups = "drop"
  ) |>
  left_join(group_sizes, by = "group") |>
  mutate(proportion_positive = n_positive / n_pairwise)

sign_audit_pooled <- tibble(
  analysis = "Covariance sign audit",
  n_groups = nrow(sign_audit),
  n_pairwise_total = sum(sign_audit$n_pairwise),
  n_negative = sum(sign_audit$n_negative),
  n_positive = sum(sign_audit$n_positive),
  proportion_positive = sum(sign_audit$n_positive) / sum(sign_audit$n_pairwise),
  minimum_correlation = min(sign_audit$minimum_correlation),
  maximum_correlation = max(sign_audit$maximum_correlation),
  mean_correlation_min = min(sign_audit$mean_correlation),
  mean_correlation_max = max(sign_audit$mean_correlation),
  n_individuals_total = sum(group_sizes$n_individuals),
  source_table = basename(upstream_correlations)
)

cat("\nCOVARIANCE SIGN AUDIT BY GROUP\n")
print(sign_audit)
cat("\nPOOLED SIGN AUDIT\n")
print(sign_audit_pooled)


# ============================================================
# 7. SAVE RESULTS TABLES
# ============================================================

outputs <- list(
  instar_matched_decomposition = trait_decomposition,
  instar_matched_joint_legs = joint_table,
  precision_planned_leg_contrasts = precision_leg_table,
  precision_weapon_leg_slopes = weapon_leg_slopes,
  precision_weapon_ear_slopes = ear_weapon_slopes,
  precision_male_weapon_ear_mde = male_min_detectable,
  covariance_sign_audit = sign_audit,
  covariance_sign_audit_pooled = sign_audit_pooled
)

iwalk(outputs, ~ readr::write_csv(.x, file.path(tables_dir, paste0(.y, ".csv"))))

cat("\nInstar-matched decomposition and precision analyses completed.\n")
cat("Tables written to:\n ", normalizePath(tables_dir), "\n")
