#!/usr/bin/env Rscript

# Reanalysis requested during manuscript review.
#
# This script adds:
#   1. group-specific size-adjusted PCA with permutation parallel analysis;
#   2. direct permutation comparisons of group covariance/correlation matrices;
#   3. one-stage raw-trait refits replacing residual predictors;
#   4. head-adjusted eye models;
#   5. equivalence tests for ninth-instar intermediacy and weapon-ear slopes;
#   6. morph-frequency and classification-sensitivity checks; and
#   7. an across-analysis multiplicity sensitivity table.

suppressPackageStartupMessages({
  library(tidyverse)
  library(emmeans)
  library(nlme)
  library(lme4)
  library(lmerTest)
})

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
data_file <- if (length(args) >= 1) args[[1]] else file.path(project_root, "data", "trait_data.csv")
output_root <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(project_root, "analysis_outputs", "weta_trait_analysis")
}

# Usage: Rscript scripts/reviewer_reanalysis.R [data_file] [output_dir]
#   [permutations] [rarefactions]
if (length(args) > 4L) stop("Expected at most four arguments; see usage.")
parallel_reps <- parse_replicates(
  if (length(args) >= 3L) args[[3L]] else 9999L, "Permutations", 999L
)
matrix_permutations <- parallel_reps
rarefaction_reps <- parse_replicates(
  if (length(args) >= 4L) args[[4L]] else 2000L, "Rarefactions"
)
if (!file.exists(data_file)) stop("Trait data not found: ", data_file)

tables_dir <- file.path(output_root, "tables")
figures_dir <- file.path(output_root, "figures")
models_dir <- file.path(output_root, "models")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260921)
options(contrasts = c("contr.treatment", "contr.poly"))

reference_pronotum <- 7.2
random_skewers_reps <- 10000L

group_levels <- c("Female", "Eighth instar", "Ninth instar", "Tenth instar")
male_levels <- c("Eighth instar", "Ninth instar", "Tenth instar")

group_colours <- c(
  "Female" = "#4D4D4D",
  "Eighth instar" = "#2C7FB8",
  "Ninth instar" = "#7FCDBB",
  "Tenth instar" = "#D95F0E"
)

theme_manuscript <- function() {
  theme_classic(base_size = 11) +
    theme(
      text = element_text(family = "Times New Roman", colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

nested_f_test <- function(reduced, full, analysis) {
  tab <- anova(reduced, full)
  tibble(
    analysis = analysis,
    numerator_df = tab$Df[[2]],
    residual_df = tab$Res.Df[[2]],
    F = tab$F[[2]],
    p_value = tab$`Pr(>F)`[[2]]
  )
}

tidy_trends <- function(model, variable) {
  emmeans::emtrends(model, ~ group, var = variable) |>
    summary(infer = c(TRUE, TRUE)) |>
    as.data.frame() |>
    as_tibble()
}

tost_summary <- function(estimate, se, df, lower_bound, upper_bound) {
  t_lower <- (estimate - lower_bound) / se
  t_upper <- (estimate - upper_bound) / se
  p_lower <- pt(t_lower, df = df, lower.tail = FALSE)
  p_upper <- pt(t_upper, df = df, lower.tail = TRUE)
  crit <- qt(0.95, df = df)

  tibble(
    estimate = estimate,
    SE = se,
    df = df,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    lower_90 = estimate - crit * se,
    upper_90 = estimate + crit * se,
    p_lower = p_lower,
    p_upper = p_upper,
    TOST_p = pmax(p_lower, p_upper),
    equivalent = p_lower < 0.05 & p_upper < 0.05
  )
}

parallel_pca <- function(data_matrix, groups, reps = parallel_reps) {
  z <- scale(data_matrix)
  pca <- prcomp(z, center = FALSE, scale. = FALSE)

  # Principal-component signs are arbitrary. Orient every component
  # deterministically so saved loadings, scores, centroids and figures cannot
  # reverse sign when the analysis is rerun with a different BLAS/LAPACK build.
  # The sum of loadings sets the sign unless it is numerically zero, in which
  # case the loading with the largest absolute value is made positive.
  orientation <- apply(
    pca$rotation,
    2,
    function(loadings) {
      anchor <- sum(loadings)
      if (abs(anchor) <= sqrt(.Machine$double.eps)) {
        anchor <- loadings[[which.max(abs(loadings))]]
      }
      if (anchor < 0) -1 else 1
    }
  )
  pca$rotation <- sweep(pca$rotation, 2, orientation, "*")
  pca$x <- sweep(pca$x, 2, orientation, "*")
  observed <- pca$sdev^2

  null_eigenvalues <- replicate(
    reps,
    {
      permuted <- apply(z, 2, sample, replace = FALSE)
      eigen(cor(permuted), symmetric = TRUE, only.values = TRUE)$values
    }
  )

  null_95 <- apply(null_eigenvalues, 1, quantile, probs = 0.95, na.rm = TRUE)
  retained <- observed > null_95

  scores <- as_tibble(pca$x)
  scores$group <- groups

  group_tests <- map_dfr(
    seq_len(ncol(pca$x)),
    function(component) {
      response <- pca$x[, component]
      fit <- lm(response ~ groups)
      a <- anova(fit)
      tibble(
        component = paste0("PC", component),
        df = a$Df[[1]],
        residual_df = a$Df[[2]],
        F = a$`F value`[[1]],
        p_value = a$`Pr(>F)`[[1]]
      )
    }
  ) |>
    mutate(
      p_BH = p.adjust(p_value, method = "BH"),
      p_Holm = p.adjust(p_value, method = "holm")
    )

  loadings <- as_tibble(pca$rotation, rownames = "trait") |>
    pivot_longer(-trait, names_to = "component", values_to = "loading")

  variance <- tibble(
    component = paste0("PC", seq_along(observed)),
    observed_eigenvalue = observed,
    permutations = reps,
    parallel_95_eigenvalue = null_95,
    variance_percent = 100 * observed / sum(observed),
    cumulative_percent = cumsum(100 * observed / sum(observed)),
    retained = retained
  )

  list(
    pca = pca,
    scores = scores,
    group_tests = group_tests,
    loadings = loadings,
    variance = variance
  )
}

adjust_traits_to_reference <- function(data, trait_map) {
  adjusted <- map_dfc(
    names(trait_map),
    function(variable) {
      fit <- lm(
        stats::as.formula(paste(variable, "~ logP_c * group")),
        data = data
      )
      reference_data <- data
      reference_data$logP_c <- 0

      tibble(
        !!variable := residuals(fit) +
          predict(fit, newdata = reference_data)
      )
    }
  )
  names(adjusted) <- unname(trait_map)
  adjusted
}

group_specific_residual_matrix <- function(data, trait_map) {
  residuals_by_trait <- map_dfc(
    names(trait_map),
    function(variable) {
      tibble(
        !!variable := residuals(
          lm(
            stats::as.formula(paste(variable, "~ logP_c * group")),
            data = data
          )
        )
      )
    }
  )
  names(residuals_by_trait) <- unname(trait_map)
  scale(as.matrix(residuals_by_trait))
}

pooled_within_covariance <- function(mat, labels) {
  groups <- split(seq_len(nrow(mat)), labels)
  numerator <- Reduce(
    `+`,
    lapply(groups, function(index) (length(index) - 1) * cov(mat[index, , drop = FALSE]))
  )
  numerator / (nrow(mat) - length(groups))
}

matrix_distance_stat <- function(mat, labels, matrix_type = c("covariance", "correlation")) {
  matrix_type <- match.arg(matrix_type)
  groups <- split(seq_len(nrow(mat)), labels)

  if (matrix_type == "covariance") {
    pooled <- pooled_within_covariance(mat, labels)
    matrices <- lapply(groups, function(index) cov(mat[index, , drop = FALSE]))
  } else {
    pooled <- cor(mat)
    matrices <- lapply(groups, function(index) cor(mat[index, , drop = FALSE]))
  }

  sum(vapply(
    names(groups),
    function(group_name) {
      weight <- if (matrix_type == "covariance") length(groups[[group_name]]) - 1 else 1
      weight * sum((matrices[[group_name]] - pooled)^2)
    },
    numeric(1)
  ))
}

omnibus_matrix_permutation <- function(mat, labels, matrix_type, reps = matrix_permutations) {
  observed <- matrix_distance_stat(mat, labels, matrix_type)
  null <- replicate(
    reps,
    matrix_distance_stat(mat, sample(labels, replace = FALSE), matrix_type)
  )
  tibble(
    matrix_type = matrix_type,
    statistic = observed,
    permutations = reps,
    p_value = (1 + sum(null >= observed)) / (reps + 1)
  )
}

pairwise_matrix_permutation <- function(mat, labels, matrix_type, reps = matrix_permutations) {
  pairs <- combn(levels(droplevels(labels)), 2, simplify = FALSE)

  map_dfr(
    pairs,
    function(pair) {
      keep <- labels %in% pair
      pair_mat <- mat[keep, , drop = FALSE]
      pair_labels <- droplevels(labels[keep])

      matrix_for <- function(level_name) {
        x <- pair_mat[pair_labels == level_name, , drop = FALSE]
        if (matrix_type == "covariance") cov(x) else cor(x)
      }

      observed <- sum((matrix_for(pair[[1]]) - matrix_for(pair[[2]]))^2)
      null <- replicate(
        reps,
        {
          shuffled <- sample(pair_labels, replace = FALSE)
          a <- pair_mat[shuffled == pair[[1]], , drop = FALSE]
          b <- pair_mat[shuffled == pair[[2]], , drop = FALSE]
          ma <- if (matrix_type == "covariance") cov(a) else cor(a)
          mb <- if (matrix_type == "covariance") cov(b) else cor(b)
          sum((ma - mb)^2)
        }
      )

      tibble(
        matrix_type = matrix_type,
        group_1 = pair[[1]],
        group_2 = pair[[2]],
        frobenius_distance = sqrt(observed),
        permutations = reps,
        p_value = (1 + sum(null >= observed)) / (reps + 1)
      )
    }
  ) |>
    group_by(matrix_type) |>
    mutate(
      p_BH = p.adjust(p_value, method = "BH"),
      p_Holm = p.adjust(p_value, method = "holm")
    ) |>
    ungroup()
}

random_skewers_similarity <- function(a, b, reps = random_skewers_reps) {
  betas <- matrix(rnorm(nrow(a) * reps), nrow = nrow(a), ncol = reps)
  betas <- sweep(betas, 2, sqrt(colSums(betas^2)), "/")
  response_a <- a %*% betas
  response_b <- b %*% betas
  cosine <- colSums(response_a * response_b) /
    sqrt(colSums(response_a^2) * colSums(response_b^2))
  mean(cosine, na.rm = TRUE)
}

matrix_similarity <- function(mat, labels, matrix_type) {
  group_names <- levels(droplevels(labels))
  matrices <- setNames(
    lapply(group_names, function(g) {
      x <- mat[labels == g, , drop = FALSE]
      if (matrix_type == "covariance") cov(x) else cor(x)
    }),
    group_names
  )
  pairs <- combn(group_names, 2, simplify = FALSE)
  upper <- upper.tri(matrices[[1]], diag = TRUE)

  map_dfr(
    pairs,
    function(pair) {
      a <- matrices[[pair[[1]]]]
      b <- matrices[[pair[[2]]]]
      tibble(
        matrix_type = matrix_type,
        group_1 = pair[[1]],
        group_2 = pair[[2]],
        element_correlation = cor(a[upper], b[upper]),
        random_skewers_similarity = random_skewers_similarity(a, b)
      )
    }
  )
}

rarefied_matrix_distances <- function(
    mat,
    labels,
    matrix_type,
    target_n = min(table(labels)),
    reps = rarefaction_reps
) {
  group_names <- levels(droplevels(labels))
  pairs <- combn(group_names, 2, simplify = FALSE)

  distances <- map_dfr(
    seq_len(reps),
    function(iteration) {
      sampled <- unlist(lapply(
        group_names,
        function(g) sample(which(labels == g), target_n, replace = FALSE)
      ))
      sampled_labels <- droplevels(labels[sampled])
      sampled_mat <- mat[sampled, , drop = FALSE]
      matrices <- setNames(
        lapply(group_names, function(g) {
          x <- sampled_mat[sampled_labels == g, , drop = FALSE]
          if (matrix_type == "covariance") cov(x) else cor(x)
        }),
        group_names
      )

      map_dfr(
        pairs,
        function(pair) {
          tibble(
            iteration = iteration,
            group_1 = pair[[1]],
            group_2 = pair[[2]],
            distance = sqrt(sum((matrices[[pair[[1]]]] - matrices[[pair[[2]]]])^2))
          )
        }
      )
    }
  )

  distances |>
    group_by(group_1, group_2) |>
    summarise(
      matrix_type = matrix_type,
      target_n = target_n,
      replicates = reps,
      median_distance = median(distance),
      lower_95 = quantile(distance, 0.025),
      upper_95 = quantile(distance, 0.975),
      .groups = "drop"
    )
}

zero_interval <- function(x) c(-0.01, 0.01)

extract_gls_joint_tests <- function(model, cov_reduce) {
  emmeans::joint_tests(
    model,
    cov.reduce = cov_reduce,
    mode = "appx-satterthwaite"
  ) |>
    as.data.frame() |>
    as_tibble() |>
    transmute(
      term = `model term`,
      `Sum Sq` = NA_real_,
      `Mean Sq` = NA_real_,
      NumDF = df1,
      DenDF = df2,
      `F value` = F.ratio,
      `Pr(>F)` = p.value
    )
}

dat_raw <- readr::read_csv(data_file, na = c("", "NA"), show_col_types = FALSE)

required <- c(
  "ID", "sex", "morph", "pronotum", "head_length", "head_width",
  "forefemur", "foretibia", "midfemur", "midtibia", "hindfemur",
  "hindtibia", "ear", "eye"
)
stopifnot(all(required %in% names(dat_raw)))

dat <- dat_raw |>
  mutate(
    group = factor(
      morph,
      levels = c("female", "eighth", "ninth", "tenth"),
      labels = group_levels
    ),
    foreleg = forefemur + foretibia,
    midleg = midfemur + midtibia,
    hindleg = hindfemur + hindtibia,
    head_size = sqrt(head_length * head_width),
    ear_linear = sqrt(ear),
    logP_c = log(pronotum) - log(reference_pronotum),
    log_head_size = log(head_size),
    log_foreleg = log(foreleg),
    log_midleg = log(midleg),
    log_hindleg = log(hindleg),
    log_foretibia = log(foretibia),
    log_ear_linear = log(ear_linear),
    log_eye = log(eye)
  )

# -----------------------------------------------------------------------------
# Morph-frequency and classification sensitivity
# -----------------------------------------------------------------------------

observed_males <- dat |>
  filter(group %in% male_levels) |>
  count(group) |>
  complete(group = factor(male_levels, levels = male_levels), fill = list(n = 0)) |>
  arrange(group)

frequency_sources <- tribble(
  ~reference, ~eighth, ~ninth, ~tenth,
  "Kelly 2026 retained assignments", 399, 263, 242,
  "Kelly and Adams 2010 field sample", 370, 264, 220
)

morph_frequency_tests <- pmap_dfr(
  frequency_sources,
  function(reference, eighth, ninth, tenth) {
    source_counts <- c(eighth, ninth, tenth)
    expected_proportions <- source_counts / sum(source_counts)
    test <- chisq.test(observed_males$n, p = expected_proportions)
    tibble(
      reference = reference,
      observed_eighth = observed_males$n[[1]],
      observed_ninth = observed_males$n[[2]],
      observed_tenth = observed_males$n[[3]],
      expected_prop_eighth = expected_proportions[[1]],
      expected_prop_ninth = expected_proportions[[2]],
      expected_prop_tenth = expected_proportions[[3]],
      chi_square = unname(test$statistic),
      df = unname(test$parameter),
      p_value = test$p.value
    )
  }
)

male_classification <- dat |>
  filter(sex == "m") |>
  mutate(
    current_label = as.character(group),
    legacy_2010_label = case_when(
      head_length < 19.04 ~ "Eighth instar",
      head_length < 24.20 ~ "Ninth instar",
      TRUE ~ "Tenth instar"
    ),
    differs_from_2010 = current_label != legacy_2010_label,
    distance_to_inferred_cutpoint = pmin(
      abs(head_length - 18.50),
      abs(head_length - 24.15)
    ),
    within_0_25_mm = distance_to_inferred_cutpoint <= 0.25,
    within_0_50_mm = distance_to_inferred_cutpoint <= 0.50
  )

classification_summary <- male_classification |>
  summarise(
    male_n = n(),
    current_eighth_max = max(head_length[current_label == "Eighth instar"]),
    current_ninth_min = min(head_length[current_label == "Ninth instar"]),
    current_ninth_max = max(head_length[current_label == "Ninth instar"]),
    current_tenth_min = min(head_length[current_label == "Tenth instar"]),
    n_disagree_legacy_2010 = sum(differs_from_2010),
    percent_disagree_legacy_2010 = 100 * mean(differs_from_2010),
    n_within_0_25_mm = sum(within_0_25_mm),
    n_within_0_50_mm = sum(within_0_50_mm),
    posterior_error_available = FALSE
  )

# -----------------------------------------------------------------------------
# Group-specific size-adjusted PCA and parallel analysis
# -----------------------------------------------------------------------------

# The five traits not used to define male morphs are the primary ordination.
# Head size is added only in a six-trait sensitivity analysis. For each trait,
# the adjusted value is the group-specific fitted value at logP_c = 0 (a
# 7.2-mm pronotum) plus the model residual. This retains adjusted group
# differences while removing each group's estimated body-size slope.
pca_traits <- c(
  log_foreleg = "Foreleg",
  log_midleg = "Midleg",
  log_hindleg = "Hindleg",
  log_ear_linear = "Linearized ear",
  log_eye = "Eye"
)
pca_six_traits <- c(
  log_head_size = "Head size",
  pca_traits
)

pca_data <- dat |>
  filter(if_all(all_of(c(names(pca_traits), "logP_c", "group")), ~ !is.na(.x)))
size_adjusted <- adjust_traits_to_reference(pca_data, pca_traits)
set.seed(2101)
pca_result <- parallel_pca(as.matrix(size_adjusted), pca_data$group)

pca_scores <- bind_cols(
  pca_data |> select(ID, group),
  as_tibble(pca_result$pca$x)
)

pca_centroids <- pca_scores |>
  group_by(group) |>
  summarise(
    n = n(),
    across(starts_with("PC"), list(mean = mean, se = ~ sd(.x) / sqrt(n()))),
    .groups = "drop"
  )

pca_loadings <- pca_result$loadings |>
  mutate(
    retained = component %in% pca_result$variance$component[pca_result$variance$retained]
  )

# Compatibility alias: the primary analysis is now already head-excluded.
pca_no_head <- pca_result
pca_no_head_summary <- pca_result$variance |>
  left_join(pca_result$group_tests, by = "component")

# Six-trait sensitivity, including head size.
pca_six_data <- dat |>
  filter(if_all(all_of(c(names(pca_six_traits), "logP_c", "group")), ~ !is.na(.x)))
pca_six_adjusted <- adjust_traits_to_reference(pca_six_data, pca_six_traits)
set.seed(2102)
pca_six_result <- parallel_pca(
  as.matrix(pca_six_adjusted),
  pca_six_data$group
)
pca_six_scores <- bind_cols(
  pca_six_data |> select(ID, group),
  as_tibble(pca_six_result$pca$x)
)
pca_six_centroids <- pca_six_scores |>
  group_by(group) |>
  summarise(
    n = n(),
    across(starts_with("PC"), list(mean = mean, se = ~ sd(.x) / sqrt(n()))),
    .groups = "drop"
  )
pca_six_loadings <- pca_six_result$loadings |>
  mutate(
    retained = component %in%
      pca_six_result$variance$component[pca_six_result$variance$retained]
  )

# Classification sensitivity for the primary five-trait PCA: exclude males
# whose archived label differs from the category implied by the published 2010
# cutoffs, then repeat the group-specific adjustment and parallel analysis.
classification_disagreement_ids <- male_classification$ID[male_classification$differs_from_2010]
pca_sensitivity_data <- dat |>
  filter(!ID %in% classification_disagreement_ids) |>
  filter(if_all(all_of(c(names(pca_traits), "logP_c", "group")), ~ !is.na(.x)))
pca_sensitivity_matrix <- adjust_traits_to_reference(
  pca_sensitivity_data,
  pca_traits
)
set.seed(2103)
pca_classification_sensitivity <- parallel_pca(
  as.matrix(pca_sensitivity_matrix),
  pca_sensitivity_data$group
)

# -----------------------------------------------------------------------------
# Direct group covariance/correlation matrix comparisons
# -----------------------------------------------------------------------------

# The primary matrix comparison uses the same five non-head traits as the
# primary PCA. Pairwise tests, matrix similarities and rarefaction are retained
# as supporting summaries; the omnibus test is the inferential focus.
within_group_matrix <- group_specific_residual_matrix(pca_data, pca_traits)
matrix_groups <- droplevels(pca_data$group)

set.seed(2201)
covariance_omnibus <- bind_rows(
  omnibus_matrix_permutation(within_group_matrix, matrix_groups, "covariance"),
  omnibus_matrix_permutation(within_group_matrix, matrix_groups, "correlation")
) |>
  mutate(
    trait_set = "five_non_head_traits",
    analysis_role = "primary_omnibus"
)

set.seed(2202)
covariance_pairwise <- bind_rows(
  pairwise_matrix_permutation(within_group_matrix, matrix_groups, "covariance"),
  pairwise_matrix_permutation(within_group_matrix, matrix_groups, "correlation")
) |>
  mutate(
    trait_set = "five_non_head_traits",
    analysis_role = "supporting_pairwise"
)

set.seed(2203)
covariance_similarity <- bind_rows(
  matrix_similarity(within_group_matrix, matrix_groups, "covariance"),
  matrix_similarity(within_group_matrix, matrix_groups, "correlation")
) |>
  mutate(
    trait_set = "five_non_head_traits",
    analysis_role = "supporting_similarity"
)

set.seed(2204)
covariance_rarefaction <- bind_rows(
  rarefied_matrix_distances(within_group_matrix, matrix_groups, "covariance"),
  rarefied_matrix_distances(within_group_matrix, matrix_groups, "correlation")
) |>
  mutate(
    trait_set = "five_non_head_traits",
    analysis_role = "supporting_rarefaction"
  )

group_correlation_long <- map_dfr(
  levels(matrix_groups),
  function(group_name) {
    m <- cor(within_group_matrix[matrix_groups == group_name, , drop = FALSE])
    as.data.frame(as.table(m), responseName = "correlation") |>
      as_tibble() |>
      mutate(group = group_name)
  }
) |>
  rename(trait_1 = Var1, trait_2 = Var2) |>
  mutate(trait_set = "five_non_head_traits")

# Compatibility alias for the earlier filename. The head-excluded result is
# now the primary five-trait analysis rather than a secondary sensitivity.
covariance_no_head_omnibus <- covariance_omnibus |>
  mutate(sensitivity = "Primary analysis: head trait excluded")

# Six-trait sensitivity: add head size, the variable that contributes to male
# morph classification, and repeat only the prespecified omnibus comparisons.
within_group_matrix_six <- group_specific_residual_matrix(
  pca_six_data,
  pca_six_traits
)
matrix_groups_six <- droplevels(pca_six_data$group)
set.seed(2205)
covariance_six_trait_omnibus <- bind_rows(
  omnibus_matrix_permutation(
    within_group_matrix_six,
    matrix_groups_six,
    "covariance"
  ),
  omnibus_matrix_permutation(
    within_group_matrix_six,
    matrix_groups_six,
    "correlation"
  )
) |>
  mutate(
    trait_set = "six_traits_including_head",
    analysis_role = "sensitivity_omnibus"
  )

# -----------------------------------------------------------------------------
# One-stage raw-trait refits
# -----------------------------------------------------------------------------

# Keep approximate model degrees of freedom independent of resampling counts.
set.seed(2301)

legs_long <- dat |>
  select(ID, group, logP_c, log_head_size, log_foreleg, log_midleg, log_hindleg) |>
  pivot_longer(
    c(log_foreleg, log_midleg, log_hindleg),
    names_to = "leg",
    values_to = "log_leg_length"
  ) |>
  mutate(
    ID = factor(ID),
    leg = factor(
      leg,
      levels = c("log_foreleg", "log_midleg", "log_hindleg"),
      labels = c("Foreleg", "Midleg", "Hindleg")
    ),
    leg_order = as.integer(leg)
  ) |>
  arrange(ID, leg_order)

if (anyDuplicated(legs_long[c("ID", "leg")])) {
  stop("Each individual must contribute at most one observation per leg pair.")
}

contrasts(legs_long$group) <- contr.sum(nlevels(legs_long$group))
contrasts(legs_long$leg) <- contr.sum(nlevels(legs_long$leg))

mod_head_leg_raw <- nlme::gls(
  log_leg_length ~
    logP_c * group * leg +
    log_head_size * group * leg,
  data = legs_long,
  correlation = nlme::corSymm(form = ~ leg_order | ID),
  weights = nlme::varIdent(form = ~ 1 | leg),
  method = "REML",
  na.action = na.omit,
  control = nlme::glsControl(opt = "optim", msMaxIter = 1000)
)

head_leg_raw_anova <- extract_gls_joint_tests(
  mod_head_leg_raw,
  cov_reduce = list(
    logP_c = zero_interval,
    log_head_size = emmeans::make.meanint(1)
  )
)
head_leg_raw_trends <- emtrends(
  mod_head_leg_raw,
  ~ group * leg,
  var = "log_head_size",
  mode = "appx-satterthwaite"
)
head_leg_raw_slopes <- head_leg_raw_trends |>
  summary(infer = c(TRUE, TRUE)) |>
  as.data.frame() |>
  as_tibble()

# Planned biological contrasts in the one-stage raw-head model. These replace
# the legacy contrasts obtained from a model that used a residual as a
# predictor. Holm adjustment is applied to the three contrasts within each leg.
head_leg_strategy_methods <- list(
  "Eighth instar - Tenth instar" = c(0, 1, 0, -1),
  "Ninth instar - midpoint(Eighth instar,Tenth instar)" =
    c(0, -0.5, 1, -0.5),
  "Tenth instar - Female" = c(-1, 0, 0, 1)
)
head_leg_raw_strategy_slopes <- contrast(
  head_leg_raw_trends,
  method = head_leg_strategy_methods,
  by = "leg",
  adjust = "holm"
) |>
  summary(infer = c(TRUE, TRUE), adjust = "holm") |>
  as.data.frame() |>
  as_tibble() |>
  mutate(Holm_family = paste0("Three planned strategy slope contrasts within ", leg))

mod_ear_tibia_raw_common <- lm(
  log_ear_linear ~ logP_c * group + log_foretibia,
  data = dat,
  na.action = na.omit
)
mod_ear_tibia_raw_full <- lm(
  log_ear_linear ~ logP_c * group + log_foretibia * group,
  data = dat,
  na.action = na.omit
)
ear_tibia_raw_test <- nested_f_test(
  mod_ear_tibia_raw_common,
  mod_ear_tibia_raw_full,
  "Ear by raw foretibia: group-specific slopes"
)
ear_tibia_raw_slopes <- tidy_trends(mod_ear_tibia_raw_full, "log_foretibia")

mod_eye_head_raw_common <- lm(
  log_eye ~ logP_c * group + log_head_size,
  data = dat,
  na.action = na.omit
)
mod_eye_head_raw_full <- lm(
  log_eye ~ logP_c * group + log_head_size * group,
  data = dat,
  na.action = na.omit
)
eye_head_raw_test <- nested_f_test(
  mod_eye_head_raw_common,
  mod_eye_head_raw_full,
  "Eye by raw head size: group-specific slopes"
)
eye_head_raw_slopes <- tidy_trends(mod_eye_head_raw_full, "log_head_size")

mod_ear_head_raw_common <- lm(
  log_ear_linear ~ logP_c * group + log_head_size,
  data = dat,
  na.action = na.omit
)
mod_ear_head_raw_full <- lm(
  log_ear_linear ~ logP_c * group + log_head_size * group,
  data = dat,
  na.action = na.omit
)
ear_head_raw_test <- nested_f_test(
  mod_ear_head_raw_common,
  mod_ear_head_raw_full,
  "Ear by raw head size: group-specific slopes"
)
ear_head_raw_slopes <- tidy_trends(mod_ear_head_raw_full, "log_head_size")

# Does group explain eye length after both pronotum and head size are included?
mod_eye_head_no_group <- lm(
  log_eye ~ logP_c + log_head_size,
  data = dat,
  na.action = na.omit
)
mod_eye_head_with_group <- lm(
  log_eye ~ logP_c + log_head_size + group,
  data = dat,
  na.action = na.omit
)
eye_group_after_head_test <- nested_f_test(
  mod_eye_head_no_group,
  mod_eye_head_with_group,
  "Group effect on eye length after pronotum and head size"
)

mod_eye_group_body <- lm(
  log_eye ~ logP_c + group,
  data = dat,
  na.action = na.omit
)
mod_eye_group_body_head <- lm(
  log_eye ~ logP_c + group + log_head_size,
  data = dat,
  na.action = na.omit
)
eye_head_increment_test <- nested_f_test(
  mod_eye_group_body,
  mod_eye_group_body_head,
  "Head-size contribution to eye length after pronotum and group"
)

head_support <- dat |>
  group_by(group) |>
  summarise(
    n = sum(!is.na(head_size) & !is.na(log_eye)),
    minimum_head_size = min(head_size, na.rm = TRUE),
    maximum_head_size = max(head_size, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    common_lower = max(minimum_head_size),
    common_upper = min(maximum_head_size),
    common_overlap = common_lower <= common_upper
  )

# Refit the head-adjusted eye model after excluding the 10 classifications that
# disagree with the published 2010 thresholds.
eye_classification_sensitivity <- bind_rows(
  eye_group_after_head_test |>
    mutate(dataset = "Current archived labels"),
  {
    sensitivity_data <- dat |> filter(!ID %in% classification_disagreement_ids)
    reduced <- lm(log_eye ~ logP_c + log_head_size, data = sensitivity_data)
    full <- lm(log_eye ~ logP_c + log_head_size + group, data = sensitivity_data)
    nested_f_test(reduced, full, "Group effect on eye length after pronotum and head size") |>
      mutate(dataset = "Excluded labels disagreeing with 2010 cutoffs")
  }
)

# -----------------------------------------------------------------------------
# Equivalence and sensitivity analyses
# -----------------------------------------------------------------------------

leg_strategy_path <- file.path(tables_dir, "joint_leg_strategy_contrasts.csv")
if (!file.exists(leg_strategy_path)) {
  stop("The existing joint-leg contrast table is required: ", leg_strategy_path)
}

leg_strategy <- readr::read_csv(leg_strategy_path, show_col_types = FALSE)
ninth_midpoint <- leg_strategy |>
  filter(str_detect(contrast, regex("Ninth instar - midpoint", ignore_case = TRUE)))

leg_intermediacy_equivalence <- map_dfr(
  c(0.025, 0.05),
  function(margin) {
    bounds <- c(log(1 - margin), log(1 + margin))
    pmap_dfr(
      ninth_midpoint,
      function(...) {
        row <- list(...)
        tost_summary(
          estimate = as.numeric(row$estimate),
          se = as.numeric(row$SE),
          df = as.numeric(row$df),
          lower_bound = bounds[[1]],
          upper_bound = bounds[[2]]
        ) |>
          mutate(
            leg = as.character(row$leg),
            equivalence_margin_percent = 100 * margin,
            estimate_ratio = exp(estimate),
            lower_90_ratio = exp(lower_90),
            upper_90_ratio = exp(upper_90)
          )
      }
    )
  }
) |>
  select(
    leg, equivalence_margin_percent, estimate_ratio,
    lower_90_ratio, upper_90_ratio, TOST_p, equivalent,
    everything()
  )

# Apply the same equivalence framework to the sensory intermediacy contrasts.
# Eye length is modelled directly on the log scale. Ear area is modelled as
# log(sqrt(area)), so area-scale equivalence bounds and ratios use a factor of 2.
sensory_intermediacy_inputs <- tribble(
  ~trait, ~filename, ~log_scale_multiplier, ~reported_scale,
  "Tympanal area", "ear_body_strategy_contrasts.csv", 2, "area ratio",
  "Eye length", "eye_body_strategy_contrasts.csv", 1, "length ratio"
) |>
  mutate(
    contrast_row = map(filename, function(filename) {
      readr::read_csv(file.path(tables_dir, filename), show_col_types = FALSE) |>
        filter(str_detect(contrast, regex("Ninth instar - midpoint", ignore_case = TRUE)))
    })
  )

sensory_intermediacy_equivalence <- pmap_dfr(
  sensory_intermediacy_inputs,
  function(trait, filename, log_scale_multiplier, reported_scale, contrast_row) {
    stopifnot(nrow(contrast_row) == 1)

    map_dfr(
      c(0.025, 0.05),
      function(margin) {
        bounds <- c(
          log(1 - margin) / log_scale_multiplier,
          log(1 + margin) / log_scale_multiplier
        )
        result <- tost_summary(
          estimate = contrast_row$estimate[[1]],
          se = contrast_row$SE[[1]],
          df = contrast_row$df[[1]],
          lower_bound = bounds[[1]],
          upper_bound = bounds[[2]]
        )

        result |>
          mutate(
            trait = trait,
            reported_scale = reported_scale,
            equivalence_margin_percent = 100 * margin,
            estimate_ratio = exp(log_scale_multiplier * estimate),
            lower_90_ratio = exp(log_scale_multiplier * lower_90),
            upper_90_ratio = exp(log_scale_multiplier * upper_90)
          )
      }
    )
  }
) |>
  select(
    trait, reported_scale, equivalence_margin_percent, estimate_ratio,
    lower_90_ratio, upper_90_ratio, TOST_p, equivalent, everything()
  )

# Standardized partial slopes allow a common small-effect equivalence region.
weapon_ear_equivalence <- map_dfr(
  group_levels,
  function(group_name) {
    x <- dat |> filter(group == group_name) |>
      drop_na(log_ear_linear, logP_c, log_head_size)
    standardized <- x |>
      transmute(
        z_ear = as.numeric(scale(log_ear_linear)),
        z_body = as.numeric(scale(logP_c)),
        z_head = as.numeric(scale(log_head_size))
      )
    fit <- lm(z_ear ~ z_body + z_head, data = standardized)
    co <- summary(fit)$coefficients["z_head", ]
    eq <- tost_summary(
      estimate = co[["Estimate"]],
      se = co[["Std. Error"]],
      df = df.residual(fit),
      lower_bound = -0.20,
      upper_bound = 0.20
    )
    robust_fit <- MASS::rlm(z_ear ~ z_body + z_head, data = standardized, maxit = 200)

    eq |>
      mutate(
        group = group_name,
        n = nrow(standardized),
        standardized_beta = estimate,
        robust_beta = coef(robust_fit)[["z_head"]],
        equivalence_bound = 0.20
      )
  }
) |>
  select(
    group, n, standardized_beta, lower_90, upper_90,
    robust_beta, equivalence_bound, TOST_p, equivalent, everything()
  )

# -----------------------------------------------------------------------------
# Multiplicity sensitivity
# -----------------------------------------------------------------------------

read_single_p <- function(filename) {
  readr::read_csv(file.path(tables_dir, filename), show_col_types = FALSE)$p_value[[1]]
}

joint_anova <- readr::read_csv(
  file.path(tables_dir, "joint_leg_model_type3_anova.csv"),
  show_col_types = FALSE
)

head_leg_term <- head_leg_raw_anova |>
  filter(term == "group:leg:log_head_size") |>
  pull(`Pr(>F)`)

focal_tests <- tribble(
  ~analysis, ~tier, ~p_value,
  "Leg allocation group x leg", "Primary", joint_anova$`Pr(>F)`[joint_anova$term == "group:leg"],
  "Leg allometry body size x group", "Secondary", joint_anova$`Pr(>F)`[joint_anova$term == "logP_c:group"],
  "Leg allometry body size x group x leg", "Secondary", joint_anova$`Pr(>F)`[joint_anova$term == "logP_c:group:leg"],
  "Raw head size x group x leg", "Primary", head_leg_term,
  "Ear body-size slope x group", "Secondary", read_single_p("ear_body_slope_model_test.csv"),
  "Eye body-size slope x group", "Secondary", read_single_p("eye_body_slope_model_test.csv"),
  "Raw ear-foretibia slope x group", "Secondary", ear_tibia_raw_test$p_value,
  "Raw eye-head slope x group", "Secondary", eye_head_raw_test$p_value,
  "Raw weapon-ear slope x group", "Secondary", ear_head_raw_test$p_value,
  "Five-trait correlation-matrix equality", "Primary", covariance_omnibus$p_value[covariance_omnibus$matrix_type == "correlation"]
) |>
  mutate(
    p_BH = p.adjust(p_value, method = "BH"),
    p_Holm = p.adjust(p_value, method = "holm")
  )

fragile_contrast <- readr::read_csv(
  file.path(tables_dir, "joint_three_way_slope_contrasts.csv"),
  show_col_types = FALSE
) |>
  slice_min(abs(p.value - 0.047646), n = 1, with_ties = FALSE) |>
  transmute(
    analysis = paste0("Exploratory contrast: ", group_pairwise, "; ", leg_pairwise),
    reported_family_adjusted_p = p.value,
    interpretation = "Not treated as confirmatory after cross-analysis multiplicity review"
  )

# -----------------------------------------------------------------------------
# Analysis hierarchy and combined adjusted-trait summary
# -----------------------------------------------------------------------------

analysis_hierarchy <- tribble(
  ~analysis, ~classification, ~recommended_location, ~rationale,
  "Female-referenced rank-one contrast-vector test", "central", "main text", "Direct test of whether adjusted morph divergence requires more than one direction",
  "Joint repeated-trait leg-allocation model", "central", "main text", "Primary test of locomotory allocation across all three leg pairs",
  "Adjusted sensory-trait group means", "central", "main text", "Tests the directional sensory prediction at a common body size",
  "One-stage raw-head weapon-leg model", "supporting", "main text and supplement", "Tests weapon-leg integration without residual predictors but is not needed to establish multivariate dimensionality",
  "Five-trait covariance/correlation omnibus tests", "central", "main text", "Direct test of morph-specific integration architecture excluding the classification trait",
  "Five-trait PCA with parallel analysis", "supporting", "main text or supplement", "Describes the dominant dimension of individual variation but is not the direct centroid-vector test",
  "Tenth-instar male-female anterior-posterior contrast", "central", "main text", "Holds terminal instar constant while contrasting weapon expression",
  "Ninth-instar midpoint equivalence tests", "central", "main text", "Provides positive evidence for the retained intermediacy claim",
  "Measurement repeatability", "supporting", "supplement", "Establishes measurement reliability without expanding the biological hypothesis set",
  "Matrix pairwise comparisons", "supporting", "supplement", "Localizes any omnibus matrix result",
  "Matrix similarities and rarefaction", "supporting", "supplement", "Describes orientation and unequal-sample-size sensitivity",
  "Six-trait PCA including head size", "sensitivity", "supplement", "Assesses dependence on the trait used in male morph classification",
  "Six-trait covariance/correlation omnibus tests", "sensitivity", "supplement", "Assesses whether including head size changes the matrix conclusion",
  "Head-adjusted eye model", "sensitivity", "main text or supplement", "Addresses structural coupling and classification circularity",
  "Classification-disagreement exclusion", "sensitivity", "supplement", "Checks dependence on archived morph labels near historical cutoffs",
  "Head-included and alternative-scaling rank-one fits", "sensitivity", "supplement", "Checks scale and trait-set dependence of vector geometry",
  "Morph-frequency comparison", "sensitivity", "supplement", "Quantifies possible ascertainment or collection bias",
  "Separate leg-pair and segment allometries", "exploratory", "supplement", "Decomposes the joint leg result but adds many weakly powered tests",
  "Local ear-foretibia and eye-head slopes", "exploratory", "supplement", "Describes local anatomical coupling without testing the central allocation contrast",
  "Weapon-ear slopes and equivalence", "exploratory", "supplement", "Relevant to decoupling but too imprecise for a strong positive claim",
  "Across-analysis multiplicity sensitivity", "exploratory", "supplement", "Post hoc robustness summary rather than a prespecified family"
)

adjusted_mean_paths <- c(
  legs = file.path(tables_dir, "joint_adjusted_leg_means.csv"),
  ear = file.path(tables_dir, "ear_body_adjusted_group_means.csv"),
  eye = file.path(tables_dir, "eye_body_adjusted_group_means.csv")
)
missing_adjusted_mean_paths <- adjusted_mean_paths[!file.exists(adjusted_mean_paths)]
if (length(missing_adjusted_mean_paths) > 0) {
  stop(
    "Combined adjusted-trait figure requires: ",
    paste(missing_adjusted_mean_paths, collapse = ", ")
  )
}

adjusted_leg_figure_data <- readr::read_csv(
  adjusted_mean_paths[["legs"]],
  show_col_types = FALSE
) |>
  transmute(
    group,
    trait = paste0(leg, " length"),
    response_scale = "mm",
    estimate = predicted,
    lower_95 = lower,
    upper_95 = upper
  )

adjusted_ear_figure_data <- readr::read_csv(
  adjusted_mean_paths[["ear"]],
  show_col_types = FALSE
) |>
  transmute(
    group,
    trait = "Tympanal area",
    response_scale = "mm^2",
    estimate = predicted_ear_area,
    lower_95 = ear_area_lower,
    upper_95 = ear_area_upper
  )

adjusted_eye_figure_data <- readr::read_csv(
  adjusted_mean_paths[["eye"]],
  show_col_types = FALSE
) |>
  transmute(
    group,
    trait = "Eye length",
    response_scale = "mm",
    estimate = predicted_eye_length,
    lower_95 = eye_lower,
    upper_95 = eye_upper
  )

adjusted_trait_means_combined <- bind_rows(
  adjusted_leg_figure_data,
  adjusted_ear_figure_data,
  adjusted_eye_figure_data
) |>
  mutate(
    group = factor(group, levels = group_levels),
    trait = factor(
      trait,
      levels = c(
        "Foreleg length", "Midleg length", "Hindleg length",
        "Tympanal area", "Eye length"
      )
    )
  )

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

p_scree <- pca_result$variance |>
  select(component, observed_eigenvalue, parallel_95_eigenvalue) |>
  pivot_longer(-component, names_to = "series", values_to = "eigenvalue") |>
  mutate(
    component = factor(component, levels = paste0("PC", seq_len(nrow(pca_result$variance)))),
    series = recode(
      series,
      observed_eigenvalue = "Observed",
      parallel_95_eigenvalue = "Parallel-analysis 95th percentile"
    )
  ) |>
  ggplot(aes(component, eigenvalue, group = series, colour = series, shape = series)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("Observed" = "black", "Parallel-analysis 95th percentile" = "#D95F0E")) +
  labs(x = "Principal component", y = "Eigenvalue", colour = NULL, shape = NULL) +
  theme_manuscript()

p_scores <- pca_scores |>
  ggplot(aes(PC1, PC2, colour = group)) +
  geom_point(alpha = 0.45, size = 1.6) +
  stat_ellipse(level = 0.68, linewidth = 0.8, show.legend = FALSE) +
  geom_point(
    data = pca_scores |> group_by(group) |> summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop"),
    size = 3.2,
    shape = 18
  ) +
  scale_colour_manual(values = group_colours) +
  labs(
    x = paste0("PC1 (", round(pca_result$variance$variance_percent[[1]], 1), "%)"),
    y = paste0("PC2 (", round(pca_result$variance$variance_percent[[2]], 1), "%)"),
    colour = "Group"
  ) +
  theme_manuscript()

p_correlations <- group_correlation_long |>
  mutate(
    group = factor(group, levels = group_levels),
    trait_1 = factor(trait_1, levels = colnames(within_group_matrix)),
    trait_2 = factor(trait_2, levels = rev(colnames(within_group_matrix)))
  ) |>
  ggplot(aes(trait_1, trait_2, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_wrap(~ group, ncol = 2) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  theme_manuscript() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_adjusted_traits <- adjusted_trait_means_combined |>
  ggplot(aes(group, estimate, colour = group)) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0.12,
    linewidth = 0.55
  ) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ trait,
    scales = "free_y",
    ncol = 3,
    labeller = as_labeller(c(
      "Foreleg length" = "Foreleg length (mm)",
      "Midleg length" = "Midleg length (mm)",
      "Hindleg length" = "Hindleg length (mm)",
      "Tympanal area" = "Tympanal area (mm²)",
      "Eye length" = "Eye length (mm)"
    ))
  ) +
  scale_colour_manual(values = group_colours, drop = FALSE) +
  labs(
    x = NULL,
    y = paste0("Adjusted mean at pronotum = ", reference_pronotum, " mm")
  ) +
  theme_manuscript() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(figures_dir, "pca_parallel_analysis.png"), p_scree, width = 6.6, height = 4.2, dpi = 300)
ggsave(file.path(figures_dir, "pca_group_scores.png"), p_scores, width = 6.6, height = 5.0, dpi = 300)
ggsave(file.path(figures_dir, "group_trait_correlation_matrices.png"), p_correlations, width = 7.2, height = 7.0, dpi = 300)
ggsave(
  file.path(figures_dir, "adjusted_trait_means_combined.png"),
  p_adjusted_traits,
  width = 8.2,
  height = 6.2,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Save tables and model objects
# -----------------------------------------------------------------------------

outputs <- list(
  reviewer_resampling_settings = tibble(
    procedure = c("pca_parallel_analysis", "covariance_matrix_permutation",
                  "covariance_rarefaction", "random_skewers"),
    replicates = c(parallel_reps, matrix_permutations,
                   rarefaction_reps, random_skewers_reps)
  ),
  morph_frequency_tests = morph_frequency_tests,
  morph_classification_summary = classification_summary,
  morph_classification_records = male_classification,
  pca_parallel_analysis = pca_result$variance |>
    mutate(
      trait_set = "five_non_head_traits",
      size_adjustment = "group_specific_at_7.2_mm"
    ),
  pca_loadings = pca_loadings |>
    mutate(trait_set = "five_non_head_traits"),
  pca_group_tests = pca_result$group_tests |>
    mutate(trait_set = "five_non_head_traits"),
  pca_group_centroids = pca_centroids |>
    mutate(trait_set = "five_non_head_traits"),
  pca_head_excluded_summary = pca_no_head_summary |>
    mutate(
      trait_set = "five_non_head_traits",
      size_adjustment = "group_specific_at_7.2_mm"
    ),
  pca_six_trait_sensitivity_parallel_analysis = pca_six_result$variance |>
    mutate(
      trait_set = "six_traits_including_head",
      size_adjustment = "group_specific_at_7.2_mm"
    ),
  pca_six_trait_sensitivity_loadings = pca_six_loadings |>
    mutate(trait_set = "six_traits_including_head"),
  pca_six_trait_sensitivity_group_tests = pca_six_result$group_tests |>
    mutate(trait_set = "six_traits_including_head"),
  pca_six_trait_sensitivity_group_centroids = pca_six_centroids |>
    mutate(trait_set = "six_traits_including_head"),
  pca_classification_sensitivity =
    pca_classification_sensitivity$group_tests |>
    mutate(trait_set = "five_non_head_traits"),
  pca_classification_sensitivity_variance =
    pca_classification_sensitivity$variance |>
    mutate(
      trait_set = "five_non_head_traits",
      size_adjustment = "group_specific_at_7.2_mm"
    ),
  covariance_matrix_omnibus_tests = covariance_omnibus,
  covariance_matrix_pairwise_tests = covariance_pairwise,
  covariance_matrix_similarity = covariance_similarity,
  covariance_matrix_rarefaction = covariance_rarefaction,
  covariance_head_excluded_omnibus = covariance_no_head_omnibus,
  covariance_six_trait_sensitivity_omnibus = covariance_six_trait_omnibus,
  group_trait_correlations = group_correlation_long,
  head_leg_raw_type3_anova = head_leg_raw_anova,
  head_leg_raw_slopes = head_leg_raw_slopes,
  head_leg_raw_strategy_slope_contrasts = head_leg_raw_strategy_slopes,
  ear_foretibia_raw_model_test = ear_tibia_raw_test,
  ear_foretibia_raw_slopes = ear_tibia_raw_slopes,
  eye_head_raw_model_test = eye_head_raw_test,
  eye_head_raw_slopes = eye_head_raw_slopes,
  ear_head_raw_model_test = ear_head_raw_test,
  ear_head_raw_slopes = ear_head_raw_slopes,
  eye_group_after_head_test = eye_group_after_head_test,
  eye_head_increment_test = eye_head_increment_test,
  eye_head_support = head_support,
  eye_classification_sensitivity = eye_classification_sensitivity,
  leg_intermediacy_equivalence = leg_intermediacy_equivalence,
  sensory_intermediacy_equivalence = sensory_intermediacy_equivalence,
  weapon_ear_equivalence = weapon_ear_equivalence,
  focal_multiplicity_sensitivity = focal_tests,
  fragile_contrast_disposition = fragile_contrast,
  adjusted_trait_means_combined = adjusted_trait_means_combined |>
    mutate(
      group = as.character(group),
      trait = as.character(trait)
    ),
  analysis_hierarchy = analysis_hierarchy
)

iwalk(outputs, ~ readr::write_csv(.x, file.path(tables_dir, paste0(.y, ".csv"))))

saveRDS(
  list(
    pca = pca_result,
    pca_no_head = pca_no_head,
    pca_six_trait_sensitivity = pca_six_result,
    pca_classification_sensitivity = pca_classification_sensitivity,
    covariance_five_trait_omnibus = covariance_omnibus,
    covariance_six_trait_sensitivity_omnibus = covariance_six_trait_omnibus,
    head_leg_raw = mod_head_leg_raw,
    head_leg_raw_trends = head_leg_raw_trends,
    head_leg_raw_strategy_slopes = head_leg_raw_strategy_slopes,
    ear_tibia_raw_common = mod_ear_tibia_raw_common,
    ear_tibia_raw_full = mod_ear_tibia_raw_full,
    eye_head_raw_common = mod_eye_head_raw_common,
    eye_head_raw_full = mod_eye_head_raw_full,
    ear_head_raw_common = mod_ear_head_raw_common,
    ear_head_raw_full = mod_ear_head_raw_full,
    eye_head_no_group = mod_eye_head_no_group,
    eye_head_with_group = mod_eye_head_with_group
  ),
  file.path(models_dir, "reviewer_reanalysis_models.rds")
)

capture.output(
  {
    cat("REVIEWER REANALYSIS SUMMARY\n\n")
    cat("Primary five-trait PCA parallel analysis\n")
    print(pca_result$variance)
    print(pca_result$group_tests)
    cat("\nSix-trait PCA sensitivity\n")
    print(pca_six_result$variance)
    cat("\nPrimary five-trait covariance/correlation equality\n")
    print(covariance_omnibus)
    print(covariance_pairwise)
    cat("\nSix-trait covariance/correlation sensitivity\n")
    print(covariance_six_trait_omnibus)
    cat("\nRaw-trait refits\n")
    print(head_leg_raw_anova)
    print(head_leg_raw_strategy_slopes)
    print(ear_tibia_raw_test)
    print(eye_head_raw_test)
    print(ear_head_raw_test)
    print(eye_group_after_head_test)
    cat("\nEquivalence analyses\n")
    print(leg_intermediacy_equivalence)
    print(sensory_intermediacy_equivalence)
    print(weapon_ear_equivalence)
    cat("\nSampling and classification checks\n")
    print(morph_frequency_tests)
    print(classification_summary)
    cat("\nMultiplicity sensitivity\n")
    print(focal_tests)
  },
  file = file.path(output_root, "reviewer_reanalysis_summary.txt")
)

capture.output(sessionInfo(), file = file.path(output_root, "reviewer_reanalysis_session_info.txt"))

cat("Reviewer reanalysis completed. Outputs written to:\n", normalizePath(output_root), "\n")
