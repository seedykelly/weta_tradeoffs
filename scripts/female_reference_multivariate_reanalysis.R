#!/usr/bin/env Rscript

# Female-referenced multivariate contrast-vector analysis
#
# This script deliberately does not fit a phenotypic trajectory model. The
# three male morphs are cross-sectional terminal groups, and females form one
# shared reference group. Repeating the same females as three independent
# starting states would pseudoreplicate that reference. Instead, all contrast
# vectors are estimated jointly from one multivariate model in which every
# specimen appears once.
#
# The primary analysis excludes head size because head length contributes to
# the archived male morph assignment. It therefore uses five logged
# linear-dimension traits (foreleg, midleg, hindleg, ear and eye) standardized
# by their pooled residual SD after body-size correction. Six-trait
# standardized and raw-log analyses are retained as sensitivities. Inference
# uses Freedman-Lane residual randomization with an add-one p-value correction.

args <- commandArgs(trailingOnly = TRUE)

data_file <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "data/trait_data.csv"
}

output_root <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "analysis_outputs/weta_trait_analysis"
}

permutations <- if (length(args) >= 3L) as.integer(args[[3L]]) else 9999L
bootstrap_reps <- if (length(args) >= 4L) as.integer(args[[4L]]) else 9999L

if (is.na(permutations) || permutations < 999L) {
  stop("Use at least 999 residual randomizations.")
}
if (is.na(bootstrap_reps) || bootstrap_reps < 999L) {
  stop("Use at least 999 stratified bootstrap replicates.")
}

tables_dir <- file.path(output_root, "tables")
models_dir <- file.path(output_root, "models")
figures_dir <- file.path(output_root, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

reference_pronotum <- 7.2
group_levels <- c("female", "eighth", "ninth", "tenth")
group_labels <- c(
  female = "Female",
  eighth = "Eighth instar",
  ninth = "Ninth instar",
  tenth = "Tenth instar"
)
vector_labels <- c(eighth = "D8", ninth = "D9", tenth = "D10")

matrix_fit <- function(X, Y) {
  qx <- qr(X)
  coefficients <- qr.coef(qx, Y)
  fitted_values <- X %*% coefficients
  residuals <- Y - fitted_values
  list(
    coefficients = coefficients,
    fitted = fitted_values,
    residuals = residuals,
    SSE = sum(residuals^2),
    rank = qx$rank,
    residual_df = nrow(X) - qx$rank
  )
}

nested_statistics <- function(Y, X_reduced, X_full) {
  reduced <- matrix_fit(X_reduced, Y)
  full <- matrix_fit(X_full, Y)
  term_df <- full$rank - reduced$rank
  term_SS <- reduced$SSE - full$SSE
  pseudo_F <- (term_SS / term_df) / (full$SSE / full$residual_df)
  total_SS <- sum(scale(Y, center = TRUE, scale = FALSE)^2)
  c(
    term_SS = term_SS,
    residual_SS = full$SSE,
    numerator_df = term_df,
    residual_df = full$residual_df,
    pseudo_F = pseudo_F,
    R2_total = term_SS / total_SS,
    R2_partial = term_SS / (term_SS + full$SSE)
  )
}

rrpp_nested <- function(Y, X_reduced, X_full, reps, seed) {
  set.seed(seed)
  observed <- nested_statistics(Y, X_reduced, X_full)
  reduced <- matrix_fit(X_reduced, Y)
  null_F <- numeric(reps)

  for (iteration in seq_len(reps)) {
    pseudo_Y <- reduced$fitted +
      reduced$residuals[sample.int(nrow(Y)), , drop = FALSE]
    null_F[[iteration]] <- nested_statistics(
      pseudo_Y,
      X_reduced,
      X_full
    )[["pseudo_F"]]
  }

  p_value <- (1 + sum(null_F >= observed[["pseudo_F"]])) / (reps + 1)
  Z_logF <- (
    log(observed[["pseudo_F"]]) - mean(log(null_F))
  ) / stats::sd(log(null_F))

  list(
    statistics = c(observed, p_value = p_value, Z_logF = Z_logF),
    null_F = null_F
  )
}

reference_mean_rows <- function(interaction = TRUE) {
  new_data <- data.frame(
    logP_c = rep(0, length(group_levels)),
    group = factor(group_levels, levels = group_levels)
  )
  if (interaction) {
    model.matrix(~ logP_c * group, data = new_data)
  } else {
    model.matrix(~ logP_c + group, data = new_data)
  }
}

attributes_from_vectors <- function(vectors) {
  magnitudes <- sqrt(rowSums(vectors^2))
  vector_pairs <- combn(rownames(vectors), 2L, simplify = FALSE)
  pair_names <- vapply(
    vector_pairs,
    paste,
    collapse = "_vs_",
    FUN.VALUE = character(1)
  )

  angles <- vapply(
    vector_pairs,
    function(pair) {
      cosine <- sum(vectors[pair[[1L]], ] * vectors[pair[[2L]], ]) /
        (magnitudes[[pair[[1L]]]] * magnitudes[[pair[[2L]]]])
      acos(max(-1, min(1, cosine))) * 180 / pi
    },
    numeric(1)
  )
  magnitude_differences <- vapply(
    vector_pairs,
    function(pair) {
      abs(magnitudes[[pair[[1L]]]] - magnitudes[[pair[[2L]]]])
    },
    numeric(1)
  )
  names(angles) <- pair_names
  names(magnitude_differences) <- pair_names

  list(
    vectors = vectors,
    magnitudes = magnitudes,
    angles = angles,
    magnitude_differences = magnitude_differences
  )
}

vector_attributes <- function(coefficients, mean_rows) {
  adjusted_means <- mean_rows %*% coefficients
  rownames(adjusted_means) <- group_levels

  female_mean <- matrix(
    adjusted_means["female", ],
    nrow = 3L,
    ncol = ncol(adjusted_means),
    byrow = TRUE
  )
  vectors <- adjusted_means[c("eighth", "ninth", "tenth"), , drop = FALSE] -
    female_mean
  rownames(vectors) <- unname(vector_labels)
  attributes <- attributes_from_vectors(vectors)
  attributes$adjusted_means <- adjusted_means
  attributes[c(
    "adjusted_means", "vectors", "magnitudes", "angles",
    "magnitude_differences"
  )]
}

shared_anchor_randomization <- function(
  Y,
  data,
  X_full,
  X_shared_anchor,
  mean_rows,
  reps,
  seed
) {
  set.seed(seed)
  full <- matrix_fit(X_full, Y)
  reduced <- matrix_fit(X_shared_anchor, Y)
  observed <- vector_attributes(full$coefficients, mean_rows)
  observed_global <- nested_statistics(Y, X_shared_anchor, X_full)
  null_F <- numeric(reps)

  for (iteration in seq_len(reps)) {
    pseudo_Y <- reduced$fitted +
      reduced$residuals[sample.int(nrow(Y)), , drop = FALSE]
    pseudo_full <- matrix_fit(X_full, pseudo_Y)
    pseudo_reduced <- matrix_fit(X_shared_anchor, pseudo_Y)
    term_SS <- pseudo_reduced$SSE - pseudo_full$SSE
    term_df <- pseudo_full$rank - pseudo_reduced$rank
    null_F[[iteration]] <- (term_SS / term_df) /
      (pseudo_full$SSE / pseudo_full$residual_df)
  }

  global_p <- (
    1 + sum(null_F >= observed_global[["pseudo_F"]])
  ) / (reps + 1)
  global_Z <- (
    log(observed_global[["pseudo_F"]]) - mean(log(null_F))
  ) / stats::sd(log(null_F))

  list(
    observed = observed,
    global_test = c(
      observed_global,
      p_value = global_p,
      Z_logF = global_Z
    ),
    null_F = null_F
  )
}

stratified_vector_bootstrap <- function(
  Y,
  data,
  full_formula,
  mean_rows,
  reps,
  seed
) {
  set.seed(seed)
  row_indices <- split(seq_len(nrow(data)), data$group)
  angle_names <- c("D8_vs_D9", "D8_vs_D10", "D9_vs_D10")
  vector_names <- c("D8", "D9", "D10")

  bootstrap_angles <- matrix(
    NA_real_,
    nrow = reps,
    ncol = 3L,
    dimnames = list(NULL, angle_names)
  )
  bootstrap_magnitudes <- matrix(
    NA_real_,
    nrow = reps,
    ncol = 3L,
    dimnames = list(NULL, vector_names)
  )
  bootstrap_magnitude_differences <- matrix(
    NA_real_,
    nrow = reps,
    ncol = 3L,
    dimnames = list(
      NULL,
      c("D8_vs_D9", "D8_vs_D10", "D9_vs_D10")
    )
  )

  for (iteration in seq_len(reps)) {
    sampled_indices <- unlist(
      lapply(
        row_indices,
        function(index) sample(index, length(index), replace = TRUE)
      ),
      use.names = FALSE
    )
    bootstrap_data <- data[sampled_indices, , drop = FALSE]
    bootstrap_Y <- Y[sampled_indices, , drop = FALSE]
    bootstrap_X <- model.matrix(full_formula, data = bootstrap_data)
    attributes <- vector_attributes(
      matrix_fit(bootstrap_X, bootstrap_Y)$coefficients,
      mean_rows
    )
    bootstrap_angles[iteration, ] <- attributes$angles
    bootstrap_magnitudes[iteration, ] <- attributes$magnitudes
    bootstrap_magnitude_differences[iteration, ] <- c(
      D8_vs_D9 = attributes$magnitudes[["D9"]] -
        attributes$magnitudes[["D8"]],
      D8_vs_D10 = attributes$magnitudes[["D10"]] -
        attributes$magnitudes[["D8"]],
      D9_vs_D10 = attributes$magnitudes[["D10"]] -
        attributes$magnitudes[["D9"]]
    )
  }

  list(
    angle_CI = apply(
      bootstrap_angles,
      2L,
      quantile,
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    ),
    magnitude_CI = apply(
      bootstrap_magnitudes,
      2L,
      quantile,
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    ),
    magnitude_difference_CI = apply(
      bootstrap_magnitude_differences,
      2L,
      quantile,
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    ),
    angles = bootstrap_angles,
    magnitudes = bootstrap_magnitudes,
    magnitude_differences = bootstrap_magnitude_differences
  )
}

rank_one_fit <- function(Y, Z, W) {
  # Z contains nuisance terms and W contains two or more female-referenced
  # intercept contrasts. Rank(B) = 1 allows the W vectors to have different
  # magnitudes but requires them to share one direction. The inferential test
  # statistic is the proportion of design-weighted group-effect SS beyond the
  # first singular dimension.
  qZ <- qr(Z)
  residual_Y <- qr.resid(qZ, Y)
  residual_W <- qr.resid(qZ, W)
  R <- chol(crossprod(residual_W))
  Q <- residual_W %*% solve(R)
  A <- crossprod(Q, residual_Y)
  decomposition <- svd(A)

  A_rank_one <- (
    decomposition$u[, 1L, drop = FALSE] * decomposition$d[[1L]]
  ) %*% t(decomposition$v[, 1L, drop = FALSE])
  B_rank_one <- solve(R, A_rank_one)
  nuisance_coefficients <- qr.coef(qZ, Y - W %*% B_rank_one)
  rank_one_fitted <- Z %*% nuisance_coefficients + W %*% B_rank_one
  rank_one_residuals <- Y - rank_one_fitted

  full_X <- cbind(Z, W)
  full_fit <- matrix_fit(full_X, Y)
  B_unrestricted <- solve(R, A)
  beyond_rank_one_SS <- sum(decomposition$d[-1L]^2)
  design_weighted_group_effect_SS <- sum(decomposition$d^2)
  proportion_beyond_rank_one <-
    beyond_rank_one_SS / design_weighted_group_effect_SS

  list(
    singular_values = decomposition$d,
    B_unrestricted = B_unrestricted,
    B_rank_one = B_rank_one,
    fitted = rank_one_fitted,
    residuals = rank_one_residuals,
    full_residual_SS = full_fit$SSE,
    beyond_rank_one_SS = beyond_rank_one_SS,
    design_weighted_group_effect_SS = design_weighted_group_effect_SS,
    proportion_beyond_rank_one = proportion_beyond_rank_one,
    rank_one_fraction = 1 - proportion_beyond_rank_one
  )
}

rank_one_randomization <- function(Y, Z, W, reps, seed) {
  set.seed(seed)
  observed <- rank_one_fit(Y, Z, W)
  observed_attributes <- attributes_from_vectors(observed$B_unrestricted)
  null_proportion_beyond_rank_one <- numeric(reps)
  null_angles <- matrix(
    NA_real_,
    nrow = reps,
    ncol = length(observed_attributes$angles),
    dimnames = list(NULL, names(observed_attributes$angles))
  )

  for (iteration in seq_len(reps)) {
    pseudo_Y <- observed$fitted +
      observed$residuals[sample.int(nrow(Y)), , drop = FALSE]
    pseudo_fit <- rank_one_fit(
      pseudo_Y,
      Z,
      W
    )
    null_proportion_beyond_rank_one[[iteration]] <-
      pseudo_fit$proportion_beyond_rank_one
    null_angles[iteration, ] <- attributes_from_vectors(
      pseudo_fit$B_unrestricted
    )$angles
  }

  p_value <- (
    1 + sum(
      null_proportion_beyond_rank_one >=
        observed$proportion_beyond_rank_one
    )
  ) / (reps + 1)
  angle_p <- vapply(
    seq_along(observed_attributes$angles),
    function(column) {
      (
        1 + sum(
          null_angles[, column] >= observed_attributes$angles[[column]]
        )
      ) / (reps + 1)
    },
    numeric(1)
  )
  names(angle_p) <- names(observed_attributes$angles)

  list(
    observed = observed,
    observed_attributes = observed_attributes,
    p_value = p_value,
    angle_p = angle_p,
    angle_p_Holm = p.adjust(angle_p, method = "holm"),
    null_proportion_beyond_rank_one =
      null_proportion_beyond_rank_one,
    null_angles = null_angles,
    null_quantiles = quantile(
      null_proportion_beyond_rank_one,
      probs = c(0.5, 0.95, 0.975, 0.99),
      names = TRUE
    )
  )
}

pairwise_rank_one_randomization <- function(Y, Z, W, reps, seed) {
  # Each test constrains only the named pair to share a direction. The third
  # female-referenced vector is added to the nuisance design and is therefore
  # estimated without a directional constraint. This produces genuine
  # pair-specific angular tests rather than local statistics from the global
  # three-vector rank-one null.
  vector_pairs <- combn(colnames(W), 2L, simplify = FALSE)
  angle_p <- numeric(length(vector_pairs))
  observed_angles <- numeric(length(vector_pairs))
  pair_results <- vector("list", length(vector_pairs))

  for (pair_index in seq_along(vector_pairs)) {
    pair <- vector_pairs[[pair_index]]
    third_vector <- setdiff(colnames(W), pair)
    pair_nuisance <- cbind(
      Z,
      W[, third_vector, drop = FALSE]
    )
    pair_contrasts <- W[, pair, drop = FALSE]
    pair_result <- rank_one_randomization(
      Y,
      pair_nuisance,
      pair_contrasts,
      reps = reps,
      seed = seed + pair_index - 1L
    )

    observed_angles[[pair_index]] <-
      pair_result$observed_attributes$angles[[1L]]
    angle_p[[pair_index]] <- pair_result$angle_p[[1L]]
    pair_results[[pair_index]] <- pair_result
  }

  pair_names <- vapply(
    vector_pairs,
    paste,
    collapse = "_vs_",
    FUN.VALUE = character(1)
  )
  names(observed_angles) <- pair_names
  names(angle_p) <- pair_names
  names(pair_results) <- pair_names

  list(
    observed_angles = observed_angles,
    angle_p = angle_p,
    angle_p_Holm = p.adjust(angle_p, method = "holm"),
    pair_results = pair_results
  )
}

dat_raw <- read.csv(data_file, check.names = FALSE)
required_columns <- c(
  "ID", "sex", "morph", "pronotum", "head_length", "head_width",
  "forefemur", "foretibia", "midfemur", "midtibia", "hindfemur",
  "hindtibia", "ear", "eye"
)
if (!all(required_columns %in% names(dat_raw))) {
  stop("The input data do not contain all required columns.")
}

dat_raw$group <- factor(dat_raw$morph, levels = group_levels)
dat_raw$sex_class <- factor(
  ifelse(dat_raw$group == "female", "female", "male"),
  levels = c("female", "male")
)
dat_raw$head_size <- sqrt(dat_raw$head_length * dat_raw$head_width)
dat_raw$foreleg <- dat_raw$forefemur + dat_raw$foretibia
dat_raw$midleg <- dat_raw$midfemur + dat_raw$midtibia
dat_raw$hindleg <- dat_raw$hindfemur + dat_raw$hindtibia
dat_raw$ear_linear <- sqrt(dat_raw$ear)
dat_raw$logP_c <- log(dat_raw$pronotum) - log(reference_pronotum)

Y_raw_log <- cbind(
  head = log(dat_raw$head_size),
  foreleg = log(dat_raw$foreleg),
  midleg = log(dat_raw$midleg),
  hindleg = log(dat_raw$hindleg),
  ear = log(dat_raw$ear_linear),
  eye = log(dat_raw$eye)
)

complete_rows <- complete.cases(
  Y_raw_log,
  dat_raw$logP_c,
  dat_raw$group
)
dat <- droplevels(dat_raw[complete_rows, , drop = FALSE])
Y_raw_log <- Y_raw_log[complete_rows, , drop = FALSE]

X_size <- model.matrix(~ logP_c, data = dat)
size_residuals <- qr.resid(qr(X_size), Y_raw_log)
response_residual_SD <- apply(size_residuals, 2L, stats::sd)
Y_standardized <- sweep(Y_raw_log, 2L, response_residual_SD, "/")

X_additive <- model.matrix(~ logP_c + group, data = dat)
X_interaction <- model.matrix(~ logP_c * group, data = dat)

# The shared-anchor null preserves all four group-specific size slopes but
# makes the three male adjusted centroids identical at logP_c = 0.
X_shared_anchor <- model.matrix(
  ~ sex_class + group:logP_c,
  data = dat
)
mean_rows_interaction <- reference_mean_rows(interaction = TRUE)

# Equivalent rank-one parameterization. The single intercept in Z is the
# female adjusted mean at logP_c = 0; W adds the three female-referenced male
# contrasts. The four columns after the intercept retain group-specific slopes.
group_indicators <- model.matrix(~ 0 + group, data = dat)
Z_rank_one <- cbind(
  intercept = 1,
  sweep(group_indicators, 1L, dat$logP_c, "*")
)
W_rank_one <- model.matrix(~ group, data = dat)[, -1L, drop = FALSE]
colnames(W_rank_one) <- unname(vector_labels)

scenario_specifications <- list(
  standardized_head_excluded = list(
    analysis_tier = "primary",
    primary_analysis = TRUE,
    analysis_reason = paste(
      "Head size is excluded because head length contributes to archived",
      "male morph assignment."
    ),
    analysis_space = "standardized_size_adjusted_SD",
    trait_set = "head_excluded_five_traits",
    Y = Y_standardized[, setdiff(colnames(Y_standardized), "head"), drop = FALSE],
    seeds = c(group = 3000L, slope = 3004L, shared = 3001L,
              bootstrap = 3003L, rank_one = 5002L,
              pairwise_rank_one = 5104L)
  ),
  standardized_six = list(
    analysis_tier = "sensitivity",
    primary_analysis = FALSE,
    analysis_reason = paste(
      "Adds head size to test whether inference depends on including the",
      "classification-linked weapon trait."
    ),
    analysis_space = "standardized_size_adjusted_SD",
    trait_set = "six_traits",
    Y = Y_standardized,
    seeds = c(group = 1001L, slope = 1002L, shared = 1005L,
              bootstrap = 1006L, rank_one = 5001L,
              pairwise_rank_one = 5101L)
  ),
  raw_log_six = list(
    analysis_tier = "sensitivity",
    primary_analysis = FALSE,
    analysis_reason = paste(
      "Uses unstandardized log traits to test dependence on residual-SD",
      "standardization while retaining all six traits."
    ),
    analysis_space = "raw_log_linear_dimensions",
    trait_set = "six_traits",
    Y = Y_raw_log,
    seeds = c(group = 2001L, slope = 2002L, shared = 2005L,
              bootstrap = 2006L, rank_one = 5003L,
              pairwise_rank_one = 5107L)
  )
)

scenario_metadata <- do.call(
  rbind,
  lapply(
    names(scenario_specifications),
    function(scenario_name) {
      specification <- scenario_specifications[[scenario_name]]
      data.frame(
        scenario = scenario_name,
        analysis_tier = specification$analysis_tier,
        primary_analysis = specification$primary_analysis,
        analysis_reason = specification$analysis_reason,
        analysis_space = specification$analysis_space,
        trait_set = specification$trait_set,
        stringsAsFactors = FALSE
      )
    }
  )
)
primary_scenarios <- scenario_metadata$scenario[
  scenario_metadata$primary_analysis
]
if (length(primary_scenarios) != 1L) {
  stop("Exactly one scenario must be marked as the primary analysis.")
}
primary_scenario <- primary_scenarios[[1L]]

run_scenario <- function(specification) {
  Y <- specification$Y
  seeds <- specification$seeds

  group_test <- rrpp_nested(
    Y,
    X_size,
    X_additive,
    reps = permutations,
    seed = seeds[["group"]]
  )
  slope_test <- rrpp_nested(
    Y,
    X_additive,
    X_interaction,
    reps = permutations,
    seed = seeds[["slope"]]
  )
  shared_anchor <- shared_anchor_randomization(
    Y,
    dat,
    X_interaction,
    X_shared_anchor,
    mean_rows_interaction,
    reps = permutations,
    seed = seeds[["shared"]]
  )
  bootstrap <- stratified_vector_bootstrap(
    Y,
    dat,
    ~ logP_c * group,
    mean_rows_interaction,
    reps = bootstrap_reps,
    seed = seeds[["bootstrap"]]
  )
  rank_one <- rank_one_randomization(
    Y,
    Z_rank_one,
    W_rank_one,
    reps = permutations,
    seed = seeds[["rank_one"]]
  )
  pairwise_rank_one <- pairwise_rank_one_randomization(
    Y,
    Z_rank_one,
    W_rank_one,
    reps = permutations,
    seed = seeds[["pairwise_rank_one"]]
  )

  # The two parameterizations must give the same unrestricted vectors.
  maximum_parameterization_difference <- max(abs(
    shared_anchor$observed$vectors - rank_one$observed$B_unrestricted
  ))
  if (maximum_parameterization_difference > 1e-8) {
    stop("Shared-anchor and rank-one parameterizations disagree.")
  }
  maximum_pairwise_angle_difference <- max(abs(
    shared_anchor$observed$angles -
      pairwise_rank_one$observed_angles[
        names(shared_anchor$observed$angles)
      ]
  ))
  if (maximum_pairwise_angle_difference > 1e-8) {
    stop("Global and pair-specific unrestricted angles disagree.")
  }

  list(
    specification = specification[c(
      "analysis_tier", "primary_analysis", "analysis_reason",
      "analysis_space", "trait_set", "seeds"
    )],
    group_test = group_test,
    slope_test = slope_test,
    shared_anchor = shared_anchor,
    bootstrap = bootstrap,
    rank_one = rank_one,
    pairwise_rank_one = pairwise_rank_one,
    maximum_parameterization_difference = maximum_parameterization_difference,
    maximum_pairwise_angle_difference = maximum_pairwise_angle_difference
  )
}

results <- lapply(scenario_specifications, run_scenario)

model_test_rows <- list()
component_rows <- list()
magnitude_rows <- list()
angle_rows <- list()
magnitude_difference_rows <- list()
rank_one_rows <- list()

for (scenario_name in names(results)) {
  result <- results[[scenario_name]]
  analysis_space <- result$specification$analysis_space
  trait_set <- result$specification$trait_set

  tests <- list(
    group_adjusted_for_size_common_slopes = result$group_test$statistics,
    group_by_size_slope_heterogeneity = result$slope_test$statistics,
    male_morph_centroids_at_reference_with_group_slopes =
      result$shared_anchor$global_test
  )
  for (test_name in names(tests)) {
    statistic <- tests[[test_name]]
    model_test_rows[[length(model_test_rows) + 1L]] <- data.frame(
      scenario = scenario_name,
      analysis_space = analysis_space,
      trait_set = trait_set,
      test = test_name,
      term_SS = unname(statistic[["term_SS"]]),
      residual_SS = unname(statistic[["residual_SS"]]),
      numerator_df = unname(statistic[["numerator_df"]]),
      residual_df = unname(statistic[["residual_df"]]),
      pseudo_F = unname(statistic[["pseudo_F"]]),
      R2_total = unname(statistic[["R2_total"]]),
      R2_partial = unname(statistic[["R2_partial"]]),
      Z_logF = unname(statistic[["Z_logF"]]),
      permutations = permutations,
      p_value = unname(statistic[["p_value"]]),
      stringsAsFactors = FALSE
    )
  }

  vectors <- result$shared_anchor$observed$vectors
  for (vector_name in rownames(vectors)) {
    for (trait_name in colnames(vectors)) {
      component_rows[[length(component_rows) + 1L]] <- data.frame(
        scenario = scenario_name,
        analysis_space = analysis_space,
        trait_set = trait_set,
        vector = vector_name,
        trait = trait_name,
        component = unname(vectors[vector_name, trait_name]),
        stringsAsFactors = FALSE
      )
    }
  }

  magnitudes <- result$shared_anchor$observed$magnitudes
  magnitude_CI <- result$bootstrap$magnitude_CI
  for (vector_name in names(magnitudes)) {
    magnitude_rows[[length(magnitude_rows) + 1L]] <- data.frame(
      scenario = scenario_name,
      analysis_space = analysis_space,
      trait_set = trait_set,
      vector = vector_name,
      magnitude = unname(magnitudes[[vector_name]]),
      bootstrap_lower_95 = magnitude_CI["2.5%", vector_name],
      bootstrap_median = magnitude_CI["50%", vector_name],
      bootstrap_upper_95 = magnitude_CI["97.5%", vector_name],
      bootstrap_reps = bootstrap_reps,
      stringsAsFactors = FALSE
    )
  }

  angles <- result$shared_anchor$observed$angles
  angle_CI <- result$bootstrap$angle_CI
  for (comparison in names(angles)) {
    angle_rows[[length(angle_rows) + 1L]] <- data.frame(
      scenario = scenario_name,
      analysis_space = analysis_space,
      trait_set = trait_set,
      comparison = comparison,
      angle_degrees = unname(angles[[comparison]]),
      bootstrap_lower_95 = angle_CI["2.5%", comparison],
      bootstrap_median = angle_CI["50%", comparison],
      bootstrap_upper_95 = angle_CI["97.5%", comparison],
      permutation_p = result$pairwise_rank_one$angle_p[[comparison]],
      Holm_p = result$pairwise_rank_one$angle_p_Holm[[comparison]],
      permutation_null =
        paste0(
          "pair_specific_rank_one_common_direction_",
          "third_vector_unrestricted"
        ),
      permutations = permutations,
      bootstrap_reps = bootstrap_reps,
      stringsAsFactors = FALSE
    )
  }

  magnitude_differences <-
    result$shared_anchor$observed$magnitude_differences
  magnitude_difference_CI <- result$bootstrap$magnitude_difference_CI
  ordered_difference_labels <- c(
    D8_vs_D9 = "D9 - D8",
    D8_vs_D10 = "D10 - D8",
    D9_vs_D10 = "D10 - D9"
  )
  for (comparison in names(magnitude_differences)) {
    signed_difference <- switch(
      comparison,
      D8_vs_D9 = result$shared_anchor$observed$magnitudes[["D9"]] -
        result$shared_anchor$observed$magnitudes[["D8"]],
      D8_vs_D10 = result$shared_anchor$observed$magnitudes[["D10"]] -
        result$shared_anchor$observed$magnitudes[["D8"]],
      D9_vs_D10 = result$shared_anchor$observed$magnitudes[["D10"]] -
        result$shared_anchor$observed$magnitudes[["D9"]]
    )
    lower_CI <- magnitude_difference_CI["2.5%", comparison]
    upper_CI <- magnitude_difference_CI["97.5%", comparison]
    magnitude_difference_rows[[length(magnitude_difference_rows) + 1L]] <-
      data.frame(
        scenario = scenario_name,
        analysis_space = analysis_space,
        trait_set = trait_set,
        comparison = comparison,
        ordered_difference = ordered_difference_labels[[comparison]],
        signed_magnitude_difference = unname(signed_difference),
        absolute_magnitude_difference =
          unname(magnitude_differences[[comparison]]),
        bootstrap_lower_95 = lower_CI,
        bootstrap_median =
          magnitude_difference_CI["50%", comparison],
        bootstrap_upper_95 = upper_CI,
        bootstrap_95_excludes_zero =
          (lower_CI > 0) || (upper_CI < 0),
        bootstrap_reps = bootstrap_reps,
        stringsAsFactors = FALSE
      )
  }

  rank_result <- result$rank_one
  singular_values <- rank_result$observed$singular_values
  singular_values <- c(singular_values, rep(NA_real_, 3L - length(singular_values)))
  rank_one_rows[[length(rank_one_rows) + 1L]] <- data.frame(
    scenario = scenario_name,
    analysis_space = analysis_space,
    trait_set = trait_set,
    singular_value_1 = singular_values[[1L]],
    singular_value_2 = singular_values[[2L]],
    singular_value_3 = singular_values[[3L]],
    design_weighted_group_effect_SS =
      rank_result$observed$design_weighted_group_effect_SS,
    beyond_rank_one_SS = rank_result$observed$beyond_rank_one_SS,
    rank_one_fraction = rank_result$observed$rank_one_fraction,
    proportion_beyond_rank_one =
      rank_result$observed$proportion_beyond_rank_one,
    null_median = rank_result$null_quantiles[["50%"]],
    null_95 = rank_result$null_quantiles[["95%"]],
    null_97_5 = rank_result$null_quantiles[["97.5%"]],
    null_99 = rank_result$null_quantiles[["99%"]],
    permutations = permutations,
    p_value = rank_result$p_value,
    stringsAsFactors = FALSE
  )
}

annotate_analysis_tier <- function(table) {
  metadata_index <- match(table$scenario, scenario_metadata$scenario)
  table$analysis_tier <- scenario_metadata$analysis_tier[metadata_index]
  table$primary_analysis <- scenario_metadata$primary_analysis[metadata_index]
  table$analysis_reason <- scenario_metadata$analysis_reason[metadata_index]
  leading_columns <- c(
    "scenario", "analysis_tier", "primary_analysis", "analysis_reason"
  )
  table[, c(leading_columns, setdiff(names(table), leading_columns)), drop = FALSE]
}

model_tests <- annotate_analysis_tier(do.call(rbind, model_test_rows))
vector_components <- annotate_analysis_tier(do.call(rbind, component_rows))
vector_magnitudes <- annotate_analysis_tier(do.call(rbind, magnitude_rows))
vector_angles <- annotate_analysis_tier(do.call(rbind, angle_rows))
magnitude_differences <- annotate_analysis_tier(
  do.call(rbind, magnitude_difference_rows)
)
rank_one_tests <- annotate_analysis_tier(do.call(rbind, rank_one_rows))

# Back-transform raw-log components for descriptive interpretation. Ear was
# linearized as sqrt(area), so its raw-log difference is doubled before the
# area ratio is calculated.
raw_components <- subset(
  vector_components,
  scenario == "raw_log_six"
)
raw_components$percent_difference <- ifelse(
  raw_components$trait == "ear",
  100 * (exp(2 * raw_components$component) - 1),
  100 * (exp(raw_components$component) - 1)
)
raw_components$reported_scale <- ifelse(
  raw_components$trait == "ear",
  "ear_area_percent",
  "linear_trait_percent"
)
raw_percent_contrasts <- raw_components[, c(
  "scenario", "analysis_tier", "primary_analysis", "analysis_reason",
  "vector", "trait", "component", "percent_difference", "reported_scale"
)]

group_counts <- as.data.frame(table(dat$group), stringsAsFactors = FALSE)
names(group_counts) <- c("group", "complete_case_n")
group_counts$group_label <- unname(group_labels[as.character(group_counts$group)])

pronotum_ranges <- do.call(
  rbind,
  lapply(
    split(dat$pronotum, dat$group),
    function(values) {
      c(n = length(values), minimum = min(values), maximum = max(values))
    }
  )
)
pronotum_ranges <- data.frame(
  group = rownames(pronotum_ranges),
  pronotum_ranges,
  row.names = NULL,
  check.names = FALSE
)

analysis_metadata <- data.frame(
  item = c(
    "input_rows", "complete_case_rows", "omitted_IDs",
    "reference_pronotum_mm", "permutations", "bootstrap_reps",
    "primary_scenario", "primary_trait_set", "primary_excludes_head",
    "RRPP_installed", "geomorph_installed", "evolqg_installed",
    paste0("response_residual_SD_", names(response_residual_SD))
  ),
  value = c(
    nrow(dat_raw),
    nrow(dat),
    paste(dat_raw$ID[!complete_rows], collapse = ";"),
    reference_pronotum,
    permutations,
    bootstrap_reps,
    primary_scenario,
    scenario_metadata$trait_set[scenario_metadata$primary_analysis],
    TRUE,
    requireNamespace("RRPP", quietly = TRUE),
    requireNamespace("geomorph", quietly = TRUE),
    requireNamespace("evolqg", quietly = TRUE),
    unname(response_residual_SD)
  ),
  stringsAsFactors = FALSE
)

output_tables <- list(
  female_reference_multivariate_model_tests = model_tests,
  female_reference_vector_components = vector_components,
  female_reference_vector_magnitudes = vector_magnitudes,
  female_reference_vector_angles = vector_angles,
  female_reference_magnitude_differences = magnitude_differences,
  female_reference_rank_one_tests = rank_one_tests,
  female_reference_raw_log_percent_contrasts = raw_percent_contrasts,
  female_reference_complete_case_counts = group_counts,
  female_reference_pronotum_ranges = pronotum_ranges,
  female_reference_analysis_metadata = analysis_metadata,
  female_reference_scenario_metadata = scenario_metadata
)

for (table_name in names(output_tables)) {
  write.csv(
    output_tables[[table_name]],
    file.path(tables_dir, paste0(table_name, ".csv")),
    row.names = FALSE,
    na = ""
  )
}

# Publication-ready two-panel summary of the primary standardized five-trait
# analysis. Panel A uses the ordinary SVD of the observed 3 x 5 contrast matrix
# (rather than the design-weighted decomposition used for the rank-one
# inferential test).
primary_vectors <- results[[primary_scenario]]$shared_anchor$observed$vectors
contrast_svd <- svd(primary_vectors)
contrast_coordinates <- primary_vectors %*% contrast_svd$v[, 1:2, drop = FALSE]

# Axis signs are arbitrary; orient axis 1 toward D10 and axis 2 toward D8 to
# make the display stable across platforms and reruns.
if (contrast_coordinates["D10", 1L] < 0) {
  contrast_svd$v[, 1L] <- -contrast_svd$v[, 1L]
  contrast_coordinates[, 1L] <- -contrast_coordinates[, 1L]
}
if (contrast_coordinates["D8", 2L] < 0) {
  contrast_svd$v[, 2L] <- -contrast_svd$v[, 2L]
  contrast_coordinates[, 2L] <- -contrast_coordinates[, 2L]
}
contrast_SS_percent <- 100 * contrast_svd$d^2 / sum(contrast_svd$d^2)

contrast_axis_table <- data.frame(
  scenario = primary_scenario,
  analysis_tier = "primary",
  trait_set = scenario_metadata$trait_set[scenario_metadata$primary_analysis],
  vector = rownames(contrast_coordinates),
  singular_axis_1 = contrast_coordinates[, 1L],
  singular_axis_2 = contrast_coordinates[, 2L],
  axis_1_percent_contrast_SS = contrast_SS_percent[[1L]],
  axis_2_percent_contrast_SS = contrast_SS_percent[[2L]],
  row.names = NULL,
  check.names = FALSE
)
write.csv(
  contrast_axis_table,
  file.path(tables_dir, "female_reference_contrast_singular_axes.csv"),
  row.names = FALSE
)

figure_file <- file.path(figures_dir, "female_reference_contrast_vectors.png")
group_colours <- c(
  D8 = "#2C7FB8",
  D9 = "#7FCDBB",
  D10 = "#D95F0E"
)
trait_labels <- c(
  head = "Head",
  foreleg = "Foreleg",
  midleg = "Midleg",
  hindleg = "Hindleg",
  ear = "Ear",
  eye = "Eye"
)

grDevices::png(
  filename = figure_file,
  width = 2400,
  height = 1200,
  units = "px",
  res = 300,
  bg = "white"
)
old_par <- par(no.readonly = TRUE)
par(
  mfrow = c(1, 2),
  mar = c(4.2, 4.3, 2.2, 0.8),
  oma = c(0.2, 0.2, 0.2, 0.2),
  family = "serif",
  las = 1,
  xaxs = "r",
  yaxs = "r"
)

coordinate_range <- range(c(0, contrast_coordinates), finite = TRUE)
coordinate_padding <- 0.10 * diff(coordinate_range)
x_limits <- range(c(0, contrast_coordinates[, 1L])) + c(-1, 1) *
  coordinate_padding
y_limits <- range(c(0, contrast_coordinates[, 2L])) + c(-1, 1) *
  coordinate_padding
plot(
  NA_real_,
  NA_real_,
  xlim = x_limits,
  ylim = y_limits,
  xlab = sprintf("Singular axis 1 (%.1f%% contrast SS)", contrast_SS_percent[[1L]]),
  ylab = sprintf("Singular axis 2 (%.1f%% contrast SS)", contrast_SS_percent[[2L]]),
  bty = "l",
  cex.axis = 0.85,
  cex.lab = 0.95
)
abline(h = 0, v = 0, col = "grey80", lty = 3)
for (vector_name in rownames(contrast_coordinates)) {
  arrows(
    0,
    0,
    contrast_coordinates[vector_name, 1L],
    contrast_coordinates[vector_name, 2L],
    length = 0.09,
    angle = 22,
    lwd = 2.2,
    col = group_colours[[vector_name]]
  )
  points(
    contrast_coordinates[vector_name, 1L],
    contrast_coordinates[vector_name, 2L],
    pch = 21,
    bg = group_colours[[vector_name]],
    col = "black",
    cex = 1.15
  )
  text(
    contrast_coordinates[vector_name, 1L],
    contrast_coordinates[vector_name, 2L],
    labels = vector_name,
    pos = if (vector_name == "D8") 2 else 4,
    offset = 0.45,
    cex = 0.9,
    font = 2
  )
}
mtext("A", side = 3, adj = 0, line = 0.5, font = 2, cex = 1.05)

component_limits <- range(c(0, primary_vectors), finite = TRUE)
component_padding <- 0.08 * diff(component_limits)
bar_positions <- barplot(
  primary_vectors,
  beside = TRUE,
  col = group_colours[rownames(primary_vectors)],
  border = "grey20",
  lwd = 0.5,
  names.arg = rep("", ncol(primary_vectors)),
  ylab = "Female-referenced component (SD units)",
  ylim = component_limits + c(-1, 1) * component_padding,
  cex.names = 0.78,
  cex.axis = 0.85,
  cex.lab = 0.95,
  las = 1,
  axes = FALSE
)
axis(2, las = 1, cex.axis = 0.85)
axis(
  1,
  at = colMeans(bar_positions),
  labels = unname(trait_labels[colnames(primary_vectors)]),
  tick = FALSE,
  cex.axis = 0.70,
  line = -0.35,
  gap.axis = -1
)
abline(h = 0, lwd = 0.9)
box(bty = "l")
legend(
  "topleft",
  legend = rownames(primary_vectors),
  fill = group_colours[rownames(primary_vectors)],
  border = "grey20",
  bty = "n",
  horiz = TRUE,
  cex = 0.82,
  inset = c(0, -0.01)
)
mtext("B", side = 3, adj = 0, line = 0.5, font = 2, cex = 1.05)

par(old_par)
invisible(dev.off())

model_object <- list(
  analysis_label = "female-referenced multivariate contrast vectors",
  primary_scenario = primary_scenario,
  scenario_metadata = scenario_metadata,
  data_file = normalizePath(data_file),
  reference_pronotum = reference_pronotum,
  permutations = permutations,
  bootstrap_reps = bootstrap_reps,
  complete_case_IDs = dat$ID,
  omitted_IDs = dat_raw$ID[!complete_rows],
  response_residual_SD = response_residual_SD,
  design_matrices = list(
    size = X_size,
    additive = X_additive,
    interaction = X_interaction,
    shared_anchor = X_shared_anchor,
    rank_one_nuisance = Z_rank_one,
    rank_one_contrasts = W_rank_one
  ),
  results = results,
  output_tables = output_tables,
  contrast_axis_table = contrast_axis_table,
  figure_file = normalizePath(figure_file),
  session_info = sessionInfo()
)

model_file <- file.path(
  models_dir,
  "female_reference_multivariate_reanalysis.rds"
)
saveRDS(model_object, model_file)

session_file <- file.path(
  output_root,
  "female_reference_multivariate_session_info.txt"
)
session_lines <- c(
  "Analysis label: female-referenced multivariate contrast vectors",
  paste("Input:", normalizePath(data_file)),
  paste("Primary scenario:", primary_scenario),
  paste(
    "Primary trait set:",
    scenario_metadata$trait_set[scenario_metadata$primary_analysis]
  ),
  "Primary-analysis rationale: head is excluded because head length contributes to male morph assignment.",
  paste("Complete cases:", nrow(dat), "of", nrow(dat_raw)),
  paste("Reference pronotum (mm):", reference_pronotum),
  paste("Residual randomizations:", permutations),
  paste("Group-stratified bootstrap replicates:", bootstrap_reps),
  paste("Contrast-vector figure:", normalizePath(figure_file)),
  paste(
    "Optional packages installed:",
    paste0(
      c("RRPP", "geomorph", "evolqg"),
      "=",
      vapply(
        c("RRPP", "geomorph", "evolqg"),
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
      ),
      collapse = ", "
    )
  ),
  "",
  capture.output(sessionInfo())
)
writeLines(session_lines, session_file)

cat("Female-referenced multivariate reanalysis complete.\n")
cat("Complete cases:", nrow(dat), "\n")
cat("Tables written to:", normalizePath(tables_dir), "\n")
cat("Model object:", normalizePath(model_file), "\n")
cat("Session information:", normalizePath(session_file), "\n\n")
cat("Figure:", normalizePath(figure_file), "\n\n")
cat("Primary standardized five-trait results (head excluded):\n")
print(subset(model_tests, scenario == primary_scenario))
print(subset(vector_magnitudes, scenario == primary_scenario))
print(subset(vector_angles, scenario == primary_scenario))
print(subset(magnitude_differences, scenario == primary_scenario))
print(subset(rank_one_tests, scenario == primary_scenario))
cat("\nSensitivity scenarios:\n")
print(subset(scenario_metadata, !primary_analysis))
