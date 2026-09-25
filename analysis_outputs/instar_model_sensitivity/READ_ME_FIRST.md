# Instar-model sensitivity analysis

This is an exploratory check prompted by predictive diagnostics. Recorded morph labels and the existing manuscript outputs are preserved.
The original per-trait sample sizes, estimates, p-values and residual degrees of freedom were reproduced against the archived decomposition.

## Models
M0_original: log trait ~ log pronotum + male terminal instar + sex; one residual SD.
M1_group_variance: same mean model; separate residual SDs for female, eighth, ninth and tenth groups.
M2_group_slopes: also permits separate body-size slopes; male means at 7.2 mm still follow equal log-scale steps.
M3_free_groups: separate group means, body-size slopes and residual SDs; equivalent to independent within-group regressions.

## Comparisons
Each candidate uses the same observations within each trait. ML likelihoods/AIC/AICc and nested likelihood-ratio tests compare specifications. REML estimates and covariance matrices provide inference.
For each model, average_male_step = (tenth - eighth)/2 on the log scale; matched_sex = tenth - female; step_minus_sex is their direct difference. All are evaluated at 7.2 mm pronotum length.
In M0-M2 the average male step is also the fitted constant per-instar change. M3 imposes no such constant gradient; adjacent ninth-minus-eighth and tenth-minus-ninth contrasts and the ninth-minus-midpoint departure are reported separately.
The sign of step_minus_sex compares two descriptive contrasts, not causal developmental and weapon contributions. Positive means the average male step is larger; negative means the matched-sex difference is larger.
All ear results are for square-root-transformed area (linearised size). Exponentiating step_minus_sex yields a ratio of multiplicative effects, not a subtraction of percentage changes.
Confidence intervals are pointwise 95% t intervals. M0 uses the original residual df; M1-M2 use emmeans Satterthwaite df including uncertainty in the variance parameters. For M3, the equivalent independent group regressions allow an analytic Welch-Satterthwaite calculation. Holm p-values control each five-trait contrast family separately within each model. Model-extension tests have separate five-trait Holm families.
All candidates are retained; the analysis does not choose whichever model favours a biological hypothesis. This sensitivity does not re-estimate the separate joint repeated-leg model.

## Diagnostics
Conditional simulation checks use 2000 replicates at observed covariates with fitted coefficients and residual SDs held fixed.
These are descriptive plug-in checks, not posterior predictive intervals, cross-validation or formal model-adequacy tests. Parameter uncertainty and multiplicity are not included. Free group intercepts reproduce group means by construction; passing those mean checks is therefore not independent validation.
Inspect residual_diagnostics.pdf for residual shapes and fitted-value trends as well as conditional_simulation_checks.csv for group means/spreads.

## Files
contrast_comparison.csv contains all effects, intervals, df and Holm p-values; contrast_comparison.png/pdf displays the direct step-minus-sex contrast.
model_comparison_ML.csv and nested_model_tests_ML.csv identify which assumptions affect fit. adjusted_means_and_slopes.csv records the group predictions and body-size slopes.
sample_manifest.csv, run_manifest.csv, session_info.txt and fitted_models.rds retain the analysis provenance.

Method references: https://rvlenth.github.io/emmeans/articles/models.html ; https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/anova.gls.html
