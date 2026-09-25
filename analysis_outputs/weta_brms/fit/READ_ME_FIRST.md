# Bayesian hypothesis analysis

Run mode: fit
Complete specimens: 311
Morph assignments are read directly from the supplied CSV; no reclassification is performed.

Models: morphology = group-specific body-size slopes for all five traits; integration adds group-specific raw-head slopes; instar = body size + male terminal instar + sex.
Each model has a shared unstructured residual covariance. Individual slopes are fixed effects with symmetric priors, not a hierarchical model pooling male morphs.
Matched frequentist models use identical complete cases and mean formulas. Their estimates can differ from the published trait-specific samples and selected common-slope sensory models.
Bayesian 95% intervals are marginal equal-tailed credible intervals. Frequentist 95% intervals are pointwise t intervals; named Holm families adjust p-values only. Do not equate their coverage or posterior tail areas.
Positive/negative probabilities are conditional on the likelihood and priors, not frequentist p-values. Joint probabilities are computed draw by draw and do not assume traits are independent.
The +/-2.5% midpoint region is illustrative, matching the earlier sensitivity analysis; justify a biologically meaningful margin before making a substantive equivalence claim.
Head-slope practical equivalence is omitted unless --slope-rope is supplied. That optional halfwidth applies to raw log-log slopes, NOT the earlier standardized +/-0.20 region.
Dominant-axis shares use the manuscript design weighting and pooled size-residual SDs. Probabilities above 80% and 90% are descriptive summaries; neither is a test of exactly rank one.
Group-specific covariance equality, detailed segment models and measurement-error models are outside this script.
Prior scale 2 doubles all coefficient and residual-SD prior scales; LKJ(2) is held constant. Both fits must be reported rather than choosing a favourable prior.
Prior predictive draws are simulated without conditioning on outcomes; inspect them alongside posterior predictive checks of group means and spreads.
All final fits require satisfactory Rhat, bulk/tail ESS, divergences, tree depth and energy diagnostics. Summary tables are not automatically interpretable when diagnostics fail.
Smoke-mode tables are software checks and must not be reported as scientific results.
Optional cross-validation holds out entire specimens, stratifies folds by morph, and compares mean predictions on the centred log scale. It uses prior scale 1 only; predictive accuracy does not establish a causal hypothesis.
No automated declaration of a winning framework or supported hypothesis is made. Review effect sizes, uncertainty, prior sensitivity, diagnostics and substantive relevance together.

Sampling finished: 2026-09-24 12:46:55 UTC
All diagnostic screens passed: FALSE
