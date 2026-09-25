# Bayesian comparison with brms

The standalone script is `scripts/bayesian_hypotheses_brms.R`. It preserves the
morph assignments supplied in `data/trait_data.csv` and writes only to a separate
`analysis_outputs/weta_brms/` directory. It does not change the manuscript,
frequentist outputs, or morph classifications.

Run these commands from the project root:

```sh
# Validate data, formulas, priors and generated Stan code without fitting.
Rscript scripts/bayesian_hypotheses_brms.R --mode=check

# Short end-to-end software check; do not interpret its posterior summaries.
Rscript scripts/bayesian_hypotheses_brms.R --mode=smoke

# Main run: all three models, four chains, 4,000 iterations (2,000 warmup),
# baseline and twice-as-wide priors, and prior/posterior predictive checks.
Rscript scripts/bayesian_hypotheses_brms.R --mode=fit
```

Models can be run separately, for example `--models=morphology`. Runs are cached;
changes in model data, priors or requested sampling settings trigger refitting.
Compilation and sampling can take substantial time. Dependencies are brms,
cmdstanr, posterior, and a working CmdStan installation; the script does not install
or rebuild them. On macOS it selects an SDK belonging to the selected compiler,
sets that choice only for the current R process, and disables stale precompiled
headers. Override with `--sdk=/path/to/MacOSX.sdk` or `--sdk=default` if necessary.

## Questions covered

- **Searching and detection:** posterior probabilities that eighth-instar males
  exceed tenth-instar males in relative leg and eye dimensions, plus a joint
  probability across those traits. Tympanic size is included as a separate trait.
- **Positive integration:** the reverse morph-size contrasts; partial head–trait
  slopes in each group; and contrasts testing whether tenth-instar slopes exceed
  eighth- or ninth-instar slopes. Positive associations do not establish functional
  compensation or exclude resource-allocation costs.
- **Coordinated divergence:** posterior uncertainty in the design-weighted share
  represented by the dominant axis, using the same weighting and trait scaling
  as the manuscript. Probabilities that this share exceeds 80% or 90% are labelled
  descriptive; they are not tests of exactly rank one.
- **Additional contrasts:** anterior–posterior leg allocation, ninth-instar
  midpoint equivalence within the illustrative +/-2.5% region, and descriptive
  per-instar versus instar-matched sex contrasts.

The models estimate a shared unstructured residual covariance across the five
traits. They do not test equality of covariance matrices among morphs. No morph
is assumed exchangeable with females, and no automatic partial pooling across
morphs is introduced.

## Comparing approaches fairly

All Bayesian models and their frequentist companions use the same 311 complete
specimens and the same mean formulas. This differs from some manuscript models,
which use trait-specific samples and common body-size slopes for selected sensory
analyses. Consequently, the comparison isolates the inferential framework more
closely, but its frequentist estimates are not all identical to the published
estimates. The matched morphology model's dominant-axis share reproduces 88.8%.

`bayesian_frequentist_comparison.csv` and the matching PDF show estimates and
95% intervals on a common coefficient scale. Bayesian intervals are marginal,
equal-tailed credible intervals. Frequentist intervals are pointwise t intervals;
Holm-adjusted p-values are supplied in named contrast tables. Neither these
intervals nor their meanings are interchangeable. Posterior sign probabilities
are not converted p-values.

Optional out-of-sample comparison:

```sh
Rscript scripts/bayesian_hypotheses_brms.R --mode=fit --cv-folds=5
```

This adds five-fold cross-validation with entire specimens held out and identical
folds for both approaches. It compares mean predictions using log-scale RMSE,
separately by trait and model, using Bayesian prior scale 1. It will fit additional
models and take longer. The scores are descriptive, not formal tests of the
biological hypotheses or evidence that one statistical philosophy is superior.

## Reading the outputs

Start with `READ_ME_FIRST.md`, `sample_counts.csv`, and
`sampler_diagnostics.csv`. Inspect the predictive-check PDFs and compare results
between prior scales 1 and 2. A fit passes the automated diagnostic screen only
when Rhat <= 1.01, bulk/tail ESS >= 400, no divergent transitions or maximum-depth
hits occur, and every chain has E-BFMI >= 0.3. Passing this screen is necessary but
not sufficient for an adequate model; predictive checks still need inspection.

`*_prior_only_*` files concern simulations before conditioning on the outcomes.
`*_morph_contrasts.csv`, `*_allocation_contrasts.csv`,
`*_head_trait_slopes.csv`, `*_head_slope_differences.csv`,
`*_dominant_axis.csv`, and `*_instar_sex.csv` address the questions above.
`*_joint_hypothesis_probabilities.csv` uses the joint posterior draws, retaining
trait dependence rather than multiplying marginal probabilities.

The +/-2.5% midpoint region is an illustrative sensitivity margin, not an
established biological threshold. Head-slope equivalence is omitted by default.
If a meaningful raw log-log slope margin is justified, supply `--slope-rope=...`;
this is not the standardized +/-0.20 region used in the existing supplement.
Do not choose a prior, region or statistical approach because it produces a
preferred biological result.

Pure calculation checks are available with:

```sh
Rscript tests/test_bayesian_hypotheses.R
```

## Software validation

Short sampling runs completed for all three model families, including contrast
tables, diagnostics, prior/posterior predictive checks and comparison figures.
The morphology model was also checked with twice-as-wide priors. Saved Stan data
confirmed that the intended prior scales and prior-only settings were used.
A two-fold instar-model check confirmed that every specimen was predicted once
per trait and excluded from its corresponding training fit. Calculation checks
reproduced the manuscript's 88.8% dominant-axis share and verified contrast
directions and the conversion from linearised ear size to area.

These were software checks, not final analyses. The short runs had no divergent
transitions but did not all meet the final Rhat/ESS requirements. Run `--mode=fit`
and inspect its diagnostics before drawing biological conclusions or deciding
whether the frameworks differ materially.

Technical references: [brms multivariate models](https://paulbuerkner.com/brms/articles/brms_multivariate.html),
[prior specifications](https://paulbuerkner.com/brms/reference/set_prior.html),
[Stan predictive checks](https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html).
