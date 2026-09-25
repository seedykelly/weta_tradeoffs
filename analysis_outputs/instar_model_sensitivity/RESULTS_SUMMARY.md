# Flexible instar-model sensitivity: results

The analysis reproduced the original frequentist decomposition before relaxing its assumptions. The most consequential change is allowing the relationship with body size to differ among the four groups.

## Effect estimates

All comparisons below are adjusted to a pronotum length of 7.2 mm, which is inside the observed range of every group. The male step is the average eighth-to-tenth change per instar. The matched-sex contrast compares tenth-instar males with females. Percentages refer to length or linearised tympanic size.

| Trait | Original male step | Original matched sex | M2 male step | M2 matched sex | M2 direct-difference Holm p | M3 direct-difference Holm p |
|---|---:|---:|---:|---:|---:|---:|
| Foreleg | 5.05% | 15.92% | 5.09% | 16.99% | 1.91e-39 | 2.03e-29 |
| Midleg | 4.30% | 10.91% | 4.62% | 12.30% | 2.10e-22 | 6.11e-15 |
| Hindleg | 4.55% | 2.22% | 5.10% | 4.28% | 0.300 | 0.158 |
| Tympanic linear size | 3.55% | 7.91% | 3.58% | 8.16% | 1.16e-05 | 8.76e-04 |
| Eye | 6.86% | 24.64% | 7.64% | 26.29% | 1.41e-36 | 4.33e-28 |

M2 permits separate body-size slopes and residual SDs, while retaining equal male instar steps at the reference body size. M3 additionally frees the ninth-instar mean; its two adjacent steps need not be equal. Both models, and the intermediate variance-only model, are reported in full rather than selecting by hypothesis-test significance.

## What changes for hindlegs?

- M0_original: direct difference 0.0225 log units; pointwise 95% CI 0.0093 to 0.0357; Holm p = 9.03e-04.
- M1_group_variance: direct difference 0.0230 log units; pointwise 95% CI 0.0097 to 0.0362; Holm p = 7.43e-04.
- M2_group_slopes: direct difference 0.0078 log units; pointwise 95% CI -0.0070 to 0.0225; Holm p = 0.300.
- M3_free_groups: direct difference 0.0127 log units; pointwise 95% CI -0.0050 to 0.0303; Holm p = 0.158.

The estimated hindleg instar step remains larger, but the more flexible models do not clearly establish that it exceeds the matched-sex contrast. This reflects both a smaller estimated difference and uncertainty around it. It is not evidence that the two effects are equivalent.

The original claim that the hindleg instar change exceeds the sex contrast should therefore be qualified as dependent on the common body-size slope assumption.

For forelegs, midlegs, tympanic size and eyes, the matched-sex contrast remains larger than the average male instar step under every fitted specification, including after Holm adjustment.

## Which changes improve fit?

| Trait | Lowest AICc among candidates | Original minus best AICc | M2 minus best AICc | M3 minus best AICc |
|---|---|---:|---:|---:|
| Foreleg | M2_group_slopes | 11.55 | 0.00 | 2.15 |
| Midleg | M2_group_slopes | 11.15 | 0.00 | 1.89 |
| Hindleg | M2_group_slopes | 15.28 | 0.00 | 1.14 |
| Tympanic linear size | M0_original | 0.00 | 9.56 | 11.28 |
| Eye | M0_original | 0.00 | 6.48 | 7.60 |

AICc values compare models fitted by maximum likelihood to the same observations within each trait. The lowest value is a relative comparison among these candidates, not proof that a model is true.

| Trait | Add residual SDs: Holm p | Add size slopes: Holm p | Free group means: Holm p |
|---|---:|---:|---:|
| Foreleg | 0.765 | 8.67e-04 | 1.000 |
| Midleg | 0.395 | 0.002 | 1.000 |
| Hindleg | 1.000 | 2.10e-05 | 1.000 |
| Tympanic linear size | 1.000 | 1.000 | 1.000 |
| Eye | 0.765 | 1.000 | 1.000 |

Each column tests a successive extension of the previous model. The evidence for different leg slopes is stronger than the evidence for different residual variances or for freeing the ninth-instar mean. The original specification has the lowest AICc for tympanic size and eyes among the candidates examined.

## Diagnostics and limits

The conditional simulations check whether observed group means and SDs fall within 95% bands from simulations at the observed body sizes. These use fitted parameters, omit their uncertainty, and are descriptive checks rather than formal adequacy tests.

| Model | Group means outside bands | Group SDs outside bands | Checks per column |
|---|---:|---:|---:|
| M0_original | 3 | 3 | 20 |
| M1_group_variance | 3 | 3 | 20 |
| M2_group_slopes | 0 | 0 | 20 |
| M3_free_groups | 0 | 0 | 20 |

The group-mean checks in M3 are satisfied by construction because its group intercepts are free. Residual plots are supplied for the original, group-slope and free-group models; passing the group summaries does not guarantee every distributional assumption.

This analysis is exploratory because it was prompted by diagnostics. Confidence intervals are pointwise; Holm p-values control five traits separately for each contrast/model, not the entire process of examining several models. The contrasts describe associations and do not isolate causal effects of development, sex or weapon investment.

With separate body-size slopes, the contrasts are specific to the stated 7.2 mm reference. This analysis concerns the secondary instar-sex decomposition and does not re-run or replace the main multivariate analyses or the joint repeated-leg model.

## Re-running

From R/RStudio, source the project script `scripts/instar_model_sensitivity.R`. It locates its input and output folders from the script location, reproduces the original baseline, and runs all four models for all five traits. Required R packages: nlme and emmeans. Default diagnostic simulations: 2,000 per fitted model, with a recorded random seed.

Detailed estimates and all confidence intervals are in `contrast_comparison.csv`; model definitions, methods and output descriptions are in `READ_ME_FIRST.md`. Fitted models and software versions are saved alongside them.
