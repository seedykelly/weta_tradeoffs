# Final analytical approach

Specification settled on 24 September 2026 for the current manuscript and data.

## Decisions

1. **Retain the primary frequentist analyses.** The repeated-leg model retains group-specific body-size slopes, separate variances by leg pair and an unstructured within-individual correlation matrix. It supplies the tests of leg allocation. The existing one-stage raw-head models supply weapon–trait associations. Sensory group comparisons retain the established common-slope models supported by their diagnostics. The multivariate group-effect, individual-level PCA and matrix analyses retain their respective trait sets and are identified as exploratory.
2. **Use the fully flexible group model for the secondary instar comparison.** Each trait has separate group intercepts, body-size slopes and residual variances. No linear sequence of male morph means is imposed. This is the M3 specification already examined in the sensitivity analysis. Separate group regressions reproduce its GLS estimates and permit an analytic Welch–Satterthwaite calculation.
3. **Define the biological comparisons explicitly.** At pronotum length 7.2 mm, the average male change per instar is half the tenth-minus-eighth log difference. The matched-sex contrast is tenth-instar males minus females. Their direct difference accounts for the shared tenth-instar estimate. These are comparisons of adult groups, not longitudinal developmental effects or causal shares of weapon investment.
4. **Retain uncertainty and multiplicity rules.** For the secondary comparison, use pointwise 95% intervals and Holm adjustment across five traits separately for the average male change, matched-sex contrast and direct difference. Retain the primary contrast families and their reported interval adjustments. Separate trait tests do not make a joint test across leg pairs; the main repeated-leg model performs the allocation tests.
5. **Report the assumption sensitivity.** Supplementary Table S20 compares the original constrained model with the final flexible model. The final specification is used consistently across all five traits, irrespective of significance. Allowing separate group means costs some precision but avoids an unnecessary developmental constraint on categorical morphs. This is a design-based decision, not selection of the most favourable p-value or the lowest AIC separately for each trait.
6. **Retain supplied morph labels and the original data.** Historical threshold exclusions remain sensitivity analyses. The Bayesian analyses remain separate material for possible reviewer queries.

## Final interpretation

The main findings remain coordinated enlargement of locomotor and sensory traits across male morphs, widespread positive within-group associations, and greater anterior-leg enlargement in tenth-instar males relative to females. The secondary comparison supports a larger matched-sex difference than average male instar change for forelegs, midlegs, tympanic size and eyes. For hindlegs, the estimates are approximately 5.25% and 3.93%, respectively, but their ordering is uncertain (Holm p = 0.158). The manuscript no longer claims a demonstrated larger hindleg instar effect.

The fully flexible secondary model need not produce identical sensory sex contrasts to the primary common-slope sensory models. Their assumptions and purposes are stated separately; the primary estimates are not replaced selectively.

## Reproduction and provenance

Run `Rscript scripts/run_all.R` to reproduce the maintained pipeline. To update the final stage from the existing verified upstream outputs, run `Rscript scripts/instar_matched_precision_analysis.R`. The final main-text table reads `analysis_outputs/weta_trait_analysis/tables/instar_matched_contrasts.csv`; the supplement reads `instar_contrast_sensitivity.csv`.

The original constrained per-trait output remains for sensitivity and historical reproducibility. The constrained joint instar model is retired. Primary repeated-leg models are unchanged.

The upstream analyses used R 4.4.2. Final secondary contrasts and document checks used R 4.3.2. Software versions and input hashes are recorded with the results; the input checksum is 2d67b2d772f6f28b3cdf77c87b35b5fc. No raw measurements or recorded morph assignments were changed.

## Completion criteria

- Final model estimates independently agree with a GLS implementation.
- Contrast definitions, variance propagation, degrees of freedom, transformations and multiplicity adjustments pass numerical checks.
- Retained primary results remain unchanged and use the same input data.
- Manuscript and supplement cite the final tables and express the hindleg uncertainty consistently.
- The maintained documents render successfully, with readable tables, figures, references and intact cross-references.

After these checks pass, the analytical work is complete for submission. Further model exploration is not part of this finalization; reopening an analysis would require a concrete data/code error, new data or a specific reviewer request.
