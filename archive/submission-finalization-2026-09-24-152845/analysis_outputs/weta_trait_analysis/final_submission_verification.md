# Final submission verification

Completed 24 September 2026. The final specification is recorded in `FINAL_ANALYTICAL_APPROACH.md`.

## Statistical verification

- The final secondary instar analysis was run successfully from the unchanged analysis input (MD5 `2d67b2d772f6f28b3cdf77c87b35b5fc`). It allows separate means, body-size slopes and residual variances in all four groups for each of five traits.
- `tests/test_statistical_reporting.R` passed all seven checks. An independent GLS fit reproduces the final group estimates. Independent calculations verify prediction variances, the three contrast definitions, their covariance propagation, Welch–Satterthwaite degrees of freedom, confidence intervals, percentage transformations and Holm adjustments. Recorded specimen IDs and morph labels agree with the input data.
- The same reporting checks reproduce the retained primary contrast tables, slope and precision reporting, original-model sensitivity, classification sensitivity and multiplicity audit.
- `tests/test_workflow.R` and `tests/test_qmd_tables.R` passed.
- All pre-existing numerical tables, saved model files and analysis figures included in the finalization workspace remained byte-identical to their starting versions. The primary analyses were retained and their reported calculations checked; the full upstream simulation pipeline was not rerun during this finalization.
- The final hindleg direct comparison has Holm-adjusted p = 0.158 and a pointwise log-scale 95% interval spanning zero. Manuscript Table 3, supplementary Table S20 and the prose consistently describe its ordering as uncertain.

## Document verification

- The revised manuscript, supplement and title page rendered successfully from their Quarto sources. The rendering helper preserves all three outputs in the submission folder.
- Every page was visually reviewed: manuscript 32 pages, supplement 33 pages, title page 2 pages. Unchanged pages were verified against previously inspected page images. Tables, figures, captions, mathematical notation and references were readable, with captions kept with their associated content.
- The main manuscript contains three figures and the supplement four; all seven have descriptive alternative text. Figure 2 labels are visible.
- Structural checks found no unresolved cross-reference markers, unevaluated inline R, missing-value p-values or tracked changes. Creator metadata is empty; author details remain on the separate title page.
- The final manuscript retains its 222-word abstract. Results and discussion use the final secondary model, while the original constrained model appears only as sensitivity evidence.

## Provenance and completion

Primary outputs retain their original R 4.4.2 provenance. The final secondary analysis and reporting checks used R 4.3.2, documented in the session record. No measurements or supplied morph classifications were changed. `final_submission_manifest.csv` records hashes for the final sources, scripts, results and Word files.

The stated completion criteria have passed. The analytical work is complete for submission; additional exploratory model comparisons are not required.
