# Statistical reporting verification

The checks below passed against the maintained analysis input and result tables.

- Eight reported contrast tables reproduce Holm tests, Bonferroni intervals and response-scale ratios.
- All 12 leg-precision rows use the main contrast estimates, SEs and Satterthwaite degrees of freedom.
- All 12 head-leg estimates are positive; 11 pointwise intervals exclude zero. Precision intervals match the slope table.
- Five instar/sex regressions reproduce coefficients, intervals, sample sizes, direct coefficient comparisons, percentage changes and three Holm families.
- Classification sensitivity and the new figure use the same design-weighted group-effect definition as the main analysis.
- The ten-test audit reproduces Holm and Benjamini–Hochberg adjustments; the matrix test is labelled exploratory.

Excluding ten category-disagreement males changed the dominant-axis share from 88.8% to 88.6%; the residual rank-one test remained p = 0.0001 (9999 randomizations).

Holm-adjusted p-values and Bonferroni-adjusted simultaneous intervals are explicitly distinguished in the documents. Precision intervals and minimum detectable effects refer to single comparisons without multiplicity adjustment.

Leg-precision estimates now reuse the main contrast/slope tables rather than substituting residual observation degrees of freedom. Main-model estimates, Holm p-values and biological conclusions are unchanged.

Morph categories are the externally derived categories supplied with the dataset; no mixture model was fitted anew. The exclusion analysis is a sensitivity to historical category definitions, not a reclassification of the main sample.

Scope: targeted verification of reported contrasts, transformations, precision calculations, instar/sex regressions and the classification sensitivity. This is not an independent rerun of every model or a verification of raw measurements.
