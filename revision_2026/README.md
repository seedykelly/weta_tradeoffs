# Synchronized 2026 revision

This directory is the canonical revised manuscript bundle. The manuscript and supplement read only the analysis outputs stored here, so they cannot silently pick up older project-level results.

## Render the documents

From this directory:

```sh
quarto render kelly_2026_weta_morph_investment_revised.qmd --to docx
quarto render kelly_2026_weta_morph_investment_supplement.qmd --to docx
```

## Re-run the focal trait analyses

The bundled `data/trait_data.csv` contains the 313 analysed specimens. From this directory:

```sh
Rscript scripts/run_all.R
```

This regenerates `analysis_outputs/weta_trait_analysis`. The archived outputs were produced with 9,999 residual randomizations and 9,999 bootstrap replicates.

The repeatability workflow uses the separate raw remeasurement archive. Run it with the parent research-project directory as its first argument; its maintained default already points to the current local project.

## Analysis hierarchy

The main inference is based on the female-referenced multivariate dimensionality tests, the unstructured-covariance repeated-leg model, the one-stage raw-head weapon–leg model, body-size-adjusted sensory models, covariance/correlation-matrix comparisons, and defined equivalence contrasts. Residual-index models, femur–tibia decompositions, and broad group-specific slope searches remain in the analysis archive and are not used as headline evidence.

## Author details still required

Before submission, replace the visible author-action notes for:

- the institutional animal-care approval/protocol or formal exemption, including confirmation that freezing at −20 °C was approved;
- whether measurements were made blind to group; and
- the NSERC grant number.

The 2016 collection year and Department of Conservation acknowledgement wording are already incorporated.
