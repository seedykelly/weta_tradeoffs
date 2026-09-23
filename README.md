# Wētā trait-compensation project

This directory is the canonical manuscript and analysis project.

## Project layout

- `kelly_2026_weta_morph_investment.qmd` — main manuscript
- `kelly_2026_weta_morph_investment_supplement.qmd` — supplementary material
- `scripts/` — data-cleaning, repeatability, and manuscript-analysis scripts
- `data/` — raw, cleaned, and analysis-ready data
- `analysis_outputs/` — outputs used by both manuscripts
- `figures/` — Figure 1 and its editable source material
- `archive/` — historical source files retained for reference

## Re-run the analyses

From the project root:

```sh
Rscript scripts/run_all.R
```

This uses `data/trait_data.csv` and regenerates
`analysis_outputs/weta_trait_analysis`. The repeatability workflow uses the
separate raw remeasurement archive:

```sh
Rscript scripts/repeatability_analysis.R
```

## Render the documents

```sh
quarto render kelly_2026_weta_morph_investment.qmd --to docx
quarto render kelly_2026_weta_morph_investment_supplement.qmd --to docx
```
