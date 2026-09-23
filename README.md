# Wētā trait-compensation project

This directory is the canonical manuscript and analysis project.

## Project layout

- `manuscript.qmd` — sole maintained main manuscript, anonymized for review
- `supplement.qmd` — sole maintained supplement, anonymized for review
- `title-page.qmd` — author information, declarations and acknowledgements; submit separately
- `submission-review/` — checked Word snapshots from the maintained QMD sources (23 September 2026); regenerate after source edits
- `kelly_2026_weta_morph_investment*.qmd` — superseded historical drafts; do not edit for submission or copy over the maintained pair
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
`analysis_outputs/weta_trait_analysis`. Default paths are resolved from the
script's project directory, so an absolute script path also works from another
working directory. Explicit relative arguments remain relative to the caller.

Optional arguments are data file, output directory, permutation count,
bootstrap count, and rarefaction count, in that order. Defaults are 9,999
permutations, 9,999 bootstraps, and 2,000 rarefactions. Permutation and bootstrap
counts must be integers of at least 999; rarefactions must be a positive integer.
For a reduced-cost verification run, use a separate output directory:

```sh
Rscript scripts/run_all.R data/trait_data.csv /tmp/weta-check 999 999 25
```

The permutation count controls both the female-reference randomizations and
the reviewer's PCA and covariance/correlation permutation analyses. The
rarefaction count controls covariance-matrix subsampling. Random skewers remain
fixed at 10,000 draws. Actual counts are saved in `analysis_run_manifest.csv`,
`tables/reviewer_resampling_settings.csv`, and the relevant result tables.
Reduced-count runs are for testing, not final reporting. Fixed seeds also cover
the simulated calculations used by approximate Satterthwaite degrees of freedom.

The repeatability workflow uses the
separate raw remeasurement archive:

```sh
Rscript scripts/repeatability_analysis.R
```

Each `data_cleaning_*.R` script accepts optional raw and clean measurement
directories. `join_cleaned_files.R` accepts the clean measurement directory and
an optional output filename. It rejects duplicate or missing IDs and missing or
invalid sex metadata before joining. Males without a measured head length retain
a missing morph assignment. Its default output remains `data/clean/trait_data.csv`;
it does not overwrite the analysis input `data/trait_data.csv`.

Run the regression checks with:

```sh
Rscript tests/test_workflow.R
Rscript tests/test_qmd_tables.R
```

These checks use temporary copies, including a relocated project, and do not
overwrite the original data or analysis results. Full-pipeline reproducibility
should additionally be checked by comparing two runs with identical counts.

### Analysis messages

The current package versions still emit advisory messages. They are not hidden:

- `contrasts dropped` is produced while R rebuilds factor levels for an
  `emmeans` reference grid; the GLS method then explicitly reapplies the fitted
  model's contrast matrices. Sum coding and the fitted models are unchanged.
- Approximate Satterthwaite notices describe the existing GLS inference method;
  fixed seeds make its simulated calculations repeatable on the same software.
- Interaction and VIF notices require interpretation of conditional effects;
  they are not a reason to remove the prespecified interactions.
- Packages built under R 4.4.3 are currently being used with R 4.4.2. Some
  diagnostic plotting code also emits `NULL`-conversion notices. These are
  dependency/environment issues, not permission to suppress statistical warnings.

## Render the documents

```sh
quarto render manuscript.qmd --to docx
quarto render supplement.qmd --to docx
quarto render title-page.qmd --to docx
```

Documents read the analysis outputs beside their own QMD file; they no longer
select a different project copy based on a hard-coded home/work computer path.
To refresh the Word snapshots in `submission-review/`, add
`--output-dir submission-review` to each rendering command above. Edit the QMD
sources, not those generated Word copies.

## Submission sources and editorial checks

Only `manuscript.qmd` and `supplement.qmd` are maintained reviewer documents.
The older named pair is explicitly superseded; historical files are not sources
for future updates. Author identity, repository links, contributions, funding,
acknowledgements and declarations belong in `title-page.qmd`, not in the reviewer files.
The Methods retain an anonymized ethics statement.

The abstract must remain below 250 words, with 4–10 keywords, and the main text
below 7,500 words under the current JEB Research Article requirements:
https://academic.oup.com/jeb/pages/author-guidelines.
Count rendered prose, not R source or inline expressions. Preserve the distinction
between the primary five-trait matrix results and the head-inclusive sensitivity
when editing the abstract, Results, Discussion or conclusions. Nonsignificant
tests are not evidence of equivalence or absence of allocation costs.

Before submission, confirm the corresponding author's full postal address and
any applicable grant or permit identifiers on the title page. None were invented
during editing. Retain the AI-assistance disclosure and also disclose this assistance
in any accompanying cover letter, as requested by the journal. Check the final
Word metadata as well as visible text for anonymity.
