# README for trait_data.csv

## Dataset overview

`trait_data.csv` contains morphological measurements for 313 adult Wellington tree wētā (*Hemideina crassidens*). It is the analysis input for the manuscript **A single dominant axis of locomotor and sensory divergence across alternative mating strategies in Wellington tree wētā**.

The dataset supports comparisons of locomotor, sensory and weapon-associated morphology among females and three adult male morphotypes. Each row represents one individual. The file contains specimen identifiers, sex, morph assignments based on the 2016 head-length cut-offs, and 11 morphological measurements.

**Data creator and contact:** Clint D. Kelly, Département des Sciences biologiques, Université du Québec à Montréal, Montréal, Canada. Email: kelly.clint@uqam.ca. ORCID: https://orcid.org/0000-0002-0693-7211.

## Collection and measurement methods

Adults were collected haphazardly at night from vegetation on Te Pākeka/Maud Island, New Zealand (41°02′ S, 173°54′ E), during March and April 2016. Specimens were euthanized by freezing at −20 °C. The head, pronotum and both legs of each pair were dissected and retained under the individual's identifier.

Dissected structures were photographed with a Canon EOS 5D Mark IV camera and Canon EF 100 mm macro lens. A ruler was photographed in the same plane as each structure to calibrate image measurements. Measurements were made in Fiji/ImageJ using consistent anatomical landmarks. Lengths were measured with the straight-line tool; tympanic membrane area was measured by tracing the membrane perimeter. Measurement landmarks are illustrated in Figure 1 of the associated manuscript.

Linear measurements are in millimetres (mm). Tympanic membrane area (`ear`) is in square millimetres (mm²). Values are stored on their measured scales, before logarithmic or square-root transformations.

## File structure and sample composition

- **Filename:** `trait_data.csv`; located at `data/trait_data.csv` in the analysis project.
- **Format:** UTF-8 comma-separated values, with a header row and a decimal point as the decimal separator.
- **Dimensions:** 313 data rows and 14 columns.
- **Observation unit:** one adult individual per row.
- **Unique key:** `ID`; all 313 identifiers are unique and nonmissing.
- **Missing-value code:** an empty CSV field.

| Sex | `morph` value | Description | Individuals |
| --- | --- | --- | ---: |
| `f` | `female` | Adult female | 118 |
| `m` | `eighth` | Adult male assigned to the eighth-instar morphotype | 78 |
| `m` | `ninth` | Adult male assigned to the ninth-instar morphotype | 79 |
| `m` | `tenth` | Adult male assigned to the tenth-instar morphotype | 38 |
| **Total** | | | **313** |

The male morph labels describe inferred terminal developmental morphotypes. All specimens in this file are adults; these rows are not observations of individuals followed through successive juvenile instars.

## Variable dictionary

Columns are listed in their order in the CSV. “Bilateral mean” means the arithmetic mean of available left- and right-side measurements, as explained below.

| Variable | Type | Units or permitted values | Definition |
| --- | --- | --- | --- |
| `ID` | Text | Uppercase alphabetic identifier | Unique specimen identifier linking measurements and metadata. Preserve as text. |
| `sex` | Categorical text | `f`, `m` | Recorded sex: female or male. |
| `morph` | Categorical text | `female`, `eighth`, `ninth`, `tenth` | Female reference group or male morphotype assigned using the 2016 head-length cut-offs. |
| `pronotum` | Numeric | mm | Pronotum length, used as the measure of structural body size. |
| `head_length` | Numeric | mm | Head length; the measurement used for male morph assignment. |
| `head_width` | Numeric | mm | Head width. |
| `forefemur` | Numeric | mm | Foreleg femur length, bilateral mean. |
| `foretibia` | Numeric | mm | Foreleg tibia length, bilateral mean. |
| `midfemur` | Numeric | mm | Midleg femur length, bilateral mean. |
| `midtibia` | Numeric | mm | Midleg tibia length, bilateral mean. |
| `hindfemur` | Numeric | mm | Hindleg femur length, bilateral mean. |
| `hindtibia` | Numeric | mm | Hindleg tibia length, bilateral mean. |
| `ear` | Numeric | mm² | Tympanic membrane area on the foreleg tibia, bilateral mean. This column stores area. |
| `eye` | Numeric | mm | Eye length, bilateral mean. |

## Male morph assignment: 2016 cut-offs

Morph assignments use head length in millimetres:

| `morph` value | Rule |
| --- | --- |
| `female` | `sex = f`, irrespective of head length |
| `eighth` | `sex = m` and `head_length < 18.50579` |
| `ninth` | `sex = m` and `18.50579 <= head_length <= 24.15225` |
| `tenth` | `sex = m` and `head_length > 24.15225` |

Both boundary values belong to the ninth-instar category. Apply the cut-offs to the unrounded `head_length` values. All 195 male assignments in this file match these rules. The recorded 2016 categorisations are retained throughout the manuscript analyses; mixture models were not refitted to this dataset.

## Bilateral averaging and data processing

Separate head, pronotum and leg measurement files were cleaned, joined by `ID`, and linked to sex metadata. For the leg segments, eyes and tympanic membranes, left and right measurements were combined as follows:

1. If both sides were available, their arithmetic mean was used.
2. If only one side was available, that measurement was used.
3. If neither side was available, the combined value was left missing.

This CSV contains the combined values. It does not retain side-specific measurements or a flag identifying individuals represented by one side only. Some values have long decimal representations arising from numerical storage and averaging; the number of printed decimal places should not be interpreted as measurement accuracy.

The preparation steps are documented in the project's `scripts/data_cleaning_*.R` files and `scripts/join_cleaned_files.R`. The latter normally writes `data/clean/trait_data.csv`; the fixed input used for manuscript analyses is `data/trait_data.csv`, described here.

## Missing observations

There are five empty measurement cells affecting two individuals:

| Variable | Missing values | Affected `ID` values |
| --- | ---: | --- |
| `midfemur` | 1 | `QUAR` |
| `midtibia` | 1 | `QUAR` |
| `hindfemur` | 1 | `QUAR` |
| `hindtibia` | 2 | `QUAR`, `TWQN` |

All other fields are complete. A total of 311 individuals have measurements for all 11 traits. All recorded morphological values are positive. Empty fields represent unavailable measurements, not zero values; the reason for each missing measurement is not encoded in this file. Missing observations are excluded for the variables required by each model, so analysis sample sizes can differ.

## Variables calculated during analysis

The following quantities are calculated from this CSV and are not additional columns in it:

| Derived quantity | Calculation | Units before logarithmic transformation |
| --- | --- | --- |
| Total foreleg length | `forefemur + foretibia` | mm |
| Total midleg length | `midfemur + midtibia` | mm |
| Total hindleg length | `hindfemur + hindtibia` | mm |
| Composite head size | `sqrt(head_length * head_width)` | mm |
| Linearised tympanic size | `sqrt(ear)` | mm |
| Centred log pronotum length | `log(pronotum) - log(7.2)` | Dimensionless |

Total leg length here is the sum of femur and tibia lengths. It does not include other leg segments. The analysis scripts use natural logarithms. Composite head size is a geometric mean of two lengths; tympanic area is square-root transformed before logarithmic analysis. The reference pronotum length for adjusted comparisons is 7.2 mm.

## Reuse and reproducibility

The CSV can be read in R without additional packages:

```r
dat <- read.csv(
  "trait_data.csv",
  na.strings = "",
  stringsAsFactors = FALSE
)
```

This example assumes the working directory contains the CSV. Within the analysis project, use the path `data/trait_data.csv`.

The project's `scripts/run_all.R` uses this file as its default input for the trait analyses. Running the analyses also requires the accompanying analysis scripts and their R packages. The reported analyses used R 4.4.2 and R 4.3.2; detailed software versions are recorded with the analysis outputs.

The repeated-measurement rounds and randomized measurement-order metadata used to assess repeatability are separate data files. They are not contained in `trait_data.csv`.

The sample represents adults encountered during vegetation sampling. Collection microhabitat, specimen-specific capture dates and search effort are not columns in this dataset. Male morph frequencies should therefore be interpreted as the composition of the sampled adults rather than as an unbiased estimate of population frequencies.

## File identification

This README was prepared on 25 September 2026 for the 313-row analysis input. The CSV was not changed during documentation.

- **File size:** 35,369 bytes.
- **SHA-256:** `6a8601014911b7b74c0821eb067d88f41a6b8f83db1e7b2383c0b9e65add4674`.
- **MD5 recorded in the analysis provenance:** `2d67b2d772f6f28b3cdf77c87b35b5fc`.
