# Example code: long-term temperature-attributable mortality spectrum

This repository provides **example R code** illustrating the principal analytical procedures used in the study:

> *Climate, Demographic and Societal Transitions Reshaped the Burden of Disease for Non-optimal Temperature-attributable Mortality over the Past Half-century*

This is a compact methodological demonstration, not a complete reproduction package. The mortality data used in the study are subject to legal restrictions and cannot be redistributed. The included datasets are fully synthetic and contain no Japanese Vital Statistics records.

## What is demonstrated

1. `code/01_temperature_mortality_model.R` fits a one-stage space-time-stratified case-crossover model with a distributed lag non-linear temperature term, identifies the minimum mortality temperature (MMT), and exports a cumulative relative-risk curve and P1/P99 relative risks.
2. `code/02_attributable_mortality.R` uses the fitted cross-basis coefficients to estimate cold-, heat-, and total temperature-attributable fractions and numbers with a backward attributable-risk formulation and Monte Carlo uncertainty intervals.
3. `code/03_shapley_decomposition.R` demonstrates exact five-factor Shapley decomposition of changes in (i) annual attributable deaths and (ii) cause-specific attributable mortality shares. The five factors are population size, age structure, cause-specific mortality rates, the temperature–death distribution, and the exposure–response function.

## Repository structure

```text
temperature_mortality_spectrum_example/
├── README.md
├── code/
│   ├── 01_temperature_mortality_model.R
│   ├── 02_attributable_mortality.R
│   └── 03_shapley_decomposition.R
├── data/
│   ├── README.md
│   ├── example_daily_data.csv
│   └── example_decomposition_data.csv
└── output/
    └── README.md
```

## System requirements

- The study analyses were conducted in R 4.4.2. The example scripts were tested successfully on both Windows and macOS, but have not been systematically tested across other operating systems or multiple R versions.
- The example scripts require the R packages `dplyr`, `readr`, `lubridate`, `tibble`, `dlnm`, and `gnm`. The base R package `splines` is also used by the temperature–mortality model.
- No non-standard hardware is required; a standard desktop or laptop capable of running R is sufficient.

## Installation

Install R from [CRAN](https://cran.r-project.org/) if it is not already available. In R, install the required packages with:

```r
install.packages(c(
  "dplyr", "readr", "lubridate", "tibble",
  "dlnm", "gnm"
))
```

On a typical desktop computer with a stable internet connection and R already installed, package installation is expected to take about 5 minutes. The actual time depends on the connection and whether the packages are already installed.

## Running the example

Start R with the repository root as the working directory, then run:

```r
source("code/01_temperature_mortality_model.R")
source("code/02_attributable_mortality.R")
source("code/03_shapley_decomposition.R")
```

The scripts use relative paths and write generated files to `output/`. The model demonstration uses one synthetic mortality outcome; the decomposition demonstration uses a separate synthetic dataset. On a typical desktop computer, the complete example is expected to take about 2 minute. Runtime depends on the computer and package setup.

## Expected output

The scripts generate model summaries and relative-risk results, cold-, heat-, and total-attributable fractions and numbers, uncertainty summaries, and five-factor Shapley contributions. These outputs demonstrate the analysis workflow and are not estimates from the study. See `output/README.md` for the output directory note.

## Using your own data

The model and attributable-mortality scripts use the structure shown in `data/example_daily_data.csv`: one row per prefecture and date, with date, prefecture identifier, daily mean temperature, daily death count, and a holiday indicator. To use other data, prepare the columns in the same format, update the input file path and any study-specific settings in the relevant scripts, and ensure that the data cover the required dates and locations. The decomposition script uses the separate factor structure documented in `data/README.md`. The synthetic files are examples only; users must ensure that their own data are legally available for their intended use.

## Relationship to the manuscript methods

The example preserves the core specifications of the main analysis where applicable:

- space-time strata defined by prefecture, year, month, and day of week;
- conditional quasi-Poisson regression;
- natural cubic spline for temperature with internal knots at the 10th, 75th, and 90th percentiles;
- natural cubic spline for lag with three internal knots equally spaced on the log scale;
- maximum lag of 21 days;
- MMT search restricted to the 50th–99th percentiles;
- cold and heat relative risks evaluated at P1 and P99 relative to the MMT;
- backward attributable-risk calculation;
- exact five-factor Shapley averaging across all 32 counterfactual states.

The synthetic example is smaller than the study dataset and uses one mortality outcome for the model demonstration to keep runtime and dependencies manageable.

## Synthetic data

`data/example_daily_data.csv` contains synthetic prefecture-day mortality counts, temperature, and a holiday indicator. `data/example_decomposition_data.csv` contains synthetic factor values for demonstrating the five-factor decomposition. See `data/README.md` for details. Both datasets were generated for this repository and contain no records from the Japanese Vital Statistics data used in the study.

## Reproducibility

Random seeds are fixed in scripts that use simulation. No absolute file paths are used. The scripts include basic input validation and decomposition reconciliation checks. The example demonstrates the methods but cannot reproduce the manuscript estimates without access to the restricted mortality data and the full study inputs.

## Data availability

The mortality data used in the study were obtained from the Ministry of Health, Labour and Welfare of Japan. Their redistribution is restricted. Meteorological and population data are available from the original public data providers described in the manuscript.

## License

The example code in this repository is released under the MIT License.
