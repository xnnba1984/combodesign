# Component contribution design for combination therapy trials

This repository contains the R package and reproducible analysis code for a
paper on component-contribution designs in fixed three-arm combination therapy
trials. The primary setting compares a combination arm, `AB`, with its two
single-component arms, `A` and `B`, and evaluates whether both components make
the required contribution to the combination.

The code supports:

- operating-characteristic calculation for the two contribution contrasts;
- joint-power and allocation calculations for full contribution claims;
- optional separate component-claim calculations using maxT, Holm, or
  Bonferroni adjustments;
- simulation outputs used in the manuscript;
- the public PDX case-study analysis and display files.

## Repository layout

```
.
├── README.md
├── data/
│   └── combination.csv
├── package/
│   └── combodesign/
├── analysis/
├── simulation/
├── case_study/
│   └── pdx/
├── figures/
├── tables/
└── reproducibility/
```

`package/combodesign` is the R package. The other folders contain the analysis
scripts, simulation scripts, public-data case-study scripts, generated tables,
generated figures, and a run-all QA driver.

## R setup

Install package dependencies first.

```r
install.packages(c(
  "mvtnorm",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "purrr",
  "ggplot2",
  "patchwork"
))
```

The figure scripts can use `ragg` if it is installed, but it is optional.

To install the local R package from a terminal:

```bash
R CMD INSTALL package/combodesign
```

Or install it from R:

```r
install.packages("remotes")
remotes::install_local("package/combodesign")
```

After the repository is on GitHub, users can install it from R with:

```r
install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USER/YOUR_REPO", subdir = "package/combodesign")
```

Replace `YOUR_GITHUB_USER/YOUR_REPO` with the actual repository path.

## Quick example

```r
library(combodesign)

delta <- c(AB_minus_A = 0.55, AB_minus_B = 0.65)

fit <- component_contribution_design(
  delta = delta,
  target_power = 0.80,
  alpha = 0.025,
  min_allocation = 0.10
)

fit$summary

oc <- contribution_operating_characteristics(
  N = 180,
  allocation = c(A = 1 / 3, B = 1 / 3, AB = 1 / 3),
  delta = delta,
  alpha = 0.025
)

oc$joint_power
```

In this package, `joint_power` is the probability that both required
component-contribution tests reject under the specified design assumptions.

## Reproduce the paper outputs

From the repository root, run the quick QA driver:

```bash
Rscript reproducibility/run_all_analysis.R --mode=quick
```

Quick mode reruns the core method checks, builds and checks the R package, and
verifies the included simulation and PDX output files.

To regenerate the simulation displays from the included main simulation output:

```bash
Rscript simulation/create_simulation_displays.R
Rscript simulation/create_combined_simulation_figure3.R
```

To rerun the PDX case study and rebuild its displays, install the Python
dependencies and run:

```bash
python3 -m pip install -r requirements-python.txt
Rscript case_study/pdx/run_pdx_case_study.R
python3 case_study/pdx/create_pdx_displays.py
```

To rerun the larger analysis chain:

```bash
Rscript reproducibility/run_all_analysis.R --mode=full --python=$(which python3)
```

The full mode reruns core checks, package checks, participant-level simulation,
simulation displays, the PDX case study, and PDX displays. It may take longer
than quick mode.

## Data note

`data/combination.csv` is the public PDX combination-treatment data file used
for the case-study analysis. Before making a public repository, confirm that
redistributing the raw data file is acceptable under the source data terms and
your institution's policy. If redistribution is not allowed, remove
`data/combination.csv` and replace it with download instructions while keeping
the same local file path for reproducibility.

## Outputs

Important generated outputs are stored in:

- `tables/simulation/`
- `figures/simulation/`
- `tables/pdx/`
- `figures/pdx/`

The `MANIFEST.csv` file lists every included file with byte size and SHA-256
checksum.

## License

The R package uses the MIT license. See `LICENSE.md`.
