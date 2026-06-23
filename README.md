# Higher-order Likelihood-based Inference for Weibull Models

This repository contains the R code used to reproduce the numerical results for the manuscript **"Higher-order Likelihood-based Inference for Weibull Models"**.

The code focuses on simulation studies and a real-data-inspired case-study experiment for Weibull reliability models. It compares confidence interval procedures for Weibull parameters and lifetime quantiles under complete samples, Type II right censoring, and accelerated life testing (ALT).

## What this repository reproduces

The repository is organized around four self-contained R scripts:

| Script | Main purpose | Main outputs |
|---|---|---|
| `Sim_complete.R` | Reproduces the complete-sample Weibull simulation study. | Coverage probabilities and average interval lengths for shape, log scale, and log quantile. |
| `Sim_typeII.R` | Reproduces the Type II right-censoring simulation study. | Coverage probabilities and average interval lengths under different numbers of observed failures. |
| `Sim_alt.R` | Reproduces the accelerated life testing simulation study. | Coverage probabilities and average interval lengths for ALT regression parameters and use-stress quantiles. |
| `Case_study.R` | Reproduces the fitted-model case-study simulation based on the insulating-fluid breakdown data. | Coverage summaries and interval-length boxplots for the use-stress log quantile. |

Each script can be run independently.

## Repository structure

A simple repository structure is:

```text
.
├── README.md
├── Sim_complete.R
├── Sim_typeII.R
├── Sim_alt.R
└── Case_study.R
```

If preferred, the scripts can also be placed in a `code/` folder:

```text
.
├── README.md
├── code/
│   ├── Sim_complete.R
│   ├── Sim_typeII.R
│   ├── Sim_alt.R
│   └── Case_study.R
├── figure/
└── output/
```

If you use this structure, please update the file paths or working directory before running the scripts.

## Software requirements

The code was written in R. The main external packages are:

```r
install.packages(c("dplyr", "ggplot2", "ggh4x"))
```

`ggh4x` is used only for some plots with independent facet scales. The core simulation code does not depend on it.

## How to run

From R or RStudio, run one script at a time:

```r
source("Sim_complete.R")
source("Sim_typeII.R")
source("Sim_alt.R")
source("Case_study.R")
```

From a terminal, the scripts can also be run with:

```bash
Rscript Sim_complete.R
Rscript Sim_typeII.R
Rscript Sim_alt.R
Rscript Case_study.R
```

Each script starts by clearing the R workspace, so it is recommended to run them separately.

## Script descriptions

### 1. Complete Weibull samples

File:

```text
Sim_complete.R
```

This script reproduces the simulation study for complete Weibull samples. It evaluates confidence interval performance for:

- the Weibull shape parameter;
- the log scale parameter;
- the log 0.1 quantile.

The simulation varies the Weibull parameter setting and sample size. It produces coverage plots and interval-length plots comparing the competing methods considered in the manuscript.

Typical output objects include:

```r
res_complete
res_length
p_cov
p_len
```

### 2. Type II right-censored samples

File:

```text
Sim_typeII.R
```

This script reproduces the Type II right-censoring simulation study. The total number of units is fixed, and the experiment stops after a specified number of failures.

The script evaluates confidence intervals for:

- the Weibull shape parameter;
- the log scale parameter;
- the log 0.1 quantile.

Typical output objects include:

```r
res_type2
res_length
p_type2_95
p_length
```

The script saves the interval-length figure as:

```text
type2_length.pdf
```

### 3. Accelerated life testing simulation

File:

```text
Sim_alt.R
```

This script reproduces the ALT simulation study under a log-Weibull regression model. It evaluates interval performance for:

- the intercept parameter;
- the stress-effect parameter;
- the Weibull shape-related parameter;
- the log 0.1 quantile at the use stress.

Typical output objects include:

```r
res_alt_all
res_len
p_alt
p_len
```

The script saves the interval-length figure as:

```text
alt_length.pdf
```

### 4. Insulating-fluid fitted-model experiment

File:

```text
Case_study.R
```

This script reproduces the fitted-model experiment based on the insulating-fluid breakdown-time data used in the manuscript.

The data are included directly in the script. The experiment fits the ALT model, treats the fitted values as the data-generating parameters, and compares interval procedures for the use-stress log 0.1 quantile.

Typical output object:

```r
sim_result
```

The script saves the interval-length boxplot as:

```text
case_IT.pdf
```

## Reproducing manuscript-scale results

Some scripts may use smaller Monte Carlo sizes for quicker testing. To reproduce the full manuscript-scale experiments, check the values of:

```r
n_rep
boot_B
```

and increase them as needed. In the manuscript experiments, large simulation sizes are used, so the full runs may take substantial time.

For a quick test, use smaller values such as:

```r
n_rep = 100
boot_B = 100
```

For final reproduction, use the larger values reported in the manuscript.

## Expected figure files

Depending on which `ggsave()` lines are active, the scripts may generate files such as:

```text
weibull_complete_coverage.pdf
complete_length.pdf
weibull_typeII_coverage_95.pdf
type2_length.pdf
weibull_alt_coverage.pdf
alt_length.pdf
case_IT.pdf
```

Some saving commands are commented out by default. Uncomment the corresponding `ggsave()` command if a figure should be written to disk.

## Notes

- The scripts are self-contained and can be run independently.
- Random seeds are set in the simulation calls to make the results reproducible.
- Long simulation runs may require substantial computation time.
- The code is intended for reproducing the numerical studies in the manuscript rather than as a general-purpose R package.

## Citation

If you use this code, please cite the accompanying manuscript:

```bibtex
@unpublished{weibull_higher_order_likelihood,
  title  = {Higher-order Likelihood-based Inference for Weibull Models},
  author = {Author names to be added},
  year   = {2026},
  note   = {Manuscript}
}
```

Please replace the placeholder author information with the final citation information before public release.

## License

No license is specified yet. If the repository is made public, please add a `LICENSE` file to clarify reuse conditions.
