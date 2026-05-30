# Higher-order Likelihood-based Inference for Weibull Models

This repository contains the R scripts used for the paper **“Higher-order Likelihood-based Inference for Weibull Models”**, submitted to **JQT**. The code reproduces the simulation studies for Weibull models under complete data, Type-II right censoring, and accelerated life testing (ALT), together with the R code used for the real-data case study.

The repository has been anonymised for double-blind review. Upon acceptance, the anonymous repository link will be replaced by a permanent public repository.

## Repository structure

```text
.
├── Bartlett_comparison.R
├── Complete_data.R
├── Type_2.R
├── ALT.R
├── Case_study.R
└── README.md
```

If the scripts are provided as `.txt` files during review, they can be renamed to `.R` before execution.

## Software requirements

The code was written in R and mainly relies on base R numerical optimisation routines, together with several plotting and data manipulation packages.

Required or commonly used packages include:

```r
install.packages(c("ggplot2", "dplyr", "tidyr"))
```

Some plotting commands also use functions from optional packages:

```r
install.packages(c("cowplot", "ggpubr", "ggh4x"))
```

These optional packages are only needed for specific plotting layouts, such as combined figures or independent facet axes.

## Description of scripts

### `Bartlett_comparison.R`

This script implements the complete-data Weibull simulation study for the shape parameter. It compares three likelihood-based inference methods:

* the signed root likelihood ratio statistic (r);
* the modified signed root statistic (r^\ast);
* a Bartlett-type correction based on the adjusted likelihood ratio statistic.

The script defines the complete-data Weibull profile log-likelihood for the shape parameter, computes the maximum likelihood estimator, constructs confidence intervals using numerical root-finding, and evaluates empirical coverage probabilities and average interval lengths.

The main simulation function is:

```r
sim_beta_cp_length()
```

By default, the simulation considers several combinations of Weibull shape and scale parameters, sample sizes (n = 5, 10, 20), and 1,000 Monte Carlo replications. The output is a data frame containing coverage probabilities and average interval lengths for the three methods. The script also produces coverage and interval-length plots.

### `Complete_data.R`

This script provides an additional complete-data comparison for the Weibull shape parameter. It compares:

* signed root likelihood inference;
* modified signed root likelihood inference;
* a profile-likelihood method based on the Technometrics critical-level approximation.

The code includes helper functions for the Technometrics critical levels, complete-data Weibull likelihood profiling, signed root and modified signed root statistics, and numerical confidence interval construction.

The main functions are:

```r
sim_length_grid()
sim_coverage_grid()
```

`sim_length_grid()` computes average confidence interval lengths, while `sim_coverage_grid()` computes empirical coverage probabilities. The simulation is run over a grid of Weibull parameter settings and sample sizes.

This script is intended as a complementary comparison for the complete-data case, especially for assessing interval length and coverage differences between higher-order likelihood corrections and profile-likelihood critical-level approximations.

### `Type_2.R`

This script implements the Type-II right-censoring simulation study for Weibull models. The data consist of the first (r) observed failures out of a total sample size (n). The likelihood accounts for both the observed failures and the remaining right-censored observations.

The script studies inference for three parameters:

* the Weibull shape parameter (\beta);
* the log-scale parameter (a = \log(\eta));
* the log quantile (\zeta_p = \log(Q_p)), with default focus on (p = 0.1).

For each parameter, the script computes and compares:

* the signed root statistic (r);
* the modified signed root statistic (r^\ast).

The main coverage simulation function is:

```r
sim_type2_grid()
```

The main interval-length simulation function is:

```r
sim_length_type2()
```

The default settings use (n = 10), observed failure counts (r = 3, 5, 8), several Weibull parameter scenarios, and Monte Carlo replications. The script produces coverage plots and average interval-length plots. It also saves the interval-length results to:

```r
type2_buchang.Rdata
```

### `ALT.R`

This script implements the accelerated life testing simulation study under a log-Weibull regression model:

```text
Y_i = log(T_i) = mu0 + mu1 * x_i + sigma * W_i,
```

where (W_i) follows a standard extreme-value distribution. The model is used to study reliability inference at a use-stress level.

The main parameter of interest is the log quantile at use stress:

```text
zeta_p(x_use) = mu0 + mu1 * x_use + sigma * log(-log(1-p)).
```

The script also studies inference for the model parameters:

* (\mu_0);
* (\mu_1);
* (\sigma);
* the log quantile (\zeta_p(x_{\text{use}})).

The code compares:

* signed root likelihood inference;
* modified signed root likelihood inference.

The main simulation functions are:

```r
sim_alt_quantile()
sim_alt_params_quantile()
sim_alt_length_styled()
```

`sim_alt_params_quantile()` computes empirical coverage probabilities for the model parameters and the use-stress log quantile. `sim_alt_length_styled()` computes average interval lengths.

The script produces coverage and interval-length plots for the ALT setting and saves the interval-length results to:

```r
alt_buchang.Rdata
```

### `Case_study.R`

This script contains the real-data case study based on a log-Weibull ALT model. The code defines the likelihood, canonical parameter calculations, numerical observed information matrix, maximum likelihood estimation, profile likelihood for the use-stress log quantile, signed root statistic, modified signed root statistic, and parametric bootstrap intervals.

The main function for the real-data analysis is:

```r
run_alt_case_study()
```

The function takes stress levels, grouped lifetime observations, a use-stress level, a grid of quantile probabilities, and a confidence level. It returns point estimates and confidence intervals for the use-stress quantile curve using:

* signed root likelihood intervals;
* modified signed root likelihood intervals;
* parametric bootstrap intervals.

The script also includes plotting functions:

```r
plot_alt_quantiles()
weibull_prob_plot()
```

The first function plots estimated quantile curves and confidence limits. The second function creates Weibull probability plots by stress level.

The case-study script further includes a simulation study calibrated from the real dataset. This additional simulation evaluates:

* mean squared error of the log-quantile estimator;
* empirical coverage probabilities of signed root, modified signed root, and bootstrap intervals;
* interval length distributions.

## How to run the code

Each script is self-contained and can be run directly in R or from the command line.

For example:

```bash
Rscript Bartlett_comparison.R
Rscript Complete_data.R
Rscript Type_2.R
Rscript ALT.R
Rscript Case_study.R
```

The scripts clear the R workspace at the beginning using:

```r
rm(list = ls())
```

Therefore, they should preferably be run in separate R sessions.

## Reproducibility

Most simulation scripts set random seeds internally, for example:

```r
seed = 123
```

or call:

```r
set.seed(123)
```

Using the default seeds and simulation settings should reproduce the reported simulation results, subject to minor numerical differences across R versions, operating systems, and optimisation routines.

For faster test runs, users may reduce the number of Monte Carlo replications, for example:

```r
n_rep = 100
```

For final reproduction of the simulation results, the larger settings used in the scripts, such as `n_rep = 1000`, should be used.

## Expected outputs

The scripts produce some or all of the following outputs:

* printed simulation result tables;
* coverage probability plots;
* average interval length plots;
* saved `.Rdata` files for selected length simulations;
* case-study quantile curves and Weibull probability plots.

Some `ggsave()` commands are included but commented out. To save figures as PDF files, uncomment the corresponding `ggsave()` lines in the relevant scripts.

## Notes for reviewers

The repository is designed to reproduce the computational results in the submitted manuscript. The scripts are intentionally kept transparent and self-contained, with likelihood functions, profiling routines, higher-order likelihood statistics, simulation loops, and plotting commands shown explicitly.

The anonymised version contains no author-identifying information. A permanent public version of the repository will be provided upon acceptance.

