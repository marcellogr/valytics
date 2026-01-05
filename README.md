
<!-- README.md is generated from README.Rmd. Please edit that file -->

# valytics <img src="man/figures/logo.png" align="right" height="139" alt="valytics logo" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/marcellogr/valytics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/marcellogr/valytics/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Statistical methods for analytical method comparison and validation
studies. The package implements Bland-Altman analysis and Passing-Bablok
regression — approaches commonly used in clinical laboratory method
validation.

## Installation

You can install the development version of valytics from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("marcellogr/valytics")
```

or from CRAN (when available)

``` r
install.packages("valytics")
```

## Overview

`valytics` provides tools for comparing two measurement methods:

- **Bland-Altman analysis**: Assess agreement through bias and limits of
  agreement
- **Passing-Bablok regression**: Non-parametric regression robust to
  outliers

Both methods produce publication-ready plots and comprehensive
statistical summaries.

## Example: Bland-Altman Analysis

``` r
library(valytics)

# Compare two creatinine measurement methods
data(creatinine_serum)

ba <- ba_analysis(
 x = creatinine_serum$enzymatic,
 y = creatinine_serum$jaffe
)

ba
#> 
#> Bland-Altman Analysis
#> ---------------------------------------- 
#> n = 80 paired observations
#> 
#> Difference type: Absolute (y - x)
#> Confidence level: 95%
#> 
#> Results:
#>   Bias (mean difference): 0.174
#>     95% CI: [0.127, 0.220]
#>   SD of differences: 0.209
#> 
#> Limits of Agreement:
#>   Lower LoA: -0.236
#>     95% CI: [-0.316, -0.156]
#>   Upper LoA: 0.584
#>     95% CI: [0.504, 0.663]
```

``` r
plot(ba)
```

<img src="man/figures/README-ba-plot-1.png" alt="Bland-Altman plot showing differences between Jaffe and enzymatic creatinine methods" width="100%" />

## Example: Passing-Bablok Regression

``` r
pb <- pb_regression(
 x = creatinine_serum$enzymatic,
 y = creatinine_serum$jaffe
)

pb
#> 
#> Passing-Bablok Regression
#> ---------------------------------------- 
#> n = 80 paired observations
#> 
#> CI method: Analytical (Passing-Bablok 1983)
#> Confidence level: 95%
#> 
#> Regression equation:
#>   creatinine_serum$jaffe = 0.234 + 0.971 * creatinine_serum$enzymatic
#> 
#> Results:
#>   Intercept: 0.234
#>     95% CI: [0.229, 0.239]
#>     (excludes 0: significant constant bias)
#> 
#>   Slope: 0.971
#>     95% CI: [0.966, 0.974]
#>     (excludes 1: significant proportional bias)
```

``` r
plot(pb)
```

<img src="man/figures/README-pb-plot-1.png" alt="Passing-Bablok regression scatter plot with regression line and confidence band" width="100%" />

## Features

- **Multiple interfaces**: Vector input or formula syntax
  (`method1 ~ method2`)
- **Flexible CI methods**: Analytical (Passing-Bablok 1983) or bootstrap
  BCa
- **Assumption checking**: CUSUM linearity test, Shapiro-Wilk normality
  test
- **Publication-ready plots**: Built on ggplot2, fully customizable
- **Tidy workflows**: Consistent API, informative error messages

## Example Datasets

The package includes three realistic clinical datasets:

| Dataset            | Description                       | n   |
|--------------------|-----------------------------------|-----|
| `glucose_methods`  | POC meter vs. laboratory analyzer | 60  |
| `creatinine_serum` | Enzymatic vs. Jaffe methods       | 80  |
| `troponin_cardiac` | Two hs-cTnI platforms             | 50  |

## References

**Bland-Altman:**

- Bland JM, Altman DG (1986). Statistical methods for assessing
  agreement between two methods of clinical measurement. *Lancet*,
  1(8476):307-310.
- Bland JM, Altman DG (1999). Measuring agreement in method comparison
  studies. *Statistical Methods in Medical Research*, 8(2):135-160.

**Passing-Bablok:**

- Passing H, Bablok W (1983). A new biometrical procedure for testing
  the equality of measurements from two different analytical methods.
  *Journal of Clinical Chemistry and Clinical Biochemistry*,
  21(11):709-720.

## License

GPL-3
