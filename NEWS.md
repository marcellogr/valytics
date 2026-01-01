# valytics 0.1.0

Initial CRAN release.

## Features

### Bland-Altman Analysis
* `ba_analysis()`: Bland-Altman analysis for method comparison
  - Vector and formula interfaces
  - Absolute and percent difference types
  - Confidence intervals for bias and limits of agreement (Bland & Altman 1999)
  - `print()`, `summary()`, and `plot()` methods
  - `autoplot()` method for ggplot2-style plotting

### Passing-Bablok Regression
* `pb_regression()`: Non-parametric regression for method comparison

  - Fast O(n log n) algorithm via `robslopes` package
  - Analytical confidence intervals (Passing & Bablok 1983)
  - Bootstrap BCa confidence intervals (optional)
  - CUSUM test for linearity assessment
  - `print()`, `summary()`, and `plot()` methods
  - Multiple plot types: scatter, residuals, CUSUM

### Example Datasets
* `glucose_methods`: POC glucose meter vs laboratory analyzer (n=60)
* `creatinine_serum`: Enzymatic vs Jaffe creatinine methods (n=80)
* `troponin_cardiac`: Two high-sensitivity troponin platforms (n=50)

## Dependencies
* `ggplot2` for publication-ready visualizations
* `robslopes` for fast Passing-Bablok point estimation
