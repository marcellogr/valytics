# valytics 0.2.0

## New features

* `deming_regression()`: Deming regression (errors-in-variables regression) for 
  method comparison. Unlike ordinary least squares, Deming regression accounts 
  for measurement error in both variables. Key features:
  
  - Customizable error ratio (lambda) for methods with different precision
  - Orthogonal regression (lambda = 1) as default
  - Jackknife confidence intervals (default) following Linnet (1990)
  - Bootstrap BCa confidence intervals as alternative
  - Perpendicular residuals for model diagnostics

* S3 methods for Deming regression: `print()`, `summary()`, `plot()`, and 
  `autoplot()` (ggplot2).

* Publication-ready visualizations for Deming regression including scatter plots 
  with confidence bands and residual plots.

## Documentation

* Added comprehensive documentation for Deming regression with references to 
  Linnet (1990, 1993) and Cornbleet & Gochman (1979).

## References

* Linnet K (1990). Estimation of the linear relationship between the 
  measurements of two methods with proportional errors. Statistics in Medicine, 
  9(12):1463-1473.
  
* Linnet K (1993). Evaluation of regression procedures for methods comparison 
  studies. Clinical Chemistry, 39(3):424-432.


# valytics 0.1.0

Initial CRAN release.
 
## New features

* `ba_analysis()`: Bland-Altman method comparison analysis with bias estimation, 
  limits of agreement, and confidence intervals. Supports both absolute and 
  percentage difference scaling.

* `pb_regression()`: Passing-Bablok regression with fast O(n log n) algorithm 
  via the robslopes package. Includes analytical confidence intervals 
  (Passing & Bablok 1983) and optional bootstrap BCa intervals. CUSUM test 
  for linearity assessment with Kolmogorov-Smirnov p-value.

* S3 methods for both analyses: `print()`, `summary()`, `plot()`, and 
  `autoplot()` (ggplot2).

* Publication-ready visualizations using ggplot2, including Bland-Altman 
  plots, regression scatter plots with confidence bands, residual plots, 
  and CUSUM plots for linearity assessment.

## Datasets
 
* `glucose_methods`: Point-of-care glucose meter vs laboratory analyzer (n=60)

* `creatinine_serum`: Enzymatic vs Jaffe creatinine methods (n=80)

* `troponin_cardiac`: Two high-sensitivity cardiac troponin I platforms (n=50)

## Documentation

* Vignette: "Method Comparison Workflow" - step-by-step analysis guide

* Vignette: "Understanding Method Comparison Statistics" - educational overview 
  of statistical concepts for method comparison studies