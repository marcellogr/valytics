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

* Vignette: "Method Comparison Workflow" — step-by-step analysis guide

* Vignette: "Understanding Method Comparison Statistics" — educational overview 
  of statistical concepts for method comparison studies
