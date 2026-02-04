# valytics 0.4.0

## New features

### Precision Experiments (EP05/EP15-aligned)

* `precision_study()`: Comprehensive variance component analysis for precision 
  experiments. Supports multiple experimental designs with automatic detection:
  - EP05-style: 20 days × 2 runs × 2 replicates
  - EP15-style: 5 days × 5 replicates (user verification)
  - Multi-site: 3+ sites × days × replicates (reproducibility)
  - Custom nested designs with any factor combination

* Variance component estimation via two methods:
  - ANOVA (Method of Moments) — default, works for balanced designs
  - REML (Restricted Maximum Likelihood) — requires lme4, better for unbalanced

* Confidence intervals for precision estimates:
  - Satterthwaite approximation (default)
  - Modified Large Sample (MLS) method
  - Bootstrap BCa method

* `verify_precision()`: Statistical verification of precision against 
  manufacturer claims using chi-square hypothesis testing. Features:
  - Accepts numeric vectors, precision_study objects, or data frames
  - Upper verification limit (UVL) calculation
  - Full hypothesis test output with p-values
  - Pass/fail determination with interpretation

* S3 methods: `print()`, `summary()`, `plot()`, `autoplot()` for both functions

* Three visualization types for precision studies:
  - `type = "variance"`: Variance component bar chart
  - `type = "cv"`: CV profile across concentration levels
  - `type = "precision"`: Forest plot with confidence intervals

## Dependencies

* `lme4` added to Suggests for REML estimation (optional)

## Documentation

* Comprehensive roxygen2 documentation for all exported functions
* Examples demonstrating typical precision study workflows

---

# valytics 0.3.0

## New features

* `ate_from_bv()`: Calculate allowable total error (ATE) specifications from 
  biological variation data using the Fraser-Petersen model. Supports three 
  performance levels (optimal, desirable, minimum) and provides allowable 
  imprecision, allowable bias, and total allowable error specifications.

* `sigma_metric()`: Calculate the Six Sigma metric for analytical performance 
  assessment. Returns sigma value with interpretation category (World Class 
  to Unacceptable) and approximate defect rates.

* `ate_assessment()`: Comprehensive evaluation of observed method performance 
  against allowable total error specifications. Provides pass/fail assessment 
  for individual components (bias, CV, total error) and overall method 
  acceptability, integrated with sigma metric calculation.

## Documentation

* New vignette: "Setting Quality Goals with Biological Variation" -- explains 
  the biological variation model for analytical performance specifications 
  with practical examples using the new ATE functions.

# valytics 0.2.0

## New features

* `deming_regression()`: Deming regression for method comparison, accounting 
  for measurement error in both variables. Supports known error ratio or 
  estimation from replicates. Includes jackknife and bootstrap BCa confidence 
  intervals.

* S3 methods for Deming regression: `print()`, `summary()`, `plot()`, and 
  `autoplot()` (ggplot2).

## Documentation

* New vignette: "Deming Regression for Method Comparison" -- comprehensive guide 
  to Deming regression theory and practical application.

* Updated vignette: "Understanding Method Comparison Statistics" -- added 
  guidance on choosing between regression methods.

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

* Vignette: "Method Comparison Workflow" -- step-by-step analysis guide

* Vignette: "Understanding Method Comparison Statistics" -- educational overview 
  of statistical concepts for method comparison studies
