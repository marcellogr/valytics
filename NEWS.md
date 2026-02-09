# valytics 0.4.0

## New Features

### Precision Experiments (EP05/EP15-aligned)

* `precision_study()`: Comprehensive variance component analysis for precision 
  experiments with nested experimental designs.
  
  - Supports multiple design types: single-site (day/run/replicate), 
    multi-site (site/day/run/replicate), and custom nested designs
  - Automatic design detection from data structure
  - Two estimation methods: ANOVA (method of moments) and REML (via lme4)
  - Three confidence interval methods: Satterthwaite (default), Modified Large 
    Sample (MLS), and bootstrap BCa
  - Multi-sample support for analyzing multiple concentration levels
  - Complete S3 methods: `print()`, `summary()`, `plot()`, `autoplot()`

* `verify_precision()`: Statistical verification of observed precision against 
  manufacturer claims using chi-square hypothesis testing.
  
  - Accepts numeric vectors, precision_study objects, or data frames
  - Calculates Upper Verification Limit (UVL)
  - Provides confidence intervals for observed precision
  - Clear pass/fail determination with detailed interpretation

* `precision_profile()`: Models the relationship between CV and concentration
  for functional sensitivity estimation.
  
  - Hyperbolic model: CV = sqrt(a² + (b/x)²) (default)
  - Linear model: CV = a + b/x
  - Calculates functional sensitivity at user-specified CV targets
  - Optional bootstrap confidence intervals for functional sensitivity
  - Integrates with precision_study objects for seamless workflow

### Precision Measures

The package now calculates and reports:
- **Repeatability**: Within-run precision
- **Between-run precision**: Additional variability between runs within a day
- **Between-day precision**: Additional variability between days
- **Within-laboratory precision**: Combined day + run + error variance
- **Between-site precision**: Additional variability between sites (multi-site only
- **Reproducibility**: Total precision including all variance components

### Visualization

* `plot.precision_study()` with three plot types:
  - `type = "variance"`: Variance component bar chart
  - `type = "cv"`: CV profile across precision measures with CIs
  - `type = "precision"`: Forest plot of precision estimates

* `plot.precision_profile()`: Publication-ready precision profile visualization
  - Fitted curve with prediction intervals
  - Functional sensitivity target lines
  - Optional logarithmic x-axis scale

## New Dataset

* `troponin_precision`: High-sensitivity cardiac troponin I precision study data
  with 6 concentration levels (5-500 ng/L), designed for demonstrating 
  `precision_study()` and `precision_profile()` workflows.

## Documentation

* New vignette: "Precision Profiles and Functional Sensitivity" - demonstrates 
  the complete workflow from raw precision data to functional sensitivity 
  estimation with clinical interpretation.

## Dependencies

* `lme4` added to Suggests for REML estimation (optional)

---

# valytics 0.3.0

Initial CRAN release.

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
