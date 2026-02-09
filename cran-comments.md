## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: macOS (aarch64-apple-darwin), R 4.4.x
* GitHub Actions:
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (devel)
  - ubuntu-latest (release)
  - ubuntu-latest (oldrel-1)
* R-hub: [platforms tested]
* win-builder: R-devel

## Downstream dependencies

This package has no reverse dependencies.

## Notes for CRAN

This is the second release of valytics to CRAN.

### Changes in this version (0.4.0)

Major new functionality for precision experiments:

* `precision_study()`: Variance component analysis for nested experimental 
  designs with ANOVA and REML estimation methods

* `verify_precision()`: Statistical verification of precision against 
  manufacturer claims using chi-square hypothesis testing

* `precision_profile()`: Modeling CV-concentration relationships for 
  functional sensitivity estimation

* New dataset: `troponin_precision` for precision study examples

* New vignette: "Precision Profiles and Functional Sensitivity"

The lme4 package has been added to Suggests for optional REML estimation.