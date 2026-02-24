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
* win-builder: R-devel

## Downstream dependencies

This package has no reverse dependencies.

## Notes for CRAN

This is a patch release (v0.4.1) that fixes a data generation bug in the 
`troponin_precision` dataset included in v0.4.0. The dataset was generated 
with constant values (zero variance) at each concentration level, which 
caused `precision_profile()` examples and the associated vignette to fail.

The dataset has been regenerated with realistic variability following a 
hyperbolic CV model appropriate for high-sensitivity immunoassays.

No changes to package functionality - only the example dataset was corrected.