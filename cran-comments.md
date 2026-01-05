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

This is a new package with no reverse dependencies.

## Notes for CRAN

This is the first submission of valytics to CRAN.

The package provides statistical methods for analytical method comparison, 
including Bland-Altman analysis and Passing-Bablok regression. These methods 
are widely used in clinical laboratory validation and biomedical research.

The package depends on robslopes (GPL-3) for the fast Passing-Bablok algorithm, 
which requires this package to also be GPL-3 licensed.
