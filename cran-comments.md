## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE is for new submission and possibly misspelled words:
- Bablok: Correct spelling (Passing-Bablok regression method)
- CLSI: Acronym for Clinical and Laboratory Standards Institute

## Test environments

* Local: macOS (aarch64-apple-darwin), R 4.4.x
* GitHub Actions:
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (devel)
  - ubuntu-latest (release)
  - ubuntu-latest (oldrel-1)

## Downstream dependencies

This is a new package with no reverse dependencies.

## Resubmission

This is a resubmission. In this version I have:
 
* Replaced Unicode characters (lambda, not-equal signs) with LaTeX-safe 
  alternatives in roxygen2 documentation to fix PDF manual generation errors.

## Notes for CRAN

The package provides statistical methods for analytical method comparison 
and validation, including Bland-Altman analysis, Passing-Bablok regression, 
Deming regression, and allowable total error utilities. These methods are 
widely used in clinical laboratory validation and biomedical research.

The package depends on robslopes (GPL-3) for the fast Passing-Bablok algorithm, 
which requires this package to also be GPL-3 licensed.
