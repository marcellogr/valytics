#' Calculate Allowable Total Error from Biological Variation
#'
#' @description
#' Calculates analytical performance specifications (allowable imprecision,
#' allowable bias, and total allowable error) based on biological variation
#' data using the hierarchical model from Fraser & Petersen (1993).
#'
#' @param cvi Numeric. Within-subject (intra-individual) biological variation
#'   coefficient of variation, expressed as a percentage.
#' @param cvg Numeric. Between-subject (inter-individual) biological variation
#'   coefficient of variation, expressed as a percentage. If `NULL` (default),
#'   only imprecision specifications are calculated.
#' @param level Character. Performance level: `"desirable"` (default),
#'   `"optimal"`, or `"minimum"`. See Details.
#' @param k Numeric. Coverage factor for total allowable error calculation
#'   (default: 1.65 for ~95% coverage assuming normal distribution).
#'
#' @return An object of class `c("ate_specs", "valytics_ate", "valytics_result")`,
#'   which is a list containing:
#'
#'   \describe{
#'     \item{specifications}{List with calculated specifications:
#'       \itemize{
#'         \item `allowable_cv`: Allowable analytical imprecision (CV_A)
#'         \item `allowable_bias`: Allowable analytical bias (NULL if cvg not provided)
#'         \item `tea`: Total allowable error (NULL if cvg not provided)
#'       }
#'     }
#'     \item{input}{List with input parameters:
#'       \itemize{
#'         \item `cvi`: Within-subject CV
#'         \item `cvg`: Between-subject CV (or NULL)
#'         \item `level`: Performance level used
#'         \item `k`: Coverage factor used
#'       }
#'     }
#'     \item{multipliers}{List with level-specific multipliers used:
#'       \itemize{
#'         \item `imprecision`: Multiplier for CV_I
#'         \item `bias`: Multiplier for sqrt(CV_I^2 + CV_G^2)
#'       }
#'     }
#'   }
#'
#' @details
#' The biological variation model for analytical performance specifications
#' was developed by Fraser, Petersen, and colleagues. The model derives
#' allowable analytical error from the inherent biological variability of
#' the measurand.
#'
#' **Formulas (Desirable level):**
#'
#' \deqn{CV_A \leq 0.50 \times CV_I}{CV_A <= 0.50 * CV_I}
#'
#' \deqn{Bias \leq 0.25 \times \sqrt{CV_I^2 + CV_G^2}}{Bias <= 0.25 * sqrt(CV_I^2 + CV_G^2)}
#'
#' \deqn{TEa \leq k \times CV_A + Bias}{TEa <= k * CV_A + Bias}
#'
#' Where:
#' \itemize{
#'   \item CV_I = within-subject biological variation
#'   \item CV_G = between-subject biological variation
#'   \item CV_A = allowable analytical imprecision
#'   \item k = coverage factor (default 1.65)
#' }
#'
#' **Performance Levels:**
#'
#' Three hierarchical performance levels are defined:
#' \itemize{
#'   \item **Optimal**: Most stringent; multipliers are 0.25x desirable
#'     (i.e., 0.125 for CV, 0.0625 for bias)
#'   \item **Desirable**: Standard target; multipliers are 0.50 for CV,
#'     0.25 for bias
#'   \item **Minimum**: Least stringent; multipliers are 1.5x desirable
#'     (i.e., 0.75 for CV, 0.375 for bias)
#' }
#'
#' @section Data Sources:
#' Biological variation data (CV_I and CV_G) should be obtained from
#' authoritative sources. The recommended current source is the
#' **EFLM Biological Variation Database**: \url{https://biologicalvariation.eu/}
#'
#' This database provides rigorously reviewed BV estimates derived from
#' published studies meeting defined quality specifications.
#'
#' @references
#' Fraser CG, Petersen PH (1993). Desirable standards for laboratory tests
#' if they are to fulfill medical needs. \emph{Clinical Chemistry},
#' 39(7):1447-1453.
#'
#' Ricos C, Alvarez V, Cava F, et al. (1999). Current databases on biological
#' variation: pros, cons and progress. \emph{Scandinavian Journal of Clinical
#' and Laboratory Investigation}, 59(7):491-500.
#'
#' Aarsand AK, Fernandez-Calle P, Webster C, et al. (2020). The EFLM
#' Biological Variation Database. \url{https://biologicalvariation.eu/}
#'
#' Westgard JO (2008). \emph{Basic Method Validation} (3rd ed.).
#' Westgard QC, Inc.
#'
#' @seealso
#' [sigma_metric()] for calculating Six Sigma metrics,
#' [ate_assessment()] for comparing observed performance to specifications
#'
#' @examples
#' # Glucose: CV_I = 5.6%, CV_G = 7.5% (example values)
#' ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)
#' ate
#'
#' # Optimal performance level (more stringent)
#' ate_optimal <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "optimal")
#' ate_optimal
#'
#' # Minimum acceptable performance
#' ate_min <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "minimum")
#' ate_min
#'
#' # When only within-subject CV is known (bias goal not calculable)
#' ate_cv_only <- ate_from_bv(cvi = 5.6)
#' ate_cv_only
#'
#' # Custom coverage factor (e.g., 2.0 for ~97.5% coverage)
#' ate_custom <- ate_from_bv(cvi = 5.6, cvg = 7.5, k = 2.0)
#'
#' # Access individual specifications
#' ate$specifications$allowable_cv
#' ate$specifications$allowable_bias
#' ate$specifications$tea
#'
#' @export
ate_from_bv <- function(cvi,
                        cvg = NULL,
                        level = c("desirable", "optimal", "minimum"),
                        k = 1.65) {


  # Input validation ----
 .validate_ate_input(cvi, cvg, k)

  level <- match.arg(level)

  # Get multipliers for the specified level ----
  multipliers <- .get_ate_multipliers(level)

  # Calculate specifications ----
  # Allowable imprecision: CV_A <= multiplier * CV_I

  allowable_cv <- multipliers$imprecision * cvi

 # Allowable bias and TEa require CV_G
  if (!is.null(cvg)) {
    # Allowable bias: Bias <= multiplier * sqrt(CV_I^2 + CV_G^2)
    total_bv <- sqrt(cvi^2 + cvg^2)
    allowable_bias <- multipliers$bias * total_bv

    # Total allowable error: TEa = k * CV_A + Bias
    tea <- k * allowable_cv + allowable_bias
  } else {
    allowable_bias <- NULL
    tea <- NULL
  }

  # Construct output object ----
  structure(
    list(
      specifications = list(
        allowable_cv = allowable_cv,
        allowable_bias = allowable_bias,
        tea = tea
      ),
      input = list(
        cvi = cvi,
        cvg = cvg,
        level = level,
        k = k
      ),
      multipliers = multipliers
    ),
    class = c("ate_specs", "valytics_ate", "valytics_result")
  )
}


# Helper Functions ----
# =============================================================================

#' Validate input for ate_from_bv
#' @noRd
#' @keywords internal
.validate_ate_input <- function(cvi, cvg, k) {

  # Check cvi
  if (length(cvi) != 1 || !is.numeric(cvi) || is.na(cvi)) {
    stop("`cvi` must be a single numeric value.", call. = FALSE)
  }
  if (cvi <= 0) {
    stop("`cvi` must be a positive number.", call. = FALSE)
  }

  # Check cvg (if provided)
  if (!is.null(cvg)) {
    if (length(cvg) != 1 || !is.numeric(cvg) || is.na(cvg)) {
      stop("`cvg` must be a single numeric value or NULL.", call. = FALSE)
    }
    if (cvg <= 0) {
      stop("`cvg` must be a positive number.", call. = FALSE)
    }
  }

  # Check k
  if (length(k) != 1 || !is.numeric(k) || is.na(k)) {
    stop("`k` must be a single numeric value.", call. = FALSE)
  }
  if (k <= 0) {
    stop("`k` must be a positive number.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Get multipliers for ATE calculation based on performance level
#' @noRd
#' @keywords internal
.get_ate_multipliers <- function(level) {

  # Desirable level multipliers (Fraser & Petersen 1993)
  # CV_A <= 0.50 * CV_I
  # Bias <= 0.25 * sqrt(CV_I^2 + CV_G^2)
  desirable_cv <- 0.50
  desirable_bias <- 0.25

  switch(level,
    "optimal" = list(
      imprecision = desirable_cv * 0.50,    # 0.25
      bias = desirable_bias * 0.50          # 0.125
    ),
    "desirable" = list(
      imprecision = desirable_cv,           # 0.50
      bias = desirable_bias                 # 0.25
    ),
    "minimum" = list(
      imprecision = desirable_cv * 1.50,    # 0.75
      bias = desirable_bias * 1.50          # 0.375
    )
  )
}
