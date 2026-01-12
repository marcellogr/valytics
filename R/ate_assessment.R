#' Assess Analytical Performance Against Allowable Total Error
#'
#' @description
#' Evaluates observed analytical performance (bias and imprecision) against
#' allowable total error specifications. Provides pass/fail assessment for
#' individual components and overall method acceptability, along with the
#' sigma metric.
#'
#' @param bias Numeric. Observed bias (systematic error), expressed as a
#'   percentage.
#' @param cv Numeric. Observed coefficient of variation (imprecision),
#'   expressed as a percentage.
#' @param tea Numeric. Total allowable error specification. Can be provided
#'   directly or will be calculated if `allowable_bias` and `allowable_cv`
#'   are provided with `k`.
#' @param allowable_bias Numeric. Allowable bias specification (optional).
#'   If provided, enables individual bias assessment.
#' @param allowable_cv Numeric. Allowable imprecision specification (optional).
#'   If provided, enables individual CV assessment.
#' @param k Numeric. Coverage factor for TEa calculation when using component
#'   specifications (default: 1.65).
#'
#' @return An object of class `c("ate_assessment", "valytics_ate", "valytics_result")`,
#'   which is a list containing:
#'
#'   \describe{
#'     \item{assessment}{List with pass/fail results:
#'       \itemize{
#'         \item `bias_acceptable`: Logical; TRUE if |bias| <= allowable_bias
#'         \item `cv_acceptable`: Logical; TRUE if cv <= allowable_cv
#'         \item `tea_acceptable`: Logical; TRUE if observed TE <= TEa
#'         \item `overall`: Logical; TRUE if method meets specifications
#'       }
#'     }
#'     \item{observed}{List with observed performance:
#'       \itemize{
#'         \item `bias`: Observed bias
#'         \item `cv`: Observed CV
#'         \item `te`: Observed total error (k * CV + |Bias|)
#'       }
#'     }
#'     \item{specifications}{List with allowable specifications:
#'       \itemize{
#'         \item `allowable_bias`: Allowable bias (or NULL)
#'         \item `allowable_cv`: Allowable CV (or NULL)
#'         \item `tea`: Total allowable error
#'       }
#'     }
#'     \item{sigma}{List with sigma metric results:
#'       \itemize{
#'         \item `value`: Sigma metric value
#'         \item `category`: Performance category
#'       }
#'     }
#'     \item{settings}{List with settings:
#'       \itemize{
#'         \item `k`: Coverage factor used
#'       }
#'     }
#'   }
#'
#' @details
#' The assessment evaluates method performance at multiple levels:
#'
#' **Component Assessment** (if specifications provided):
#' \itemize{
#'   \item Bias: Pass if |observed bias| <= allowable bias
#'   \item CV: Pass if observed CV <= allowable CV
#' }
#'
#' **Total Error Assessment**:
#' \itemize{
#'   \item Observed TE = k * CV + |Bias| (linear model)
#'   \item Pass if observed TE <= TEa
#' }
#'
#' **Sigma Metric**:
#' \itemize{
#'   \item Sigma = (TEa - |Bias|) / CV
#'   \item Provides quality rating from "World Class" to "Unacceptable"
#' }
#'
#' @section Overall Assessment:
#' The overall assessment is determined as follows:
#' \itemize{
#'   \item If only TEa is provided: based on total error assessment
#'   \item If component specs provided: all components must pass
#'   \item Sigma >= 3 is generally considered minimum acceptable
#' }
#'
#' @references
#' Westgard JO (2008). \emph{Basic Method Validation} (3rd ed.).
#' Westgard QC, Inc.
#'
#' Fraser CG (2001). \emph{Biological Variation: From Principles to Practice}.
#' AACC Press.
#'
#' @seealso
#' [ate_from_bv()] for calculating specifications from biological variation,
#' [sigma_metric()] for sigma calculation details
#'
#' @examples
#' # Basic assessment with TEa only
#' assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)
#' assess
#'
#' # Assessment with all component specifications
#' assess_full <- ate_assessment(
#'   bias = 1.5,
#'   cv = 2.5,
#'   tea = 10,
#'   allowable_bias = 3.0,
#'   allowable_cv = 4.0
#' )
#' assess_full
#'
#' # Using specifications from ate_from_bv()
#' specs <- ate_from_bv(cvi = 5.6, cvg = 7.5)
#' assess <- ate_assessment(
#'   bias = 1.5,
#'   cv = 2.5,
#'   tea = specs$specifications$tea,
#'   allowable_bias = specs$specifications$allowable_bias,
#'   allowable_cv = specs$specifications$allowable_cv
#' )
#' summary(assess)
#'
#' # Check if method passes
#' assess$assessment$overall
#'
#' @export
ate_assessment <- function(bias,
                           cv,
                           tea,
                           allowable_bias = NULL,
                           allowable_cv = NULL,
                           k = 1.65) {

  # Input validation ----
  .validate_assessment_input(bias, cv, tea, allowable_bias, allowable_cv, k)

  # Calculate observed total error ----
  observed_te <- k * cv + abs(bias)

  # Perform assessments ----

  # Bias assessment (if specification provided)
  if (!is.null(allowable_bias)) {
    bias_acceptable <- abs(bias) <= allowable_bias
  } else {
    bias_acceptable <- NA
  }

  # CV assessment (if specification provided)
  if (!is.null(allowable_cv)) {
    cv_acceptable <- cv <= allowable_cv
  } else {
    cv_acceptable <- NA
  }

  # Total error assessment
  tea_acceptable <- observed_te <= tea

  # Overall assessment
  # If component specs provided, all must pass
 # If only TEa, use TEa assessment
  if (!is.null(allowable_bias) && !is.null(allowable_cv)) {
    overall <- bias_acceptable && cv_acceptable && tea_acceptable
  } else {
    overall <- tea_acceptable
  }

  # Calculate sigma metric ----
  sigma_value <- (tea - abs(bias)) / cv
  sigma_interp <- .interpret_sigma(sigma_value)

  # Construct output object ----
  structure(
    list(
      assessment = list(
        bias_acceptable = bias_acceptable,
        cv_acceptable = cv_acceptable,
        tea_acceptable = tea_acceptable,
        overall = overall
      ),
      observed = list(
        bias = bias,
        cv = cv,
        te = observed_te
      ),
      specifications = list(
        allowable_bias = allowable_bias,
        allowable_cv = allowable_cv,
        tea = tea
      ),
      sigma = list(
        value = sigma_value,
        category = sigma_interp$category
      ),
      settings = list(
        k = k
      )
    ),
    class = c("ate_assessment", "valytics_ate", "valytics_result")
  )
}


# Helper Functions ----
# =============================================================================

#' Validate input for ate_assessment
#' @noRd
#' @keywords internal
.validate_assessment_input <- function(bias, cv, tea,
                                        allowable_bias, allowable_cv, k) {

  # Check bias (can be any numeric, including negative)
  if (length(bias) != 1 || !is.numeric(bias) || is.na(bias)) {
    stop("`bias` must be a single numeric value.", call. = FALSE)
  }

  # Check cv (must be positive)
  if (length(cv) != 1 || !is.numeric(cv) || is.na(cv)) {
    stop("`cv` must be a single numeric value.", call. = FALSE)
  }
  if (cv <= 0) {
    stop("`cv` must be a positive number.", call. = FALSE)
  }

  # Check tea (must be positive)
  if (length(tea) != 1 || !is.numeric(tea) || is.na(tea)) {
    stop("`tea` must be a single numeric value.", call. = FALSE)
  }
  if (tea <= 0) {
    stop("`tea` must be a positive number.", call. = FALSE)
  }

  # Check allowable_bias (if provided)
  if (!is.null(allowable_bias)) {
    if (length(allowable_bias) != 1 || !is.numeric(allowable_bias) ||
        is.na(allowable_bias)) {
      stop("`allowable_bias` must be a single numeric value or NULL.",
           call. = FALSE)
    }
    if (allowable_bias <= 0) {
      stop("`allowable_bias` must be a positive number.", call. = FALSE)
    }
  }

  # Check allowable_cv (if provided)
  if (!is.null(allowable_cv)) {
    if (length(allowable_cv) != 1 || !is.numeric(allowable_cv) ||
        is.na(allowable_cv)) {
      stop("`allowable_cv` must be a single numeric value or NULL.",
           call. = FALSE)
    }
    if (allowable_cv <= 0) {
      stop("`allowable_cv` must be a positive number.", call. = FALSE)
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


#' Print method for ate_assessment objects
#'
#' @description
#' Displays a concise summary of the performance assessment.
#'
#' @param x An object of class `ate_assessment`.
#' @param digits Number of decimal places to display (default: 2).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)
#' print(assess)
#'
#' @export
print.ate_assessment <- function(x, digits = 2, ...) {

  cat("\n")
  cat("Analytical Performance Assessment\n")
  cat(strrep("-", 50), "\n\n")

  # Overall result banner
  if (x$assessment$overall) {
    cat("  >>> METHOD ACCEPTABLE <<<\n\n")
  } else {
    cat("  >>> METHOD NOT ACCEPTABLE <<<\n\n")
  }

  # Observed vs Allowable comparison
  cat("Performance Summary:\n")
  cat(sprintf("  %-20s %10s %10s %10s\n",
              "Parameter", "Observed", "Allowable", "Status"))
  cat(sprintf("  %-20s %10s %10s %10s\n",
              strrep("-", 20), strrep("-", 10), strrep("-", 10), strrep("-", 10)))

  # Bias
  if (!is.na(x$assessment$bias_acceptable)) {
    status_bias <- if (x$assessment$bias_acceptable) "PASS" else "FAIL"
    cat(sprintf("  %-20s %9.*f%% %9.*f%% %10s\n",
                "Bias",
                digits, abs(x$observed$bias),
                digits, x$specifications$allowable_bias,
                status_bias))
  } else {
    cat(sprintf("  %-20s %9.*f%% %10s %10s\n",
                "Bias",
                digits, abs(x$observed$bias),
                "---",
                "---"))
  }

  # CV
  if (!is.na(x$assessment$cv_acceptable)) {
    status_cv <- if (x$assessment$cv_acceptable) "PASS" else "FAIL"
    cat(sprintf("  %-20s %9.*f%% %9.*f%% %10s\n",
                "CV (Imprecision)",
                digits, x$observed$cv,
                digits, x$specifications$allowable_cv,
                status_cv))
  } else {
    cat(sprintf("  %-20s %9.*f%% %10s %10s\n",
                "CV (Imprecision)",
                digits, x$observed$cv,
                "---",
                "---"))
  }

  # Total Error
  status_te <- if (x$assessment$tea_acceptable) "PASS" else "FAIL"
  cat(sprintf("  %-20s %9.*f%% %9.*f%% %10s\n",
              "Total Error",
              digits, x$observed$te,
              digits, x$specifications$tea,
              status_te))

  cat("\n")

  # Sigma metric
  cat(sprintf("Sigma Metric: %.*f (%s)\n",
              digits, x$sigma$value, x$sigma$category))

  cat("\n")

  invisible(x)
}


#' Summary method for ate_assessment objects
#'
#' @description
#' Provides a detailed summary of the performance assessment,
#' including calculations and interpretation guidance.
#'
#' @param object An object of class `ate_assessment`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the object.
#'
#' @examples
#' assess <- ate_assessment(
#'   bias = 1.5, cv = 2.5, tea = 10,
#'   allowable_bias = 3.0, allowable_cv = 4.0
#' )
#' summary(assess)
#'
#' @export
summary.ate_assessment <- function(object, ...) {

  x <- object

  cat("\n")
  cat("Analytical Performance Assessment - Detailed Summary\n")
  cat(strrep("=", 60), "\n\n")

  # Overall result
  cat("Overall Result: ")
  if (x$assessment$overall) {
    cat("METHOD ACCEPTABLE\n\n")
  } else {
    cat("METHOD NOT ACCEPTABLE\n\n")
  }

  # Observed performance
  cat("Observed Performance:\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("  Bias: %.2f%%\n", x$observed$bias))
  cat(sprintf("  CV (Imprecision): %.2f%%\n", x$observed$cv))
  cat(sprintf("  Total Error (k=%.2f): %.2f%%\n",
              x$settings$k, x$observed$te))
  cat(sprintf("    [TE = %.2f x %.2f + |%.2f| = %.2f]\n\n",
              x$settings$k, x$observed$cv, x$observed$bias, x$observed$te))

  # Specifications
  cat("Allowable Specifications:\n")
  cat(strrep("-", 60), "\n")
  if (!is.null(x$specifications$allowable_bias)) {
    cat(sprintf("  Allowable Bias: %.2f%%\n", x$specifications$allowable_bias))
  } else {
    cat("  Allowable Bias: not specified\n")
  }
  if (!is.null(x$specifications$allowable_cv)) {
    cat(sprintf("  Allowable CV: %.2f%%\n", x$specifications$allowable_cv))
  } else {
    cat("  Allowable CV: not specified\n")
  }
  cat(sprintf("  Total Allowable Error (TEa): %.2f%%\n\n", x$specifications$tea))

  # Component assessments
  cat("Component Assessment:\n")
  cat(strrep("-", 60), "\n")

  # Bias
  if (!is.na(x$assessment$bias_acceptable)) {
    bias_result <- if (x$assessment$bias_acceptable) "PASS" else "FAIL"
    bias_margin <- x$specifications$allowable_bias - abs(x$observed$bias)
    cat(sprintf("  Bias: %s (margin: %+.2f%%)\n", bias_result, bias_margin))
    cat(sprintf("    |%.2f| %s %.2f\n",
                x$observed$bias,
                if (x$assessment$bias_acceptable) "<=" else ">",
                x$specifications$allowable_bias))
  } else {
    cat("  Bias: not assessed (no specification provided)\n")
  }

  # CV
  if (!is.na(x$assessment$cv_acceptable)) {
    cv_result <- if (x$assessment$cv_acceptable) "PASS" else "FAIL"
    cv_margin <- x$specifications$allowable_cv - x$observed$cv
    cat(sprintf("  CV: %s (margin: %+.2f%%)\n", cv_result, cv_margin))
    cat(sprintf("    %.2f %s %.2f\n",
                x$observed$cv,
                if (x$assessment$cv_acceptable) "<=" else ">",
                x$specifications$allowable_cv))
  } else {
    cat("  CV: not assessed (no specification provided)\n")
  }

  # Total Error
  te_result <- if (x$assessment$tea_acceptable) "PASS" else "FAIL"
  te_margin <- x$specifications$tea - x$observed$te
  cat(sprintf("  Total Error: %s (margin: %+.2f%%)\n", te_result, te_margin))
  cat(sprintf("    %.2f %s %.2f\n\n",
              x$observed$te,
              if (x$assessment$tea_acceptable) "<=" else ">",
              x$specifications$tea))

  # Sigma metric
  cat("Sigma Metric:\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("  Sigma = (TEa - |Bias|) / CV\n"))
  cat(sprintf("  Sigma = (%.2f - %.2f) / %.2f = %.2f\n",
              x$specifications$tea, abs(x$observed$bias),
              x$observed$cv, x$sigma$value))
  cat(sprintf("  Category: %s\n\n", x$sigma$category))

  # Sigma scale
  cat("  Sigma Scale:\n")
  cat("    >= 6: World Class | >= 5: Excellent | >= 4: Good\n")
  cat("    >= 3: Marginal    | >= 2: Poor      | < 2: Unacceptable\n")

  cat("\n")

  # Recommendations
  cat("Interpretation:\n")
  cat(strrep("-", 60), "\n")

  if (x$assessment$overall && x$sigma$value >= 4) {
    cat("  Method performance is acceptable with good quality margin.\n")
  } else if (x$assessment$overall && x$sigma$value >= 3) {
    cat("  Method meets minimum specifications but has limited margin.\n")
    cat("  Consider implementing stringent QC procedures.\n")
  } else if (x$assessment$overall) {
    cat("  Method technically passes but sigma < 3 indicates high risk.\n")
    cat("  Strongly recommend method improvement or enhanced QC.\n")
  } else {
    cat("  Method does not meet specifications.\n")
    if (!is.na(x$assessment$bias_acceptable) && !x$assessment$bias_acceptable) {
      cat("  - Bias exceeds allowable limit: consider recalibration.\n")
    }
    if (!is.na(x$assessment$cv_acceptable) && !x$assessment$cv_acceptable) {
      cat("  - Imprecision exceeds allowable limit: investigate sources.\n")
    }
    if (!x$assessment$tea_acceptable) {
      cat("  - Total error exceeds TEa: method requires improvement.\n")
    }
  }

  cat("\n")

  invisible(x)
}
