#' Calculate Six Sigma Metric for Analytical Performance
#'
#' @description
#' Calculates the sigma metric, which quantifies analytical performance
#' in terms of the number of standard deviations between observed performance
#' and the allowable total error limit. Higher sigma values indicate better
#' performance and lower defect rates.
#'
#' @param bias Numeric. Observed bias (systematic error), expressed as a
#'   percentage or in the same units as `tea`.
#' @param cv Numeric. Observed coefficient of variation (imprecision),
#'   expressed as a percentage.
#' @param tea Numeric. Total allowable error specification, expressed as a
#'   percentage or in the same units as `bias`.
#'
#' @return An object of class `c("sigma_metric", "valytics_ate", "valytics_result")`,
#'   which is a list containing:
#'
#'   \describe{
#'     \item{sigma}{Numeric. The calculated sigma metric value.}
#'     \item{input}{List with input parameters:
#'       \itemize{
#'         \item `bias`: Observed bias
#'         \item `cv`: Observed CV
#'         \item `tea`: Total allowable error
#'       }
#'     }
#'     \item{interpretation}{List with performance interpretation:
#'       \itemize{
#'         \item `category`: Performance category (e.g., "World Class", "Good")
#'         \item `defect_rate`: Approximate defect rate per million
#'       }
#'     }
#'   }
#'
#' @details
#' The sigma metric is calculated as:
#'
#' \deqn{\sigma = \frac{TEa - |Bias|}{CV}}{Sigma = (TEa - |Bias|) / CV}
#'
#' Where:
#' \itemize{
#'   \item TEa = Total allowable error (quality specification)
#'   \item Bias = Observed systematic error (absolute value used)
#'   \item CV = Observed coefficient of variation
#' }
#'
#' **Interpretation Guidelines:**
#'
#' The sigma metric provides a standardized way to assess method performance:
#' \itemize{
#'   \item **>= 6 sigma**: World class performance (<3.4 defects per million)
#'   \item **>= 5 sigma**: Excellent performance (~230 defects per million)
#'   \item **>= 4 sigma**: Good performance (~6,200 defects per million)
#'   \item **>= 3 sigma**: Marginal performance (~66,800 defects per million)
#'   \item **< 3 sigma**: Poor performance (unacceptable defect rates)
#' }
#'
#' Note: These defect rates assume a 1.5 sigma shift (industry standard for
#' long-term process variation).
#'
#' @section Clinical Laboratory Context:
#' In clinical laboratories, a sigma metric of 4 or higher is generally
#' considered acceptable for routine testing, while 6 sigma is the gold
#' standard. Methods with sigma < 3 require stringent QC procedures and
#' may not be suitable for clinical use without improvement.
#'
#' @references
#' Westgard JO, Westgard SA (2006). The quality of laboratory testing today:
#' an assessment of sigma metrics for analytic quality using performance data
#' from proficiency testing surveys and the CLIA criteria for acceptable
#' performance. \emph{American Journal of Clinical Pathology}, 125(3):343-354.
#'
#' Westgard JO (2008). \emph{Basic Method Validation} (3rd ed.).
#' Westgard QC, Inc.
#'
#' @seealso
#' [ate_from_bv()] for calculating TEa from biological variation,
#' [ate_assessment()] for comprehensive performance assessment
#'
#' @examples
#' # Basic sigma calculation
#' sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)
#' sm
#'
#' # World-class performance example
#' sm_excellent <- sigma_metric(bias = 0.5, cv = 1.0, tea = 8)
#' sm_excellent
#'
#' # Marginal performance example
#' sm_marginal <- sigma_metric(bias = 3.0, cv = 3.0, tea = 12)
#' sm_marginal
#'
#' # Using with ate_from_bv() for glucose
#' ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)
#' # Assume observed bias = 1.5%, CV = 2.5%
#' sm <- sigma_metric(bias = 1.5, cv = 2.5, tea = ate$specifications$tea)
#' sm
#'
#' # Access the sigma value directly
#' sm$sigma
#'
#' @export
sigma_metric <- function(bias, cv, tea) {

 # Input validation ----
 .validate_sigma_input(bias, cv, tea)

  # Calculate sigma ----
  # Sigma = (TEa - |Bias|) / CV
  sigma <- (tea - abs(bias)) / cv

  # Determine interpretation ----
  interpretation <- .interpret_sigma(sigma)

  # Construct output object ----
  structure(
    list(
      sigma = sigma,
      input = list(
        bias = bias,
        cv = cv,
        tea = tea
      ),
      interpretation = interpretation
    ),
    class = c("sigma_metric", "valytics_ate", "valytics_result")
  )
}


# Helper Functions ----

#' Validate input for sigma_metric
#' @noRd
#' @keywords internal
.validate_sigma_input <- function(bias, cv, tea) {

  # Check bias (can be negative, zero, or positive)
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

  invisible(TRUE)
}


#' Interpret sigma metric value
#' @noRd
#' @keywords internal
.interpret_sigma <- function(sigma) {

  # Defect rates based on 1.5 sigma shift (industry standard)
  # These are approximate values for the long-term process
  if (sigma >= 6) {
    category <- "World Class"
    defect_rate <- 3.4
  } else if (sigma >= 5) {
    category <- "Excellent"
    defect_rate <- 230
  } else if (sigma >= 4) {
    category <- "Good"
    defect_rate <- 6210
  } else if (sigma >= 3) {
    category <- "Marginal"
    defect_rate <- 66800
  } else if (sigma >= 2) {
    category <- "Poor"
    defect_rate <- 308500
  } else if (sigma >= 1) {
    category <- "Unacceptable"
    defect_rate <- 690000
  } else {
    category <- "Unacceptable"
    defect_rate <- NA_real_  # Beyond typical tables
  }

  list(
    category = category,
    defect_rate = defect_rate
  )
}


#' Print method for sigma_metric objects
#'
#' @description
#' Displays a concise summary of the sigma metric calculation.
#'
#' @param x An object of class `sigma_metric`.
#' @param digits Number of decimal places to display (default: 2).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)
#' print(sm)
#'
#' @export
print.sigma_metric <- function(x, digits = 2, ...) {

  cat("\n")
  cat("Six Sigma Metric\n")
  cat(strrep("-", 40), "\n\n")

  # Input
  cat("Input:\n")
  cat(sprintf("  Observed bias: %.*f%%\n", digits, x$input$bias))
  cat(sprintf("  Observed CV: %.*f%%\n", digits, x$input$cv))
  cat(sprintf("  Total allowable error (TEa): %.*f%%\n\n", digits, x$input$tea))

  # Result
  cat("Result:\n")
  cat(sprintf("  Sigma: %.*f\n", digits, x$sigma))
  cat(sprintf("  Performance: %s\n", x$interpretation$category))

  if (!is.na(x$interpretation$defect_rate)) {
    if (x$interpretation$defect_rate < 100) {
      cat(sprintf("  Defect rate: ~%.1f per million\n",
                  x$interpretation$defect_rate))
    } else {
      cat(sprintf("  Defect rate: ~%s per million\n",
                  format(round(x$interpretation$defect_rate), big.mark = ",")))
    }
  }

  cat("\n")

  invisible(x)
}


#' Summary method for sigma_metric objects
#'
#' @description
#' Provides a detailed summary of the sigma metric calculation,
#' including the formula and interpretation scale.
#'
#' @param object An object of class `sigma_metric`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the object.
#'
#' @examples
#' sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)
#' summary(sm)
#'
#' @export
summary.sigma_metric <- function(object, ...) {

  x <- object

  cat("\n")
  cat("Six Sigma Metric - Detailed Summary\n")
  cat(strrep("=", 50), "\n\n")

  # Formula
  cat("Formula:\n")
  cat(strrep("-", 50), "\n")
  cat("  Sigma = (TEa - |Bias|) / CV\n")
  cat(sprintf("  Sigma = (%.2f - |%.2f|) / %.2f\n",
              x$input$tea, x$input$bias, x$input$cv))
  cat(sprintf("  Sigma = (%.2f - %.2f) / %.2f\n",
              x$input$tea, abs(x$input$bias), x$input$cv))
  cat(sprintf("  Sigma = %.2f / %.2f\n",
              x$input$tea - abs(x$input$bias), x$input$cv))
  cat(sprintf("  Sigma = %.2f\n\n", x$sigma))

  # Result and interpretation
  cat("Result:\n")
  cat(strrep("-", 50), "\n")
  cat(sprintf("  Sigma metric: %.2f\n", x$sigma))
  cat(sprintf("  Performance category: %s\n", x$interpretation$category))
  if (!is.na(x$interpretation$defect_rate)) {
    cat(sprintf("  Expected defect rate: ~%s per million\n\n",
                format(round(x$interpretation$defect_rate), big.mark = ",")))
  } else {
    cat("  Expected defect rate: Beyond standard tables\n\n")
  }

  # Interpretation scale
  cat("Sigma Scale Reference:\n")
  cat(strrep("-", 50), "\n")
  cat("  Sigma    Category        Defects/Million\n")
  cat("  ------   -------------   ---------------\n")

  scale_data <- data.frame(
    sigma = c(">= 6", ">= 5", ">= 4", ">= 3", ">= 2", "< 2"),
    category = c("World Class", "Excellent", "Good",
                 "Marginal", "Poor", "Unacceptable"),
    defects = c("3.4", "230", "6,210", "66,800", "308,500", "> 690,000"),
    stringsAsFactors = FALSE
  )

  # Mark current level
  current_idx <- which(c(x$sigma >= 6, x$sigma >= 5 & x$sigma < 6,
                         x$sigma >= 4 & x$sigma < 5, x$sigma >= 3 & x$sigma < 4,
                         x$sigma >= 2 & x$sigma < 3, x$sigma < 2))

  for (i in 1:6) {
    marker <- if (i == current_idx) " *" else "  "
    cat(sprintf("  %-6s   %-13s   %15s%s\n",
                scale_data$sigma[i],
                scale_data$category[i],
                scale_data$defects[i],
                marker))
  }
  cat("\n  * Current performance level\n")

  cat("\n")
  cat("Note: Defect rates assume 1.5 sigma long-term shift.\n")
  cat("\n")

  invisible(x)
}
