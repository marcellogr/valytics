#' Print method for deming_regression objects
#'
#' @description
#' Displays a concise summary of Deming regression results, including
#' slope and intercept estimates with confidence intervals.
#'
#' @param x An object of class `deming_regression`.
#' @param digits Number of significant digits to display (default: 3).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' set.seed(42)
#' true_vals <- rnorm(50, 100, 20)
#' method_a <- true_vals + rnorm(50, sd = 5)
#' method_b <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)
#' dm <- deming_regression(method_a, method_b)
#' print(dm)
#'
#' @seealso [summary.deming_regression()] for detailed output
#' @export
print.deming_regression <- function(x, digits = 3, ...) {

  cat("\nDeming Regression\n")
  cat(strrep("-", 40), "\n")

  # Sample size
  cat(sprintf("n = %d paired observations", x$input$n))
  if (x$input$n_excluded > 0) {
    cat(sprintf(" (%d excluded)", x$input$n_excluded))
  }
  cat("\n\n")

  # Settings
  ci_pct <- paste0(x$settings$conf_level * 100, "%")
  ci_method_str <- if (x$settings$ci_method == "jackknife") {
    "Jackknife"
  } else {
    sprintf("Bootstrap BCa (n = %d)", x$settings$boot_n)
  }

  cat(sprintf("Error ratio (lambda): %.*f\n", digits, x$settings$error_ratio))
  cat(sprintf("CI method: %s\n", ci_method_str))
  cat(sprintf("Confidence level: %s\n\n", ci_pct))

  # Regression equation
  cat("Regression equation:\n")
  cat(sprintf("  %s = %.3f + %.3f * %s\n\n",
              x$input$var_names["y"],
              x$results$intercept,
              x$results$slope,
              x$input$var_names["x"]))

  # Results
  cat("Results:\n")

  # Intercept
  cat(sprintf("  Intercept: %.*f (SE = %.*f)\n",
              digits, x$results$intercept,
              digits, x$results$intercept_se))
  if (!any(is.na(x$results$intercept_ci))) {
    cat(sprintf("    %s CI: [%.*f, %.*f]\n",
                ci_pct, digits, x$results$intercept_ci["lower"],
                digits, x$results$intercept_ci["upper"]))

    # Check if 0 is in CI
    if (x$results$intercept_ci["lower"] <= 0 &&
        x$results$intercept_ci["upper"] >= 0) {
      cat("    (includes 0: no significant constant bias)\n")
    } else {
      cat("    (excludes 0: significant constant bias)\n")
    }
  }

  cat("\n")

  # Slope
  cat(sprintf("  Slope: %.*f (SE = %.*f)\n",
              digits, x$results$slope,
              digits, x$results$slope_se))
  if (!any(is.na(x$results$slope_ci))) {
    cat(sprintf("    %s CI: [%.*f, %.*f]\n",
                ci_pct, digits, x$results$slope_ci["lower"],
                digits, x$results$slope_ci["upper"]))

    # Check if 1 is in CI
    if (x$results$slope_ci["lower"] <= 1 &&
        x$results$slope_ci["upper"] >= 1) {
      cat("    (includes 1: no significant proportional bias)\n")
    } else {
      cat("    (excludes 1: significant proportional bias)\n")
    }
  }

  cat("\n")

  invisible(x)
}


#' Summary method for deming_regression objects
#'
#' @description
#' Provides a detailed summary of Deming regression results, including
#' regression coefficients, confidence intervals, standard errors,
#' and interpretation guidance.
#'
#' @param object An object of class `deming_regression`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns a list with summary statistics.
#'
#' @details
#' The summary includes:
#' \itemize{
#'   \item Regression coefficients with standard errors and confidence intervals
#'   \item Interpretation of slope and intercept CIs
#'   \item Method agreement conclusion
#'   \item Residual summary statistics
#' }
#'
#' @examples
#' set.seed(42)
#' true_vals <- rnorm(50, 100, 20)
#' method_a <- true_vals + rnorm(50, sd = 5)
#' method_b <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)
#' dm <- deming_regression(method_a, method_b)
#' summary(dm)
#'
#' @seealso [print.deming_regression()] for concise output
#' @export
summary.deming_regression <- function(object, ...) {

  x <- object
  ci_pct <- paste0(x$settings$conf_level * 100, "%")

  cat("\n")
  cat("Deming Regression - Detailed Summary\n")
  cat(strrep("=", 50), "\n\n")

  # Input summary ----
  cat("Data:\n")
  cat(sprintf("  X variable: %s\n", x$input$var_names["x"]))
  cat(sprintf("  Y variable: %s\n", x$input$var_names["y"]))
  cat(sprintf("  Sample size: %d\n", x$input$n))
  if (x$input$n_excluded > 0) {
    cat(sprintf("  Excluded (NA): %d\n", x$input$n_excluded))
  }
  cat("\n")

  # Settings ----
  cat("Settings:\n")
  cat(sprintf("  Error ratio (lambda): %.4f\n", x$settings$error_ratio))
  if (x$settings$error_ratio == 1) {
    cat("    (orthogonal regression - equal error variances assumed)\n")
  } else if (x$settings$error_ratio > 1) {
    cat("    (Y has higher error variance than X)\n")
  } else {
    cat("    (X has higher error variance than Y)\n")
  }
  cat(sprintf("  Confidence level: %s\n", ci_pct))
  ci_method_str <- if (x$settings$ci_method == "jackknife") {
    "Jackknife"
  } else {
    sprintf("Bootstrap BCa (n = %d)", x$settings$boot_n)
  }
  cat(sprintf("  CI method: %s\n", ci_method_str))
  cat("\n")

  # Regression coefficients ----
  cat("Regression Coefficients:\n")
  cat(strrep("-", 50), "\n")

  # Create coefficient table
  coef_table <- data.frame(
    Estimate = c(x$results$intercept, x$results$slope),
    SE = c(x$results$intercept_se, x$results$slope_se),
    Lower = c(x$results$intercept_ci["lower"], x$results$slope_ci["lower"]),
    Upper = c(x$results$intercept_ci["upper"], x$results$slope_ci["upper"]),
    row.names = c("Intercept", "Slope")
  )
  names(coef_table) <- c("Estimate", "Std. Error",
                         paste0(ci_pct, " Lower"),
                         paste0(ci_pct, " Upper"))

  print(round(coef_table, 4))
  cat("\n")

  # Regression equation
  cat("Regression equation:\n")
  cat(sprintf("  %s = %.4f + %.4f * %s\n\n",
              x$input$var_names["y"],
              x$results$intercept,
              x$results$slope,
              x$input$var_names["x"]))

  # Interpretation ----
  cat("Interpretation:\n")
  cat(strrep("-", 50), "\n")

  # Intercept interpretation
  intercept_in_ci <- !any(is.na(x$results$intercept_ci)) &&
    x$results$intercept_ci["lower"] <= 0 &&
    x$results$intercept_ci["upper"] >= 0

  if (is.na(x$results$intercept_ci["lower"])) {
    cat("  Intercept: CI not available\n")
  } else if (intercept_in_ci) {
    cat("  Intercept: CI includes 0\n")
    cat("    -> No significant constant (additive) bias\n")
  } else {
    direction <- if (x$results$intercept > 0) "positive" else "negative"
    cat(sprintf("  Intercept: CI excludes 0 (%.3f to %.3f)\n",
                x$results$intercept_ci["lower"],
                x$results$intercept_ci["upper"]))
    cat(sprintf("    -> Significant %s constant bias of %.3f\n",
                direction, x$results$intercept))
  }

  # Slope interpretation
  slope_in_ci <- !any(is.na(x$results$slope_ci)) &&
    x$results$slope_ci["lower"] <= 1 &&
    x$results$slope_ci["upper"] >= 1

  if (is.na(x$results$slope_ci["lower"])) {
    cat("  Slope: CI not available\n")
  } else if (slope_in_ci) {
    cat("  Slope: CI includes 1\n")
    cat("    -> No significant proportional (multiplicative) bias\n")
  } else {
    pct_diff <- (x$results$slope - 1) * 100
    cat(sprintf("  Slope: CI excludes 1 (%.3f to %.3f)\n",
                x$results$slope_ci["lower"],
                x$results$slope_ci["upper"]))
    cat(sprintf("    -> Significant proportional bias of %.1f%%\n", pct_diff))
  }

  cat("\n")

  # Overall conclusion ----
  cat("Conclusion:\n")
  cat(strrep("-", 50), "\n")

  if (intercept_in_ci && slope_in_ci) {
    cat("  The two methods are EQUIVALENT within the measured range.\n")
    cat("  No systematic differences detected.\n")
  } else {
    cat("  The two methods show SYSTEMATIC DIFFERENCES:\n")

    if (!is.na(x$results$intercept_ci["lower"]) && !intercept_in_ci) {
      cat(sprintf("    - Constant bias: %.3f %s\n",
                  abs(x$results$intercept),
                  x$input$var_names["y"]))
    }
    if (!is.na(x$results$slope_ci["lower"]) && !slope_in_ci) {
      cat(sprintf("    - Proportional bias: %.1f%%\n",
                  (x$results$slope - 1) * 100))
    }
  }
  cat("\n")

  # Residual summary ----
  cat("Residuals (perpendicular):\n")
  cat(strrep("-", 50), "\n")
  resid_summary <- summary(x$results$residuals)
  print(resid_summary)
  cat("\n")

  # Comparison note ----
  cat("Note on error ratio:\n")
  cat(strrep("-", 50), "\n")
  if (x$settings$error_ratio == 1) {
    cat("  Using lambda = 1 (orthogonal regression).\n")
    cat("  This assumes both methods have equal measurement error variance.\n")
    cat("  If this assumption is violated, consider specifying 'error_ratio'\n")
    cat("  based on replicate measurements or known precision data.\n")
  } else {
    cat(sprintf("  Using lambda = %.3f (specified error ratio).\n",
                x$settings$error_ratio))
    cat("  The Y method is assumed to have variance ratio lambda times\n")
    cat("  the X method's error variance.\n")
  }
  cat("\n")

  # Return summary statistics invisibly
  invisible(list(
    coefficients = coef_table,
    intercept_includes_zero = intercept_in_ci,
    slope_includes_one = slope_in_ci,
    methods_equivalent = intercept_in_ci && slope_in_ci
  ))
}
