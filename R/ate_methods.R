#' Print method for ate_specs objects
#'
#' @description
#' Displays a concise summary of allowable total error specifications
#' calculated from biological variation.
#'
#' @param x An object of class `ate_specs`.
#' @param digits Number of decimal places to display (default: 2).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)
#' print(ate)
#'
#' @export
print.ate_specs <- function(x, digits = 2, ...) {

  cat("\n")
  cat("Analytical Performance Specifications from Biological Variation\n")
  cat(strrep("-", 60), "\n\n")

  # Input parameters
  cat("Input:\n")
  cat(sprintf("  Within-subject CV (CV_I): %.*f%%\n", digits, x$input$cvi))
  if (!is.null(x$input$cvg)) {
    cat(sprintf("  Between-subject CV (CV_G): %.*f%%\n", digits, x$input$cvg))
  } else {
    cat("  Between-subject CV (CV_G): not provided\n")
  }
  cat(sprintf("  Performance level: %s\n", x$input$level))
  cat(sprintf("  Coverage factor (k): %.*f\n\n", 2, x$input$k))

  # Specifications
  cat("Specifications:\n")
  cat(sprintf("  Allowable imprecision (CV_A): %.*f%%\n",
              digits, x$specifications$allowable_cv))

  if (!is.null(x$specifications$allowable_bias)) {
    cat(sprintf("  Allowable bias: %.*f%%\n",
                digits, x$specifications$allowable_bias))
    cat(sprintf("  Total allowable error (TEa): %.*f%%\n",
                digits, x$specifications$tea))
  } else {
    cat("  Allowable bias: requires CV_G\n")
    cat("  Total allowable error (TEa): requires CV_G\n")
  }

  cat("\n")

  invisible(x)
}


#' Summary method for ate_specs objects
#'
#' @description
#' Provides a detailed summary of allowable total error specifications,
#' including the formulas used and all three performance tiers for comparison.
#'
#' @param object An object of class `ate_specs`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class `summary.ate_specs` containing detailed
#'   specification information, printed as a side effect.
#'
#' @examples
#' ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)
#' summary(ate)
#'
#' @export
summary.ate_specs <- function(object, ...) {

  x <- object

  cat("\n")
  cat("Analytical Performance Specifications - Detailed Summary\n")
  cat(strrep("=", 60), "\n\n")

  # Input section
  cat("Biological Variation Data:\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("  Within-subject CV (CV_I):  %.2f%%\n", x$input$cvi))
  if (!is.null(x$input$cvg)) {
    cat(sprintf("  Between-subject CV (CV_G): %.2f%%\n", x$input$cvg))
    total_bv <- sqrt(x$input$cvi^2 + x$input$cvg^2)
    cat(sprintf("  Total BV [sqrt(CV_I^2 + CV_G^2)]: %.2f%%\n", total_bv))
  } else {
    cat("  Between-subject CV (CV_G): not provided\n")
  }
  cat("\n")

  # Settings
  cat("Settings:\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("  Selected performance level: %s\n", x$input$level))
  cat(sprintf("  Coverage factor (k): %.2f\n", x$input$k))
  cat("\n")

  # Formulas
  cat("Formulas (Fraser & Petersen 1993):\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("  CV_A  = %.2f x CV_I\n", x$multipliers$imprecision))
  if (!is.null(x$input$cvg)) {
    cat(sprintf("  Bias  = %.3f x sqrt(CV_I^2 + CV_G^2)\n", x$multipliers$bias))
    cat(sprintf("  TEa   = k x CV_A + Bias\n"))
  }
  cat("\n")

  # Results for selected level
  cat(sprintf("Specifications (%s level):\n", x$input$level))
  cat(strrep("-", 60), "\n")
  cat(sprintf("  Allowable imprecision (CV_A): %.2f%%\n",
              x$specifications$allowable_cv))
  if (!is.null(x$specifications$allowable_bias)) {
    cat(sprintf("  Allowable bias:               %.2f%%\n",
                x$specifications$allowable_bias))
    cat(sprintf("  Total allowable error (TEa):  %.2f%%\n",
                x$specifications$tea))
  }
  cat("\n")

  # Comparison across all levels (if cvg provided)
  if (!is.null(x$input$cvg)) {
    cat("Comparison Across Performance Levels:\n")
    cat(strrep("-", 60), "\n")

    # Calculate for all levels
    levels <- c("optimal", "desirable", "minimum")
    comparison <- data.frame(
      Level = levels,
      CV_A = numeric(3),
      Bias = numeric(3),
      TEa = numeric(3),
      stringsAsFactors = FALSE
    )

    for (i in seq_along(levels)) {
      mult <- .get_ate_multipliers(levels[i])
      comparison$CV_A[i] <- mult$imprecision * x$input$cvi
      total_bv <- sqrt(x$input$cvi^2 + x$input$cvg^2)
      comparison$Bias[i] <- mult$bias * total_bv
      comparison$TEa[i] <- x$input$k * comparison$CV_A[i] + comparison$Bias[i]
    }

    # Mark selected level
    comparison$Level <- ifelse(
      comparison$Level == x$input$level,
      paste0(comparison$Level, " *"),
      comparison$Level
    )

    # Format for display
    cat(sprintf("  %-12s %8s %8s %8s\n", "Level", "CV_A", "Bias", "TEa"))
    cat(sprintf("  %-12s %8s %8s %8s\n", "-----", "----", "----", "---"))
    for (i in 1:3) {
      cat(sprintf("  %-12s %7.2f%% %7.2f%% %7.2f%%\n",
                  comparison$Level[i],
                  comparison$CV_A[i],
                  comparison$Bias[i],
                  comparison$TEa[i]))
    }
    cat("\n  * Selected level\n")
  }

  cat("\n")

  # Data source note
  cat("Data Source:\n")
  cat(strrep("-", 60), "\n")
  cat("  Biological variation values should be obtained from the EFLM\n")
  cat("  Biological Variation Database: https://biologicalvariation.eu/\n")
  cat("\n")

  # Return summary object invisibly
  invisible(structure(
    list(
      specifications = x$specifications,
      input = x$input,
      multipliers = x$multipliers
    ),
    class = "summary.ate_specs"
  ))
}
