#' Print method for linearity_study objects
#'
#' @description
#' Displays a concise summary of linearity evaluation results, including
#' the best-fitting model, maximum deviation, and ADL assessment.
#'
#' @param x An object of class `linearity_study`.
#' @param digits Number of significant digits to display (default: 3).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' set.seed(42)
#' target <- rep(c(10, 50, 100, 200, 300, 400, 500), each = 4)
#' measured <- target * rnorm(28, 1.0, 0.02) + rnorm(28, 0, 2)
#' lin <- linearity_study(target, measured)
#' print(lin)
#'
#' @export
print.linearity_study <- function(x, digits = 3, ...) {
  
  cat("\n")
  cat("Linearity Evaluation\n")
  cat(strrep("-", 40), "\n")
  
  # Sample info
  cat(sprintf("n = %d observations across %d levels",
              x$input$n, x$input$n_levels))
  if (x$input$n_excluded > 0) {
    cat(sprintf(" (%d excluded due to NAs)", x$input$n_excluded))
  }
  cat("\n")
  
  # Dilution info (if applicable)
  if (!is.null(x$input$dilution)) {
    cat(sprintf("Design: Dilution-based (low pool mean = %s, high pool mean = %s)\n",
                format(x$input$endpoint_means["low"], digits = digits),
                format(x$input$endpoint_means["high"], digits = digits)))
  }
  
  # Concentration range
  target_range <- range(x$input$level_stats$target)
  cat(sprintf("Range: %s to %s\n",
              format(target_range[1], digits = digits),
              format(target_range[2], digits = digits)))
  cat("\n")
  
  # Weights
  if (x$settings$weights != "equal") {
    cat(sprintf("Weighting: %s\n", x$settings$weights))
  }
  
  # Best model
  cat(sprintf("Best-fitting model: %s\n", x$results$best_model))
  
  # Lack-of-fit summary
  lof <- x$results$lack_of_fit
  
  # Pure error LOF test (if available)
  if (!is.null(lof$pure_error)) {
    pe_lin <- lof$pure_error[lof$pure_error$model == "linear", ]
    if (nrow(pe_lin) > 0 && !is.na(pe_lin$f_statistic)) {
      cat("Lack-of-fit (pure error):\n")
      cat(sprintf("  Linear model:    F = %s, p = %s %s\n",
                  format(round(pe_lin$f_statistic, digits), nsmall = digits),
                  format(pe_lin$p_value, digits = digits, scientific = TRUE),
                  if (!is.na(pe_lin$significant) && pe_lin$significant)
                    "*" else ""))
      pe_quad <- lof$pure_error[lof$pure_error$model == "quadratic", ]
      if (nrow(pe_quad) > 0 && !is.na(pe_quad$f_statistic)) {
        cat(sprintf("  Quadratic model: F = %s, p = %s %s\n",
                    format(round(pe_quad$f_statistic, digits), nsmall = digits),
                    format(pe_quad$p_value, digits = digits, scientific = TRUE),
                    if (!is.na(pe_quad$significant) && pe_quad$significant)
                      "*" else ""))
      }
      cat("\n")
    }
  }
  
  # Nested polynomial tests
  cat("Polynomial comparison:\n")
  lq <- lof$polynomial[lof$polynomial$comparison == "Linear vs Quadratic", ]
  if (nrow(lq) > 0) {
    cat(sprintf("  Linear vs Quadratic: F = %s, p = %s %s\n",
                format(round(lq$f_statistic, digits), nsmall = digits),
                format(lq$p_value, digits = digits, scientific = TRUE),
                if (!is.na(lq$significant) && lq$significant) "*" else ""))
  }
  if (x$settings$max_poly >= 3) {
    qc <- lof$polynomial[lof$polynomial$comparison == "Quadratic vs Cubic", ]
    if (nrow(qc) > 0) {
      cat(sprintf("  Quadratic vs Cubic: F = %s, p = %s %s\n",
                  format(round(qc$f_statistic, digits), nsmall = digits),
                  format(qc$p_value, digits = digits, scientific = TRUE),
                  if (!is.na(qc$significant) && qc$significant) "*" else ""))
    }
  }
  cat("\n")
  
  # Maximum deviation
  md <- x$results$max_deviation
  cat("Maximum deviation from linearity:\n")
  cat(sprintf("  Absolute: %s (at concentration %s)\n",
              format(round(md$max_abs, digits), nsmall = digits),
              format(md$max_abs_at, digits = digits)))
  cat(sprintf("  Percent:  %s%% (at concentration %s)\n",
              format(round(md$max_pct, digits), nsmall = digits),
              format(md$max_pct_at, digits = digits)))
  
  # ADL assessment
  if (!is.null(x$settings$adl)) {
    cat("\n")
    adl_label <- if (x$settings$adl_type == "absolute") {
      format(x$settings$adl, digits = digits)
    } else {
      paste0(format(x$settings$adl, digits = digits), "%")
    }
    
    if (isTRUE(x$results$linear)) {
      cat(sprintf("ADL assessment: LINEAR (all deviations within +/-%s)\n",
                  adl_label))
    } else if (isFALSE(x$results$linear)) {
      cat(sprintf("ADL assessment: NON-LINEAR (deviations exceed +/-%s)\n",
                  adl_label))
      if (!is.null(x$results$linear_range)) {
        cat(sprintf("  Linear range: %s to %s\n",
                    format(x$results$linear_range[1], digits = digits),
                    format(x$results$linear_range[2], digits = digits)))
      }
    }
  }
  
  cat("\n")
  invisible(x)
}


#' Summary method for linearity_study objects
#'
#' @description
#' Provides a detailed summary of linearity evaluation results, including
#' per-level statistics, polynomial coefficients, deviations, and recovery.
#'
#' @param object An object of class `linearity_study`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class `summary.linearity_study` containing all
#'   summary components for formatted printing.
#'
#' @examples
#' set.seed(42)
#' target <- rep(c(10, 50, 100, 200, 300, 400, 500), each = 4)
#' measured <- target * rnorm(28, 1.0, 0.02) + rnorm(28, 0, 2)
#' lin <- linearity_study(target, measured, adl = 5, adl_type = "percent")
#' summary(lin)
#'
#' @export
summary.linearity_study <- function(object, ...) {
  
  structure(
    list(
      call = object$call,
      n = object$input$n,
      n_excluded = object$input$n_excluded,
      n_levels = object$input$n_levels,
      var_names = object$input$var_names,
      level_stats = object$input$level_stats,
      best_model = object$results$best_model,
      lack_of_fit = object$results$lack_of_fit,
      coefficients = object$results$coefficients,
      deviations = object$results$deviations,
      max_deviation = object$results$max_deviation,
      recovery = object$results$recovery,
      linear = object$results$linear,
      linear_range = object$results$linear_range,
      settings = object$settings
    ),
    class = "summary.linearity_study"
  )
}


#' Print method for summary.linearity_study objects
#'
#' @param x An object of class `summary.linearity_study`.
#' @param digits Number of significant digits to display (default: 4).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @export
print.summary.linearity_study <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat("Linearity Evaluation - Detailed Summary\n")
  cat(strrep("=", 50), "\n\n")
  
  # Call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  # Data summary
  cat(sprintf("Sample size: n = %d observations, %d levels",
              x$n, x$n_levels))
  if (x$n_excluded > 0) {
    cat(sprintf(" (%d excluded)", x$n_excluded))
  }
  cat("\n")
  cat(sprintf("Variables: target = '%s', measured = '%s'\n",
              x$var_names["x"], x$var_names["y"]))
  if (x$settings$weights != "equal") {
    cat(sprintf("Weighting: %s\n", x$settings$weights))
  }
  cat("\n")
  
  # Per-level statistics
  cat(strrep("-", 50), "\n")
  cat("Per-Level Statistics:\n")
  cat(strrep("-", 50), "\n")
  ls_print <- x$level_stats
  ls_print$mean <- round(ls_print$mean, digits)
  ls_print$sd <- round(ls_print$sd, digits)
  ls_print$cv <- round(ls_print$cv, digits - 1)
  print(ls_print, row.names = FALSE)
  cat("\n")
  
  # Polynomial coefficients
  cat(strrep("-", 50), "\n")
  cat("Polynomial Coefficients:\n")
  cat(strrep("-", 50), "\n")
  coef_print <- x$coefficients
  coef_print$estimate <- format(coef_print$estimate, digits = digits,
                                scientific = TRUE)
  print(coef_print, row.names = FALSE)
  cat("\n")
  
  # Lack-of-fit tests
  cat(strrep("-", 50), "\n")
  cat("Lack-of-Fit Tests:\n")
  cat(strrep("-", 50), "\n")
  
  # Pure error LOF test
  if (!is.null(x$lack_of_fit$pure_error)) {
    cat("\nPure Error Lack-of-Fit:\n")
    pe_print <- x$lack_of_fit$pure_error[,
                                         c("model", "ss_lof", "df_lof", "ss_pe", "df_pe",
                                           "f_statistic", "p_value", "significant")]
    pe_print$ss_lof <- round(pe_print$ss_lof, digits)
    pe_print$ss_pe <- round(pe_print$ss_pe, digits)
    pe_print$f_statistic <- round(pe_print$f_statistic, digits)
    pe_print$p_value <- format(pe_print$p_value, digits = digits,
                               scientific = TRUE)
    print(pe_print, row.names = FALSE)
    cat("\n")
  }
  
  # Nested polynomial tests
  cat("Nested Polynomial Comparisons:\n")
  poly_print <- x$lack_of_fit$polynomial[,
                                         c("comparison", "f_statistic", "p_value", "significant")]
  poly_print$f_statistic <- round(poly_print$f_statistic, digits)
  poly_print$p_value <- format(poly_print$p_value, digits = digits,
                               scientific = TRUE)
  print(poly_print, row.names = FALSE)
  cat(sprintf("\nBest-fitting model: %s\n\n", x$best_model))
  
  # Deviations from linearity
  cat(strrep("-", 50), "\n")
  cat("Deviations from Linearity:\n")
  cat(strrep("-", 50), "\n")
  dev_print <- x$deviations[, c("target", "observed_mean", "linear_predicted",
                                "deviation_abs", "deviation_pct")]
  dev_print$observed_mean <- round(dev_print$observed_mean, digits)
  dev_print$linear_predicted <- round(dev_print$linear_predicted, digits)
  dev_print$deviation_abs <- round(dev_print$deviation_abs, digits)
  dev_print$deviation_pct <- round(dev_print$deviation_pct, digits - 1)
  print(dev_print, row.names = FALSE)
  
  cat(sprintf("\nMax absolute deviation: %s at concentration %s\n",
              format(round(x$max_deviation$max_abs, digits), nsmall = digits),
              format(x$max_deviation$max_abs_at, digits = digits)))
  cat(sprintf("Max percent deviation:  %s%% at concentration %s\n",
              format(round(x$max_deviation$max_pct, digits), nsmall = digits),
              format(x$max_deviation$max_pct_at, digits = digits)))
  cat("\n")
  
  # Recovery
  cat(strrep("-", 50), "\n")
  cat("Recovery:\n")
  cat(strrep("-", 50), "\n")
  rec_print <- x$recovery
  rec_print$observed_mean <- round(rec_print$observed_mean, digits)
  rec_print$recovery_pct <- round(rec_print$recovery_pct, digits - 1)
  print(rec_print, row.names = FALSE)
  cat("\n")
  
  # ADL assessment
  if (!is.null(x$settings$adl)) {
    cat(strrep("-", 50), "\n")
    cat("ADL Assessment:\n")
    cat(strrep("-", 50), "\n")
    
    adl_label <- if (x$settings$adl_type == "absolute") {
      format(x$settings$adl, digits = digits)
    } else {
      paste0(format(x$settings$adl, digits = digits), "%")
    }
    
    cat(sprintf("Allowable Deviation from Linearity (ADL): +/-%s (%s)\n",
                adl_label, x$settings$adl_type))
    
    if (isTRUE(x$linear)) {
      cat("Result: LINEAR\n")
      cat("  All deviations are within the ADL.\n")
    } else if (isFALSE(x$linear)) {
      cat("Result: NON-LINEAR\n")
      cat("  One or more deviations exceed the ADL.\n")
      
      # Identify failing levels
      if (x$settings$adl_type == "absolute") {
        failing <- x$deviations$target[
          abs(x$deviations$deviation_abs) > x$settings$adl
        ]
      } else {
        failing <- x$deviations$target[
          abs(x$deviations$deviation_pct) > x$settings$adl
        ]
      }
      if (length(failing) > 0) {
        cat(sprintf("  Failing levels: %s\n",
                    paste(format(failing, digits = digits), collapse = ", ")))
      }
      
      if (!is.null(x$linear_range)) {
        cat(sprintf("  Estimated linear range: %s to %s\n",
                    format(x$linear_range[1], digits = digits),
                    format(x$linear_range[2], digits = digits)))
      }
    }
    cat("\n")
  }
  
  # Conclusion
  cat(strrep("-", 50), "\n")
  cat("Conclusion:\n")
  cat(strrep("-", 50), "\n")
  
  if (x$best_model == "linear") {
    cat("  No significant deviation from linearity detected.\n")
    cat("  The measurement procedure is linear across the tested range.\n")
  } else {
    cat(sprintf("  Significant %s component detected (p < 0.05).\n",
                x$best_model))
    cat("  The measurement procedure shows non-linear behavior.\n")
    cat("  Evaluate whether deviations are clinically acceptable\n")
    cat("  for your intended application.\n")
  }
  cat("\n")
  
  invisible(x)
}