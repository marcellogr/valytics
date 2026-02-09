#' Print method for precision_profile objects
#'
#' @description
#' Displays a concise summary of precision profile results, including model
#' fit and functional sensitivity estimates.
#'
#' @param x An object of class `precision_profile`.
#' @param digits Number of significant digits to display (default: 3).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' # See ?precision_profile for examples
#'
#' @export
print.precision_profile <- function(x, digits = 3, ...) {
  
  cat("\n")
  cat("Precision Profile Analysis\n")
  cat(strrep("-", 40), "\n")
  
  # Sample info
  cat(sprintf("n = %d concentration levels\n", x$input$n_levels))
  cat(sprintf("Concentration range: %.3g to %.3g (%.2f-fold)\n",
              x$input$conc_range["min"],
              x$input$conc_range["max"],
              x$input$conc_span))
  cat("\n")
  
  # Model info
  cat(sprintf("Model: %s\n", tools::toTitleCase(x$model$type)))
  cat(sprintf("  %s\n", x$model$equation))
  cat("\n")
  
  # Model parameters
  cat("Parameters:\n")
  for (i in seq_along(x$model$parameters)) {
    param_name <- names(x$model$parameters)[i]
    param_val <- x$model$parameters[i]
    cat(sprintf("  %s = %.*f\n", param_name, digits, param_val))
  }
  cat("\n")
  
  # Fit quality
  cat("Fit Quality:\n")
  cat(sprintf("  R-squared = %.*f\n", digits, x$fit_quality$r_squared))
  cat(sprintf("  RMSE = %.*f\n", digits, x$fit_quality$rmse))
  cat("\n")
  
  # Functional sensitivity
  cat("Functional Sensitivity:\n")
  
  for (i in seq_len(nrow(x$functional_sensitivity))) {
    fs <- x$functional_sensitivity[i, ]
    target <- fs$cv_target
    conc <- fs$concentration
    achievable <- fs$achievable
    
    if (achievable && !is.na(conc)) {
      ci_available <- !is.na(fs$ci_lower) && !is.na(fs$ci_upper)
      
      if (ci_available) {
        cat(sprintf("  CV = %g%%: concentration = %.*f (%.*f to %.*f)\n",
                    target, digits, conc,
                    digits, fs$ci_lower,
                    digits, fs$ci_upper))
      } else {
        cat(sprintf("  CV = %g%%: concentration = %.*f\n",
                    target, digits, conc))
      }
    } else {
      cat(sprintf("  CV = %g%%: not achievable within observed range\n", target))
    }
  }
  
  cat("\n")
  
  invisible(x)
}


#' Summary method for precision_profile objects
#'
#' @description
#' Provides a detailed summary of precision profile results, including fitted
#' values, residuals, fit statistics, and functional sensitivity estimates.
#'
#' @param object An object of class `precision_profile`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class `summary.precision_profile` containing summary
#'   statistics.
#'
#' @examples
#' # See ?precision_profile for examples
#'
#' @export
summary.precision_profile <- function(object, ...) {
  
  x <- object
  
  cat("\n")
  cat("Precision Profile Analysis - Detailed Summary\n")
  cat(strrep("=", 50), "\n\n")
  
  # Input summary ----
  cat("Data:\n")
  cat(sprintf("  Number of concentration levels: %d\n", x$input$n_levels))
  cat(sprintf("  Concentration range: %.3g to %.3g\n",
              x$input$conc_range["min"],
              x$input$conc_range["max"]))
  cat(sprintf("  Concentration span: %.2f-fold\n", x$input$conc_span))
  cat("\n")
  
  # Model summary ----
  cat("Model:\n")
  cat(sprintf("  Type: %s\n", tools::toTitleCase(x$model$type)))
  cat(sprintf("  Equation: %s\n", x$model$equation))
  cat("\n")
  
  # Parameters
  cat("Parameters:\n")
  cat(strrep("-", 50), "\n")
  param_df <- data.frame(
    Parameter = names(x$model$parameters),
    Estimate = x$model$parameters,
    row.names = NULL
  )
  print(param_df, digits = 4)
  cat("\n")
  
  # Fit quality ----
  cat("Fit Quality:\n")
  cat(strrep("-", 50), "\n")
  fit_df <- data.frame(
    Statistic = c("R-squared", "Adjusted R-squared", "RMSE", "MAE"),
    Value = c(
      x$fit_quality$r_squared,
      x$fit_quality$adj_r_squared,
      x$fit_quality$rmse,
      x$fit_quality$mae
    ),
    row.names = NULL
  )
  print(fit_df, digits = 4)
  cat("\n")
  
  # Fitted values and residuals ----
  cat("Fitted Values:\n")
  cat(strrep("-", 50), "\n")
  fitted_display <- x$fitted[, c("concentration", "cv_observed", "cv_fitted", "residual")]
  print(fitted_display, digits = 3, row.names = FALSE)
  cat("\n")
  
  # Residual summary
  cat("Residual Summary:\n")
  cat(strrep("-", 50), "\n")
  residual_summary <- summary(x$fitted$residual)
  print(residual_summary)
  cat("\n")
  
  # Functional sensitivity ----
  cat("Functional Sensitivity:\n")
  cat(strrep("-", 50), "\n")
  
  fs_display <- x$functional_sensitivity
  
  # Add interpretation column
  fs_display$interpretation <- ifelse(
    fs_display$achievable,
    "Achievable",
    "Not achievable"
  )
  
  # Reorder columns
  fs_cols <- c("cv_target", "concentration", "ci_lower", "ci_upper", "interpretation")
  fs_display <- fs_display[, fs_cols]
  
  print(fs_display, digits = 3, row.names = FALSE)
  cat("\n")
  
  # Interpretation ----
  cat("Interpretation:\n")
  cat(strrep("-", 50), "\n")
  
  # Check achievability of targets
  achievable_targets <- x$functional_sensitivity$cv_target[x$functional_sensitivity$achievable]
  not_achievable_targets <- x$functional_sensitivity$cv_target[!x$functional_sensitivity$achievable]
  
  if (length(achievable_targets) > 0) {
    cat("Achievable CV targets:\n")
    for (target in achievable_targets) {
      fs_row <- x$functional_sensitivity[x$functional_sensitivity$cv_target == target, ]
      cat(sprintf("  - %g%% CV at concentration %.3g\n",
                  target, fs_row$concentration))
    }
  }
  
  if (length(not_achievable_targets) > 0) {
    cat("\n")
    cat("CV targets not achievable within observed range:\n")
    for (target in not_achievable_targets) {
      cat(sprintf("  - %g%%\n", target))
    }
    
    # Suggest what's achievable
    min_cv_achievable <- min(x$fitted$cv_fitted)
    cat(sprintf("\nMinimum achievable CV (at high concentration): %.2f%%\n",
                min_cv_achievable))
  }
  
  # Model interpretation
  cat("\n")
  if (x$model$type == "hyperbolic") {
    a <- x$model$parameters["a"]
    cat(sprintf("Asymptotic CV at high concentration: %.2f%%\n", a))
  } else {
    a <- x$model$parameters["a"]
    if (a > 0) {
      cat(sprintf("Baseline CV (intercept): %.2f%%\n", a))
    }
  }
  
  cat("\n")
  
  # Return summary object invisibly
  invisible(list(
    input = x$input,
    model = x$model,
    fitted = x$fitted,
    fit_quality = x$fit_quality,
    functional_sensitivity = x$functional_sensitivity
  ))
}


#' Print method for summary.precision_profile objects
#'
#' @param x An object of class `summary.precision_profile`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @export
print.summary.precision_profile <- function(x, ...) {
  # The summary method already prints everything
  # This is just to maintain S3 consistency
  invisible(x)
}