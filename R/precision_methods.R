#' Print method for precision_study objects
#'
#' @description
#' Displays a concise summary of precision study results, including
#' variance components and key precision estimates.
#'
#' @param x An object of class `precision_study`.
#' @param digits Number of significant digits to display (default: 3).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' data <- data.frame(
#'   day = rep(1:5, each = 4),
#'   value = rnorm(20, mean = 100, sd = 5)
#' )
#' data$value <- data$value + rep(rnorm(5, 0, 3), each = 4)
#'
#' prec <- precision_study(data, value = "value", day = "day")
#' print(prec)
#'
#' @seealso [summary.precision_study()] for detailed output
#' @export
print.precision_study <- function(x, digits = 3, ...) {
  
  cat("\n")
  cat("Precision Study Analysis\n")
  cat(strrep("-", 45), "\n")
  
  # Sample info
  cat(sprintf("n = %d observations", x$input$n))
  if (x$input$n_excluded > 0) {
    cat(sprintf(" (%d excluded due to NAs)", x$input$n_excluded))
  }
  
  cat("\n")
  
  # Design info
  cat(sprintf("Design: %s", x$design$structure))
  if (!x$design$balanced) {
    cat(" (unbalanced)")
  }
  cat("\n")
  
  # Method info
  method_str <- if (x$settings$method == "anova") "ANOVA (Method of Moments)" else "REML"
  cat(sprintf("Estimation: %s\n", method_str))
  
  ci_str <- switch(x$settings$ci_method,
                   satterthwaite = "Satterthwaite",
                   mls = "Modified Large Sample",
                   bootstrap = sprintf("Bootstrap (n = %d)", x$settings$boot_n))
  cat(sprintf("CI method: %s, %s%% CI\n", ci_str, x$settings$conf_level * 100))
  cat("\n")
  
  # Multi-sample info
  if (!is.null(x$by_sample) && length(x$by_sample) > 1) {
    cat(sprintf("Samples: %d concentration levels\n", length(x$by_sample)))
    cat("(Showing results for first sample; use $by_sample for all)\n\n")
  }
  
  # Precision estimates
  cat("Precision Estimates:\n")
  cat(strrep("-", 45), "\n")
  
  prec <- x$precision
  ci_pct <- sprintf("%g%%", x$settings$conf_level * 100)
  
  # Format precision table
  for (i in seq_len(nrow(prec))) {
    measure <- prec$measure[i]
    sd_val <- format(round(prec$sd[i], digits), nsmall = digits)
    cv_val <- format(round(prec$cv_pct[i], 2), nsmall = 2)
    
    # Format CI
    if (!is.na(prec$ci_lower[i]) && !is.na(prec$ci_upper[i])) {
      ci_str <- sprintf("[%s, %s]",
                        format(round(prec$ci_lower[i], digits), nsmall = digits),
                        format(round(prec$ci_upper[i], digits), nsmall = digits))
    } else {
      ci_str <- "[NA, NA]"
    }
    
    cat(sprintf("  %-20s SD = %s  (CV = %s%%)\n", 
                paste0(measure, ":"), sd_val, cv_val))
    cat(sprintf("  %-20s %s CI: %s\n", "", ci_pct, ci_str))
  }
  
  cat("\n")
  
  invisible(x)
}


#' Summary method for precision_study objects
#'
#' @description
#' Provides a detailed summary of precision study results, including
#' variance components, ANOVA table (for ANOVA method), precision estimates
#' with confidence intervals, and design information.
#'
#' @param object An object of class `precision_study`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class `summary.precision_study` containing:
#'   \describe{
#'     \item{call}{The original function call.}
#'     \item{n}{Number of observations.}
#'     \item{n_excluded}{Number of observations excluded due to NAs.}
#'     \item{design}{Design information list.}
#'     \item{settings}{Analysis settings.}
#'     \item{variance_components}{Data frame of variance components.}
#'     \item{precision}{Data frame of precision estimates.}
#'     \item{anova_table}{ANOVA table (if method = "anova").}
#'     \item{by_sample}{Results by sample (if multiple samples).}
#'   }
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' data <- data.frame(
#'   day = rep(1:5, each = 4),
#'   value = rnorm(20, mean = 100, sd = 5)
#' )
#' data$value <- data$value + rep(rnorm(5, 0, 3), each = 4)
#'
#' prec <- precision_study(data, value = "value", day = "day")
#' summary(prec)
#'
#' @seealso [print.precision_study()] for concise output
#' @export
summary.precision_study <- function(object, ...) {
  
  structure(
    list(
      call = object$call,
      n = object$input$n,
      n_excluded = object$input$n_excluded,
      factors = object$input$factors,
      design = object$design,
      settings = object$settings,
      variance_components = object$variance_components,
      precision = object$precision,
      anova_table = object$anova_table,
      by_sample = object$by_sample,
      sample_means = object$sample_means
    ),
    class = "summary.precision_study"
  )
}


#' Print method for summary.precision_study objects
#'
#' @param x An object of class `summary.precision_study`.
#' @param digits Number of significant digits to display (default: 4).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @export
print.summary.precision_study <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat("Precision Study Analysis - Detailed Summary\n")
  cat(strrep("=", 55), "\n\n")
  
  # Call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  # Design Information
  cat(strrep("-", 55), "\n")
  cat("Design Information:\n")
  cat(strrep("-", 55), "\n")
  
  cat(sprintf("  Structure: %s\n", x$design$structure))
  cat(sprintf("  Type: %s\n", x$design$type))
  cat(sprintf("  Balanced: %s\n", if (x$design$balanced) "Yes" else "No"))
  
  # Factor levels
  cat("  Factor levels:\n")
  for (fname in names(x$design$levels)) {
    cat(sprintf("    %s: %d\n", fname, x$design$levels[[fname]]))
  }
  cat("\n")
  
  # Sample info
  cat(sprintf("  Total observations: %d\n", x$n))
  if (x$n_excluded > 0) {
    cat(sprintf("  Excluded (NA): %d\n", x$n_excluded))
  }
  cat("\n")
  
  # Settings
  cat(strrep("-", 55), "\n")
  cat("Analysis Settings:\n")
  cat(strrep("-", 55), "\n")
  
  method_str <- if (x$settings$method == "anova") {
    "ANOVA (Method of Moments)"
  } else {
    "REML (Restricted Maximum Likelihood)"
  }
  cat(sprintf("  Estimation method: %s\n", method_str))
  
  ci_str <- switch(x$settings$ci_method,
                   satterthwaite = "Satterthwaite approximation",
                   mls = "Modified Large Sample (MLS)",
                   bootstrap = sprintf("Bootstrap BCa (n = %d)", x$settings$boot_n))
  cat(sprintf("  CI method: %s\n", ci_str))
  cat(sprintf("  Confidence level: %g%%\n", x$settings$conf_level * 100))
  cat("\n")
  
  # Variance Components
  cat(strrep("-", 55), "\n")
  cat("Variance Components:\n")
  cat(strrep("-", 55), "\n")
  
  vc <- x$variance_components
  vc_display <- data.frame(
    Component = vc$component,
    Variance = round(vc$variance, digits),
    SD = round(vc$sd, digits),
    `Pct Total` = round(vc$pct_total, 1),
    df = vc$df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  print(vc_display, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # ANOVA Table (if available)
  if (!is.null(x$anova_table)) {
    cat(strrep("-", 55), "\n")
    cat("ANOVA Table:\n")
    cat(strrep("-", 55), "\n")
    
    aov_tbl <- x$anova_table
    aov_display <- data.frame(
      Source = aov_tbl$source,
      df = aov_tbl$df,
      SS = round(aov_tbl$ss, digits),
      MS = round(aov_tbl$ms, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    print(aov_display, row.names = FALSE, right = FALSE)
    cat("\n")
  }
  
  # Precision Estimates
  cat(strrep("-", 55), "\n")
  cat("Precision Estimates:\n")
  cat(strrep("-", 55), "\n")
  
  prec <- x$precision
  ci_pct <- sprintf("%g%%", x$settings$conf_level * 100)
  
  prec_display <- data.frame(
    Measure = prec$measure,
    SD = round(prec$sd, digits),
    `CV (%)` = round(prec$cv_pct, 2),
    `CI Lower` = round(prec$ci_lower, digits),
    `CI Upper` = round(prec$ci_upper, digits),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(prec_display)[4:5] <- c(paste0(ci_pct, " Lower"), paste0(ci_pct, " Upper"))
  print(prec_display, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # Multi-sample results
  if (!is.null(x$by_sample) && length(x$by_sample) > 1) {
    cat(strrep("-", 55), "\n")
    cat("Results by Sample:\n")
    cat(strrep("-", 55), "\n")
    
    sample_names <- names(x$by_sample)
    for (i in seq_along(x$by_sample)) {
      samp <- x$by_sample[[i]]
      samp_name <- sample_names[i]
      samp_mean <- if (!is.null(x$sample_means)) {
        round(x$sample_means[samp_name], digits)
      } else {
        NA
      }
      
      cat(sprintf("\n  Sample: %s", samp_name))
      if (!is.na(samp_mean)) {
        cat(sprintf(" (mean = %s)", samp_mean))
      }
      cat("\n")
      
      # Show key precision metrics for each sample
      sp <- samp$precision
      repeat_idx <- which(sp$measure == "Repeatability")
      if (length(repeat_idx) > 0) {
        cat(sprintf("    Repeatability:  SD = %s, CV = %s%%\n",
                    round(sp$sd[repeat_idx], digits),
                    round(sp$cv_pct[repeat_idx], 2)))
      }
      
      # Intermediate precision (look for various names)
      inter_idx <- which(grepl("Intermediate|Within-laboratory", sp$measure, 
                               ignore.case = TRUE))
      if (length(inter_idx) > 0) {
        cat(sprintf("    Intermediate:   SD = %s, CV = %s%%\n",
                    round(sp$sd[inter_idx[1]], digits),
                    round(sp$cv_pct[inter_idx[1]], 2)))
      }
    }
    cat("\n")
  }
  
  # Interpretation guidance
  cat(strrep("-", 55), "\n")
  cat("Interpretation:\n")
  cat(strrep("-", 55), "\n")
  
  # Get repeatability and intermediate precision
  prec <- x$precision
  repeat_cv <- prec$cv_pct[prec$measure == "Repeatability"]
  inter_idx <- which(grepl("Intermediate|Within-laboratory", prec$measure, 
                           ignore.case = TRUE))
  inter_cv <- if (length(inter_idx) > 0) prec$cv_pct[inter_idx[1]] else NA
  
  cat(sprintf("  Repeatability CV: %s%%\n", round(repeat_cv, 2)))
  if (!is.na(inter_cv)) {
    cat(sprintf("  Intermediate precision CV: %s%%\n", round(inter_cv, 2)))
    
    # Ratio interpretation
    if (repeat_cv > 0) {
      ratio <- inter_cv / repeat_cv
      cat(sprintf("  Ratio (intermediate/repeatability): %.2f\n", ratio))
      if (ratio > 1.5) {
        cat("  Note: Substantial between-day/run variation detected.\n")
      }
    }
  }
  
  cat("\n")
  
  invisible(x)
}