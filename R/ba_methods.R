#' Print method for ba_analysis objects
#'
#' @description
#' Displays a concise summary of Bland-Altman analysis results.
#'
#' @param x An object of class `ba_analysis`.
#' @param digits Number of significant digits to display (default: 3).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' set.seed(42)
#' method_a <- rnorm(50, mean = 100, sd = 15)
#' method_b <- method_a + rnorm(50, mean = 2, sd = 5)
#' ba <- ba_analysis(method_a, method_b)
#' print(ba)
#'
#' @export
print.ba_analysis <- function(x, digits = 3, ...) {

  cat("\n")
  cat("Bland-Altman Analysis\n")
  cat(strrep("-", 40), "\n")

  # Sample info
  cat(sprintf("n = %d paired observations", x$input$n))
  if (x$input$n_excluded > 0) {
    cat(sprintf(" (%d excluded due to NAs)", x$input$n_excluded))
  }
  cat("\n\n")

  # Type of analysis
  diff_type <- if (x$settings$type == "absolute") "Absolute" else "Percent"
  cat(sprintf("Difference type: %s (y - x)\n", diff_type))
  cat(sprintf("Confidence level: %s%%\n\n",
              format(x$settings$conf_level * 100, nsmall = 0)))

  # Results
  ci_pct <- format(x$settings$conf_level * 100, nsmall = 0)

  cat("Results:\n")
  cat(sprintf("  Bias (mean difference): %s\n",
              format(round(x$results$bias, digits), nsmall = digits)))
  cat(sprintf("    %s%% CI: [%s, %s]\n",
              ci_pct,
              format(round(x$results$bias_ci["lower"], digits), nsmall = digits),
              format(round(x$results$bias_ci["upper"], digits), nsmall = digits)))

  cat(sprintf("  SD of differences: %s\n\n",
              format(round(x$results$sd_diff, digits), nsmall = digits)))

  cat("Limits of Agreement:\n")
  cat(sprintf("  Lower LoA: %s\n",
              format(round(x$results$loa_lower, digits), nsmall = digits)))
  cat(sprintf("    %s%% CI: [%s, %s]\n",
              ci_pct,
              format(round(x$results$loa_lower_ci["lower"], digits), nsmall = digits),
              format(round(x$results$loa_lower_ci["upper"], digits), nsmall = digits)))
  cat(sprintf("  Upper LoA: %s\n",
              format(round(x$results$loa_upper, digits), nsmall = digits)))
  cat(sprintf("    %s%% CI: [%s, %s]\n",
              ci_pct,
              format(round(x$results$loa_upper_ci["lower"], digits), nsmall = digits),
              format(round(x$results$loa_upper_ci["upper"], digits), nsmall = digits)))

  cat("\n")

  invisible(x)
}


#' Summary method for ba_analysis objects
#'
#' @description
#' Provides a detailed summary of Bland-Altman analysis results, including
#' additional diagnostics and descriptive statistics.
#'
#' @param object An object of class `ba_analysis`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class `summary.ba_analysis` containing:
#'   \describe{
#'     \item{call}{The original function call.}
#'     \item{n}{Number of paired observations.}
#'     \item{n_excluded}{Number of pairs excluded due to NAs.}
#'     \item{var_names}{Variable names for x and y.}
#'     \item{type}{Type of difference calculation.}
#'     \item{conf_level}{Confidence level used.}
#'     \item{descriptives}{Data frame with descriptive statistics.}
#'     \item{agreement}{Data frame with agreement statistics.}
#'     \item{normality_test}{Shapiro-Wilk test result for differences.}
#'   }
#'
#' @examples
#' set.seed(42)
#' method_a <- rnorm(50, mean = 100, sd = 15)
#' method_b <- method_a + rnorm(50, mean = 2, sd = 5)
#' ba <- ba_analysis(method_a, method_b)
#' summary(ba)
#'
#' @export
summary.ba_analysis <- function(object, ...) {

  # Descriptive statistics for both methods
  x <- object$input$x
  y <- object$input$y

  descriptives <- data.frame(
    Variable = c(object$input$var_names["x"], object$input$var_names["y"]),
    N = c(length(x), length(y)),
    Mean = c(mean(x), mean(y)),
    SD = c(stats::sd(x), stats::sd(y)),
    Median = c(stats::median(x), stats::median(y)),
    Min = c(min(x), min(y)),
    Max = c(max(x), max(y)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Agreement statistics
  ci_level <- paste0(object$settings$conf_level * 100, "%")

  agreement <- data.frame(
    Statistic = c("Bias", "Lower LoA", "Upper LoA"),
    Estimate = c(object$results$bias,
                 object$results$loa_lower,
                 object$results$loa_upper),
    CI_Lower = c(object$results$bias_ci["lower"],
                 object$results$loa_lower_ci["lower"],
                 object$results$loa_upper_ci["lower"]),
    CI_Upper = c(object$results$bias_ci["upper"],
                 object$results$loa_lower_ci["upper"],
                 object$results$loa_upper_ci["upper"]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  names(agreement)[3:4] <- c(paste0("CI_Lower_", ci_level),
                              paste0("CI_Upper_", ci_level))

  # Normality test for differences (Shapiro-Wilk)
  # Only run if n is in valid range for shapiro.test (3-5000)
  n <- length(object$results$differences)
  if (n >= 3 && n <= 5000) {
    normality_test <- stats::shapiro.test(object$results$differences)
  } else {
    normality_test <- NULL
  }

  # Construct summary object
  structure(
    list(
      call = object$call,
      n = object$input$n,
      n_excluded = object$input$n_excluded,
      var_names = object$input$var_names,
      type = object$settings$type,
      conf_level = object$settings$conf_level,
      descriptives = descriptives,
      agreement = agreement,
      normality_test = normality_test,
      sd_diff = object$results$sd_diff
    ),
    class = "summary.ba_analysis"
  )
}


#' Print method for summary.ba_analysis objects
#'
#' @param x An object of class `summary.ba_analysis`.
#' @param digits Number of significant digits to display (default: 4).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @export
print.summary.ba_analysis <- function(x, digits = 4, ...) {

  cat("\n")
  cat("Bland-Altman Analysis - Detailed Summary\n")
  cat(strrep("=", 50), "\n\n")

  # Call
  cat("Call:\n")
  print(x$call)
  cat("\n")

  # Sample info
  cat(sprintf("Sample size: n = %d", x$n))
  if (x$n_excluded > 0) {
    cat(sprintf(" (%d pairs excluded due to missing values)", x$n_excluded))
  }
  cat("\n")
  cat(sprintf("Variables: x = '%s', y = '%s'\n", x$var_names["x"], x$var_names["y"]))
  cat(sprintf("Difference type: %s\n",
              if (x$type == "absolute") "Absolute (y - x)" else "Percent"))
  cat(sprintf("Confidence level: %s%%\n", x$conf_level * 100))
  cat("\n")

  # Descriptive statistics
  cat(strrep("-", 50), "\n")
  cat("Descriptive Statistics:\n")
  cat(strrep("-", 50), "\n")
  print(x$descriptives, digits = digits, row.names = FALSE)
  cat("\n")

  # Agreement statistics
  cat(strrep("-", 50), "\n")
  cat("Agreement Statistics:\n")
  cat(strrep("-", 50), "\n")
  print(x$agreement, digits = digits, row.names = FALSE)
  cat(sprintf("\nSD of differences: %s\n",
              format(round(x$sd_diff, digits), nsmall = digits)))
  cat("\n")

  # Normality test
  if (!is.null(x$normality_test)) {
    cat(strrep("-", 50), "\n")
    cat("Normality of Differences (Shapiro-Wilk test):\n")
    cat(strrep("-", 50), "\n")
    cat(sprintf("W = %s, p-value = %s\n",
                format(round(x$normality_test$statistic, digits), nsmall = digits),
                format(x$normality_test$p.value, digits = digits, scientific = TRUE)))

    if (x$normality_test$p.value < 0.05) {
      cat("Note: p < 0.05 suggests differences may not be normally distributed.\n")
      cat("      Consider inspecting the Bland-Altman plot for patterns.\n")
    } else {
      cat("Note: No evidence against normality (p >= 0.05).\n")
    }
  }

  cat("\n")

  invisible(x)
}
