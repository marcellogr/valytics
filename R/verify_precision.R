#' Precision Verification Against Manufacturer Claims
#'
#' @description
#' Compares observed precision (CV or SD) to manufacturer's claimed performance
#' using statistical hypothesis testing. This function implements verification
#' protocols for validating that an analytical method meets specified precision
#' goals.
#'
#' @param x Either a numeric vector of measurements, a `precision_study` object,
#'   or a data frame containing precision study data.
#' @param claimed_cv Manufacturer's claimed coefficient of variation (as percent).
#'   Either `claimed_cv` or `claimed_sd` must be provided.
#' @param claimed_sd Manufacturer's claimed standard deviation. Either `claimed_cv`
#'   or `claimed_sd` must be provided.
#' @param mean_value Mean concentration of the sample. Required when `x` is a
#'   numeric vector and `claimed_cv` is used, or when `claimed_sd` is used and
#'   CV-based comparison is desired.
#' @param alpha Significance level for the hypothesis test (default: 0.05).
#' @param alternative Type of alternative hypothesis:
#'   `"less"` (default) tests if observed is not worse than claimed,
#'   `"two.sided"` tests for any difference,
#'   `"greater"` tests if observed is worse than claimed.
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param value Character string specifying the column name containing
#'   measurement values when `x` is a data frame. Default is `"value"`.
#' @param day Character string specifying the column name for day identifier
#'   when `x` is a data frame. Default is `"day"`.
#' @param run Character string specifying the column name for run identifier
#'   when `x` is a data frame. Default is `NULL`.
#' @param ... Additional arguments passed to `precision_study()` when `x` is
#'   a data frame.
#'
#' @return An object of class `c("verify_precision", "valytics_precision", "valytics_result")`,
#'   which is a list containing:
#'
#'   \describe{
#'     \item{input}{List with input data and metadata:
#'       \itemize{
#'         \item `n`: Number of observations
#'         \item `df`: Degrees of freedom for the test
#'         \item `mean_value`: Mean of measurements
#'         \item `source`: Description of input source
#'       }
#'     }
#'     \item{observed}{List with observed precision:
#'       \itemize{
#'         \item `sd`: Observed standard deviation
#'         \item `cv_pct`: Observed CV (percent)
#'         \item `variance`: Observed variance
#'       }
#'     }
#'     \item{claimed}{List with manufacturer's claimed precision:
#'       \itemize{
#'         \item `sd`: Claimed SD
#'         \item `cv_pct`: Claimed CV (percent)
#'         \item `variance`: Claimed variance
#'       }
#'     }
#'     \item{test}{List with hypothesis test results:
#'       \itemize{
#'         \item `statistic`: Chi-square test statistic
#'         \item `df`: Degrees of freedom
#'         \item `p_value`: P-value
#'         \item `alternative`: Alternative hypothesis used
#'         \item `method`: Test method description
#'       }
#'     }
#'     \item{verification}{List with verification outcome:
#'       \itemize{
#'         \item `verified`: Logical; TRUE if precision is verified
#'         \item `ratio`: Ratio of observed variance to claimed variance
#'         \item `cv_ratio`: Ratio of observed CV to claimed CV
#'         \item `upper_verification_limit`: Upper limit for verification
#'       }
#'     }
#'     \item{ci}{List with confidence intervals:
#'       \itemize{
#'         \item `sd_ci`: CI for standard deviation
#'         \item `cv_ci`: CI for CV (percent)
#'         \item `variance_ci`: CI for variance
#'       }
#'     }
#'     \item{settings}{List with analysis settings}
#'     \item{call}{The matched function call}
#'   }
#'
#' @details
#' **Statistical Test:**
#'
#' The verification uses a chi-square test comparing observed variance to
#' claimed variance:
#'
#' \deqn{\chi^2 = \frac{(n-1) \cdot s^2}{\sigma^2_{claimed}}}
#'
#' where \eqn{s^2} is the observed sample variance and \eqn{\sigma^2_{claimed}}
#' is the manufacturer's claimed variance.
#'
#' **Hypothesis Testing:**
#'
#' For `alternative = "less"` (default, recommended for verification):
#' \itemize{
#'   \item H0: True precision is worse than or equal to claimed
#'   \item H1: True precision is better than or equal to claimed
#'   \item Verification passes if observed precision is not significantly worse
#' }
#'
#' For typical verification studies, the observed CV should not exceed the
#' manufacturer's claimed CV by more than expected from sampling variability.
#'
#' **Verification Limit:**
#'
#' The upper verification limit (UVL) represents the maximum observed CV that
#' would still be consistent with the claimed CV at the given significance level:
#'
#' \deqn{UVL = CV_{claimed} \cdot \sqrt{\frac{\chi^2_{1-\alpha, df}}{df}}}
#'
#' If observed CV <= UVL, precision is verified.
#'
#' @section Input Options:
#' The function accepts three types of input:
#' \itemize{
#'   \item **Numeric vector**: Raw measurements (simplest case)
#'   \item **precision_study object**: Uses within-laboratory precision from
#'     a previous analysis
#'   \item **Data frame**: Runs `precision_study()` internally with specified
#'     factors
#' }
#'
#' @references
#' Chesher D (2008). Evaluating assay precision. \emph{Clinical Biochemist
#' Reviews}, 29(Suppl 1):S23-S26.
#'
#' ISO 5725-6:1994. Accuracy (trueness and precision) of measurement methods
#' and results - Part 6: Use in practice of accuracy values.
#'
#' @seealso
#' [precision_study()] for full precision analysis,
#' [ate_assessment()] for total error assessment
#'
#' @examples
#' # Example 1: Verify precision from raw measurements
#' set.seed(42)
#' measurements <- rnorm(25, mean = 100, sd = 3.5)
#' 
#' # Manufacturer claims CV = 4%
#' result <- verify_precision(measurements, claimed_cv = 4, mean_value = 100)
#' print(result)
#'
#' # Example 2: Verify precision from a precision_study object
#' prec_data <- data.frame(
#'   day = rep(1:5, each = 5),
#'   value = rnorm(25, mean = 100, sd = 3)
#' )
#' prec_data$value <- prec_data$value + rep(rnorm(5, 0, 1.5), each = 5)
#'
#' prec <- precision_study(prec_data, value = "value", day = "day")
#' result <- verify_precision(prec, claimed_cv = 5)
#' print(result)
#'
#' # Example 3: Verify precision directly from data frame
#' result <- verify_precision(
#'   prec_data,
#'   claimed_cv = 5,
#'   value = "value",
#'   day = "day"
#' )
#' print(result)
#'
#' @export
verify_precision <- function(x,
                             claimed_cv = NULL,
                             claimed_sd = NULL,
                             mean_value = NULL,
                             alpha = 0.05,
                             alternative = c("less", "two.sided", "greater"),
                             conf_level = 0.95,
                             value = "value",
                             day = "day",
                             run = NULL,
                             ...) {
  
  # Capture call
  call <- match.call()
  
  # Match arguments
  alternative <- match.arg(alternative)
  
  # Input validation ----
  .validate_verify_precision_input(
    claimed_cv = claimed_cv,
    claimed_sd = claimed_sd,
    alpha = alpha,
    conf_level = conf_level
  )
  
  # Extract precision information based on input type ----
  parsed <- .parse_verify_precision_input(
    x = x,
    mean_value = mean_value,
    value = value,
    day = day,
    run = run,
    conf_level = conf_level,
    ...
  )
  
  observed_sd <- parsed$sd
  observed_mean <- parsed$mean
  n <- parsed$n
  df <- parsed$df
  source_desc <- parsed$source
  
  # Handle claimed values ----
  claimed <- .resolve_claimed_precision(
    claimed_cv = claimed_cv,
    claimed_sd = claimed_sd,
    mean_value = observed_mean
  )
  
  claimed_sd_val <- claimed$sd
  claimed_cv_val <- claimed$cv_pct
  
  # Calculate observed statistics ----
  observed_variance <- observed_sd^2
  observed_cv <- 100 * observed_sd / observed_mean
  claimed_variance <- claimed_sd_val^2
  
  # Chi-square test ----
  test_result <- .chi_square_variance_test(
    observed_variance = observed_variance,
    claimed_variance = claimed_variance,
    df = df,
    alpha = alpha,
    alternative = alternative
  )
  
  # Calculate verification limit ----
  # Upper verification limit: maximum observed CV consistent with claimed
  chi_sq_critical <- stats::qchisq(1 - alpha, df = df)
  uvl_variance <- claimed_variance * chi_sq_critical / df
  uvl_sd <- sqrt(uvl_variance)
  uvl_cv <- 100 * uvl_sd / observed_mean
  
  # Verification outcome ----
  variance_ratio <- observed_variance / claimed_variance
  cv_ratio <- observed_cv / claimed_cv_val
  
  # For "less" alternative: verified if observed is not significantly worse
  verified <- switch(
    alternative,
    "less" = test_result$p_value > alpha,  # Can't reject that observed <= claimed
    "greater" = test_result$p_value <= alpha,  # Observed significantly better
    "two.sided" = test_result$p_value > alpha  # No significant difference
  )
  
  # Actually for "less" we want: can we show observed is acceptable?
  # The test is: H0: sigma^2 >= sigma0^2 vs H1: sigma^2 < sigma0^2
  # We verify if we can't reject that observed is worse (i.e., p > alpha)
  # OR more practically: observed CV <= upper verification limit
  verified <- observed_cv <= uvl_cv
  
  # Confidence intervals for observed precision ----
  ci <- .calculate_variance_ci(
    observed_variance = observed_variance,
    observed_sd = observed_sd,
    observed_cv = observed_cv,
    df = df,
    conf_level = conf_level
  )
  
  # Build result ----
  result <- list(
    input = list(
      n = n,
      df = df,
      mean_value = observed_mean,
      source = source_desc
    ),
    observed = list(
      sd = observed_sd,
      cv_pct = observed_cv,
      variance = observed_variance
    ),
    claimed = list(
      sd = claimed_sd_val,
      cv_pct = claimed_cv_val,
      variance = claimed_variance
    ),
    test = list(
      statistic = test_result$statistic,
      df = df,
      p_value = test_result$p_value,
      alternative = alternative,
      method = "Chi-square test for variance"
    ),
    verification = list(
      verified = verified,
      ratio = variance_ratio,
      cv_ratio = cv_ratio,
      upper_verification_limit = uvl_cv
    ),
    ci = ci,
    settings = list(
      alpha = alpha,
      conf_level = conf_level,
      alternative = alternative
    ),
    call = call
  )
  
  class(result) <- c("verify_precision", "valytics_precision", "valytics_result")
  result
}


# Input Validation ----

#' Validate verify_precision input arguments
#' @noRd
#' @keywords internal
.validate_verify_precision_input <- function(claimed_cv, claimed_sd, alpha,
                                             conf_level) {
  
  # Must have at least one of claimed_cv or claimed_sd
  
  if (is.null(claimed_cv) && is.null(claimed_sd)) {
    stop("Either `claimed_cv` or `claimed_sd` must be provided.", call. = FALSE)
  }
  
  # Validate claimed_cv
  if (!is.null(claimed_cv)) {
    if (!is.numeric(claimed_cv) || length(claimed_cv) != 1 || claimed_cv <= 0) {
      stop("`claimed_cv` must be a single positive number.", call. = FALSE)
    }
  }
  
  # Validate claimed_sd
  if (!is.null(claimed_sd)) {
    if (!is.numeric(claimed_sd) || length(claimed_sd) != 1 || claimed_sd <= 0) {
      stop("`claimed_sd` must be a single positive number.", call. = FALSE)
    }
  }
  
  # Validate alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  
  # Validate conf_level
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number between 0 and 1.", call. = FALSE)
  }
  
  invisible(TRUE)
}


# Input Parsing ----

#' Parse verify_precision input to extract precision statistics
#' @noRd
#' @keywords internal
.parse_verify_precision_input <- function(x, mean_value, value, day, run,
                                          conf_level, ...) {
  
  # Case 1: precision_study object
  if (inherits(x, "precision_study")) {
    return(.parse_precision_study_input(x))
  }
  
  # Case 2: Numeric vector
  if (is.numeric(x) && is.vector(x)) {
    return(.parse_numeric_vector_input(x, mean_value))
  }
  
  # Case 3: Data frame - run precision_study
  if (is.data.frame(x)) {
    return(.parse_dataframe_input(x, value, day, run, conf_level, ...))
  }
  
  # Unknown input type
  stop("`x` must be a numeric vector, precision_study object, or data frame.",
       call. = FALSE)
}


#' Parse precision_study object input
#' @noRd
#' @keywords internal
.parse_precision_study_input <- function(x) {
  
  # Find the within-laboratory precision (intermediate precision) row
  prec <- x$precision
  
  # Look for "Within-laboratory precision" or "Intermediate"
  inter_idx <- which(grepl("Within-laboratory|Intermediate", prec$measure,
                           ignore.case = TRUE))
  
  if (length(inter_idx) == 0) {
    # Fall back to repeatability if no intermediate precision
    inter_idx <- which(prec$measure == "Repeatability")
  }
  
  if (length(inter_idx) == 0) {
    stop("Could not find precision estimate in precision_study object.",
         call. = FALSE)
  }
  
  # Use the first match (within-laboratory precision preferred)
  idx <- inter_idx[1]
  observed_sd <- prec$sd[idx]
  
  # Get mean from the precision study
  if (!is.null(x$by_sample) && length(x$by_sample) > 0) {
    # Multi-sample: use first sample's mean
    first_sample <- x$by_sample[[1]]
    observed_mean <- first_sample$grand_mean
  } else {
    # Single sample: calculate from input data
    observed_mean <- mean(x$input$data[[x$input$value_col]], na.rm = TRUE)
  }
  
  # Degrees of freedom: use total observations minus number of groups
  # This is an approximation - for more precision, use Satterthwaite df
  n <- x$input$n
  n_days <- x$design$levels$day
  
  # For within-laboratory precision, df depends on the design
  # Conservative estimate: n - n_days
  df <- n - n_days
  
  list(
    sd = observed_sd,
    mean = observed_mean,
    n = n,
    df = df,
    source = "precision_study object"
  )
}


#' Parse numeric vector input
#' @noRd
#' @keywords internal
.parse_numeric_vector_input <- function(x, mean_value) {
  
  # Remove NAs
  x <- x[!is.na(x)]
  n <- length(x)
  
  if (n < 3) {
    stop("At least 3 observations are required.", call. = FALSE)
  }
  
  observed_sd <- stats::sd(x)
  observed_mean <- if (!is.null(mean_value)) mean_value else mean(x)
  df <- n - 1
  
  list(
    sd = observed_sd,
    mean = observed_mean,
    n = n,
    df = df,
    source = "numeric vector"
  )
}


#' Parse data frame input - runs precision_study internally
#' @noRd
#' @keywords internal
.parse_dataframe_input <- function(x, value, day, run, conf_level, ...) {
  
  # Run precision_study
  prec_result <- precision_study(
    data = x,
    value = value,
    day = day,
    run = run,
    conf_level = conf_level,
    ...
  )
  
  # Now parse the precision_study result
  parsed <- .parse_precision_study_input(prec_result)
  parsed$source <- "data frame (via precision_study)"
  
  parsed
}


# Claimed Precision Resolution ----

#' Resolve claimed CV and SD from input
#' @noRd
#' @keywords internal
.resolve_claimed_precision <- function(claimed_cv, claimed_sd, mean_value) {
  
  if (!is.null(claimed_cv) && !is.null(claimed_sd)) {
    # Both provided - use as given
    return(list(
      cv_pct = claimed_cv,
      sd = claimed_sd
    ))
  }
  
  if (!is.null(claimed_cv)) {
    # CV provided - calculate SD from mean
    if (is.null(mean_value) || mean_value <= 0) {
      stop("`mean_value` is required to convert CV to SD.", call. = FALSE)
    }
    claimed_sd <- (claimed_cv / 100) * mean_value
    return(list(
      cv_pct = claimed_cv,
      sd = claimed_sd
    ))
  }
  
  if (!is.null(claimed_sd)) {
    # SD provided - calculate CV from mean
    if (is.null(mean_value) || mean_value <= 0) {
      stop("`mean_value` is required to convert SD to CV.", call. = FALSE)
    }
    claimed_cv <- 100 * claimed_sd / mean_value
    return(list(
      cv_pct = claimed_cv,
      sd = claimed_sd
    ))
  }
  
  stop("Either `claimed_cv` or `claimed_sd` must be provided.", call. = FALSE)
}


# Statistical Tests ----

#' Chi-square test for variance
#' @noRd
#' @keywords internal
.chi_square_variance_test <- function(observed_variance, claimed_variance,
                                      df, alpha, alternative) {
  
  # Test statistic: chi_sq = (n-1) * s^2 / sigma0^2 = df * s^2 / sigma0^2
  chi_sq <- df * observed_variance / claimed_variance
  
  # P-value depends on alternative
  p_value <- switch(
    alternative,
    "less" = stats::pchisq(chi_sq, df = df, lower.tail = TRUE),
    "greater" = stats::pchisq(chi_sq, df = df, lower.tail = FALSE),
    "two.sided" = {
      # Two-sided: 2 * min(P(X <= x), P(X >= x))
      p_lower <- stats::pchisq(chi_sq, df = df, lower.tail = TRUE)
      p_upper <- stats::pchisq(chi_sq, df = df, lower.tail = FALSE)
      2 * min(p_lower, p_upper)
    }
  )
  
  list(
    statistic = chi_sq,
    p_value = p_value
  )
}


#' Calculate confidence interval for variance
#' @noRd
#' @keywords internal
.calculate_variance_ci <- function(observed_variance, observed_sd, observed_cv,
                                   df, conf_level) {
  
  alpha <- 1 - conf_level
  
  # Chi-square quantiles
  chi_lower <- stats::qchisq(alpha / 2, df = df)
  chi_upper <- stats::qchisq(1 - alpha / 2, df = df)
  
  # CI for variance: df * s^2 / chi_upper, df * s^2 / chi_lower
  var_lower <- df * observed_variance / chi_upper
  var_upper <- df * observed_variance / chi_lower
  
  # CI for SD
  sd_lower <- sqrt(var_lower)
  sd_upper <- sqrt(var_upper)
  
  # CI for CV (as percent)
  # CV = 100 * SD / mean, so CV scales linearly with SD
  cv_lower <- observed_cv * (sd_lower / observed_sd)
  cv_upper <- observed_cv * (sd_upper / observed_sd)
  
  list(
    variance_ci = c(lower = var_lower, upper = var_upper),
    sd_ci = c(lower = sd_lower, upper = sd_upper),
    cv_ci = c(lower = cv_lower, upper = cv_upper)
  )
}


# S3 Methods ----

#' Print method for verify_precision objects
#'
#' @description
#' Displays verification results including observed vs. claimed precision,
#' test statistics, and verification outcome.
#'
#' @param x An object of class `verify_precision`.
#' @param digits Number of significant digits to display (default: 3).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @examples
#' set.seed(42)
#' measurements <- rnorm(25, mean = 100, sd = 3.5)
#' result <- verify_precision(measurements, claimed_cv = 4, mean_value = 100)
#' print(result)
#'
#' @seealso [summary.verify_precision()] for detailed output
#' @export
print.verify_precision <- function(x, digits = 3, ...) {
  
  cat("\n")
  cat("Precision Verification\n")
  cat(strrep("-", 45), "\n")
  
  # Input info
  cat(sprintf("n = %d observations (df = %d)\n", x$input$n, x$input$df))
  cat(sprintf("Mean: %s\n", format(round(x$input$mean_value, digits), nsmall = digits)))
  cat(sprintf("Source: %s\n\n", x$input$source))
  
  # Comparison table
  cat("Precision Comparison:\n")
  cat(strrep("-", 45), "\n")
  
  cat(sprintf("  %-15s %10s %10s\n", "", "Observed", "Claimed"))
  cat(sprintf("  %-15s %10s %10s\n", "CV (%):",
              format(round(x$observed$cv_pct, 2), nsmall = 2),
              format(round(x$claimed$cv_pct, 2), nsmall = 2)))
  cat(sprintf("  %-15s %10s %10s\n", "SD:",
              format(round(x$observed$sd, digits), nsmall = digits),
              format(round(x$claimed$sd, digits), nsmall = digits)))
  cat("\n")
  
  # Ratios
  cat(sprintf("  CV ratio (observed/claimed): %.2f\n", x$verification$cv_ratio))
  cat(sprintf("  Upper verification limit (CV): %s%%\n",
              format(round(x$verification$upper_verification_limit, 2), nsmall = 2)))
  cat("\n")
  
  # Test results
  cat("Statistical Test:\n")
  cat(strrep("-", 45), "\n")
  cat(sprintf("  %s\n", x$test$method))
  cat(sprintf("  Chi-square statistic: %.3f (df = %d)\n",
              x$test$statistic, x$test$df))
  cat(sprintf("  p-value: %s\n", format(x$test$p_value, digits = 4)))
  cat(sprintf("  alpha: %s\n", x$settings$alpha))
  cat("\n")
  
  # Verification outcome
  cat("Verification Result:\n")
  cat(strrep("-", 45), "\n")
  
  if (x$verification$verified) {
    cat("  VERIFIED: Observed precision meets claimed specification.\n")
    cat(sprintf("  (Observed CV %.2f%% <= Upper limit %.2f%%)\n",
                x$observed$cv_pct, x$verification$upper_verification_limit))
  } else {
    cat("  NOT VERIFIED: Observed precision exceeds claimed specification.\n")
    cat(sprintf("  (Observed CV %.2f%% > Upper limit %.2f%%)\n",
                x$observed$cv_pct, x$verification$upper_verification_limit))
  }
  
  cat("\n")
  
  invisible(x)
}


#' Summary method for verify_precision objects
#'
#' @description
#' Provides a detailed summary of precision verification results, including
#' confidence intervals, test details, and interpretation guidance.
#'
#' @param object An object of class `verify_precision`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class `summary.verify_precision` containing:
#'   \describe{
#'     \item{call}{The original function call.}
#'     \item{input}{Input information.}
#'     \item{observed}{Observed precision statistics.}
#'     \item{claimed}{Claimed precision statistics.}
#'     \item{test}{Hypothesis test results.}
#'     \item{verification}{Verification outcome.}
#'     \item{ci}{Confidence intervals.}
#'     \item{settings}{Analysis settings.}
#'   }
#'
#' @examples
#' set.seed(42)
#' measurements <- rnorm(25, mean = 100, sd = 3.5)
#' result <- verify_precision(measurements, claimed_cv = 4, mean_value = 100)
#' summary(result)
#'
#' @seealso [print.verify_precision()] for concise output
#' @export
summary.verify_precision <- function(object, ...) {
  
  structure(
    list(
      call = object$call,
      input = object$input,
      observed = object$observed,
      claimed = object$claimed,
      test = object$test,
      verification = object$verification,
      ci = object$ci,
      settings = object$settings
    ),
    class = "summary.verify_precision"
  )
}


#' Print method for summary.verify_precision objects
#'
#' @param x An object of class `summary.verify_precision`.
#' @param digits Number of significant digits to display (default: 4).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object `x`.
#'
#' @export
print.summary.verify_precision <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat("Precision Verification - Detailed Summary\n")
  cat(strrep("=", 55), "\n\n")
  
  # Call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  # Input Information
  cat(strrep("-", 55), "\n")
  cat("Input Information:\n")
  cat(strrep("-", 55), "\n")
  cat(sprintf("  Source: %s\n", x$input$source))
  cat(sprintf("  Number of observations: %d\n", x$input$n))
  cat(sprintf("  Degrees of freedom: %d\n", x$input$df))
  cat(sprintf("  Mean value: %s\n",
              format(round(x$input$mean_value, digits), nsmall = digits)))
  cat("\n")
  
  # Settings
  cat(strrep("-", 55), "\n")
  cat("Analysis Settings:\n")
  cat(strrep("-", 55), "\n")
  cat(sprintf("  Significance level (alpha): %g\n", x$settings$alpha))
  cat(sprintf("  Confidence level: %g%%\n", x$settings$conf_level * 100))
  cat(sprintf("  Alternative hypothesis: %s\n", x$settings$alternative))
  cat("\n")
  
  # Precision Comparison
  cat(strrep("-", 55), "\n")
  cat("Precision Comparison:\n")
  cat(strrep("-", 55), "\n")
  
  ci_pct <- sprintf("%g%%", x$settings$conf_level * 100)
  
  # Build comparison table
  cat(sprintf("  %-20s %12s %12s\n", "", "Observed", "Claimed"))
  cat(sprintf("  %-20s %12s %12s\n",
              "CV (%):",
              format(round(x$observed$cv_pct, 2), nsmall = 2),
              format(round(x$claimed$cv_pct, 2), nsmall = 2)))
  cat(sprintf("  %-20s %12s %12s\n",
              "SD:",
              format(round(x$observed$sd, digits), nsmall = digits),
              format(round(x$claimed$sd, digits), nsmall = digits)))
  cat(sprintf("  %-20s %12s %12s\n",
              "Variance:",
              format(round(x$observed$variance, digits), nsmall = digits),
              format(round(x$claimed$variance, digits), nsmall = digits)))
  cat("\n")
  
  # Confidence Intervals
  cat(sprintf("  %s Confidence Intervals (Observed):\n", ci_pct))
  cat(sprintf("    CV:       [%s, %s]%%\n",
              format(round(x$ci$cv_ci["lower"], 2), nsmall = 2),
              format(round(x$ci$cv_ci["upper"], 2), nsmall = 2)))
  cat(sprintf("    SD:       [%s, %s]\n",
              format(round(x$ci$sd_ci["lower"], digits), nsmall = digits),
              format(round(x$ci$sd_ci["upper"], digits), nsmall = digits)))
  cat(sprintf("    Variance: [%s, %s]\n",
              format(round(x$ci$variance_ci["lower"], digits), nsmall = digits),
              format(round(x$ci$variance_ci["upper"], digits), nsmall = digits)))
  cat("\n")
  
  # Ratios
  cat(strrep("-", 55), "\n")
  cat("Ratios:\n")
  cat(strrep("-", 55), "\n")
  cat(sprintf("  CV ratio (observed/claimed): %.3f\n", x$verification$cv_ratio))
  cat(sprintf("  Variance ratio (observed/claimed): %.3f\n", x$verification$ratio))
  cat(sprintf("  Upper verification limit (CV): %s%%\n",
              format(round(x$verification$upper_verification_limit, 2), nsmall = 2)))
  cat("\n")
  
  # Statistical Test
  cat(strrep("-", 55), "\n")
  cat("Statistical Test:\n")
  cat(strrep("-", 55), "\n")
  cat(sprintf("  Method: %s\n", x$test$method))
  cat(sprintf("  Test statistic (chi-square): %.4f\n", x$test$statistic))
  cat(sprintf("  Degrees of freedom: %d\n", x$test$df))
  cat(sprintf("  p-value: %s\n", format(x$test$p_value, digits = 4)))
  cat(sprintf("  Alternative hypothesis: %s\n", x$test$alternative))
  cat("\n")
  
  # Verification Outcome
  cat(strrep("-", 55), "\n")
  cat("Verification Outcome:\n")
  cat(strrep("-", 55), "\n")
  
  if (x$verification$verified) {
    cat("  Result: VERIFIED\n")
    cat("\n")
    cat("  Interpretation:\n")
    cat("    The observed precision meets the manufacturer's claimed\n")
    cat("    specification. The observed CV does not exceed the upper\n")
    cat("    verification limit derived from the claimed CV.\n")
  } else {
    cat("  Result: NOT VERIFIED\n")
    cat("\n")
    cat("  Interpretation:\n")
    cat("    The observed precision exceeds the manufacturer's claimed\n")
    cat("    specification. The observed CV is higher than the upper\n")
    cat("    verification limit, suggesting performance may not meet\n")
    cat("    the claimed precision.\n")
    cat("\n")
    cat("  Possible Actions:\n")
    cat("    - Review experimental conditions and sample handling\n")
    cat("    - Check for outliers or systematic errors\n")
    cat("    - Repeat the study with additional samples\n")
    cat("    - Contact manufacturer if results are reproducible\n")
  }
  
  cat("\n")
  
  invisible(x)
}