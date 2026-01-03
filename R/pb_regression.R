#' Passing-Bablok Regression for Method Comparison
#'
#' @description
#' Performs Passing-Bablok regression to assess agreement between two measurement
#' methods. This non-parametric regression method is robust to outliers and does
#' not assume normally distributed errors. The implementation uses a fast
#' O(n log n) algorithm from the `robslopes` package for point estimation.
#'
#' @param x Numeric vector of measurements from method 1 (reference method),
#'   or a formula of the form `method1 ~ method2`.
#' @param y Numeric vector of measurements from method 2 (test method).
#'   Ignored if `x` is a formula.
#' @param data Optional data frame containing the variables specified in
#'   `x` and `y` (or in the formula).
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param ci_method Method for calculating confidence intervals:
#'   `"analytical"` (default) uses the original Passing-Bablok (1983) method,
#'   `"bootstrap"` uses BCa bootstrap resampling.
#' @param boot_n Number of bootstrap resamples when `ci_method = "bootstrap"`
#'   (default: 1999).
#' @param na_action How to handle missing values: `"omit"` (default) removes
#'   pairs with any NA, `"fail"` stops with an error.
#'
#' @return An object of class `c("pb_regression", "valytics_comparison", "valytics_result")`,
#'   which is a list containing:
#'
#'   \describe{
#'     \item{input}{List with original data and metadata:
#'       \itemize{
#'         \item `x`, `y`: Numeric vectors (after NA handling)
#'         \item `n`: Number of paired observations
#'         \item `n_excluded`: Number of pairs excluded due to NAs
#'         \item `var_names`: Named character vector with variable names
#'       }
#'     }
#'     \item{results}{List with statistical results:
#'       \itemize{
#'         \item `intercept`: Intercept point estimate
#'         \item `slope`: Slope point estimate
#'         \item `intercept_ci`: Named numeric vector with lower and upper CI
#'         \item `slope_ci`: Named numeric vector with lower and upper CI
#'         \item `residuals`: Perpendicular residuals
#'         \item `fitted_x`: Fitted x values
#'         \item `fitted_y`: Fitted y values
#'       }
#'     }
#'     \item{cusum}{List with CUSUM linearity test results (if calculable):
#'       \itemize{
#'         \item `statistic`: CUSUM test statistic
#'         \item `critical_value`: Critical value at alpha = 0.05
#'         \item `p_value`: Approximate p-value
#'         \item `linear`: Logical; TRUE if linearity assumption holds
#'       }
#'     }
#'     \item{settings}{List with analysis settings:
#'       \itemize{
#'         \item `conf_level`: Confidence level used
#'         \item `ci_method`: CI method used
#'         \item `boot_n`: Number of bootstrap samples (if applicable
#'       }
#'     }
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' Passing-Bablok regression is a non-parametric method for fitting a linear
#' relationship between two measurement methods. Unlike ordinary least squares,
#' it:
#'
#' \itemize{
#'   \item Makes no assumptions about error distribution
#'   \item Accounts for measurement error in both variables
#'   \item Is robust to outliers
#'   \item Produces results independent of which variable is assigned to X or Y
#'     (when using the equivariant form)
#' }
#'
#' The slope is estimated as the median of all pairwise slopes (in absolute
#' value for the equivariant version), and the intercept is the median of
#' `y - slope * x`.
#'
#' @section Interpretation:
#' \itemize{
#'   \item **Slope = 1**: No proportional difference between methods
#'   \item **Slope != 1**: Proportional (multiplicative) difference exists
#'   \item **Intercept = 0**: No constant difference between methods
#'   \item **Intercept != 0**: Constant (additive) difference exists
#' }
#'
#' Use the confidence intervals to test these hypotheses: if 1 is within the
#' slope CI and 0 is within the intercept CI, the methods are considered
#' equivalent.
#'
#' @section Assumptions:
#' \itemize{
#'   \item Linear relationship between X and Y (test with CUSUM)
#'   \item Measurement range covers the intended clinical range
#'   \item Data are continuously distributed
#' }
#'
#' @section CUSUM Test for Linearity:
#' The CUSUM (cumulative sum) test checks the linearity assumption. A
#' significant result (p < 0.05) suggests non-linearity, and Passing-Bablok
#' regression may not be appropriate.
#'
#' @references
#' Passing H, Bablok W (1983). A new biometrical procedure for testing the
#' equality of measurements from two different analytical methods. Application
#' of linear regression procedures for method comparison studies in clinical
#' chemistry, Part I. \emph{Journal of Clinical Chemistry and Clinical
#' Biochemistry}, 21(11):709-720. \doi{10.1515/cclm.1983.21.11.709}
#'
#' Passing H, Bablok W (1984). Comparison of several regression procedures for
#' method comparison studies and determination of sample sizes. Application of
#' linear regression procedures for method comparison studies in Clinical
#' Chemistry, Part II. \emph{Journal of Clinical Chemistry and Clinical
#' Biochemistry}, 22(6):431-445. \doi{10.1515/cclm.1984.22.6.431}
#'
#' Bablok W, Passing H, Bender R, Schneider B (1988). A general regression
#' procedure for method transformation. Application of linear regression
#' procedures for method comparison studies in clinical chemistry, Part III.
#' \emph{Journal of Clinical Chemistry and Clinical Biochemistry},
#' 26(11):783-790. \doi{10.1515/cclm.1988.26.11.783}
#'
#' Raymaekers J, Dufey F (2022). Equivariant Passing-Bablok regression in
#' quasilinear time. \emph{arXiv preprint}. \doi{10.48550/arXiv.2202.08060}
#'
#' @seealso
#' [plot.pb_regression()] for visualization,
#' [summary.pb_regression()] for detailed summary,
#' [ba_analysis()] for Bland-Altman analysis
#'
#' @examples
#' # Simulated method comparison data
#' set.seed(42)
#' method_a <- rnorm(50, mean = 100, sd = 15)
#' method_b <- 1.05 * method_a + 3 + rnorm(50, sd = 5)  # slope=1.05, intercept=3
#'
#' # Basic analysis
#' pb <- pb_regression(method_a, method_b)
#' pb
#'
#' # Using formula interface with data frame
#' df <- data.frame(reference = method_a, test = method_b)
#' pb <- pb_regression(reference ~ test, data = df)
#'
#' # With bootstrap confidence intervals
#' pb_boot <- pb_regression(method_a, method_b, ci_method = "bootstrap")
#'
#' # Using package example data
#' data(glucose_methods)
#' pb <- pb_regression(reference ~ poc_meter, data = glucose_methods)
#' summary(pb)
#' plot(pb)
#'
#' @export
pb_regression <- function(x,
                          y = NULL,
                          data = NULL,
                          conf_level = 0.95,
                          ci_method = c("analytical", "bootstrap"),
                          boot_n = 1999,
                          na_action = c("omit", "fail")) {

  # Capture the call
  call <- match.call()

  # Match arguments
  ci_method <- match.arg(ci_method)
  na_action <- match.arg(na_action)

  # Input parsing ----
  parsed <- .parse_pb_input(x, y, data)
  x_vec <- parsed$x
  y_vec <- parsed$y
  var_names <- parsed$var_names

  # Handle missing values ----
  complete <- complete.cases(x_vec, y_vec)
  n_total <- length(x_vec)
  n_excluded <- sum(!complete)

  if (na_action == "fail" && n_excluded > 0) {
    stop("Missing values detected. Use na_action = 'omit' to exclude them.",
         call. = FALSE)
  }

  x_vec <- x_vec[complete]
  y_vec <- y_vec[complete]
  n <- length(x_vec)

  # Validate sample size
  if (n < 10) {
    stop("At least 10 complete paired observations are required for ",
         "Passing-Bablok regression. Found: ", n, call. = FALSE)
  }

  # Validate inputs ----
  .validate_pb_inputs(x_vec, y_vec, conf_level, boot_n)

  # Point estimation ----
  pb_fit <- robslopes::PassingBablok(x_vec, y_vec, alpha = NULL, verbose = FALSE)

  slope <- pb_fit$slope
  intercept <- pb_fit$intercept

  # Confidence intervals ----
  if (ci_method == "analytical") {
    ci_result <- .pb_analytical_ci(x_vec, y_vec, slope, intercept, conf_level)
  } else {
    ci_result <- .pb_bootstrap_ci(x_vec, y_vec, conf_level, boot_n)
  }

  # Residuals and fitted values ----
  fitted_y <- intercept + slope * x_vec
  residuals_y <- y_vec - fitted_y

  # Perpendicular projection onto the regression line
  # Point (x, y) projects to (x_hat, y_hat) on line y = a + b*x
  # x_hat = (x + b*(y - a)) / (1 + b^2)
  fitted_x <- (x_vec + slope * (y_vec - intercept)) / (1 + slope^2)
  fitted_y_perp <- intercept + slope * fitted_x

  # Perpendicular distance (signed)
  perp_residuals <- sign(y_vec - fitted_y) *
    sqrt((x_vec - fitted_x)^2 + (y_vec - fitted_y_perp)^2)

  # CUSUM test for linearity ----
  cusum_result <- .pb_cusum_test(x_vec, y_vec, intercept, slope)

  # Assemble output ----
  result <- list(
    input = list(
      x = x_vec,
      y = y_vec,
      n = n,
      n_excluded = n_excluded,
      var_names = var_names
    ),
    results = list(
      intercept = intercept,
      slope = slope,
      intercept_ci = ci_result$intercept_ci,
      slope_ci = ci_result$slope_ci,
      residuals = perp_residuals,
      fitted_x = fitted_x,
      fitted_y = fitted_y_perp
    ),
    cusum = cusum_result,
    settings = list(
      conf_level = conf_level,
      ci_method = ci_method,
      boot_n = if (ci_method == "bootstrap") boot_n else NA
    ),
    call = call
  )

  class(result) <- c("pb_regression", "valytics_comparison", "valytics_result")
  result
}


# Helper Functions ----

#' Parse input for pb_regression (formula or vectors)
#' @noRd
.parse_pb_input <- function(x, y, data) {

  # Formula interface
  if (inherits(x, "formula")) {
    if (length(x) != 3) {
      stop("Formula must be of the form 'y ~ x'.", call. = FALSE)
    }

    # Extract variable names from formula
    vars <- all.vars(x)
    if (length(vars) != 2) {
      stop("Formula must contain exactly two variables.", call. = FALSE)
    }

    # Get data
    if (is.null(data)) {
      # Try to get from parent frame
      x_vec <- eval(x[[3]], envir = parent.frame(2))
      y_vec <- eval(x[[2]], envir = parent.frame(2))
    } else {
      if (!is.data.frame(data)) {
        stop("'data' must be a data frame.", call. = FALSE)
      }
      x_vec <- eval(x[[3]], envir = data)
      y_vec <- eval(x[[2]], envir = data)
    }

    var_names <- c(x = vars[2], y = vars[1])

  } else {
    # Vector interface
    if (is.null(y)) {
      stop("'y' must be provided when 'x' is not a formula.", call. = FALSE)
    }

    x_vec <- x
    y_vec <- y

    # Try to get variable names from call
    x_name <- deparse(substitute(x, env = parent.frame()))
    y_name <- deparse(substitute(y, env = parent.frame()))

    var_names <- c(x = x_name, y = y_name)
  }

  # Validate vectors
  if (!is.numeric(x_vec) || !is.numeric(y_vec)) {
    stop("Both 'x' and 'y' must be numeric vectors.", call. = FALSE)
  }

  if (length(x_vec) != length(y_vec)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }

  list(x = as.numeric(x_vec), y = as.numeric(y_vec), var_names = var_names)
}


#' Validate Passing-Bablok inputs
#' @noRd
.validate_pb_inputs <- function(x, y, conf_level, boot_n) {

  # Check for non-negative values (PB requirement)
  if (any(x < 0, na.rm = TRUE) || any(y < 0, na.rm = TRUE)) {
    warning("Passing-Bablok regression is designed for non-negative values. ",
            "Results may be unreliable with negative data.", call. = FALSE)
  }

  # Check confidence level
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number between 0 and 1.", call. = FALSE)
  }

  # Check bootstrap sample count
  if (!is.numeric(boot_n) || length(boot_n) != 1 ||
      boot_n < 100 || boot_n != floor(boot_n)) {
    stop("'boot_n' must be an integer >= 100.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Analytical confidence intervals (Passing-Bablok 1983)
#' @noRd
.pb_analytical_ci <- function(x, y, slope, intercept, conf_level) {

  n <- length(x)
  alpha <- 1 - conf_level

  # Calculate all pairwise slopes
  # Using vectorized approach for efficiency
  slopes <- numeric(n * (n - 1) / 2)
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      dx <- x[i] - x[j]
      dy <- y[i] - y[j]

      if (abs(dx) < .Machine$double.eps) {
        # Vertical line: slope = Inf or -Inf
        slopes[k] <- sign(dy) * Inf
      } else {
        slopes[k] <- dy / dx
      }
    }
  }

  # For equivariant PB, we use absolute values of slopes
  slopes_abs <- abs(slopes[is.finite(slopes)])
  N <- length(slopes_abs)

  if (N < 10) {
    warning("Too few valid pairwise slopes for reliable CI calculation.",
            call. = FALSE)
    return(list(
      slope_ci = c(lower = NA_real_, upper = NA_real_),
      intercept_ci = c(lower = NA_real_, upper = NA_real_)
    ))
  }

  # Sort slopes for quantile extraction
  slopes_sorted <- sort(slopes_abs)

  # CI for slope (Passing & Bablok 1983, equations 16-18)
  # The CI is based on order statistics of the slope distribution
  # Using normal approximation for the confidence bounds

  # Variance of the median position
  # C_alpha = z_{1-alpha/2} * sqrt(N) / 2
  z_alpha <- qnorm(1 - alpha / 2)
  C_alpha <- z_alpha * sqrt(N) / 2

  # Lower and upper indices for slope CI
  M <- (N + 1) / 2  # median position
  lower_idx <- max(1, floor(M - C_alpha))
  upper_idx <- min(N, ceiling(M + C_alpha))

  slope_lower <- slopes_sorted[lower_idx]
  slope_upper <- slopes_sorted[upper_idx]

  # CI for intercept
  # Intercept CI bounds are calculated using slope CI bounds
  # B0_lower = median(y - slope_upper * x)
  # B0_upper = median(y - slope_lower * x)
  intercept_lower <- median(y - slope_upper * x)
  intercept_upper <- median(y - slope_lower * x)

  # Ensure correct ordering
  if (intercept_lower > intercept_upper) {
    tmp <- intercept_lower
    intercept_lower <- intercept_upper
    intercept_upper <- tmp
  }

  list(
    slope_ci = c(lower = slope_lower, upper = slope_upper),
    intercept_ci = c(lower = intercept_lower, upper = intercept_upper)
  )
}


#' Bootstrap confidence intervals (BCa method)
#' @noRd
.pb_bootstrap_ci <- function(x, y, conf_level, boot_n) {

  n <- length(x)
  alpha <- 1 - conf_level

  # Storage for bootstrap estimates
  boot_slopes <- numeric(boot_n)
  boot_intercepts <- numeric(boot_n)

  # Generate bootstrap samples
  for (b in seq_len(boot_n)) {
    # Sample with replacement
    idx <- sample.int(n, n, replace = TRUE)
    x_boot <- x[idx]
    y_boot <- y[idx]

    # Fit using fast algorithm
    fit_boot <- robslopes::PassingBablok(x_boot, y_boot,
                                         alpha = NULL, verbose = FALSE)
    boot_slopes[b] <- fit_boot$slope
    boot_intercepts[b] <- fit_boot$intercept
  }

  # BCa (Bias-Corrected and Accelerated) confidence intervals
  # Calculate bias correction factor (z0)
  orig_fit <- robslopes::PassingBablok(x, y, alpha = NULL, verbose = FALSE)
  slope_orig <- orig_fit$slope
  intercept_orig <- orig_fit$intercept

  # Slope BCa CI
  slope_ci <- .bca_ci(boot_slopes, slope_orig, x, y, "slope",
                      conf_level, boot_n)

  # Intercept BCa CI
  intercept_ci <- .bca_ci(boot_intercepts, intercept_orig, x, y, "intercept",
                          conf_level, boot_n)

  list(
    slope_ci = slope_ci,
    intercept_ci = intercept_ci
  )
}


#' Calculate BCa confidence interval
#' @noRd
.bca_ci <- function(boot_stat, orig_stat, x, y, param, conf_level, boot_n) {

  alpha <- 1 - conf_level
  n <- length(x)

  # Remove any NA/Inf values
  boot_stat <- boot_stat[is.finite(boot_stat)]
  if (length(boot_stat) < boot_n * 0.9) {
    warning("Many bootstrap samples failed. CI may be unreliable.", call. = FALSE)
  }

  # Bias correction factor z0
  prop_less <- mean(boot_stat < orig_stat)
  # Handle edge cases
  prop_less <- max(0.001, min(0.999, prop_less))
  z0 <- qnorm(prop_less)

  # Acceleration factor (a) using jackknife
  jack_stat <- numeric(n)
  for (i in seq_len(n)) {
    x_jack <- x[-i]
    y_jack <- y[-i]
    fit_jack <- robslopes::PassingBablok(x_jack, y_jack,
                                         alpha = NULL, verbose = FALSE)
    jack_stat[i] <- if (param == "slope") fit_jack$slope else fit_jack$intercept
  }

  jack_mean <- mean(jack_stat)
  jack_diff <- jack_mean - jack_stat

  # Acceleration factor
  a <- sum(jack_diff^3) / (6 * (sum(jack_diff^2))^1.5)

  # Handle edge cases
  if (!is.finite(a)) a <- 0

  # Adjusted percentiles
  z_alpha_lower <- qnorm(alpha / 2)
  z_alpha_upper <- qnorm(1 - alpha / 2)

  # BCa adjusted percentiles
  alpha1 <- pnorm(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)))
  alpha2 <- pnorm(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)))

  # Bound the percentiles
  alpha1 <- max(0.001, min(0.999, alpha1))
  alpha2 <- max(0.001, min(0.999, alpha2))

  # Get quantiles
  ci <- quantile(boot_stat, probs = c(alpha1, alpha2), na.rm = TRUE)
  names(ci) <- c("lower", "upper")

  ci
}


#' CUSUM test for linearity (Passing & Bablok 1983)
#' @noRd
.pb_cusum_test <- function(x, y, intercept, slope) {

  n <- length(x)

  # Calculate residual signs
  fitted_y <- intercept + slope * x
  residuals <- y - fitted_y

  # Count positive and negative residuals
  n_pos <- sum(residuals > 0)
  n_neg <- sum(residuals < 0)
  n_zero <- sum(residuals == 0)

  # Handle ties (residuals == 0)
  if (n_zero > 0) {
    # Distribute zeros proportionally
    n_pos <- n_pos + n_zero / 2
    n_neg <- n_neg + n_zero / 2
  }

  # Assign scores
  # If y > fitted: score = sqrt(n_neg / n_pos)
  # If y < fitted: score = -sqrt(n_pos / n_neg)
  if (n_pos == 0 || n_neg == 0) {
    # All residuals same sign - cannot compute CUSUM
    return(list(
      statistic = NA_real_,
      critical_value = 1.36,
      p_value = NA_real_,
      linear = NA
    ))
  }

  score_pos <- sqrt(n_neg / n_pos)
  score_neg <- -sqrt(n_pos / n_neg)

  scores <- ifelse(residuals > 0, score_pos,
                   ifelse(residuals < 0, score_neg, 0))

  # Sort by x values
  ord <- order(x)
  scores_sorted <- scores[ord]

  # Calculate cumulative sums
  cusum <- cumsum(scores_sorted)

  # CUSUM statistic is max |cusum| / sqrt(n)
  H <- max(abs(cusum)) / sqrt(n)

  # Critical value at alpha = 0.05 is 1.36 (Kolmogorov-Smirnov)
  # At alpha = 0.01 is 1.63
  critical_05 <- 1.36

  # Approximate p-value using Kolmogorov-Smirnov distribution
  # P(H > h) ~ 2 * sum_{k=1}^{inf} (-1)^{k+1} * exp(-2 * k^2 * h^2)
  p_value <- .ks_p_value(H)

  list(
    statistic = H,
    critical_value = critical_05,
    p_value = p_value,
    linear = (H <= critical_05)
  )
}


#' Approximate p-value from Kolmogorov-Smirnov distribution
#' @noRd
.ks_p_value <- function(h) {
  if (h <= 0) return(1)
  if (h > 3) return(0)

  # Use series approximation
  # P(H > h) = 2 * sum_{k=1}^{inf} (-1)^{k+1} * exp(-2 * k^2 * h^2)
  p <- 0
  for (k in 1:100) {
    term <- ((-1)^(k + 1)) * exp(-2 * k^2 * h^2)
    p <- p + term
    if (abs(term) < 1e-10) break
  }
  p <- 2 * p
  max(0, min(1, p))
}
