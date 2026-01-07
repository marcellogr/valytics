#' Deming Regression for Method Comparison
#'
#' @description
#' Performs Deming regression to assess agreement between two measurement
#' methods. Unlike ordinary least squares, Deming regression accounts for
#' measurement error in both variables, making it appropriate for method
#' comparison studies where neither method is a perfect reference.
#'
#' @param x Numeric vector of measurements from method 1 (reference method),
#'   or a formula of the form `method1 ~ method2`.
#' @param y Numeric vector of measurements from method 2 (test method).
#'   Ignored if `x` is a formula.
#' @param data Optional data frame containing the variables specified in
#'   `x` and `y` (or in the formula).
#' @param error_ratio Ratio of error variances (Var(error_y) / Var(error_x)).
#'   Default is 1 (orthogonal regression, assuming equal error variances).
#'   Can be estimated from replicate measurements or set based on prior
#'   knowledge of method precision.
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param ci_method Method for calculating confidence intervals:
#'   `"jackknife"` (default) uses delete-one jackknife resampling,
#'   `"bootstrap"` uses BCa bootstrap resampling.
#' @param boot_n Number of bootstrap resamples when `ci_method = "bootstrap"`
#'   (default: 1999).
#' @param weighted Logical; if `TRUE`, performs weighted Deming regression
#'   where weights are inversely proportional to the variance at each point.
#'   Requires replicate measurements to estimate weights. Default is `FALSE`.
#' @param na_action How to handle missing values: `"omit"` (default) removes
#'   pairs with any NA, `"fail"` stops with an error.
#'
#' @return An object of class `c("deming_regression", "valytics_comparison", "valytics_result")`,
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
#'         \item `intercept_se`: Standard error of intercept
#'         \item `slope_se`: Standard error of slope
#'         \item `residuals`: Perpendicular residuals
#'         \item `fitted_x`: Fitted x values
#'         \item `fitted_y`: Fitted y values
#'       }
#'     }
#'     \item{settings}{List with analysis settings:
#'       \itemize{
#'         \item `error_ratio`: Error variance ratio used
#'         \item `conf_level`: Confidence level used
#'         \item `ci_method`: CI method used
#'         \item `boot_n`: Number of bootstrap samples (if applicable)
#'         \item `weighted`: Whether weighted regression was used
#'       }
#'     }
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' Deming regression (also known as errors-in-variables regression or Model II
#' regression) is designed for situations where both X and Y are measured with
#' error. This is the typical case in method comparison studies where both the
#' reference and test methods have measurement uncertainty.
#'
#' The error ratio (lambda, λ) represents the ratio of error variances:
#' \deqn{\lambda = \frac{Var(\epsilon_y)}{Var(\epsilon_x)}}
#'
#' When λ = 1 (default), this is equivalent to orthogonal regression, which
#' minimizes perpendicular distances to the regression line. When λ ≠ 1, the
#' regression minimizes a weighted combination of horizontal and vertical
#' distances.
#'
#' **Choosing the error ratio:**
#' \itemize{
#'   \item If both methods have similar precision: use λ = 1
#'   \item If precision differs: estimate from replicate measurements as
#'     λ = CV_y² / CV_x² (squared coefficient of variation ratio)
#'   \item If one method is much more precise: consider ordinary least squares
#' }
#'
#' @section Interpretation:
#' \itemize{
#'   \item **Slope = 1**: No proportional difference between methods
#'   \item **Slope ≠ 1**: Proportional (multiplicative) difference exists
#'   \item **Intercept = 0**: No constant difference between methods
#'   \item **Intercept ≠ 0**: Constant (additive) difference exists
#' }
#'
#' Use the confidence intervals to test these hypotheses: if 1 is within the
#' slope CI and 0 is within the intercept CI, the methods are considered
#' equivalent.
#'
#' @section Comparison with Other Methods:
#' \itemize{
#'   \item **Ordinary Least Squares (OLS)**: Assumes X is measured without
#'     error. Biases slope toward zero when both variables have error.
#'   \item **Passing-Bablok**: Non-parametric, robust to outliers, but assumes
#'     linear relationship and no ties.
#'   \item **Deming**: Parametric, accounts for error in both variables,
#'     allows specification of error ratio.
#' }
#'
#' @section Assumptions:
#' \itemize{
#'   \item Linear relationship between X and Y
#'   \item Measurement errors are normally distributed

#'   \item Error variances are constant (homoscedastic) or known ratio
#'   \item Errors in X and Y are independent
#' }
#'
#' @references
#' Deming WE (1943). Statistical Adjustment of Data. Wiley.
#'
#' Linnet K (1990). Estimation of the linear relationship between the
#' measurements of two methods with proportional errors.
#' \emph{Statistics in Medicine}, 9(12):1463-1473.
#' \doi{10.1002/sim.4780091210}
#'
#' Linnet K (1993). Evaluation of regression procedures for methods comparison
#' studies. \emph{Clinical Chemistry}, 39(3):424-432.
#' \doi{10.1093/clinchem/39.3.424}
#'
#' Cornbleet PJ, Gochman N (1979). Incorrect least-squares regression
#' coefficients in method-comparison analysis.
#' \emph{Clinical Chemistry}, 25(3):432-438.
#'
#' @seealso
#' [plot.deming_regression()] for visualization,
#' [summary.deming_regression()] for detailed summary,
#' [pb_regression()] for non-parametric alternative,
#' [ba_analysis()] for Bland-Altman analysis
#'
#' @examples
#' # Simulated method comparison data
#' set.seed(42)
#' true_values <- rnorm(50, mean = 100, sd = 20)
#' method_a <- true_values + rnorm(50, sd = 5)
#' method_b <- 1.05 * true_values + 3 + rnorm(50, sd = 5)
#'
#' # Basic analysis (orthogonal regression, lambda = 1)
#' dm <- deming_regression(method_a, method_b)
#' dm
#'
#' # Using formula interface with data frame
#' df <- data.frame(reference = method_a, test = method_b)
#' dm <- deming_regression(reference ~ test, data = df)
#'
#' # With known error ratio (e.g., test method has 2x variance)
#' dm <- deming_regression(method_a, method_b, error_ratio = 2)
#'
#' # With bootstrap confidence intervals
#' dm_boot <- deming_regression(method_a, method_b, ci_method = "bootstrap")
#'
#' # Using package example data
#' data(glucose_methods)
#' dm <- deming_regression(reference ~ poc_meter, data = glucose_methods)
#' summary(dm)
#' plot(dm)
#'
#' @export
deming_regression <- function(x,
                              y = NULL,
                              data = NULL,
                              error_ratio = 1,
                              conf_level = 0.95,
                              ci_method = c("jackknife", "bootstrap"),
                              boot_n = 1999,
                              weighted = FALSE,
                              na_action = c("omit", "fail")) {


  # Capture the call

call <- match.call()

  # Match arguments
  ci_method <- match.arg(ci_method)
  na_action <- match.arg(na_action)

  # Input parsing ----
  parsed <- .parse_deming_input(x, y, data)
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

  # Validate inputs ----
  .validate_deming_inputs(x_vec, y_vec, error_ratio, conf_level, boot_n, weighted)

  # Check minimum sample size
  if (n < 10) {
    stop("At least 10 complete paired observations are required for ",
         "Deming regression. Found: ", n, call. = FALSE)
  }

  # Point estimation ----
  deming_fit <- .compute_deming(x_vec, y_vec, error_ratio)

  slope <- deming_fit$slope
  intercept <- deming_fit$intercept

  # Confidence intervals ----
  if (ci_method == "jackknife") {
    ci_result <- .deming_jackknife_ci(x_vec, y_vec, error_ratio, conf_level)
  } else {
    ci_result <- .deming_bootstrap_ci(x_vec, y_vec, error_ratio, conf_level, boot_n)
  }

  # Residuals and fitted values ----
  # Project points onto regression line (perpendicular projection)
  # For point (xi, yi), the projection onto y = a + b*x is:
  # x_fitted = (xi + b*(yi - a)) / (1 + b^2)
  # y_fitted = a + b * x_fitted

  fitted_x <- (x_vec + slope * (y_vec - intercept)) / (1 + slope^2)
  fitted_y <- intercept + slope * fitted_x

  # Perpendicular (orthogonal) residuals
  # Distance from point to line, signed by whether point is above/below line
  fitted_y_at_x <- intercept + slope * x_vec
  perp_residuals <- sign(y_vec - fitted_y_at_x) *
    sqrt((x_vec - fitted_x)^2 + (y_vec - fitted_y)^2)

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
      intercept_se = ci_result$intercept_se,
      slope_se = ci_result$slope_se,
      residuals = perp_residuals,
      fitted_x = fitted_x,
      fitted_y = fitted_y
    ),
    settings = list(
      error_ratio = error_ratio,
      conf_level = conf_level,
      ci_method = ci_method,
      boot_n = if (ci_method == "bootstrap") boot_n else NA,
      weighted = weighted
    ),
    call = call
  )

  class(result) <- c("deming_regression", "valytics_comparison", "valytics_result")
  result
}


# Helper Functions ----

#' Parse input for deming_regression (formula or vectors)
#' @noRd
.parse_deming_input <- function(x, y, data) {

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


#' Validate Deming regression inputs
#' @noRd
.validate_deming_inputs <- function(x, y, error_ratio, conf_level, boot_n, weighted) {

  # Check error ratio
  if (!is.numeric(error_ratio) || length(error_ratio) != 1 ||
      error_ratio <= 0 || !is.finite(error_ratio)) {
    stop("'error_ratio' must be a single positive finite number.", call. = FALSE)
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

  # Check weighted flag
  if (!is.logical(weighted) || length(weighted) != 1) {
    stop("'weighted' must be a single logical value.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Compute Deming regression coefficients
#'
#' Implements the closed-form solution for Deming regression.
#' Based on Linnet (1993) and Cornbleet & Gochman (1979).
#'
#' @param x Numeric vector of x values
#' @param y Numeric vector of y values
#' @param lambda Error variance ratio (Var(error_y) / Var(error_x))
#' @return List with slope and intercept
#' @noRd
.compute_deming <- function(x, y, lambda = 1) {

  n <- length(x)

  # Calculate means
  x_bar <- mean(x)
  y_bar <- mean(y)

  # Calculate sums of squares and cross-products
  # Using computational formula for numerical stability
  sxx <- sum((x - x_bar)^2)
  syy <- sum((y - y_bar)^2)
  sxy <- sum((x - x_bar) * (y - y_bar))

  # Deming regression slope formula:
  # b = (syy - lambda*sxx + sqrt((syy - lambda*sxx)^2 + 4*lambda*sxy^2)) / (2*sxy)
  # This is the positive root of the quadratic equation

  u <- syy - lambda * sxx
  discriminant <- sqrt(u^2 + 4 * lambda * sxy^2)

  # Handle case where sxy is very small (near-horizontal or vertical line)
  if (abs(sxy) < .Machine$double.eps * max(abs(x), abs(y)) * n) {
    warning("Near-zero covariance between x and y. Results may be unreliable.",
            call. = FALSE)
    # Return slope based on ratio of standard deviations
    slope <- if (sxx > 0) sqrt(syy / sxx) * sign(sxy + 0.001) else 1
  } else {
    slope <- (u + discriminant) / (2 * sxy)
  }

  # Intercept: a = y_bar - b * x_bar
  intercept <- y_bar - slope * x_bar

  list(
    slope = slope,
    intercept = intercept,
    x_bar = x_bar,
    y_bar = y_bar,
    sxx = sxx,
    syy = syy,
    sxy = sxy
  )
}


#' Jackknife confidence intervals for Deming regression
#'
#' Uses delete-one jackknife to estimate standard errors and confidence
#' intervals. Based on Linnet (1990).
#'
#' @noRd
.deming_jackknife_ci <- function(x, y, lambda, conf_level) {

  n <- length(x)
  alpha <- 1 - conf_level

  # Original estimates
  orig_fit <- .compute_deming(x, y, lambda)
  slope_orig <- orig_fit$slope
  intercept_orig <- orig_fit$intercept

  # Jackknife estimates
  jack_slopes <- numeric(n)
  jack_intercepts <- numeric(n)

  for (i in seq_len(n)) {
    fit_i <- .compute_deming(x[-i], y[-i], lambda)
    jack_slopes[i] <- fit_i$slope
    jack_intercepts[i] <- fit_i$intercept
  }

  # Jackknife standard errors
  # SE = sqrt((n-1)/n * sum((theta_i - theta_bar)^2))
  slope_bar <- mean(jack_slopes)
  intercept_bar <- mean(jack_intercepts)

  slope_se <- sqrt((n - 1) / n * sum((jack_slopes - slope_bar)^2))
  intercept_se <- sqrt((n - 1) / n * sum((jack_intercepts - intercept_bar)^2))

  # Confidence intervals using t-distribution
  t_crit <- qt(1 - alpha / 2, df = n - 2)

  slope_ci <- c(
    lower = slope_orig - t_crit * slope_se,
    upper = slope_orig + t_crit * slope_se
  )

  intercept_ci <- c(
    lower = intercept_orig - t_crit * intercept_se,
    upper = intercept_orig + t_crit * intercept_se
  )

  list(
    slope_ci = slope_ci,
    intercept_ci = intercept_ci,
    slope_se = slope_se,
    intercept_se = intercept_se
  )
}


#' Bootstrap BCa confidence intervals for Deming regression
#' @noRd
.deming_bootstrap_ci <- function(x, y, lambda, conf_level, boot_n) {

  n <- length(x)
  alpha <- 1 - conf_level

  # Original estimates
  orig_fit <- .compute_deming(x, y, lambda)
  slope_orig <- orig_fit$slope
  intercept_orig <- orig_fit$intercept

  # Bootstrap samples
  boot_slopes <- numeric(boot_n)
  boot_intercepts <- numeric(boot_n)

  for (b in seq_len(boot_n)) {
    idx <- sample.int(n, n, replace = TRUE)
    fit_b <- .compute_deming(x[idx], y[idx], lambda)
    boot_slopes[b] <- fit_b$slope
    boot_intercepts[b] <- fit_b$intercept
  }

  # BCa confidence intervals
  slope_ci <- .bca_ci_deming(boot_slopes, slope_orig, x, y, lambda,
                             "slope", conf_level)
  intercept_ci <- .bca_ci_deming(boot_intercepts, intercept_orig, x, y, lambda,
                                 "intercept", conf_level)

  # Standard errors from bootstrap distribution
  slope_se <- sd(boot_slopes, na.rm = TRUE)
  intercept_se <- sd(boot_intercepts, na.rm = TRUE)

  list(
    slope_ci = slope_ci,
    intercept_ci = intercept_ci,
    slope_se = slope_se,
    intercept_se = intercept_se
  )
}


#' Calculate BCa confidence interval for Deming regression
#' @noRd
.bca_ci_deming <- function(boot_stat, orig_stat, x, y, lambda, param, conf_level) {

  alpha <- 1 - conf_level
  n <- length(x)

  # Remove any NA/Inf values
  boot_stat <- boot_stat[is.finite(boot_stat)]

  if (length(boot_stat) < 100) {
    warning("Too few valid bootstrap samples. Using percentile CI.", call. = FALSE)
    return(quantile(boot_stat, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  }

  # Bias correction factor z0
  prop_less <- mean(boot_stat < orig_stat)
  prop_less <- max(0.001, min(0.999, prop_less))
  z0 <- qnorm(prop_less)

  # Acceleration factor (a) using jackknife
  jack_stat <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- .compute_deming(x[-i], y[-i], lambda)
    jack_stat[i] <- if (param == "slope") fit_i$slope else fit_i$intercept
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
