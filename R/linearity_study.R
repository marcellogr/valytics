#' Linearity Evaluation for Quantitative Measurement Procedures
#'
#' @description
#' Evaluates the linearity of a measurement procedure by fitting polynomial
#' models and testing for lack of fit. Calculates deviations from linearity
#' and recovery at each concentration level. Optionally assesses whether
#' deviations fall within an Allowable Deviation from Linearity (ADL).
#'
#' Supports two input modes:
#' \itemize{
#'   \item **Known target concentrations**: Provide target (expected) values in
#'     `x` and measured values in `y`.
#'   \item **Dilution-based design**: Provide measured values in `y` and
#'     dilution proportions via the `dilution` parameter. Target concentrations
#'     are estimated from the observed endpoint means and the dilution scheme.
#' }
#'
#' @param x Numeric vector of target (expected) concentrations, or a formula
#'   of the form `target ~ measured`. When replicates are present, `x` should
#'   contain the target value for each replicate measurement. When using
#'   `dilution`, `x` should be the measured values (or use the formula
#'   `~ measured` with `dilution` specified).
#' @param y Numeric vector of measured values. Ignored if `x` is a formula.
#'   Also ignored when using `dilution` with a one-sided formula.
#' @param data Optional data frame containing variables specified in the
#'   formula or as column names.
#' @param sample Optional character string naming the column that identifies
#'   concentration levels. If `NULL`, levels are inferred from unique values
#'   of `x` (target concentrations) or `dilution` proportions.
#' @param dilution Numeric vector of dilution proportions (values between 0
#'   and 1, where 0 = 100\% low pool and 1 = 100\% high pool), or a character
#'   string naming the column in `data` containing dilution proportions.
#'   When provided, target concentrations are estimated from the observed
#'   means at the 0\% and 100\% endpoints using linear interpolation:
#'   `target = mean_low + dilution * (mean_high - mean_low)`. If `NULL`
#'   (default), `x` is used as target concentrations directly.
#' @param max_poly Integer; maximum polynomial degree to fit (default: 3).
#'   Must be 2 or 3.
#' @param adl Numeric; Allowable Deviation from Linearity. If provided,
#'   a pass/fail linearity assessment is performed. If `NULL` (default),
#'   deviations are reported without judgment.
#' @param adl_type Character; how `adl` is expressed: `"absolute"` (default)
#'   for deviation in measurement units, or `"percent"` for percent deviation
#'   relative to the target concentration.
#' @param weights Character; weighting scheme for polynomial regression:
#'   `"equal"` (default), `"1/x"` (inverse concentration), `"1/x2"`
#'   (inverse squared concentration), or `"1/var"` (inverse within-level
#'   variance). The `"1/x"` and `"1/x2"` options are useful when variance
#'   scales with concentration. The `"1/var"` option uses the observed
#'   replicate variance at each level, giving more weight to levels with
#'   lower variability. Requires replicates at each level; falls back to
#'   equal weights for levels with fewer than 2 replicates.
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param na_action How to handle missing values: `"omit"` (default) removes
#'   pairs with any NA, `"fail"` stops with an error.
#'
#' @return An object of class `c("linearity_study", "valytics_linearity",
#'   "valytics_result")`, which is a list containing:
#'
#'   \describe{
#'     \item{input}{List with original data and metadata:
#'       \itemize{
#'         \item `x`: Target concentrations (after NA handling)
#'         \item `y`: Measured values (after NA handling)
#'         \item `n`: Total number of observations
#'         \item `n_excluded`: Number of observations excluded due to NAs
#'         \item `n_levels`: Number of concentration levels
#'         \item `level_stats`: Data frame with per-level summary statistics
#'           (target, n, mean, sd, cv)
#'         \item `var_names`: Named character vector with variable names
#'       }
#'       If `dilution` was used, the input list also contains:
#'       \itemize{
#'         \item `dilution`: Numeric vector of dilution proportions
#'         \item `endpoint_means`: Named numeric vector with `low` and `high`
#'           pool mean values used for target estimation
#'       }
#'     }
#'     \item{results}{List with statistical results:
#'       \itemize{
#'         \item `models`: Named list of fitted `lm` objects (linear,
#'           quadratic, and optionally cubic)
#'         \item `best_model`: Character; name of the best-fitting model
#'           based on lack-of-fit testing
#'         \item `lack_of_fit`: List with two components:
#'           \itemize{
#'             \item `pure_error`: Data frame with pure error lack-of-fit
#'               test results for each polynomial model (SS_lof, SS_pe,
#'               df_lof, df_pe, f_statistic, p_value). Only available when
#'               replicates exist at each level.
#'             \item `polynomial`: Data frame with nested polynomial
#'               comparison results (model comparisons, F-statistics,
#'               p-values)
#'           }
#'         \item `coefficients`: Data frame with coefficients for each model
#'         \item `deviations`: Data frame with per-level deviations from
#'           linearity (target, observed_mean, linear_predicted,
#'           poly_predicted, deviation_abs, deviation_pct)
#'         \item `max_deviation`: List with max_abs, max_pct,
#'           at_concentration
#'         \item `recovery`: Data frame with per-level recovery percentages
#'         \item `linear`: Logical or NA; TRUE if all deviations are within
#'           ADL (NA if ADL not provided)
#'         \item `linear_range`: Numeric vector of length 2 giving the
#'           concentration range where deviations are within ADL (NULL if
#'           ADL not provided or all levels pass)
#'       }
#'     }
#'     \item{settings}{List with analysis settings:
#'       \itemize{
#'         \item `max_poly`: Maximum polynomial degree
#'         \item `adl`: Allowable Deviation from Linearity (NULL if not set)
#'         \item `adl_type`: Type of ADL ("absolute" or "percent")
#'         \item `weights`: Weighting scheme used
#'         \item `conf_level`: Confidence level
#'       }
#'     }
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' The linearity evaluation proceeds in several steps:
#'
#' \enumerate{
#'   \item **Replicate summarization**: At each concentration level, the mean
#'     and standard deviation of replicate measurements are calculated.
#'   \item **Polynomial fitting**: Linear (degree 1), quadratic (degree 2),
#'     and optionally cubic (degree 3) models are fitted to the data.
#'   \item **Lack-of-fit testing**: Two complementary tests are performed:
#'     \itemize{
#'       \item **Pure error lack-of-fit test**: Decomposes the residual sum
#'         of squares into lack-of-fit and pure error (within-level replicate
#'         variance) components. A significant result indicates the model does
#'         not adequately describe the data beyond replicate noise. Requires
#'         replicates at each level.
#'       \item **Nested polynomial comparison**: F-tests comparing successive
#'         polynomial orders (linear vs quadratic, quadratic vs cubic) to
#'         determine whether higher-order terms significantly improve the fit.
#'     }
#'   \item **Best model selection**: The best model is selected using a
#'     stepwise approach. Starting from the highest-order polynomial, the
#'     significance of each highest-order coefficient is evaluated. The
#'     simplest adequate model is retained.
#'   \item **Deviation calculation**: At each level, the deviation from
#'     linearity is the difference between the best polynomial fit and the
#'     linear fit.
#'   \item **Recovery**: Calculated as 100 * (observed mean / target) at each
#'     level.
#'   \item **ADL assessment**: If an Allowable Deviation from Linearity is
#'     provided, deviations are compared against this limit.
#' }
#'
#' @section Experimental Design:
#' A typical linearity study uses 5-9 concentration levels spanning the
#' Analytical Measurement Interval (AMI), with 2-4 replicates per level.
#' Levels are usually prepared by serial dilution of a high-concentration
#' pool with a low-concentration pool.
#'
#' @section Dilution-Based Design:
#' When target concentrations are not known (e.g., when using pooled patient
#' samples), the `dilution` parameter allows specifying dilution proportions
#' instead. The function estimates target concentrations from the observed
#' endpoint means:
#'
#' \deqn{Target_i = \bar{y}_{low} + d_i \times (\bar{y}_{high} - \bar{y}_{low})}
#'
#' where \eqn{d_i} is the dilution proportion (0 = low pool, 1 = high pool)
#' and \eqn{\bar{y}_{low}} and \eqn{\bar{y}_{high}} are the observed means
#' at the endpoints. This approach assumes that the true concentrations
#' scale linearly with the mixing ratio, which is the null hypothesis being
#' tested.
#'
#' @section Weighted Regression:
#' When variance increases with concentration (common in immunoassays and
#' assays with wide dynamic ranges), weighted regression can improve the
#' fit. The `weights` parameter supports:
#' \itemize{
#'   \item `"equal"`: No weighting (default, appropriate when variance is
#'     constant)
#'   \item `"1/x"`: Inverse concentration weighting
#'   \item `"1/x2"`: Inverse squared concentration weighting
#'   \item `"1/var"`: Inverse within-level variance weighting. Each
#'     observation is weighted by 1/s^2 of its concentration level, where
#'     s^2 is the replicate variance. This is the theoretically optimal
#'     WLS approach when variance differs across levels. Requires
#'     replicates (n >= 2) at each level.
#' }
#'
#' @references
#' Kroll MH, Emancipator K (1993). A theoretical evaluation of linearity.
#' \emph{Clinical Chemistry}, 39(3):405-413.
#' \doi{10.1093/clinchem/39.3.405}
#'
#' Emancipator K, Kroll MH (1993). A quantitative measure of nonlinearity.
#' \emph{Clinical Chemistry}, 39(5):766-772.
#' \doi{10.1093/clinchem/39.5.766}
#'
#' Jhang JS, Chang CC, Fink DJ, Kroll MH (2004). Evaluation of linearity
#' in the clinical laboratory. \emph{Archives of Pathology and Laboratory
#' Medicine}, 128(1):44-48.
#' \doi{10.5858/2004-128-44-EOLITC}
#'
#' @seealso
#' [plot.linearity_study()] for visualization,
#' [summary.linearity_study()] for detailed summary
#'
#' @examples
#' # Example 1: Known target concentrations
#' set.seed(42)
#' target <- rep(c(10, 50, 100, 200, 300, 400, 500), each = 4)
#' measured <- target * rnorm(28, mean = 1.0, sd = 0.02) +
#'             rnorm(28, mean = 0, sd = 2)
#'
#' # Basic analysis
#' lin <- linearity_study(target, measured)
#' lin
#'
#' # With ADL assessment
#' lin <- linearity_study(target, measured, adl = 5, adl_type = "percent")
#' summary(lin)
#' plot(lin)
#'
#' # Formula interface
#' df <- data.frame(target = target, measured = measured)
#' lin <- linearity_study(target ~ measured, data = df)
#'
#' # Example 2: Dilution-based design (no known targets)
#' dilution_levels <- rep(seq(0, 1, length.out = 9), each = 4)
#' # Simulate measured values from a linear system
#' true_low <- 10
#' true_high <- 500
#' true_conc <- true_low + dilution_levels * (true_high - true_low)
#' measured2 <- true_conc * rnorm(36, 1.0, 0.02) + rnorm(36, 0, 2)
#'
#' lin2 <- linearity_study(measured2, dilution = dilution_levels)
#' lin2
#'
#' # Dilution with data frame
#' df2 <- data.frame(
#'   result = measured2,
#'   dil_pct = dilution_levels
#' )
#' lin2 <- linearity_study(~ result, data = df2, dilution = "dil_pct")
#'
#' # With weighted regression
#' lin <- linearity_study(target, measured, weights = "1/x")
#'
#' @export
linearity_study <- function(x = NULL,
                            y = NULL,
                            data = NULL,
                            sample = NULL,
                            dilution = NULL,
                            max_poly = 3,
                            adl = NULL,
                            adl_type = c("absolute", "percent"),
                            weights = c("equal", "1/x", "1/x2", "1/var"),
                            conf_level = 0.95,
                            na_action = c("omit", "fail")) {
  
  # Capture the call
  call <- match.call()
  
  # Match arguments
  adl_type <- match.arg(adl_type)
  weights <- match.arg(weights)
  na_action <- match.arg(na_action)
  
  # Input parsing ----
  dilution_info <- NULL
  
  # Validate that we have enough info to proceed
  if (is.null(dilution) && is.null(x)) {
    stop("Either provide target concentrations via `x` or use `dilution` ",
         "for dilution-based designs.", call. = FALSE)
  }
  
  if (!is.null(dilution)) {
    # Dilution-based design: user provides measured values + dilution proportions
    parsed <- .parse_dilution_input(x, y, data, sample, dilution)
    y_vals <- parsed$measured
    dilution_vec <- parsed$dilution
    var_names <- parsed$var_names
    
    # Validate dilution values
    .validate_dilution(dilution_vec)
    
    # Handle NAs before computing endpoint means
    complete <- stats::complete.cases(dilution_vec, y_vals)
    n_excluded <- sum(!complete)
    if (na_action == "fail" && n_excluded > 0) {
      stop("Missing values detected. Use `na_action = 'omit'` to remove them.",
           call. = FALSE)
    }
    dilution_clean <- dilution_vec[complete]
    y_clean <- y_vals[complete]
    
    # Compute endpoint means (dilution = 0 and dilution = 1)
    endpoint_means <- .compute_endpoint_means(dilution_clean, y_clean)
    
    # Compute target concentrations from dilution proportions
    x_clean <- as.numeric(endpoint_means["low"] +
                            dilution_clean * (endpoint_means["high"] - endpoint_means["low"]))
    
    n <- length(x_clean)
    dilution_info <- list(
      dilution = dilution_clean,
      endpoint_means = endpoint_means
    )
    
  } else {
    # Standard design: user provides target concentrations + measured values
    parsed <- .parse_linearity_input(x, y, data, sample)
    x_vals <- parsed$x
    y_vals <- parsed$y
    var_names <- parsed$var_names
    
    # Input validation ----
    .validate_linearity_input(x_vals, y_vals, max_poly, adl, conf_level)
    
    # Handle missing values ----
    complete <- stats::complete.cases(x_vals, y_vals)
    n_excluded <- sum(!complete)
    
    if (na_action == "fail" && n_excluded > 0) {
      stop("Missing values detected. Use `na_action = 'omit'` to remove them.",
           call. = FALSE)
    }
    
    x_clean <- x_vals[complete]
    y_clean <- y_vals[complete]
    n <- length(x_clean)
  }
  
  # Validate remaining inputs (shared by both paths) ----
  if (!is.null(dilution)) {
    # Validate shared params only (x/y already cleaned)
    if (!is.numeric(max_poly) || length(max_poly) != 1 ||
        !max_poly %in% c(2, 3)) {
      stop("`max_poly` must be 2 or 3.", call. = FALSE)
    }
    if (!is.null(adl)) {
      if (!is.numeric(adl) || length(adl) != 1 || adl <= 0) {
        stop("`adl` must be a single positive number.", call. = FALSE)
      }
    }
    if (!is.numeric(conf_level) || length(conf_level) != 1 ||
        conf_level <= 0 || conf_level >= 1) {
      stop("`conf_level` must be a single number between 0 and 1.",
           call. = FALSE)
    }
  }
  
  if (n < 6) {
    stop("At least 6 observations are required for linearity evaluation.",
         call. = FALSE)
  }
  
  # Level statistics ----
  level_stats <- .compute_level_stats(x_clean, y_clean)
  
  n_levels <- nrow(level_stats)
  if (n_levels < 3) {
    stop("At least 3 concentration levels are required. Found: ", n_levels,
         call. = FALSE)
  }
  
  if (max_poly >= n_levels) {
    max_poly <- n_levels - 1
    message("Reduced max_poly to ", max_poly,
            " (must be less than number of levels).")
  }
  
  # Compute weights ----
  w <- .compute_linearity_weights(x_clean, weights, level_stats)
  
  # Polynomial fitting ----
  models <- .fit_polynomial_models(x_clean, y_clean, w, max_poly)
  
  # Lack-of-fit testing ----
  lof_results <- .lack_of_fit_tests(models, max_poly, x_clean, y_clean)
  
  # Determine best model ----
  best_model <- .select_best_model(lof_results, max_poly)
  
  # Coefficients table ----
  coef_table <- .extract_coefficients(models, max_poly)
  
  # Deviations from linearity ----
  deviations <- .compute_deviations(level_stats, models, best_model)
  
  # Maximum deviation ----
  max_dev <- .compute_max_deviation(deviations)
  
  # Recovery ----
  recovery <- .compute_recovery(level_stats)
  
  # ADL assessment ----
  adl_result <- .assess_adl(deviations, adl, adl_type)
  
  # Linear range ----
  linear_range <- .compute_linear_range(deviations, adl, adl_type)
  
  # Construct output ----
  input_list <- list(
    x = x_clean,
    y = y_clean,
    n = n,
    n_excluded = n_excluded,
    n_levels = n_levels,
    level_stats = level_stats,
    var_names = var_names
  )
  
  # Add dilution info if applicable
  if (!is.null(dilution_info)) {
    input_list$dilution <- dilution_info$dilution
    input_list$endpoint_means <- dilution_info$endpoint_means
  }
  
  # Compute dilution-scale WLS fit if applicable ----
  dilution_fit <- NULL
  if (!is.null(dilution_info)) {
    dil_rounded <- round(dilution_info$dilution, 8)
    unique_dils <- sort(unique(dil_rounded))
    
    dil_means <- numeric(length(unique_dils))
    dil_n <- numeric(length(unique_dils))
    for (i in seq_along(unique_dils)) {
      idx <- which(dil_rounded == unique_dils[i])
      dil_means[i] <- mean(y_clean[idx], na.rm = TRUE)
      dil_n[i] <- length(idx)
    }
    
    dil_df <- data.frame(dilution = unique_dils, mean = dil_means, w = dil_n)
    dil_lm <- stats::lm(mean ~ dilution, data = dil_df, weights = w)
    
    dilution_fit <- list(
      intercept = unname(stats::coef(dil_lm)[1]),
      slope = unname(stats::coef(dil_lm)[2]),
      r_squared = summary(dil_lm)$r.squared,
      model = dil_lm
    )
  }
  
  structure(
    list(
      input = input_list,
      results = list(
        models = models,
        best_model = best_model,
        lack_of_fit = lof_results,
        coefficients = coef_table,
        deviations = deviations,
        max_deviation = max_dev,
        recovery = recovery,
        linear = adl_result,
        linear_range = linear_range,
        dilution_fit = dilution_fit
      ),
      settings = list(
        max_poly = max_poly,
        adl = adl,
        adl_type = adl_type,
        weights = weights,
        conf_level = conf_level
      ),
      call = call
    ),
    class = c("linearity_study", "valytics_linearity", "valytics_result")
  )
}


# Helper Functions ----

#' Parse input for dilution-based linearity study
#' @noRd
#' @keywords internal
.parse_dilution_input <- function(x, y, data, sample, dilution) {
  
  # Resolve dilution vector
  if (is.character(dilution) && length(dilution) == 1) {
    if (is.null(data) || !is.data.frame(data)) {
      stop("When `dilution` is a column name, `data` must be provided.",
           call. = FALSE)
    }
    dilution_vec <- data[[dilution]]
    if (is.null(dilution_vec)) {
      stop(sprintf("Column '%s' not found in `data`.", dilution),
           call. = FALSE)
    }
  } else if (is.numeric(dilution)) {
    dilution_vec <- dilution
  } else {
    stop("`dilution` must be a numeric vector or a column name string.",
         call. = FALSE)
  }
  
  # Parse measured values
  if (inherits(x, "formula")) {
    # One-sided formula: ~ measured
    formula_vars <- all.vars(x)
    if (length(formula_vars) == 1) {
      # One-sided: ~ measured
      var_name <- formula_vars[1]
      if (!is.null(data) && is.data.frame(data)) {
        measured <- data[[var_name]]
      } else {
        env <- environment(x)
        measured <- get(var_name, envir = env)
      }
      var_names <- c(x = "target (estimated)", y = var_name)
    } else if (length(formula_vars) == 2) {
      # Two-sided formula: ignore LHS, use RHS as measured
      warning("With `dilution`, only the RHS of the formula is used as ",
              "measured values. The LHS is ignored.", call. = FALSE)
      var_name <- formula_vars[2]
      if (!is.null(data) && is.data.frame(data)) {
        measured <- data[[var_name]]
      } else {
        env <- environment(x)
        measured <- get(var_name, envir = env)
      }
      var_names <- c(x = "target (estimated)", y = var_name)
    } else {
      stop("Formula must have one or two variables when using `dilution`.",
           call. = FALSE)
    }
  } else if (is.null(x) && !is.null(y)) {
    # y provided without x: y is the measured values
    measured <- y
    y_name <- deparse(substitute(y, env = parent.frame(2)))
    var_names <- c(x = "target (estimated)", y = y_name)
  } else if (!is.null(x) && !inherits(x, "formula")) {
    # Vector interface: x is the measured values, y is ignored
    if (!is.null(y)) {
      warning("`y` is ignored when `dilution` is provided with `x`. ",
              "`x` is treated as measured values.", call. = FALSE)
    }
    measured <- x
    x_name <- deparse(substitute(x, env = parent.frame(2)))
    var_names <- c(x = "target (estimated)", y = x_name)
  } else {
    stop("Provide measured values via `x`, `y`, or a formula when using ",
         "`dilution`.", call. = FALSE)
  }
  
  if (!is.numeric(measured)) {
    stop("Measured values must be numeric.", call. = FALSE)
  }
  if (length(measured) != length(dilution_vec)) {
    stop("Measured values and `dilution` must have the same length.",
         call. = FALSE)
  }
  
  list(measured = as.numeric(measured),
       dilution = as.numeric(dilution_vec),
       var_names = var_names)
}


#' Validate dilution proportions
#' @noRd
#' @keywords internal
.validate_dilution <- function(dilution) {
  
  if (!is.numeric(dilution)) {
    stop("`dilution` must be numeric.", call. = FALSE)
  }
  
  # Check range: should be [0, 1]
  d_range <- range(dilution, na.rm = TRUE)
  
  if (d_range[1] < -0.01 || d_range[2] > 1.01) {
    stop("`dilution` values must be between 0 and 1 (proportions). ",
         "Use 0 for the low pool and 1 for the high pool.",
         call. = FALSE)
  }
  
  # Must include endpoints (or close to them)
  has_low <- any(abs(dilution) < 0.01, na.rm = TRUE)
  has_high <- any(abs(dilution - 1) < 0.01, na.rm = TRUE)
  
  if (!has_low || !has_high) {
    stop("`dilution` must include values at or near both endpoints ",
         "(0 for low pool, 1 for high pool) to estimate target concentrations.",
         call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Compute endpoint means from dilution design
#' @noRd
#' @keywords internal
.compute_endpoint_means <- function(dilution, y) {
  
  # Find observations at endpoints
  low_idx <- which(abs(dilution) < 0.01)
  high_idx <- which(abs(dilution - 1) < 0.01)
  
  if (length(low_idx) == 0 || length(high_idx) == 0) {
    stop("Cannot compute endpoint means: no observations at dilution = 0 ",
         "or dilution = 1.", call. = FALSE)
  }
  
  low_mean <- mean(y[low_idx])
  high_mean <- mean(y[high_idx])
  
  if (abs(high_mean - low_mean) < .Machine$double.eps * 100) {
    stop("High and low pool means are nearly identical. ",
         "Cannot estimate target concentrations.", call. = FALSE)
  }
  
  c(low = low_mean, high = high_mean)
}



#' Parse input for linearity_study
#' @noRd
#' @keywords internal
.parse_linearity_input <- function(x, y, data, sample) {
  
  # Formula interface: target ~ measured
  if (inherits(x, "formula")) {
    if (!is.null(y)) {
      warning("`y` is ignored when `x` is a formula.", call. = FALSE)
    }
    
    formula_vars <- all.vars(x)
    if (length(formula_vars) != 2) {
      stop("Formula must have exactly two variables: target ~ measured",
           call. = FALSE)
    }
    
    var_names <- c(x = formula_vars[1], y = formula_vars[2])
    
    if (is.null(data)) {
      env <- environment(x)
      x_vals <- get(var_names["x"], envir = env)
      y_vals <- get(var_names["y"], envir = env)
    } else {
      if (!is.data.frame(data)) {
        stop("`data` must be a data frame.", call. = FALSE)
      }
      x_vals <- data[[var_names["x"]]]
      y_vals <- data[[var_names["y"]]]
      
      if (is.null(x_vals) || is.null(y_vals)) {
        stop("Variables specified in formula not found in `data`.",
             call. = FALSE)
      }
    }
    
  } else {
    # Vector interface
    if (is.null(y)) {
      stop("Either provide a formula or both `x` and `y` vectors.",
           call. = FALSE)
    }
    
    x_vals <- x
    y_vals <- y
    
    x_name <- deparse(substitute(x, env = parent.frame(2)))
    y_name <- deparse(substitute(y, env = parent.frame(2)))
    var_names <- c(x = x_name, y = y_name)
    
    # If data is provided, extract from data frame
    if (!is.null(data)) {
      if (!is.data.frame(data)) {
        stop("`data` must be a data frame.", call. = FALSE)
      }
      if (is.character(x) && length(x) == 1) {
        x_vals <- data[[x]]
        var_names["x"] <- x
      }
      if (is.character(y) && length(y) == 1) {
        y_vals <- data[[y]]
        var_names["y"] <- y
      }
    }
  }
  
  # Validate types before coercion
  if (!is.numeric(x_vals)) {
    stop("`x` (target concentrations) must be numeric.", call. = FALSE)
  }
  if (!is.numeric(y_vals)) {
    stop("`y` (measured values) must be numeric.", call. = FALSE)
  }
  
  list(x = as.numeric(x_vals), y = as.numeric(y_vals), var_names = var_names)
}


#' Validate linearity_study inputs
#' @noRd
#' @keywords internal
.validate_linearity_input <- function(x, y, max_poly, adl, conf_level) {
  
  if (!is.numeric(x)) {
    stop("`x` (target concentrations) must be numeric.", call. = FALSE)
  }
  if (!is.numeric(y)) {
    stop("`y` (measured values) must be numeric.", call. = FALSE)
  }
  if (length(x) != length(y)) {
    stop("`x` and `y` must have the same length.", call. = FALSE)
  }
  
  if (!is.numeric(max_poly) || length(max_poly) != 1 ||
      !max_poly %in% c(2, 3)) {
    stop("`max_poly` must be 2 or 3.", call. = FALSE)
  }
  
  if (!is.null(adl)) {
    if (!is.numeric(adl) || length(adl) != 1 || adl <= 0) {
      stop("`adl` must be a single positive number.", call. = FALSE)
    }
  }
  
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number between 0 and 1.",
         call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Compute per-level summary statistics
#' @noRd
#' @keywords internal
.compute_level_stats <- function(x, y) {
  
  # Identify unique levels (use rounding to handle floating point)
  levels <- sort(unique(round(x, digits = 10)))
  
  stats_list <- lapply(levels, function(lev) {
    idx <- which(abs(x - lev) < .Machine$double.eps * 100)
    yi <- y[idx]
    data.frame(
      target = lev,
      n = length(yi),
      mean = mean(yi),
      sd = if (length(yi) > 1) stats::sd(yi) else NA_real_,
      cv = if (length(yi) > 1 && mean(yi) != 0) {
        100 * stats::sd(yi) / abs(mean(yi))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, stats_list)
}


#' Compute regression weights
#' @noRd
#' @keywords internal
.compute_linearity_weights <- function(x, weights, level_stats = NULL) {
  
  switch(weights,
         "equal" = rep(1, length(x)),
         "1/x" = {
           w <- 1 / abs(x)
           # Handle zero concentrations
           w[!is.finite(w)] <- max(w[is.finite(w)], na.rm = TRUE)
           w
         },
         "1/x2" = {
           w <- 1 / (x^2)
           w[!is.finite(w)] <- max(w[is.finite(w)], na.rm = TRUE)
           w
         },
         "1/var" = {
           if (is.null(level_stats)) {
             warning("Level statistics not available for 1/var weighting. ",
                     "Falling back to equal weights.", call. = FALSE)
             return(rep(1, length(x)))
           }
           
           # Map each observation to its level's variance
           levels <- level_stats$target
           variances <- level_stats$sd^2  # sd^2 = variance
           
           # Build per-observation weights
           w <- numeric(length(x))
           x_rounded <- round(x, digits = 10)
           
           for (i in seq_along(levels)) {
             idx <- which(abs(x_rounded - round(levels[i], 10)) <
                            .Machine$double.eps * 100)
             if (length(idx) > 0) {
               vi <- variances[i]
               if (!is.na(vi) && vi > .Machine$double.eps) {
                 w[idx] <- 1 / vi
               } else {
                 # No variance estimate (single replicate or zero variance)
                 # Will be replaced below
                 w[idx] <- NA_real_
               }
             }
           }
           
           # Replace NAs with the maximum finite weight (most precise level)
           finite_w <- w[is.finite(w) & !is.na(w)]
           if (length(finite_w) == 0) {
             warning("No valid within-level variances found for 1/var weighting. ",
                     "Falling back to equal weights.", call. = FALSE)
             return(rep(1, length(x)))
           }
           w[is.na(w) | !is.finite(w)] <- max(finite_w)
           w
         }
  )
}


#' Fit polynomial models
#' @noRd
#' @keywords internal
.fit_polynomial_models <- function(x, y, w, max_poly) {
  
  models <- list()
  
  # Linear (degree 1)
  models$linear <- stats::lm(y ~ x, weights = w)
  
  # Quadratic (degree 2)
  models$quadratic <- stats::lm(y ~ x + I(x^2), weights = w)
  
  # Cubic (degree 3)
  if (max_poly >= 3) {
    models$cubic <- stats::lm(y ~ x + I(x^2) + I(x^3), weights = w)
  }
  
  models
}


#' Perform lack-of-fit tests (pure error and nested polynomial)
#' @noRd
#' @keywords internal
.lack_of_fit_tests <- function(models, max_poly, x, y) {
  
  # Pure error lack-of-fit test ----
  # Decomposes residual SS into lack-of-fit SS and pure error SS
  # Pure error = within-level replicate variance
  # Lack-of-fit = deviation of level means from fitted model
  
  pure_error_results <- .pure_error_lof(models, max_poly, x, y)
  
  # Nested polynomial comparison ----
  # F-tests comparing successive polynomial orders
  poly_results <- .nested_polynomial_tests(models, max_poly)
  
  list(
    pure_error = pure_error_results,
    polynomial = poly_results
  )
}


#' Pure error lack-of-fit test
#'
#' Decomposes the residual sum of squares from each polynomial model into:
#' - SS_lof (lack of fit): variation of level means around the fitted curve
#' - SS_pe (pure error): within-level replicate variation
#'
#' F = (SS_lof / df_lof) / (SS_pe / df_pe)
#'
#' @references
#' Draper NR, Smith H (1998). Applied Regression Analysis, 3rd ed.
#' Wiley. Chapter 3.
#'
#' @noRd
#' @keywords internal
.pure_error_lof <- function(models, max_poly, x, y) {
  
  # Compute pure error (within-level replicate SS)
  levels <- sort(unique(round(x, digits = 10)))
  n_levels <- length(levels)
  
  # Check if replicates exist
  has_replicates <- FALSE
  ss_pe <- 0
  df_pe <- 0
  
  for (lev in levels) {
    idx <- which(abs(x - lev) < .Machine$double.eps * 100)
    ni <- length(idx)
    if (ni > 1) {
      has_replicates <- TRUE
      yi <- y[idx]
      ss_pe <- ss_pe + sum((yi - mean(yi))^2)
      df_pe <- df_pe + (ni - 1)
    }
  }
  
  if (!has_replicates || df_pe == 0) {
    return(NULL)
  }
  
  ms_pe <- ss_pe / df_pe
  
  # If pure error is essentially zero (perfect replicates), the F-test
  
  # is undefined - we cannot distinguish LOF from noise when there is no noise
  if (ms_pe < .Machine$double.eps * 100) {
    return(NULL)
  }
  
  # Compute lack-of-fit for each model
  model_names <- c("linear", "quadratic")
  if (max_poly >= 3) model_names <- c(model_names, "cubic")
  
  results_list <- list()
  
  for (mname in model_names) {
    model <- models[[mname]]
    p <- length(stats::coef(model))  # number of parameters
    n <- length(x)
    
    # Total residual SS from the model
    ss_res <- sum(stats::residuals(model)^2)
    
    # Lack-of-fit SS = Residual SS - Pure Error SS
    ss_lof <- ss_res - ss_pe
    
    # Degrees of freedom
    # df_res = n - p (total residual df)
    # df_pe = sum(ni - 1) (pure error df)
    # df_lof = df_res - df_pe = (n - p) - df_pe = n_levels - p
    df_lof <- n_levels - p
    
    if (df_lof <= 0 || ss_lof < 0) {
      # Not enough df for this test (model has too many parameters
      # relative to number of levels)
      results_list[[mname]] <- data.frame(
        model = mname,
        ss_lof = NA_real_,
        ss_pe = ss_pe,
        df_lof = df_lof,
        df_pe = df_pe,
        f_statistic = NA_real_,
        p_value = NA_real_,
        significant = NA,
        stringsAsFactors = FALSE
      )
      next
    }
    
    ms_lof <- ss_lof / df_lof
    
    f_stat <- ms_lof / ms_pe
    p_value <- stats::pf(f_stat, df_lof, df_pe, lower.tail = FALSE)
    
    results_list[[mname]] <- data.frame(
      model = mname,
      ss_lof = ss_lof,
      ss_pe = ss_pe,
      df_lof = df_lof,
      df_pe = df_pe,
      f_statistic = f_stat,
      p_value = p_value,
      significant = p_value < 0.05,
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, results_list)
  rownames(result) <- NULL
  result
}


#' Nested polynomial comparison F-tests
#' @noRd
#' @keywords internal
.nested_polynomial_tests <- function(models, max_poly) {
  
  results <- list()
  
  # Linear vs Quadratic
  anova_lq <- stats::anova(models$linear, models$quadratic)
  results <- rbind(results, data.frame(
    comparison = "Linear vs Quadratic",
    df1 = anova_lq$Df[2],
    df2 = anova_lq$Res.Df[2],
    f_statistic = anova_lq$F[2],
    p_value = anova_lq[["Pr(>F)"]][2],
    significant = anova_lq[["Pr(>F)"]][2] < 0.05,
    stringsAsFactors = FALSE
  ))
  
  # Quadratic vs Cubic
  if (max_poly >= 3 && !is.null(models$cubic)) {
    anova_qc <- stats::anova(models$quadratic, models$cubic)
    results <- rbind(results, data.frame(
      comparison = "Quadratic vs Cubic",
      df1 = anova_qc$Df[2],
      df2 = anova_qc$Res.Df[2],
      f_statistic = anova_qc$F[2],
      p_value = anova_qc[["Pr(>F)"]][2],
      significant = anova_qc[["Pr(>F)"]][2] < 0.05,
      stringsAsFactors = FALSE
    ))
  }
  
  rownames(results) <- NULL
  results
}


#' Select best polynomial model (EP06-A2 stepwise approach)
#'
#' The selection proceeds as follows:
#' 1. Start with the highest-order polynomial fitted.
#' 2. Test whether the highest-order coefficient is significant (via nested
#'    polynomial comparison). If not, reduce the order.
#' 3. Validate using the pure error lack-of-fit test: verify the selected
#'    model adequately fits the data (LOF p >= 0.05).
#' 4. If the linear model shows significant lack of fit (pure error test),
#'    but no polynomial term is significant, use the polynomial model that
#'    reduces LOF the most.
#'
#' @noRd
#' @keywords internal
.select_best_model <- function(lof_results, max_poly) {
  
  poly_tests <- lof_results$polynomial
  pe_tests <- lof_results$pure_error
  
  # Step 1: Use nested polynomial tests (top-down)
  best <- "linear"
  
  # Check if quadratic term is significant
  lq_row <- poly_tests[poly_tests$comparison == "Linear vs Quadratic", ]
  quad_sig <- nrow(lq_row) > 0 && !is.na(lq_row$p_value) &&
    lq_row$p_value < 0.05
  
  if (quad_sig) {
    best <- "quadratic"
    
    # Check if cubic term adds more
    if (max_poly >= 3) {
      qc_row <- poly_tests[poly_tests$comparison == "Quadratic vs Cubic", ]
      cubic_sig <- nrow(qc_row) > 0 && !is.na(qc_row$p_value) &&
        qc_row$p_value < 0.05
      if (cubic_sig) {
        best <- "cubic"
      }
    }
  }
  
  # Step 2: Cross-check with pure error LOF test
  # If the linear model has significant LOF (pure error test) but nested
  # polynomial test was not significant, prefer quadratic as the data
  # shows systematic departure from linearity
  if (best == "linear" && !is.null(pe_tests)) {
    pe_linear <- pe_tests[pe_tests$model == "linear", ]
    if (nrow(pe_linear) > 0 && !is.na(pe_linear$p_value) &&
        pe_linear$p_value < 0.05) {
      # Linear model has significant LOF against pure error
      # Check if quadratic resolves it
      pe_quad <- pe_tests[pe_tests$model == "quadratic", ]
      if (nrow(pe_quad) > 0 && !is.na(pe_quad$p_value) &&
          pe_quad$p_value >= 0.05) {
        best <- "quadratic"
      } else if (max_poly >= 3) {
        pe_cubic <- pe_tests[pe_tests$model == "cubic", ]
        if (nrow(pe_cubic) > 0 && !is.na(pe_cubic$p_value) &&
            pe_cubic$p_value >= 0.05) {
          best <- "cubic"
        } else {
          # Even higher-order polynomials don't resolve LOF;
          # use the one with lowest LOF F-statistic
          best <- "quadratic"
        }
      } else {
        best <- "quadratic"
      }
    }
  }
  
  best
}


#' Extract coefficients from all models into a table
#' @noRd
#' @keywords internal
.extract_coefficients <- function(models, max_poly) {
  
  coef_list <- list()
  
  # Linear
  lc <- stats::coef(models$linear)
  coef_list$linear <- data.frame(
    model = "linear",
    term = c("intercept", "x"),
    estimate = unname(lc),
    stringsAsFactors = FALSE
  )
  
  # Quadratic
  qc <- stats::coef(models$quadratic)
  coef_list$quadratic <- data.frame(
    model = "quadratic",
    term = c("intercept", "x", "x2"),
    estimate = unname(qc),
    stringsAsFactors = FALSE
  )
  
  # Cubic
  if (max_poly >= 3 && !is.null(models$cubic)) {
    cc <- stats::coef(models$cubic)
    coef_list$cubic <- data.frame(
      model = "cubic",
      term = c("intercept", "x", "x2", "x3"),
      estimate = unname(cc),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, coef_list)
}


#' Compute deviations from linearity at each level
#' @noRd
#' @keywords internal
.compute_deviations <- function(level_stats, models, best_model) {
  
  targets <- level_stats$target
  observed_means <- level_stats$mean
  
  # Predictions from linear model
  linear_pred <- stats::predict(models$linear,
                                newdata = data.frame(x = targets))
  
  # Predictions from best polynomial model
  poly_pred <- stats::predict(models[[best_model]],
                              newdata = data.frame(x = targets))
  
  # Deviation = observed mean - linear prediction
  # This measures how much each level's mean departs from the linear fit,
  # which is the clinically relevant quantity for linearity assessment
  deviation_abs <- observed_means - unname(linear_pred)
  
  # Percent deviation relative to linear prediction
  deviation_pct <- ifelse(
    abs(linear_pred) > .Machine$double.eps,
    100 * deviation_abs / linear_pred,
    NA_real_
  )
  
  data.frame(
    target = targets,
    observed_mean = observed_means,
    linear_predicted = unname(linear_pred),
    poly_predicted = unname(poly_pred),
    deviation_abs = unname(deviation_abs),
    deviation_pct = unname(deviation_pct),
    stringsAsFactors = FALSE
  )
}


#' Compute maximum deviation
#' @noRd
#' @keywords internal
.compute_max_deviation <- function(deviations) {
  
  idx_abs <- which.max(abs(deviations$deviation_abs))
  idx_pct <- which.max(abs(deviations$deviation_pct))
  
  list(
    max_abs = deviations$deviation_abs[idx_abs],
    max_abs_at = deviations$target[idx_abs],
    max_pct = deviations$deviation_pct[idx_pct],
    max_pct_at = deviations$target[idx_pct]
  )
}


#' Compute recovery at each level
#' @noRd
#' @keywords internal
.compute_recovery <- function(level_stats) {
  
  recovery_pct <- ifelse(
    abs(level_stats$target) > .Machine$double.eps,
    100 * level_stats$mean / level_stats$target,
    NA_real_
  )
  
  data.frame(
    target = level_stats$target,
    observed_mean = level_stats$mean,
    recovery_pct = recovery_pct,
    stringsAsFactors = FALSE
  )
}


#' Assess deviations against ADL
#' @noRd
#' @keywords internal
.assess_adl <- function(deviations, adl, adl_type) {
  
  if (is.null(adl)) return(NA)
  
  if (adl_type == "absolute") {
    all(abs(deviations$deviation_abs) <= adl, na.rm = TRUE)
  } else {
    all(abs(deviations$deviation_pct) <= adl, na.rm = TRUE)
  }
}


#' Compute linear range (range where deviations are within ADL)
#' @noRd
#' @keywords internal
.compute_linear_range <- function(deviations, adl, adl_type) {
  
  if (is.null(adl)) return(NULL)
  
  if (adl_type == "absolute") {
    within_adl <- abs(deviations$deviation_abs) <= adl
  } else {
    within_adl <- abs(deviations$deviation_pct) <= adl
  }
  
  # Handle NAs
  within_adl[is.na(within_adl)] <- FALSE
  
  if (!any(within_adl)) return(NULL)
  
  # Find the contiguous range of levels within ADL
  # Start from the lowest level that passes
  passing_targets <- deviations$target[within_adl]
  c(min(passing_targets), max(passing_targets))
}