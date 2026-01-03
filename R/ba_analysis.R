#' Bland-Altman Analysis for Method Comparison
#'
#' @description
#' Performs Bland-Altman analysis to assess agreement between two measurement
#' methods. Calculates bias (mean difference), limits of agreement, and
#' confidence intervals following the approach of Bland & Altman (1986, 1999).
#'
#' @param x Numeric vector of measurements from method 1 (reference method),
#'   or a formula of the form `method1 ~ method2`.
#' @param y Numeric vector of measurements from method 2 (test method).
#'   Ignored if `x` is a formula.
#' @param data Optional data frame containing the variables specified in
#'   `x` and `y` (or in the formula).
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param type Type of difference calculation: `"absolute"` (default) for
#'   `y - x`, or `"percent"` for `100 * (y - x) / mean`.
#' @param na_action How to handle missing values: `"omit"` (default) removes
#'   pairs with any NA, `"fail"` stops with an error.
#'
#' @return An object of class `c("ba_analysis", "valytics_comparison", "valytics_result")`,
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
#'         \item `differences`: Numeric vector of differences (y - x or percent)
#'         \item `averages`: Numeric vector of means ((x + y) / 2)
#'         \item `bias`: Mean difference (point estimate)
#'         \item `bias_se`: Standard error of the bias
#'         \item `bias_ci`: Named numeric vector with lower and upper CI for bias
#'         \item `sd_diff`: Standard deviation of differences
#'         \item `loa_lower`: Lower limit of agreement (bias - 1.96 * SD)
#'         \item `loa_upper`: Upper limit of agreement (bias + 1.96 * SD)
#'         \item `loa_lower_ci`: Named numeric vector with CI for lower LoA
#'         \item `loa_upper_ci`: Named numeric vector with CI for upper LoA
#'       }
#'     }
#'     \item{settings}{List with analysis settings:
#'       \itemize{
#'         \item `conf_level`: Confidence level used
#'         \item `type`: Type of difference calculation
#'         \item `multiplier`: Multiplier for LoA (default: 1.96 for 95\%)
#'       }
#'     }
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' The Bland-Altman method assesses agreement between two quantitative
#' measurements by analyzing the differences against the averages. The key
#' outputs are:
#'
#' \itemize{
#'   \item **Bias**: The mean difference between methods, indicating systematic
#'     difference. A bias significantly different from zero suggests one method
#'     consistently measures higher or lower than the other.
#'   \item **Limits of Agreement (LoA)**: The interval within which 95\% of
#'     differences are expected to lie (bias +/- 1.96 x SD). These define the
#'     range of disagreement between methods.
#'   \item **Confidence Intervals**: CIs for bias and LoA quantify the
#'     uncertainty in these estimates due to sampling variability.
#' }
#'
#' The confidence intervals for limits of agreement are calculated using the
#' exact method from Bland & Altman (1999), which accounts for the uncertainty
#' in both the mean and standard deviation.
#'
#' @section Assumptions:
#' The standard Bland-Altman analysis assumes:
#' \itemize{
#'   \item Differences are approximately normally distributed
#'   \item No proportional bias (constant bias across the measurement range)
#'   \item No repeated measurements per subject
#' }
#'
#' @references
#' Bland JM, Altman DG (1986). Statistical methods for assessing agreement
#' between two methods of clinical measurement. \emph{Lancet}, 1(8476):307-310.
#' \doi{10.1016/S0140-6736(86)90837-8}
#'
#' Bland JM, Altman DG (1999). Measuring agreement in method comparison studies.
#' \emph{Statistical Methods in Medical Research}, 8(2):135-160.
#' \doi{10.1177/096228029900800204}
#'
#' @seealso
#' [plot.ba_analysis()] for visualization,
#' [summary.ba_analysis()] for detailed summary
#'
#' @examples
#' # Simulated method comparison data
#' set.seed(42)
#' method_a <- rnorm(50, mean = 100, sd = 15)
#' method_b <- method_a + rnorm(50, mean = 2, sd = 5)  # Method B has +2 bias
#'
#' # Basic analysis
#' ba <- ba_analysis(method_a, method_b)
#' ba
#'
#' # Using formula interface with data frame
#' df <- data.frame(reference = method_a, test = method_b)
#' ba <- ba_analysis(reference ~ test, data = df)
#'
#' # Percentage differences
#' ba_pct <- ba_analysis(method_a, method_b, type = "percent")
#'
#' @export
ba_analysis <- function(x,
                        y = NULL,
                        data = NULL,
                        conf_level = 0.95,
                        type = c("absolute", "percent"),
                        na_action = c("omit", "fail")) {

  # Capture the call for reproducibility
  call <- match.call()

  # Match arguments
  type <- match.arg(type)
  na_action <- match.arg(na_action)

  # Input parsing ----
  parsed <- .parse_ba_input(x, y, data)
  x_vals <- parsed$x
  y_vals <- parsed$y
  var_names <- parsed$var_names

  # Input validation ----
  .validate_ba_input(x_vals, y_vals, conf_level)

  # Handle missing values ----
  complete_cases <- stats::complete.cases(x_vals, y_vals)
  n_excluded <- sum(!complete_cases)

  if (na_action == "fail" && n_excluded > 0) {
    stop("Missing values detected. Use `na_action = 'omit'` to remove them.",
         call. = FALSE)
  }

  x_clean <- x_vals[complete_cases]
  y_clean <- y_vals[complete_cases]
  n <- length(x_clean)

  if (n < 3) {
    stop("At least 3 complete paired observations are required.", call. = FALSE)
  }

  # Core calculations ----
  results <- .compute_ba_statistics(
    x = x_clean,
    y = y_clean,
    conf_level = conf_level,
    type = type
  )

  # Construct output object ----
  structure(
    list(
      input = list(
        x = x_clean,
        y = y_clean,
        n = n,
        n_excluded = n_excluded,
        var_names = var_names
      ),
      results = results,
      settings = list(
        conf_level = conf_level,
        type = type,
        multiplier = stats::qnorm(1 - (1 - conf_level) / 2)
      ),
      call = call
    ),
    class = c("ba_analysis", "valytics_comparison", "valytics_result")
  )
}


# Helper Functions ----

#' Parse input for ba_analysis
#' @noRd
#' @keywords internal
.parse_ba_input <- function(x, y, data) {

  # Formula interface: x ~ y
  if (inherits(x, "formula")) {
    if (!is.null(y)) {
      warning("`y` is ignored when `x` is a formula.", call. = FALSE)
    }

    # Extract variable names from formula
    formula_vars <- all.vars(x)
    if (length(formula_vars) != 2) {
      stop("Formula must have exactly two variables: method1 ~ method2",
           call. = FALSE)
    }

    var_names <- c(x = formula_vars[1], y = formula_vars[2])

    # Get data from environment or data frame
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

    # Try to get variable names from call
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

  list(x = x_vals, y = y_vals, var_names = var_names)
}


#' Validate input for ba_analysis
#' @noRd
#' @keywords internal
.validate_ba_input <- function(x, y, conf_level) {

  # Check numeric
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(y)) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }

  # Check equal length
  if (length(x) != length(y)) {
    stop("`x` and `y` must have the same length.", call. = FALSE)
  }


  # Check confidence level
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number between 0 and 1.",
         call. = FALSE)
  }

  invisible(TRUE)
}


#' Compute Bland-Altman statistics
#' @noRd
#' @keywords internal
.compute_ba_statistics <- function(x, y, conf_level, type) {

  n <- length(x)

  # Calculate averages (always the same regardless of type)
  averages <- (x + y) / 2

  # Calculate differences based on type
  if (type == "absolute") {
    differences <- y - x
  } else {
    # Percent difference: 100 * (y - x) / average
    differences <- 100 * (y - x) / averages
  }

  # Bias (mean difference)
  bias <- mean(differences)

  # Standard deviation of differences
  sd_diff <- stats::sd(differences)

  # Standard error of bias
  bias_se <- sd_diff / sqrt(n)

  # Multiplier for confidence intervals (e.g., 1.96 for 95%)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  # CI for bias (based on t-distribution)
  t_crit <- stats::qt(1 - (1 - conf_level) / 2, df = n - 1)
  bias_ci <- c(
    lower = bias - t_crit * bias_se,
    upper = bias + t_crit * bias_se
  )

  # Limits of agreement
  loa_lower <- bias - z * sd_diff
  loa_upper <- bias + z * sd_diff

  # CI for limits of agreement (Bland & Altman 1999)
  # Variance of LoA = Var(mean) + z^2 * Var(SD)
  # Var(SD) approximately SD^2 / (2 * (n - 1)) for normal data
  # SE of LoA = SD * sqrt(1/n + z^2 / (2 * (n - 1)))
  se_loa <- sd_diff * sqrt(1/n + z^2 / (2 * (n - 1)))

  loa_lower_ci <- c(
    lower = loa_lower - t_crit * se_loa,
    upper = loa_lower + t_crit * se_loa
  )

  loa_upper_ci <- c(
    lower = loa_upper - t_crit * se_loa,
    upper = loa_upper + t_crit * se_loa
  )

  list(
    differences = differences,
    averages = averages,
    bias = bias,
    bias_se = bias_se,
    bias_ci = bias_ci,
    sd_diff = sd_diff,
    loa_lower = loa_lower,
    loa_upper = loa_upper,
    loa_lower_ci = loa_lower_ci,
    loa_upper_ci = loa_upper_ci
  )
}
