#' Precision Profile Analysis
#'
#' @description
#' Constructs a precision profile (CV vs concentration relationship) from
#' precision study results and estimates functional sensitivity. The precision
#' profile characterizes how measurement imprecision changes across the
#' analytical measurement interval.
#'
#' @param x An object of class `precision_study` with multiple concentration
#'   levels, OR a data frame with columns for concentration and CV values.
#' @param concentration Character string specifying the column name for
#'   concentration values (only used if `x` is a data frame). Default is
#'   `"concentration"`.
#' @param cv Character string specifying the column name for CV values (only
#'   used if `x` is a data frame). Default is `"cv_pct"`.
#' @param model Regression model for CV-concentration relationship:
#'   `"hyperbolic"` (default) fits CV = sqrt(a^2 + (b/x)^2),
#'   `"linear"` fits CV = a + b/x.
#' @param cv_targets Numeric vector of target CV percentages for functional
#'   sensitivity estimation. Default is `c(10, 20)`.
#' @param conf_level Confidence level for prediction intervals (default: 0.95).
#' @param bootstrap Logical; if `TRUE`, uses bootstrap resampling for
#'   confidence intervals on functional sensitivity estimates. Default is `FALSE`.
#' @param boot_n Number of bootstrap resamples when `bootstrap = TRUE`
#'   (default: 1999).
#'
#' @return An object of class `c("precision_profile", "valytics_precision", "valytics_result")`,
#'   which is a list containing:
#'
#'   \describe{
#'     \item{input}{List with original data:
#'       \itemize{
#'         \item `concentration`: Numeric vector of concentrations
#'         \item `cv`: Numeric vector of CV values (percent)
#'         \item `n_levels`: Number of concentration levels
#'         \item `conc_range`: Concentration range (min, max)
#'         \item `conc_span`: Fold-difference (max/min)
#'       }
#'     }
#'     \item{model}{List with fitted model information:
#'       \itemize{
#'         \item `type`: Model type ("hyperbolic" or "linear")
#'         \item `parameters`: Named vector of fitted parameters
#'         \item `equation`: Character string describing the fitted equation
#'       }
#'     }
#'     \item{fitted}{Data frame with fitted values:
#'       \itemize{
#'         \item `concentration`: Concentration values
#'         \item `cv_observed`: Observed CV values
#'         \item `cv_fitted`: Model-fitted CV values
#'         \item `residual`: Residuals (observed - fitted)
#'         \item `ci_lower`: Lower prediction interval
#'         \item `ci_upper`: Upper prediction interval
#'       }
#'     }
#'     \item{fit_quality}{List with goodness-of-fit statistics:
#'       \itemize{
#'         \item `r_squared`: Coefficient of determination
#'         \item `adj_r_squared`: Adjusted R-squared
#'         \item `rmse`: Root mean squared error
#'         \item `mae`: Mean absolute error
#'       }
#'     }
#'     \item{functional_sensitivity}{Data frame with functional sensitivity estimates:
#'       \itemize{
#'         \item `cv_target`: Target CV percentage
#'         \item `concentration`: Estimated concentration at target CV
#'         \item `ci_lower`: Lower confidence limit (if bootstrap)
#'         \item `ci_upper`: Upper confidence limit (if bootstrap)
#'         \item `achievable`: Logical; TRUE if target CV is achievable
#'       }
#'     }
#'     \item{settings}{List with analysis settings}
#'     \item{call}{The matched function call}
#'   }
#'
#' @details
#' **Precision Profile:**
#'
#' The precision profile describes how analytical imprecision (CV) varies
#' across the analytical measurement interval. Typically, CV decreases as
#' concentration increases, following a hyperbolic relationship.
#'
#' **Hyperbolic Model:**
#'
#' The hyperbolic model is:
#' \deqn{CV = \sqrt{a^2 + (b/x)^2}}
#'
#' where:
#' - `a` represents the asymptotic CV at high concentrations
#' - `b` represents the concentration-dependent component
#' - `x` is the analyte concentration
#'
#' This model captures the characteristic behavior where CV approaches a
#' constant value at high concentrations and increases hyperbolically at
#' low concentrations.
#'
#' **Linear Model:**
#'
#' The linear model is:
#' \deqn{CV = a + b/x}
#'
#' This is a simpler alternative that may be appropriate when the relationship
#' is approximately linear when plotted as CV vs 1/concentration.
#'
#' **Functional Sensitivity:**
#'
#' Functional sensitivity is defined as the lowest concentration at which a
#' measurement procedure achieves a specified level of precision (CV). Common
#' thresholds are:
#' - **10% CV**: Modern standard for high-sensitivity assays (e.g., cardiac troponin)
#' - **20% CV**: Traditional standard (originally defined for TSH assays)
#'
#' The functional sensitivity is calculated by solving the fitted model equation
#' for the concentration that yields the target CV.
#'
#' @section Minimum Requirements:
#' - At least 4 concentration levels
#' - Concentration span of at least 2-fold (warning if less)
#' - Valid CV estimates at each level (from precision study)
#'
#' @references
#' Armbruster DA, Pry T (2008). Limit of blank, limit of detection and limit of
#' quantitation. \emph{Clinical Biochemist Reviews}, 29(Suppl 1):S49-S52.
#'
#' CLSI EP17-A2 (2012). Evaluation of Detection Capability for Clinical
#' Laboratory Measurement Procedures; Approved Guideline - Second Edition.
#' Clinical and Laboratory Standards Institute, Wayne, PA.
#'
#' Kroll MH, Emancipator K (1993). A theoretical evaluation of linearity.
#' \emph{Clinical Chemistry}, 39(3):405-413.
#'
#' @seealso
#' [precision_study()] for variance component analysis,
#' [plot.precision_profile()] for visualization
#'
#' @examples
#' # Example with simulated multi-level precision data
#' set.seed(42)
#'
#' # Generate data for 6 concentration levels
#' conc_levels <- c(5, 10, 25, 50, 100, 200)
#' n_levels <- length(conc_levels)
#'
#' prec_data <- data.frame()
#' for (i in seq_along(conc_levels)) {
#'   level_data <- expand.grid(
#'     level = conc_levels[i],
#'     day = 1:5,
#'     replicate = 1:5
#'   )
#'   
#'   # Simulate CV that decreases with concentration
#'   true_cv <- sqrt(3^2 + (20/conc_levels[i])^2)
#'   level_data$value <- conc_levels[i] * rnorm(
#'     nrow(level_data),
#'     mean = 1,
#'     sd = true_cv/100
#'   )
#'   
#'   prec_data <- rbind(prec_data, level_data)
#' }
#'
#' # Run precision study
#' prec <- precision_study(
#'   data = prec_data,
#'   value = "value",
#'   sample = "level",
#'   day = "day"
#' )
#'
#' # Generate precision profile
#' profile <- precision_profile(prec)
#' print(profile)
#' summary(profile)
#'
#' # Hyperbolic model with bootstrap CIs
#' profile_boot <- precision_profile(
#'   prec,
#'   model = "hyperbolic",
#'   cv_targets = c(10, 20),
#'   bootstrap = TRUE,
#'   boot_n = 499
#' )
#'
#' # Linear model
#' profile_linear <- precision_profile(prec, model = "linear")
#'
#' @export
precision_profile <- function(x,
                              concentration = "concentration",
                              cv = "cv_pct",
                              model = c("hyperbolic", "linear"),
                              cv_targets = c(10, 20),
                              conf_level = 0.95,
                              bootstrap = FALSE,
                              boot_n = 1999) {
  
  # Capture the call
  call <- match.call()
  
  # Match arguments
  model <- match.arg(model)
  
  # Input parsing ----
  parsed <- .parse_profile_input(
    x = x,
    concentration = concentration,
    cv = cv
  )
  
  conc_vec <- parsed$concentration
  cv_vec <- parsed$cv
  
  # Input validation ----
  .validate_profile_input(
    concentration = conc_vec,
    cv = cv_vec,
    cv_targets = cv_targets,
    conf_level = conf_level,
    boot_n = boot_n
  )
  
  # Calculate input statistics ----
  n_levels <- length(conc_vec)
  conc_range <- c(min = min(conc_vec), max = max(conc_vec))
  conc_span <- conc_range["max"] / conc_range["min"]
  
  # Warning for narrow concentration span
  if (conc_span < 2) {
    warning(
      sprintf(
        "Concentration span is only %.2f-fold. A span of at least 2-fold is recommended for reliable precision profile estimation.",
        conc_span
      ),
      call. = FALSE
    )
  }
  
  # Fit model ----
  fit_result <- .fit_precision_model(
    concentration = conc_vec,
    cv = cv_vec,
    model = model
  )
  
  # Calculate prediction intervals ----
  pred_intervals <- .calculate_prediction_intervals(
    concentration = conc_vec,
    cv = cv_vec,
    fitted = fit_result$fitted,
    residual_se = fit_result$residual_se,
    conf_level = conf_level
  )
  
  # Assemble fitted data frame
  fitted_df <- data.frame(
    concentration = conc_vec,
    cv_observed = cv_vec,
    cv_fitted = fit_result$fitted,
    residual = cv_vec - fit_result$fitted,
    ci_lower = pred_intervals$lower,
    ci_upper = pred_intervals$upper
  )
  
  # Calculate functional sensitivity ----
  func_sens <- .calculate_functional_sensitivity(
    parameters = fit_result$parameters,
    model = model,
    cv_targets = cv_targets,
    concentration = conc_vec,
    cv = cv_vec,
    bootstrap = bootstrap,
    boot_n = boot_n,
    conf_level = conf_level
  )
  
  # Assemble output ----
  result <- structure(
    list(
      input = list(
        concentration = conc_vec,
        cv = cv_vec,
        n_levels = n_levels,
        conc_range = conc_range,
        conc_span = conc_span
      ),
      model = list(
        type = model,
        parameters = fit_result$parameters,
        equation = fit_result$equation
      ),
      fitted = fitted_df,
      fit_quality = fit_result$fit_quality,
      functional_sensitivity = func_sens,
      settings = list(
        conf_level = conf_level,
        bootstrap = bootstrap,
        boot_n = if (bootstrap) boot_n else NA
      ),
      call = call
    ),
    class = c("precision_profile", "valytics_precision", "valytics_result")
  )
  
  result
}


# Helper Functions ----

#' Parse input for precision_profile
#' @noRd
.parse_profile_input <- function(x, concentration, cv) {
  
  # Check if x is a precision_study object
  if (inherits(x, "precision_study")) {
    
    # Verify that the study has multiple samples
    if (is.null(x$by_sample) || length(x$by_sample) < 4) {
      stop(
        "precision_profile() requires a precision_study object with at least 4 concentration levels. ",
        "Found: ", if (is.null(x$by_sample)) 1 else length(x$by_sample),
        call. = FALSE
      )
    }
    
    # Extract concentration from sample_means (stored at top level)
    if (is.null(x$sample_means)) {
      stop(
        "precision_study object does not contain sample_means. ",
        "This should not happen - please report as a bug.",
        call. = FALSE
      )
    }
    
    conc_vec <- as.numeric(x$sample_means)
    cv_vec <- numeric(length(x$by_sample))
    
    for (i in seq_along(x$by_sample)) {
      sample_result <- x$by_sample[[i]]
      
      # Get within-laboratory precision CV (or repeatability if not available)
      prec_summary <- sample_result$precision
      
      # Try to get within-laboratory precision, fall back to repeatability
      # Note: measure names are title case with hyphens
      if ("Within-laboratory precision" %in% prec_summary$measure) {
        cv_vec[i] <- prec_summary$cv_pct[prec_summary$measure == "Within-laboratory precision"]
      } else if ("Repeatability" %in% prec_summary$measure) {
        cv_vec[i] <- prec_summary$cv_pct[prec_summary$measure == "Repeatability"]
      } else {
        # Use the first available CV
        cv_vec[i] <- prec_summary$cv_pct[1]
      }
    }
    
  } else if (is.data.frame(x)) {
    
    # Data frame interface
    if (!concentration %in% names(x)) {
      stop(
        "Column '", concentration, "' not found in data frame.",
        call. = FALSE
      )
    }
    
    if (!cv %in% names(x)) {
      stop(
        "Column '", cv, "' not found in data frame.",
        call. = FALSE
      )
    }
    
    conc_vec <- x[[concentration]]
    cv_vec <- x[[cv]]
    
  } else {
    stop(
      "x must be either a precision_study object with multiple samples or a data frame.",
      call. = FALSE
    )
  }
  
  list(
    concentration = conc_vec,
    cv = cv_vec
  )
}


#' Validate input for precision_profile
#' @noRd
.validate_profile_input <- function(concentration, cv, cv_targets,
                                    conf_level, boot_n) {
  
  # Check that vectors are numeric
  if (!is.numeric(concentration)) {
    stop("concentration must be numeric.", call. = FALSE)
  }
  
  if (!is.numeric(cv)) {
    stop("cv must be numeric.", call. = FALSE)
  }
  
  # Check equal length
  if (length(concentration) != length(cv)) {
    stop("concentration and cv must have the same length.", call. = FALSE)
  }
  
  # Check minimum number of levels
  if (length(concentration) < 4) {
    stop(
      "At least 4 concentration levels are required. Found: ",
      length(concentration),
      call. = FALSE
    )
  }
  
  # Check for NAs
  if (any(is.na(concentration)) || any(is.na(cv))) {
    stop("concentration and cv must not contain NA values.", call. = FALSE)
  }
  
  # Check for non-positive values
  if (any(concentration <= 0)) {
    stop("concentration values must be positive.", call. = FALSE)
  }
  
  if (any(cv <= 0)) {
    stop("cv values must be positive.", call. = FALSE)
  }
  
  # Check CV targets
  if (!is.numeric(cv_targets) || any(cv_targets <= 0) || any(cv_targets > 100)) {
    stop("cv_targets must be positive numbers between 0 and 100.", call. = FALSE)
  }
  
  # Check conf_level
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number between 0 and 1.", call. = FALSE)
  }
  
  # Check boot_n
  if (!is.numeric(boot_n) || length(boot_n) != 1 ||
      boot_n < 100 || boot_n != floor(boot_n)) {
    stop("boot_n must be an integer >= 100.", call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Fit precision profile model
#' @noRd
.fit_precision_model <- function(concentration, cv, model) {
  
  if (model == "hyperbolic") {
    fit <- .fit_hyperbolic_model(concentration, cv)
  } else {
    fit <- .fit_linear_model(concentration, cv)
  }
  
  fit
}


#' Fit hyperbolic model: CV = sqrt(a^2 + (b/x)^2)
#' @noRd
.fit_hyperbolic_model <- function(concentration, cv) {
  
  # Transform to linearized form for initial parameter estimates
  # CV^2 = a^2 + (b/x)^2
  # CV^2 = a^2 + b^2 * (1/x^2)
  
  cv_sq <- cv^2
  inv_conc_sq <- 1 / (concentration^2)
  
  # Linear regression: CV^2 ~ 1 + 1/x^2
  lm_init <- lm(cv_sq ~ inv_conc_sq)
  
  # Initial parameter estimates
  a_init <- sqrt(max(0, coef(lm_init)[1]))
  b_init <- sqrt(max(0, coef(lm_init)[2]))
  
  # Ensure reasonable starting values
  if (a_init < 0.01) a_init <- min(cv) * 0.5
  if (b_init < 0.01) b_init <- min(cv) * min(concentration) * 0.5
  
  # Non-linear least squares for refined estimates
  # Use try-catch in case nls fails
  fit_nls <- tryCatch(
    {
      nls(
        cv ~ sqrt(a^2 + (b/concentration)^2),
        start = list(a = a_init, b = b_init),
        control = nls.control(maxiter = 200, warnOnly = TRUE)
      )
    },
    error = function(e) {
      # If NLS fails, return NULL to use linearized estimates
      NULL
    },
    warning = function(w) {
      # Suppress convergence warnings, will use linearized fallback
      NULL
    }
  )
  
  if (!is.null(fit_nls) && !any(is.na(coef(fit_nls)))) {
    # Use NLS estimates if successful
    params <- coef(fit_nls)
    a <- as.numeric(params["a"])
    b <- as.numeric(params["b"])
    fitted <- fitted(fit_nls)
  } else {
    # Fall back to linearized estimates
    a <- a_init
    b <- b_init
    fitted <- sqrt(a^2 + (b/concentration)^2)
  }
  
  # Residuals and fit statistics
  residuals <- cv - fitted
  ss_res <- sum(residuals^2)
  ss_tot <- sum((cv - mean(cv))^2)
  
  r_squared <- 1 - (ss_res / ss_tot)
  n <- length(cv)
  p <- 2  # Number of parameters
  adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))
  
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  residual_se <- sqrt(ss_res / (n - p))
  
  # Equation string
  equation <- sprintf("CV = sqrt(%.3f^2 + (%.3f/x)^2)", a, b)
  
  list(
    parameters = c(a = a, b = b),
    fitted = fitted,
    equation = equation,
    residual_se = residual_se,
    fit_quality = list(
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      rmse = rmse,
      mae = mae
    )
  )
}


#' Fit linear model: CV = a + b/x
#' @noRd
.fit_linear_model <- function(concentration, cv) {
  
  # Transform: CV = a + b * (1/x)
  inv_conc <- 1 / concentration
  
  # Linear regression
  lm_fit <- lm(cv ~ inv_conc)
  
  params <- coef(lm_fit)
  a <- as.numeric(params[1])
  b <- as.numeric(params[2])
  
  fitted <- fitted(lm_fit)
  residuals <- residuals(lm_fit)
  
  # Fit statistics
  ss_res <- sum(residuals^2)
  ss_tot <- sum((cv - mean(cv))^2)
  
  r_squared <- 1 - (ss_res / ss_tot)
  n <- length(cv)
  p <- 2
  adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))
  
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  residual_se <- summary(lm_fit)$sigma
  
  # Equation string
  equation <- sprintf("CV = %.3f + %.3f/x", a, b)
  
  list(
    parameters = c(a = a, b = b),
    fitted = fitted,
    equation = equation,
    residual_se = residual_se,
    fit_quality = list(
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      rmse = rmse,
      mae = mae
    )
  )
}


#' Calculate prediction intervals for fitted CV values
#' @noRd
.calculate_prediction_intervals <- function(concentration, cv, fitted,
                                            residual_se, conf_level) {
  
  n <- length(cv)
  
  # t-critical value
  t_crit <- qt(1 - (1 - conf_level) / 2, df = n - 2)
  
  # Prediction interval
  # SE of prediction = residual_se * sqrt(1 + 1/n + (x - mean(x))^2 / sum((x - mean(x))^2))
  
  # For simplicity, use constant prediction interval width
  # (Proper prediction intervals would require the model matrix)
  pred_se <- residual_se * sqrt(1 + 1/n)
  
  lower <- fitted - t_crit * pred_se
  upper <- fitted + t_crit * pred_se
  
  # Ensure non-negative
  lower <- pmax(lower, 0)
  
  list(lower = lower, upper = upper)
}


#' Calculate functional sensitivity
#' @noRd
.calculate_functional_sensitivity <- function(parameters, model, cv_targets,
                                              concentration, cv, bootstrap,
                                              boot_n, conf_level) {
  
  # Calculate point estimates
  func_sens_df <- data.frame(
    cv_target = cv_targets,
    concentration = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    achievable = FALSE
  )
  
  for (i in seq_along(cv_targets)) {
    target <- cv_targets[i]
    
    conc_est <- .solve_for_concentration(
      target_cv = target,
      parameters = parameters,
      model = model
    )
    
    func_sens_df$concentration[i] <- conc_est$concentration
    func_sens_df$achievable[i] <- conc_est$achievable
  }
  
  # Bootstrap confidence intervals if requested
  if (bootstrap && any(func_sens_df$achievable)) {
    
    boot_results <- .bootstrap_functional_sensitivity(
      concentration = concentration,
      cv = cv,
      model = model,
      cv_targets = cv_targets,
      boot_n = boot_n,
      conf_level = conf_level
    )
    
    # Merge bootstrap CIs
    for (i in seq_along(cv_targets)) {
      if (func_sens_df$achievable[i]) {
        func_sens_df$ci_lower[i] <- boot_results$ci_lower[i]
        func_sens_df$ci_upper[i] <- boot_results$ci_upper[i]
      }
    }
  }
  
  func_sens_df
}


#' Solve fitted model for concentration at target CV
#' @noRd
.solve_for_concentration <- function(target_cv, parameters, model) {
  
  a <- parameters["a"]
  b <- parameters["b"]
  
  # Check for NA parameters
  if (is.na(a) || is.na(b)) {
    return(list(concentration = NA_real_, achievable = FALSE))
  }
  
  if (model == "hyperbolic") {
    # CV_target = sqrt(a^2 + (b/x)^2)
    # CV_target^2 = a^2 + (b/x)^2
    # (b/x)^2 = CV_target^2 - a^2
    
    if (target_cv^2 < a^2) {
      # Target CV is below the asymptotic minimum
      return(list(concentration = NA_real_, achievable = FALSE))
    }
    
    # x = b / sqrt(CV_target^2 - a^2)
    conc <- abs(b) / sqrt(target_cv^2 - a^2)
    
  } else {
    # Linear: CV_target = a + b/x
    # b/x = CV_target - a
    # x = b / (CV_target - a)
    
    if (b / (target_cv - a) <= 0) {
      # Solution not feasible
      return(list(concentration = NA_real_, achievable = FALSE))
    }
    
    conc <- b / (target_cv - a)
  }
  
  list(concentration = conc, achievable = TRUE)
}


#' Bootstrap functional sensitivity estimates
#' @noRd
.bootstrap_functional_sensitivity <- function(concentration, cv, model,
                                              cv_targets, boot_n, conf_level) {
  
  n <- length(concentration)
  
  # Storage for bootstrap estimates
  boot_conc <- matrix(NA, nrow = boot_n, ncol = length(cv_targets))
  
  for (b in seq_len(boot_n)) {
    # Bootstrap sample
    idx <- sample.int(n, n, replace = TRUE)
    conc_boot <- concentration[idx]
    cv_boot <- cv[idx]
    
    # Fit model
    fit_boot <- tryCatch(
      {
        .fit_precision_model(conc_boot, cv_boot, model)
      },
      error = function(e) NULL
    )
    
    if (!is.null(fit_boot)) {
      # Solve for each target
      for (i in seq_along(cv_targets)) {
        result <- .solve_for_concentration(
          target_cv = cv_targets[i],
          parameters = fit_boot$parameters,
          model = model
        )
        
        if (result$achievable) {
          boot_conc[b, i] <- result$concentration
        }
      }
    }
  }
  
  # Calculate BCa confidence intervals
  ci_lower <- numeric(length(cv_targets))
  ci_upper <- numeric(length(cv_targets))
  
  for (i in seq_along(cv_targets)) {
    boot_vals <- boot_conc[, i]
    boot_vals <- boot_vals[!is.na(boot_vals)]
    
    if (length(boot_vals) >= boot_n * 0.9) {
      # Sufficient bootstrap samples succeeded
      ci <- quantile(boot_vals, probs = c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2))
      ci_lower[i] <- ci[1]
      ci_upper[i] <- ci[2]
    } else {
      ci_lower[i] <- NA_real_
      ci_upper[i] <- NA_real_
    }
  }
  
  list(ci_lower = ci_lower, ci_upper = ci_upper)
}