#' Precision Study Analysis
#'
#' @description
#' Performs variance component analysis for precision experiments following
#' established methodology for clinical laboratory method validation. Estimates
#' repeatability, intermediate precision, and reproducibility from nested
#' experimental designs.
#'
#' @param data A data frame containing the precision experiment data.
#' @param value Character string specifying the column name containing
#'   measurement values. Default is `"value"`.
#' @param sample Character string specifying the column name for sample/level
#'   identifier. Use when multiple concentration levels are tested. Default
#'   is `NULL` (single sample).
#' @param site Character string specifying the column name for site/device
#'   identifier. Use for multi-site reproducibility studies. Default is `NULL`
#'   (single site).
#' @param day Character string specifying the column name for day identifier.
#'   Default is `"day"`.
#' @param run Character string specifying the column name for run identifier
#'   (within day). Default is `NULL` (assumes single run per day).
#' @param replicate Character string specifying the column name for replicate
#'   identifier. If `NULL` (default), replicates are inferred from the data
#'   structure.
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param ci_method Method for calculating confidence intervals:
#'   `"satterthwaite"` (default) uses the Satterthwaite approximation,
#'   `"mls"` uses the Modified Large Sample method,
#'   `"bootstrap"` uses BCa bootstrap resampling.
#' @param boot_n Number of bootstrap resamples when `ci_method = "bootstrap"`
#'   (default: 1999).
#' @param method Estimation method for variance components:
#'   `"anova"` (default) uses ANOVA-based method of moments,
#'   `"reml"` uses Restricted Maximum Likelihood (requires lme4 package).
#'
#' @return An object of class `c("precision_study", "valytics_precision", "valytics_result")`,
#'   which is a list containing:
#'
#'   \describe{
#'     \item{input}{List with original data and metadata:
#'       \itemize{
#'         \item `data`: The input data frame (after validation)
#'         \item `n`: Total number of observations
#'         \item `n_excluded`: Number of observations excluded due to NAs
#'         \item `factors`: Named list of factor column names used
#'         \item `value_col`: Name of the value column
#'       }
#'     }
#'     \item{design}{List describing the experimental design:
#'       \itemize{
#'         \item `type`: Design type (e.g., "single_site", "multi_site")
#'         \item `structure`: Character string describing nesting (e.g., "day/run")
#'         \item `levels`: Named list with number of levels for each factor
#'         \item `balanced`: Logical; TRUE if design is balanced
#'         \item `n_samples`: Number of distinct samples/concentration levels
#'       }
#'     }
#'     \item{variance_components}{Data frame with variance component estimates:
#'       \itemize{
#'         \item `component`: Name of variance component
#'         \item `variance`: Estimated variance
#'         \item `sd`: Standard deviation (sqrt of variance
#'         \item `pct_total`: Percentage of total variance
#'         \item `df`: Degrees of freedom
#'       }
#'     }
#'     \item{precision}{Data frame with precision estimates:
#'       \itemize{
#'         \item `measure`: Precision measure name (repeatability, intermediate, etc.)
#'         \item `sd`: Standard deviation
#'         \item `cv_pct`: Coefficient of variation (percent)
#'         \item `ci_lower`: Lower confidence limit
#'         \item `ci_upper`: Upper confidence limit
#'       }
#'     }
#'     \item{anova_table}{ANOVA table with SS, MS, DF for each source of variation}
#'     \item{by_sample}{If multiple samples: list of results per sample}
#'     \item{settings}{List with analysis settings}
#'     \item{call}{The matched function call}
#'   }
#'
#' @details
#' This function implements variance component analysis for nested experimental
#' designs commonly used in clinical laboratory precision studies. The analysis
#' follows methodology consistent with international standards.
#'
#' **Supported Experimental Designs:**
#'
#' \itemize{
#'   \item **Single-site, day/run/replicate**: Classic 20 x 2 x 2 design
#'     (20 days, 2 runs per day, 2 replicates per run)
#'   \item **Single-site, day/replicate**: Simplified design without run factor
#'     (e.g., 5 days x 5 replicates for verification)
#'   \item **Multi-site**: 3 sites x 5 days x 5 replicates for reproducibility
#'   \item **Custom designs**: Any fully-nested combination of factors
#' }
#'
#' **Variance Components:**
#'
#' For a design with site/day/run/replicate, the model is:
#' \deqn{y_{ijkl} = \mu + S_i + D_{j(i)} + R_{k(ij)} + \epsilon_{l(ijk)}}
#'
#' where S = site, D = day (nested in site), R = run (nested in day),
#' and epsilon = residual error.
#'
#' **Precision Measures:**
#'
#' \itemize{
#'   \item **Repeatability**: Within-run precision (sqrt of error variance)
#'   \item **Between-run precision**: Additional variability between runs
#'   \item **Between-day precision**: Additional variability between days
#'   \item **Intermediate precision**: Within-laboratory precision
#'     (combines day, run, and error variance)
#'   \item **Reproducibility**: Total precision including between-site
#'     variability (for multi-site designs)
#' }
#'
#' @section Confidence Intervals:
#' Three methods are available for confidence interval estimation:
#'
#' \itemize{
#'   \item **Satterthwaite** (default): Uses Satterthwaite's approximation
#'     for degrees of freedom of linear combinations of variance components.
#'   \item **MLS**: Modified Large Sample method, which can provide better
#'     coverage when variance components may be estimated as negative.
#'   \item **Bootstrap**: BCa bootstrap resampling. Most robust but
#'     computationally intensive.
#' }
#'
#' @section ANOVA vs REML:
#' \itemize{
#'   \item **ANOVA** (default): Method of moments estimation. Works well for
#'     balanced designs. May produce negative variance estimates for small
#'     variance components (set to zero by default).
#'   \item **REML**: Restricted Maximum Likelihood. Preferred for unbalanced
#'     designs. Requires the lme4 package. Always produces non-negative
#'     estimates.
#' }
#'
#' @references
#' Chesher D (2008). Evaluating assay precision. \emph{Clinical Biochemist
#' Reviews}, 29(Suppl 1):S23-S26.
#'
#' ISO 5725-2:2019. Accuracy (trueness and precision) of measurement methods
#' and results - Part 2: Basic method for the determination of repeatability
#' and reproducibility of a standard measurement method.
#'
#' Searle SR, Casella G, McCulloch CE (1992). \emph{Variance Components}.
#' Wiley, New York.
#'
#' Satterthwaite FE (1946). An approximate distribution of estimates of
#' variance components. \emph{Biometrics Bulletin}, 2:110-114.
#'
#' @seealso
#' [verify_precision()] for comparing results to manufacturer claims,
#' [plot.precision_study()] for visualization,
#' [summary.precision_study()] for detailed summary
#'
#' @examples
#' # Example with simulated precision data
#' set.seed(42)
#'
#' # Generate study design: 20 days x 2 runs x 2 replicates
#' n_days <- 20
#' n_runs <- 2
#' n_reps <- 2
#'
#' prec_data <- expand.grid(
#'   day = 1:n_days,
#'   run = 1:n_runs,
#'   replicate = 1:n_reps
#' )
#'
#' # Add realistic variance components
#' day_effect <- rep(rnorm(n_days, 0, 1.5), each = n_runs * n_reps)
#' run_effect <- rep(rnorm(n_days * n_runs, 0, 1.0), each = n_reps)
#' error <- rnorm(nrow(prec_data), 0, 2.0)
#'
#' prec_data$value <- 100 + day_effect + run_effect + error
#'
#' # Run precision study
#' prec <- precision_study(
#'   data = prec_data,
#'   value = "value",
#'   day = "day",
#'   run = "run"
#' )
#'
#' print(prec)
#' summary(prec)
#'
#' @export
precision_study <- function(data,
                            value = "value",
                            sample = NULL,
                            site = NULL,
                            day = "day",
                            run = NULL,
                            replicate = NULL,
                            conf_level = 0.95,
                            ci_method = c("satterthwaite", "mls", "bootstrap"),
                            boot_n = 1999,
                            method = c("anova", "reml")) {
  
  
  # Capture the call
  
  call <- match.call()
  
  # Match arguments
  
  ci_method <- match.arg(ci_method)
  method <- match.arg(method)
  
  # Check REML availability
  
  if (method == "reml") {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package 'lme4' is required for REML estimation. ",
           "Install it with install.packages('lme4') or use method = 'anova'.",
           call. = FALSE)
    }
  }
  
  # Input validation ----
  validated <- .validate_precision_input(
    data = data,
    value = value,
    sample = sample,
    site = site,
    day = day,
    run = run,
    replicate = replicate,
    conf_level = conf_level,
    boot_n = boot_n
  )
  
  data_clean <- validated$data
  factors <- validated$factors
  n_excluded <- validated$n_excluded
  
  # Design detection ----
  design <- .detect_precision_design(data_clean, factors)
  
  # Convert factors to proper factor type for analysis
  data_clean <- .prepare_factors(data_clean, factors)
  
  # Check for multiple samples ----
  if (!is.null(factors$sample)) {
    # Analyze each sample separately, then combine
    samples <- unique(data_clean[[factors$sample]])
    n_samples <- length(samples)
    
    results_by_sample <- lapply(samples, function(s) {
      sample_data <- data_clean[data_clean[[factors$sample]] == s, ]
      sample_factors <- factors
      sample_factors$sample <- NULL  # Remove sample from nested structure
      
      sample_design <- .detect_precision_design(sample_data, sample_factors)
      
      .compute_precision_single(
        data = sample_data,
        factors = sample_factors,
        design = sample_design,
        conf_level = conf_level,
        ci_method = ci_method,
        boot_n = boot_n,
        method = method,
        value_col = value
      )
    })
    names(results_by_sample) <- as.character(samples)
    
    # Get overall mean for CV calculation
    overall_means <- sapply(samples, function(s) {
      mean(data_clean[[value]][data_clean[[factors$sample]] == s], na.rm = TRUE)
    })
    names(overall_means) <- as.character(samples)
    
    # Use first sample's structure for main results (typical case)
    # Users can access by_sample for individual results
    main_result <- results_by_sample[[1]]
    
    result <- list(
      input = list(
        data = data_clean,
        n = nrow(data_clean),
        n_excluded = n_excluded,
        factors = factors,
        value_col = value
      ),
      design = design,
      variance_components = main_result$variance_components,
      precision = main_result$precision,
      anova_table = main_result$anova_table,
      by_sample = results_by_sample,
      sample_means = overall_means,
      settings = list(
        conf_level = conf_level,
        ci_method = ci_method,
        boot_n = if (ci_method == "bootstrap") boot_n else NA,
        method = method
      ),
      call = call
    )
    
  } else {
    # Single sample analysis
    computed <- .compute_precision_single(
      data = data_clean,
      factors = factors,
      design = design,
      conf_level = conf_level,
      ci_method = ci_method,
      boot_n = boot_n,
      method = method,
      value_col = value
    )
    
    result <- list(
      input = list(
        data = data_clean,
        n = nrow(data_clean),
        n_excluded = n_excluded,
        factors = factors,
        value_col = value
      ),
      design = design,
      variance_components = computed$variance_components,
      precision = computed$precision,
      anova_table = computed$anova_table,
      by_sample = NULL,
      sample_means = NULL,
      settings = list(
        conf_level = conf_level,
        ci_method = ci_method,
        boot_n = if (ci_method == "bootstrap") boot_n else NA,
        method = method
      ),
      call = call
    )
  }
  
  class(result) <- c("precision_study", "valytics_precision", "valytics_result")
  result
}


# Input Validation ----

#' Validate precision study input
#' @noRd
#' @keywords internal
.validate_precision_input <- function(data, value, sample, site, day, run,
                                      replicate, conf_level, boot_n) {
  
  # Check data is a data frame
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  
  if (nrow(data) == 0) {
    stop("`data` cannot be empty.", call. = FALSE)
  }
  
  # Check value column exists and is numeric
  if (!value %in% names(data)) {
    stop(sprintf("Column '%s' not found in data.", value), call. = FALSE)
  }
  
  if (!is.numeric(data[[value]])) {
    stop(sprintf("Column '%s' must be numeric.", value), call. = FALSE)
  }
  
  # Build factors list (only non-NULL factors)
  factors <- list()
  
  # Check each factor column
  factor_specs <- list(
    sample = sample,
    site = site,
    day = day,
    run = run,
    replicate = replicate
  )
  
  for (fname in names(factor_specs)) {
    fcol <- factor_specs[[fname]]
    if (!is.null(fcol)) {
      if (!fcol %in% names(data)) {
        stop(sprintf("Column '%s' (specified for %s) not found in data.",
                     fcol, fname), call. = FALSE)
      }
      factors[[fname]] <- fcol
    }
  }
  
  # Must have at least day factor
  if (is.null(factors$day)) {
    stop("At least 'day' factor must be specified.", call. = FALSE)
  }
  
  # Validate conf_level
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number between 0 and 1.", call. = FALSE)
  }
  
  # Validate boot_n
  if (!is.numeric(boot_n) || length(boot_n) != 1 ||
      boot_n < 100 || boot_n != floor(boot_n)) {
    stop("`boot_n` must be an integer >= 100.", call. = FALSE)
  }
  
  # Handle missing values
  # Identify columns to check for NA
  cols_to_check <- c(value, unlist(factors))
  complete_rows <- complete.cases(data[, cols_to_check, drop = FALSE])
  n_excluded <- sum(!complete_rows)
  
  if (n_excluded > 0) {
    message(sprintf("Note: %d observations excluded due to missing values.",
                    n_excluded))
  }
  
  data_clean <- data[complete_rows, , drop = FALSE]
  
  if (nrow(data_clean) < 3) {
    stop("At least 3 complete observations are required.", call. = FALSE)
  }
  
  list(
    data = data_clean,
    factors = factors,
    n_excluded = n_excluded
  )
}

# Design Detection ----

#' Detect experimental design from data structure
#' @noRd
#' @keywords internal
.detect_precision_design <- function(data, factors) {
  
  # Determine design type based on factors present
  has_site <- !is.null(factors$site)
  has_day <- !is.null(factors$day)
  has_run <- !is.null(factors$run)
  has_replicate <- !is.null(factors$replicate)
  
  # Design type
  if (has_site) {
    design_type <- "multi_site"
  } else {
    design_type <- "single_site"
  }
  
  # Build structure string (nesting hierarchy)
  structure_parts <- c()
  if (has_site) structure_parts <- c(structure_parts, "site")
  if (has_day) structure_parts <- c(structure_parts, "day")
  if (has_run) structure_parts <- c(structure_parts, "run")
  if (has_replicate) structure_parts <- c(structure_parts, "replicate")
  
  # If no explicit replicate column, replicates are inferred
  if (!has_replicate) {
    structure_parts <- c(structure_parts, "replicate (inferred)")
  }
  
  structure_string <- paste(structure_parts, collapse = "/")
  
  # Count levels for each factor
  levels_list <- list()
  for (fname in names(factors)) {
    if (fname != "sample") {  # Sample is not part of nesting
      fcol <- factors[[fname]]
      levels_list[[fname]] <- length(unique(data[[fcol]]))
    }
  }
  
  # Determine number of replicates (inferred if not explicit)
  if (!has_replicate) {
    # Count observations per lowest grouping
    grouping_cols <- c()
    if (has_site) grouping_cols <- c(grouping_cols, factors$site)
    if (has_day) grouping_cols <- c(grouping_cols, factors$day)
    if (has_run) grouping_cols <- c(grouping_cols, factors$run)
    
    if (length(grouping_cols) > 0) {
      counts <- aggregate(
        rep(1, nrow(data)),
        by = data[, grouping_cols, drop = FALSE],
        FUN = length
      )
      n_reps <- unique(counts$x)
      levels_list$replicate <- if (length(n_reps) == 1) n_reps else "varies"
    } else {
      levels_list$replicate <- nrow(data)
    }
  }
  
  # Check if design is balanced
  balanced <- .check_balance(data, factors)
  
  # Number of samples (concentration levels)
  if (!is.null(factors$sample)) {
    n_samples <- length(unique(data[[factors$sample]]))
  } else {
    n_samples <- 1
  }
  
  # Describe design in human-readable format
  design_desc <- .describe_design(levels_list, design_type)
  
  list(
    type = design_type,
    structure = structure_string,
    levels = levels_list,
    balanced = balanced,
    n_samples = n_samples,
    description = design_desc
  )
}


#' Check if design is balanced
#' @noRd
#' @keywords internal
.check_balance <- function(data, factors) {
  
  # Get nesting factors (exclude sample)
  nesting_factors <- factors[!names(factors) %in% c("sample", "replicate")]
  
  if (length(nesting_factors) == 0) {
    return(TRUE)  # No nesting = balanced by default
  }
  
  # Check counts at each level of nesting
  # For a balanced design, all groups should have equal counts
  
  # Start from the highest level and check down
  factor_names <- names(nesting_factors)
  
  for (i in seq_along(factor_names)) {
    # Group by factors up to this level
    group_cols <- unlist(nesting_factors[1:i])
    
    counts <- aggregate(
      rep(1, nrow(data)),
      by = data[, group_cols, drop = FALSE],
      FUN = length
    )
    
    # Check if all counts are equal
    if (length(unique(counts$x)) > 1) {
      return(FALSE)
    }
  }
  
  TRUE
}


#' Generate human-readable design description
#' @noRd
#' @keywords internal
.describe_design <- function(levels_list, design_type) {
  
  parts <- c()
  
  if ("site" %in% names(levels_list)) {
    parts <- c(parts, sprintf("%d sites", levels_list$site))
  }
  
  if ("day" %in% names(levels_list)) {
    parts <- c(parts, sprintf("%d days", levels_list$day))
  }
  
  if ("run" %in% names(levels_list)) {
    parts <- c(parts, sprintf("%d runs/day", levels_list$run))
  }
  
  if ("replicate" %in% names(levels_list)) {
    rep_val <- levels_list$replicate
    if (is.numeric(rep_val)) {
      parts <- c(parts, sprintf("%d replicates", rep_val))
    } else {
      parts <- c(parts, "varying replicates")
    }
  }
  
  if (length(parts) == 0) {
    return("Unknown design")
  }
  
  paste(parts, collapse = " x ")
}


#' Prepare factors for analysis (convert to proper factor type)
#' @noRd
#' @keywords internal
.prepare_factors <- function(data, factors) {
  
  for (fname in names(factors)) {
    fcol <- factors[[fname]]
    if (!is.factor(data[[fcol]])) {
      data[[fcol]] <- as.factor(data[[fcol]])
    }
  }
  
  data
}


# Core Computation ----

#' Compute precision for a single sample
#' @noRd
#' @keywords internal
.compute_precision_single <- function(data, factors, design, conf_level,
                                      ci_method, boot_n, method, value_col) {
  
  # Get the mean for CV calculation
  grand_mean <- mean(data[[value_col]], na.rm = TRUE)
  
  # Dispatch to appropriate estimation method
  if (method == "anova") {
    vc_result <- .estimate_vc_anova(data, factors, value_col)
  } else {
    vc_result <- .estimate_vc_reml(data, factors, value_col)
  }
  
  # Calculate confidence intervals
  ci_result <- .calculate_precision_ci(
    vc_result = vc_result,
    conf_level = conf_level,
    ci_method = ci_method,
    boot_n = boot_n,
    data = data,
    factors = factors,
    value_col = value_col,
    method = method
  )
  
  # Build precision summary
  precision <- .build_precision_summary(
    vc_result = vc_result,
    ci_result = ci_result,
    grand_mean = grand_mean,
    factors = factors
  )
  
  list(
    variance_components = vc_result$variance_components,
    precision = precision,
    anova_table = vc_result$anova_table,
    grand_mean = grand_mean
  )
}


# ANOVA Variance Component Estimation ----

#' Estimate variance components using ANOVA method
#'
#' Uses nested ANOVA (Type I SS) to estimate variance components via
#' method of moments. Supports hierarchies: site/day/run/replicate.
#'
#' @noRd
#' @keywords internal
.estimate_vc_anova <- function(data, factors, value_col) {
  
  n <- nrow(data)
  grand_mean <- mean(data[[value_col]], na.rm = TRUE)
  
  # Determine which factors are present
  
  has_site <- !is.null(factors$site)
  has_day <- !is.null(factors$day)
  
  has_run <- !is.null(factors$run)
  
  # Build the appropriate ANOVA based on available factors
  if (has_site && has_day && has_run) {
    # Full model: site/day/run/replicate
    result <- .anova_site_day_run(data, factors, value_col)
  } else if (has_site && has_day && !has_run) {
    # site/day/replicate (no run)
    result <- .anova_site_day(data, factors, value_col)
  } else if (!has_site && has_day && has_run) {
    # day/run/replicate (single site)
    result <- .anova_day_run(data, factors, value_col)
  } else if (!has_site && has_day && !has_run) {
    # day/replicate only (simplest case)
    result <- .anova_day_only(data, factors, value_col)
  } else {
    stop("Unsupported factor combination.", call. = FALSE)
  }
  
  result$grand_mean <- grand_mean
  result
}


#' ANOVA for day-only design (day/replicate)
#'
#' Model: y_ij = mu + D_i + e_ij
#' @noRd
#' @keywords internal
.anova_day_only <- function(data, factors, value_col) {
  
  y <- data[[value_col]]
  day <- data[[factors$day]]
  
  n <- length(y)
  n_days <- length(unique(day))
  
  # Calculate group sizes
  n_per_day <- as.numeric(table(day))
  
  # Grand mean
  grand_mean <- mean(y)
  
  # Day means
  day_means <- tapply(y, day, mean)
  
  # Sum of Squares
  # SS_total = sum((y - grand_mean)^2)
  ss_total <- sum((y - grand_mean)^2)
  
  # SS_day = sum(n_i * (day_mean_i - grand_mean)^2)
  ss_day <- sum(n_per_day * (day_means - grand_mean)^2)
  
  # SS_error = SS_total - SS_day
  ss_error <- ss_total - ss_day
  
  # Degrees of freedom
  df_day <- n_days - 1
  df_error <- n - n_days
  df_total <- n - 1
  
  # Mean Squares
  ms_day <- ss_day / df_day
  ms_error <- ss_error / df_error
  
  # Variance components (Method of Moments)
  # E[MS_day] = sigma^2_error + n_0 * sigma^2_day
  # E[MS_error] = sigma^2_error
  # where n_0 is the harmonic-like mean of group sizes for unbalanced designs
  
  # For balanced: n_0 = n_per_day (all equal)
  # For unbalanced: n_0 = (n - sum(n_i^2)/n) / (a - 1)
  if (length(unique(n_per_day)) == 1) {
    n_0 <- n_per_day[1]
  } else {
    n_0 <- (n - sum(n_per_day^2) / n) / (n_days - 1)
  }
  
  var_error <- ms_error
  var_day <- (ms_day - ms_error) / n_0
  
  # Handle negative variance estimates (set to 0)
  var_day <- max(0, var_day)
  
  var_total <- var_day + var_error
  
  # Build variance components table
  variance_components <- data.frame(
    component = c("between_day", "error", "total"),
    variance = c(var_day, var_error, var_total),
    sd = c(sqrt(var_day), sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_day / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_day, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  # Build ANOVA table
  anova_table <- data.frame(
    source = c("day", "error", "total"),
    df = c(df_day, df_error, df_total),
    ss = c(ss_day, ss_error, ss_total),
    ms = c(ms_day, ms_error, NA_real_),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = anova_table,
    n_0 = n_0
  )
}


#' ANOVA for day/run design (day/run/replicate)
#'
#' Model: y_ijk = mu + D_i + R_j(i) + e_ijk
#' Run is nested within day.
#' @noRd
#' @keywords internal
.anova_day_run <- function(data, factors, value_col) {
  
  y <- data[[value_col]]
  day <- data[[factors$day]]
  run <- data[[factors$run]]
  
  n <- length(y)
  n_days <- length(unique(day))
  
  # Create day:run interaction for nested structure
  day_run <- interaction(day, run, drop = TRUE)
  n_cells <- length(unique(day_run))
  
  # Grand mean
  grand_mean <- mean(y)
  
  # Day means
  day_means <- tapply(y, day, mean)
  n_per_day <- as.numeric(table(day))
  
  # Cell (day:run) means
  cell_means <- tapply(y, day_run, mean)
  n_per_cell <- as.numeric(table(day_run))
  
  # Sum of Squares
  ss_total <- sum((y - grand_mean)^2)
  
  # SS_day
  ss_day <- sum(n_per_day * (day_means - grand_mean)^2)
  
  # SS_run(day) = SS_cells - SS_day
  # where SS_cells = sum(n_ij * (cell_mean_ij - grand_mean)^2)
  ss_cells <- sum(n_per_cell * (cell_means - grand_mean)^2)
  ss_run <- ss_cells - ss_day
  
  # SS_error
  ss_error <- ss_total - ss_cells
  
  # Degrees of freedom
  df_day <- n_days - 1
  df_run <- n_cells - n_days  # runs nested in days
  
  df_error <- n - n_cells
  df_total <- n - 1
  
  # Mean Squares
  ms_day <- ss_day / df_day
  ms_run <- if (df_run > 0) ss_run / df_run else 0
  ms_error <- if (df_error > 0) ss_error / df_error else 0
  
  # Expected Mean Squares coefficients for unbalanced designs
  # For balanced: straightforward
  # For unbalanced: use synthesis coefficients
  
  # Calculate n_0 coefficients
  # These depend on the design balance
  design_info <- .get_design_coefficients_day_run(data, factors)
  
  # Variance components (Method of Moments)
  # E[MS_error] = sigma^2_e
  # E[MS_run] = sigma^2_e + n_r * sigma^2_run
  # E[MS_day] = sigma^2_e + n_r * sigma^2_run + n_d * sigma^2_day
  
  var_error <- ms_error
  var_run <- if (design_info$n_r > 0) (ms_run - ms_error) / design_info$n_r else 0
  var_day <- if (design_info$n_d > 0) {
    (ms_day - ms_run) / design_info$n_d
  } else {
    0
  }
  
  # Handle negative variance estimates
  var_run <- max(0, var_run)
  var_day <- max(0, var_day)
  
  var_total <- var_day + var_run + var_error
  
  # Build variance components table
  variance_components <- data.frame(
    component = c("between_day", "between_run", "error", "total"),
    variance = c(var_day, var_run, var_error, var_total),
    sd = c(sqrt(var_day), sqrt(var_run), sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_day / var_total,
      100 * var_run / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_day, df_run, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  # Build ANOVA table
  anova_table <- data.frame(
    source = c("day", "run(day)", "error", "total"),
    df = c(df_day, df_run, df_error, df_total),
    ss = c(ss_day, ss_run, ss_error, ss_total),
    ms = c(ms_day, ms_run, ms_error, NA_real_),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = anova_table,
    coefficients = design_info
  )
}


#' Calculate EMS coefficients for day/run design
#' @noRd
#' @keywords internal
.get_design_coefficients_day_run <- function(data, factors) {
  
  day <- data[[factors$day]]
  run <- data[[factors$run]]
  day_run <- interaction(day, run, drop = TRUE)
  
  n <- nrow(data)
  n_days <- length(unique(day))
  n_cells <- length(unique(day_run))
  
  # Group sizes
  n_per_day <- as.numeric(table(day))
  n_per_cell <- as.numeric(table(day_run))
  
  # Runs per day
  runs_per_day <- tapply(run, day, function(x) length(unique(x)))
  
  # For balanced design: n_r = replicates per cell, n_d = reps * runs
  if (length(unique(n_per_cell)) == 1 && length(unique(runs_per_day)) == 1) {
    # Balanced
    n_r <- n_per_cell[1]  # replicates per run
    n_d <- n_r * runs_per_day[1]  # replicates per day
  } else {
    # Unbalanced - use harmonic-like means
    # n_r for run effect
    n_r <- (n - sum(n_per_cell^2) / n) / (n_cells - n_days)
    if (!is.finite(n_r) || n_r <= 0) n_r <- mean(n_per_cell)
    
    # n_d for day effect
    # More complex for unbalanced - approximate
    n_d <- mean(n_per_day)
  }
  
  list(n_r = n_r, n_d = n_d)
}


#' ANOVA for site/day design (no run factor)
#'
#' Model: y_ijk = mu + S_i + D_j(i) + e_ijk
#' Day is nested within site.
#' @noRd
#' @keywords internal
.anova_site_day <- function(data, factors, value_col) {
  
  y <- data[[value_col]]
  site <- data[[factors$site]]
  day <- data[[factors$day]]
  
  n <- length(y)
  n_sites <- length(unique(site))
  
  # Create site:day interaction for nested structure
  site_day <- interaction(site, day, drop = TRUE)
  n_cells <- length(unique(site_day))
  
  # Grand mean
  grand_mean <- mean(y)
  
  # Site means
  site_means <- tapply(y, site, mean)
  n_per_site <- as.numeric(table(site))
  
  # Cell (site:day) means
  cell_means <- tapply(y, site_day, mean)
  n_per_cell <- as.numeric(table(site_day))
  
  # Sum of Squares
  ss_total <- sum((y - grand_mean)^2)
  ss_site <- sum(n_per_site * (site_means - grand_mean)^2)
  ss_cells <- sum(n_per_cell * (cell_means - grand_mean)^2)
  ss_day <- ss_cells - ss_site
  ss_error <- ss_total - ss_cells
  
  # Degrees of freedom
  df_site <- n_sites - 1
  df_day <- n_cells - n_sites
  df_error <- n - n_cells
  df_total <- n - 1
  
  # Mean Squares
  ms_site <- ss_site / df_site
  ms_day <- if (df_day > 0) ss_day / df_day else 0
  ms_error <- if (df_error > 0) ss_error / df_error else 0
  
  # Design coefficients
  design_info <- .get_design_coefficients_site_day(data, factors)
  
  # Variance components
  var_error <- ms_error
  var_day <- if (design_info$n_d > 0) (ms_day - ms_error) / design_info$n_d else 0
  var_site <- if (design_info$n_s > 0) (ms_site - ms_day) / design_info$n_s else 0
  
  var_day <- max(0, var_day)
  var_site <- max(0, var_site)
  
  var_total <- var_site + var_day + var_error
  
  variance_components <- data.frame(
    component = c("between_site", "between_day", "error", "total"),
    variance = c(var_site, var_day, var_error, var_total),
    sd = c(sqrt(var_site), sqrt(var_day), sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_site / var_total,
      100 * var_day / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_site, df_day, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  anova_table <- data.frame(
    source = c("site", "day(site)", "error", "total"),
    df = c(df_site, df_day, df_error, df_total),
    ss = c(ss_site, ss_day, ss_error, ss_total),
    ms = c(ms_site, ms_day, ms_error, NA_real_),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = anova_table,
    coefficients = design_info
  )
}


#' Calculate EMS coefficients for site/day design
#' @noRd
#' @keywords internal
.get_design_coefficients_site_day <- function(data, factors) {
  
  site <- data[[factors$site]]
  day <- data[[factors$day]]
  site_day <- interaction(site, day, drop = TRUE)
  
  n <- nrow(data)
  n_sites <- length(unique(site))
  n_cells <- length(unique(site_day))
  
  n_per_site <- as.numeric(table(site))
  n_per_cell <- as.numeric(table(site_day))
  days_per_site <- tapply(day, site, function(x) length(unique(x)))
  
  if (length(unique(n_per_cell)) == 1 && length(unique(days_per_site)) == 1) {
    n_d <- n_per_cell[1]
    n_s <- n_d * days_per_site[1]
  } else {
    n_d <- (n - sum(n_per_cell^2) / n) / (n_cells - n_sites)
    if (!is.finite(n_d) || n_d <= 0) n_d <- mean(n_per_cell)
    n_s <- mean(n_per_site)
  }
  
  list(n_d = n_d, n_s = n_s)
}


#' ANOVA for full site/day/run design
#'
#' Model: y_ijkl = mu + S_i + D_j(i) + R_k(ij) + e_l(ijk)
#' @noRd
#' @keywords internal
.anova_site_day_run <- function(data, factors, value_col) {
  
  y <- data[[value_col]]
  site <- data[[factors$site]]
  day <- data[[factors$day]]
  run <- data[[factors$run]]
  
  n <- length(y)
  n_sites <- length(unique(site))
  
  # Create nested interaction terms
  site_day <- interaction(site, day, drop = TRUE)
  site_day_run <- interaction(site, day, run, drop = TRUE)
  
  n_site_days <- length(unique(site_day))
  n_cells <- length(unique(site_day_run))
  
  # Grand mean
  grand_mean <- mean(y)
  
  # Means at each level
  site_means <- tapply(y, site, mean)
  site_day_means <- tapply(y, site_day, mean)
  cell_means <- tapply(y, site_day_run, mean)
  
  # Group sizes
  n_per_site <- as.numeric(table(site))
  n_per_site_day <- as.numeric(table(site_day))
  n_per_cell <- as.numeric(table(site_day_run))
  
  # Sum of Squares
  ss_total <- sum((y - grand_mean)^2)
  ss_site <- sum(n_per_site * (site_means - grand_mean)^2)
  ss_site_day <- sum(n_per_site_day * (site_day_means - grand_mean)^2)
  ss_cells <- sum(n_per_cell * (cell_means - grand_mean)^2)
  
  ss_day <- ss_site_day - ss_site
  ss_run <- ss_cells - ss_site_day
  ss_error <- ss_total - ss_cells
  
  # Degrees of freedom
  df_site <- n_sites - 1
  df_day <- n_site_days - n_sites
  df_run <- n_cells - n_site_days
  df_error <- n - n_cells
  df_total <- n - 1
  
  # Mean Squares
  ms_site <- ss_site / df_site
  ms_day <- if (df_day > 0) ss_day / df_day else 0
  ms_run <- if (df_run > 0) ss_run / df_run else 0
  ms_error <- if (df_error > 0) ss_error / df_error else 0
  
  # Design coefficients (simplified for balanced designs)
  design_info <- .get_design_coefficients_full(data, factors)
  
  # Variance components
  var_error <- ms_error
  var_run <- if (design_info$n_r > 0) (ms_run - ms_error) / design_info$n_r else 0
  var_day <- if (design_info$n_d > 0) (ms_day - ms_run) / design_info$n_d else 0
  var_site <- if (design_info$n_s > 0) (ms_site - ms_day) / design_info$n_s else 0
  
  var_run <- max(0, var_run)
  var_day <- max(0, var_day)
  var_site <- max(0, var_site)
  
  var_total <- var_site + var_day + var_run + var_error
  
  variance_components <- data.frame(
    component = c("between_site", "between_day", "between_run", "error", "total"),
    variance = c(var_site, var_day, var_run, var_error, var_total),
    sd = c(sqrt(var_site), sqrt(var_day), sqrt(var_run),
           sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_site / var_total,
      100 * var_day / var_total,
      100 * var_run / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_site, df_day, df_run, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  anova_table <- data.frame(
    source = c("site", "day(site)", "run(site:day)", "error", "total"),
    df = c(df_site, df_day, df_run, df_error, df_total),
    ss = c(ss_site, ss_day, ss_run, ss_error, ss_total),
    ms = c(ms_site, ms_day, ms_run, ms_error, NA_real_),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = anova_table,
    coefficients = design_info
  )
}


#' Calculate EMS coefficients for full site/day/run design
#' @noRd
#' @keywords internal
.get_design_coefficients_full <- function(data, factors) {
  
  site <- data[[factors$site]]
  day <- data[[factors$day]]
  run <- data[[factors$run]]
  
  site_day <- interaction(site, day, drop = TRUE)
  site_day_run <- interaction(site, day, run, drop = TRUE)
  
  n <- nrow(data)
  n_sites <- length(unique(site))
  n_site_days <- length(unique(site_day))
  n_cells <- length(unique(site_day_run))
  
  n_per_site <- as.numeric(table(site))
  n_per_site_day <- as.numeric(table(site_day))
  n_per_cell <- as.numeric(table(site_day_run))
  
  # Check if balanced
  balanced <- length(unique(n_per_cell)) == 1 &&
    length(unique(n_per_site_day)) == 1 &&
    length(unique(n_per_site)) == 1
  
  if (balanced) {
    n_r <- n_per_cell[1]
    runs_per_day <- n_cells / n_site_days
    n_d <- n_r * runs_per_day
    days_per_site <- n_site_days / n_sites
    n_s <- n_d * days_per_site
  } else {
    # Approximate for unbalanced
    n_r <- mean(n_per_cell)
    n_d <- mean(n_per_site_day)
    n_s <- mean(n_per_site)
  }
  
  list(n_r = n_r, n_d = n_d, n_s = n_s)
}


# REML Estimation ----

#' Estimate variance components using REML
#' @noRd
#' @keywords internal
.estimate_vc_reml <- function(data, factors, value_col) {
  
  # REML estimation using lme4::lmer()
  # Requires lme4 package (checked in main function)
  
  n <- nrow(data)
  grand_mean <- mean(data[[value_col]], na.rm = TRUE)
  
  # Determine which factors are present
  has_site <- !is.null(factors$site)
  has_day <- !is.null(factors$day)
  has_run <- !is.null(factors$run)
  
  # Prepare data with factors
  model_data <- data.frame(
    y = data[[value_col]]
  )
  
  # Add factors as proper factor variables
  if (has_day) {
    model_data$day <- factor(data[[factors$day]])
  }
  if (has_run) {
    # For nested run within day, create unique run identifier
    model_data$run <- factor(interaction(data[[factors$day]], 
                                         data[[factors$run]], 
                                         drop = TRUE))
  }
  if (has_site) {
    model_data$site <- factor(data[[factors$site]])
  }
  
  # Build formula and fit model based on available factors
  if (has_site && has_day && has_run) {
    result <- .reml_site_day_run(model_data, n)
  } else if (has_site && has_day && !has_run) {
    result <- .reml_site_day(model_data, n)
  } else if (!has_site && has_day && has_run) {
    result <- .reml_day_run(model_data, n)
  } else if (!has_site && has_day && !has_run) {
    result <- .reml_day_only(model_data, n)
  } else {
    stop("Unsupported factor combination for REML.", call. = FALSE)
  }
  
  result$grand_mean <- grand_mean
  result
}


#' REML estimation for day-only design
#' @noRd
#' @keywords internal
.reml_day_only <- function(model_data, n) {
  
  # Model: y ~ 1 + (1|day)
  fit <- lme4::lmer(y ~ 1 + (1 | day), data = model_data, REML = TRUE)
  
  # Extract variance components
  vc <- lme4::VarCorr(fit)
  var_day <- as.numeric(vc$day)
  var_error <- attr(vc, "sc")^2
  var_total <- var_day + var_error
  
  # Degrees of freedom (approximate for REML)
  n_days <- length(unique(model_data$day))
  df_day <- n_days - 1
  df_error <- n - n_days
  df_total <- n - 1
  
  # Build variance components table
  variance_components <- data.frame(
    component = c("between_day", "error", "total"),
    variance = c(var_day, var_error, var_total),
    sd = c(sqrt(var_day), sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_day / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_day, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = NULL,  # REML doesn't produce traditional ANOVA table
    model = fit,
    method = "reml"
  )
}


#' REML estimation for day/run design
#' @noRd
#' @keywords internal
.reml_day_run <- function(model_data, n) {
  
  # Model: y ~ 1 + (1|day) + (1|run)
  # run is already coded as day:run interaction (unique identifier)
  fit <- lme4::lmer(y ~ 1 + (1 | day) + (1 | run), data = model_data, REML = TRUE)
  
  # Extract variance components
  vc <- lme4::VarCorr(fit)
  var_day <- as.numeric(vc$day)
  var_run <- as.numeric(vc$run)
  var_error <- attr(vc, "sc")^2
  var_total <- var_day + var_run + var_error
  
  # Degrees of freedom (approximate for REML)
  n_days <- length(unique(model_data$day))
  n_runs <- length(unique(model_data$run))
  df_day <- n_days - 1
  df_run <- n_runs - n_days
  df_error <- n - n_runs
  df_total <- n - 1
  
  # Build variance components table
  variance_components <- data.frame(
    component = c("between_day", "between_run", "error", "total"),
    variance = c(var_day, var_run, var_error, var_total),
    sd = c(sqrt(var_day), sqrt(var_run), sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_day / var_total,
      100 * var_run / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_day, df_run, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = NULL,
    model = fit,
    method = "reml"
  )
}


#' REML estimation for site/day design
#' @noRd
#' @keywords internal
.reml_site_day <- function(model_data, n) {
  
  # Model: y ~ 1 + (1|site) + (1|site:day)
  # Create nested day within site
  model_data$site_day <- interaction(model_data$site, model_data$day, drop = TRUE)
  
  fit <- lme4::lmer(y ~ 1 + (1 | site) + (1 | site_day), data = model_data, REML = TRUE)
  
  # Extract variance components
  vc <- lme4::VarCorr(fit)
  var_site <- as.numeric(vc$site)
  var_day <- as.numeric(vc$site_day)
  var_error <- attr(vc, "sc")^2
  var_total <- var_site + var_day + var_error
  
  # Degrees of freedom (approximate)
  n_sites <- length(unique(model_data$site))
  n_site_days <- length(unique(model_data$site_day))
  df_site <- n_sites - 1
  df_day <- n_site_days - n_sites
  df_error <- n - n_site_days
  df_total <- n - 1
  
  # Build variance components table
  variance_components <- data.frame(
    component = c("between_site", "between_day", "error", "total"),
    variance = c(var_site, var_day, var_error, var_total),
    sd = c(sqrt(var_site), sqrt(var_day), sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_site / var_total,
      100 * var_day / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_site, df_day, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = NULL,
    model = fit,
    method = "reml"
  )
}


#' REML estimation for site/day/run design
#' @noRd
#' @keywords internal
.reml_site_day_run <- function(model_data, n) {
  
  # Model: y ~ 1 + (1|site) + (1|site:day) + (1|site:day:run)
  # Create nested factors
  model_data$site_day <- interaction(model_data$site, model_data$day, drop = TRUE)
  # run is already site:day:run (created in parent function)
  
  fit <- lme4::lmer(y ~ 1 + (1 | site) + (1 | site_day) + (1 | run), 
                    data = model_data, REML = TRUE)
  
  # Extract variance components
  vc <- lme4::VarCorr(fit)
  var_site <- as.numeric(vc$site)
  var_day <- as.numeric(vc$site_day)
  var_run <- as.numeric(vc$run)
  var_error <- attr(vc, "sc")^2
  var_total <- var_site + var_day + var_run + var_error
  
  # Degrees of freedom (approximate)
  n_sites <- length(unique(model_data$site))
  n_site_days <- length(unique(model_data$site_day))
  n_runs <- length(unique(model_data$run))
  df_site <- n_sites - 1
  df_day <- n_site_days - n_sites
  df_run <- n_runs - n_site_days
  df_error <- n - n_runs
  df_total <- n - 1
  
  # Build variance components table
  variance_components <- data.frame(
    component = c("between_site", "between_day", "between_run", "error", "total"),
    variance = c(var_site, var_day, var_run, var_error, var_total),
    sd = c(sqrt(var_site), sqrt(var_day), sqrt(var_run), sqrt(var_error), sqrt(var_total)),
    pct_total = c(
      100 * var_site / var_total,
      100 * var_day / var_total,
      100 * var_run / var_total,
      100 * var_error / var_total,
      100
    ),
    df = c(df_site, df_day, df_run, df_error, df_total),
    stringsAsFactors = FALSE
  )
  
  list(
    variance_components = variance_components,
    anova_table = NULL,
    model = fit,
    method = "reml"
  )
}


# Confidence Intervals ----

#' Calculate confidence intervals for precision estimates
#' @noRd
#' @keywords internal
#' @param factors Factor column mapping
#' @param value_col Name of value column
#' @param method Estimation method ("anova" or "reml")
#'
#' @return List with CI for each variance component and precision measure
#' @noRd
#' @keywords internal
.calculate_precision_ci <- function(vc_result, conf_level, ci_method, boot_n,
                                    data, factors, value_col, method) {
  
  if (ci_method == "satterthwaite") {
    ci_result <- .ci_satterthwaite(vc_result, conf_level)
  } else if (ci_method == "mls") {
    ci_result <- .ci_mls(vc_result, conf_level)
  } else if (ci_method == "bootstrap") {
    ci_result <- .ci_bootstrap(vc_result, conf_level, boot_n,
                               data, factors, value_col, method)
  } else {
    stop("Unknown CI method: ", ci_method, call. = FALSE)
  }
  
  ci_result
}


# Satterthwaite CI ----

#' Confidence intervals using Satterthwaite approximation
#'
#' For a single variance component sigma^2 with df degrees of freedom:
#'   CI = [df * sigma^2 / chi^2_{alpha/2}, df * sigma^2 / chi^2_{1-alpha/2}]
#'
#' For a linear combination L = sum(sigma^2_i):
#'   df_L = L^2 / sum(sigma^4_i / df_i)
#'   Then use chi-square CI with df_L
#'
#' @noRd
#' @keywords internal
.ci_satterthwaite <- function(vc_result, conf_level) {
  
  alpha <- 1 - conf_level
  vc <- vc_result$variance_components
  
  # Identify which components are present
  has_site <- "between_site" %in% vc$component
  has_run <- "between_run" %in% vc$component
  has_day <- "between_day" %in% vc$component
  
  # Extract variance estimates and degrees of freedom
  var_error <- vc$variance[vc$component == "error"]
  df_error <- vc$df[vc$component == "error"]
  
  var_day <- if (has_day) vc$variance[vc$component == "between_day"] else 0
  df_day <- if (has_day) vc$df[vc$component == "between_day"] else 0
  
  var_run <- if (has_run) vc$variance[vc$component == "between_run"] else 0
  df_run <- if (has_run) vc$df[vc$component == "between_run"] else 0
  
  var_site <- if (has_site) vc$variance[vc$component == "between_site"] else 0
  df_site <- if (has_site) vc$df[vc$component == "between_site"] else 0
  
  # CI for repeatability (error variance only)
  repeatability_ci <- .ci_single_variance(var_error, df_error, alpha)
  
  # CI for between-day variance
  between_day_ci <- if (has_day && df_day > 0) {
    .ci_single_variance(var_day, df_day, alpha)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  
  # CI for between-run variance
  between_run_ci <- if (has_run && df_run > 0) {
    .ci_single_variance(var_run, df_run, alpha)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  
  # CI for between-site variance
  between_site_ci <- if (has_site && df_site > 0) {
    .ci_single_variance(var_site, df_site, alpha)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  
  # CI for intermediate precision (within-lab)
  if (has_run) {
    intermediate_components <- c(var_day, var_run, var_error)
    intermediate_dfs <- c(df_day, df_run, df_error)
  } else {
    intermediate_components <- c(var_day, var_error)
    intermediate_dfs <- c(df_day, df_error)
  }
  intermediate_ci <- .ci_variance_sum(intermediate_components, intermediate_dfs, alpha)
  
  # CI for reproducibility (total, including site)
  if (has_site) {
    repro_components <- c(var_site, var_day, var_run, var_error)
    repro_dfs <- c(df_site, df_day, df_run, df_error)
    reproducibility_ci <- .ci_variance_sum(repro_components, repro_dfs, alpha)
  } else {
    reproducibility_ci <- intermediate_ci
  }
  
  list(
    repeatability = repeatability_ci,
    between_day = between_day_ci,
    between_run = between_run_ci,
    between_site = between_site_ci,
    intermediate = intermediate_ci,
    reproducibility = reproducibility_ci,
    method = "satterthwaite"
  )
}


#' CI for a single variance component using chi-square
#'
#' For a variance component with point estimate sigma^2 and df degrees of
#' freedom, the CI is: [df * sigma^2 / chi^2_{1-alpha/2}, df * sigma^2 / chi^2_{alpha/2}]
#'
#' When variance estimate is 0 (constrained from negative), we return
#' [0, 0] since the point estimate is on the boundary.
#'
#' @noRd
#' @keywords internal
.ci_single_variance <- function(variance, df, alpha) {
  
  if (df <= 0 || !is.finite(df) || !is.finite(variance)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  # When variance is 0 (constrained), return [0, 0]
  # This is a boundary estimate - the true variance could be 0 or small positive
  if (variance <= 0) {
    return(c(lower = 0, upper = 0))
  }
  
  chi_lower <- stats::qchisq(1 - alpha / 2, df)
  chi_upper <- stats::qchisq(alpha / 2, df)
  
  ci_lower <- df * variance / chi_lower
  ci_upper <- df * variance / chi_upper
  
  c(lower = ci_lower, upper = ci_upper)
}


#' CI for a sum of variance components using Satterthwaite approximation
#'
#' @param variances Vector of variance component estimates
#' @param dfs Vector of degrees of freedom for each component
#' @param alpha Significance level (1 - conf_level)
#'
#' @noRd
#' @keywords internal
.ci_variance_sum <- function(variances, dfs, alpha) {
  
  valid <- dfs > 0 & is.finite(variances) & is.finite(dfs)
  
  if (sum(valid) == 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  variances <- variances[valid]
  dfs <- dfs[valid]
  
  L <- sum(variances)
  
  # When sum of variances is 0 (constrained), return [0, 0]
  if (L <= 0) {
    return(c(lower = 0, upper = 0))
  }
  
  # Satterthwaite approximation: df_L = L^2 / sum(sigma^4_i / df_i)
  denominator <- sum(variances^2 / dfs)
  
  if (denominator <= 0 || !is.finite(denominator)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  df_satt <- L^2 / denominator
  
  chi_lower <- stats::qchisq(1 - alpha / 2, df_satt)
  chi_upper <- stats::qchisq(alpha / 2, df_satt)
  
  ci_lower <- df_satt * L / chi_lower
  ci_upper <- df_satt * L / chi_upper
  
  c(lower = ci_lower, upper = ci_upper)
}


# MLS CI ----

#' Confidence intervals using Modified Large Sample method
#'
#' MLS provides better coverage when variance components may be near zero.
#'
#' @noRd
#' @keywords internal
.ci_mls <- function(vc_result, conf_level) {
  
  alpha <- 1 - conf_level
  vc <- vc_result$variance_components
  anova <- vc_result$anova_table
  
  if (is.null(anova)) {
    return(.ci_satterthwaite(vc_result, conf_level))
  }
  
  has_site <- "between_site" %in% vc$component
  has_run <- "between_run" %in% vc$component
  has_day <- "between_day" %in% vc$component
  
  var_error <- vc$variance[vc$component == "error"]
  df_error <- vc$df[vc$component == "error"]
  
  var_day <- if (has_day) vc$variance[vc$component == "between_day"] else 0
  df_day <- if (has_day) vc$df[vc$component == "between_day"] else 0
  
  var_run <- if (has_run) vc$variance[vc$component == "between_run"] else 0
  df_run <- if (has_run) vc$df[vc$component == "between_run"] else 0
  
  var_site <- if (has_site) vc$variance[vc$component == "between_site"] else 0
  df_site <- if (has_site) vc$df[vc$component == "between_site"] else 0
  
  repeatability_ci <- .ci_single_variance(var_error, df_error, alpha)
  
  between_day_ci <- if (has_day && df_day > 0) {
    .ci_mls_single(var_day, df_day, alpha)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  
  between_run_ci <- if (has_run && df_run > 0) {
    .ci_mls_single(var_run, df_run, alpha)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  
  between_site_ci <- if (has_site && df_site > 0) {
    .ci_mls_single(var_site, df_site, alpha)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  
  if (has_run) {
    intermediate_components <- c(var_day, var_run, var_error)
    intermediate_dfs <- c(df_day, df_run, df_error)
  } else {
    intermediate_components <- c(var_day, var_error)
    intermediate_dfs <- c(df_day, df_error)
  }
  intermediate_ci <- .ci_mls_sum(intermediate_components, intermediate_dfs, alpha)
  
  if (has_site) {
    repro_components <- c(var_site, var_day, var_run, var_error)
    repro_dfs <- c(df_site, df_day, df_run, df_error)
    reproducibility_ci <- .ci_mls_sum(repro_components, repro_dfs, alpha)
  } else {
    reproducibility_ci <- intermediate_ci
  }
  
  list(
    repeatability = repeatability_ci,
    between_day = between_day_ci,
    between_run = between_run_ci,
    between_site = between_site_ci,
    intermediate = intermediate_ci,
    reproducibility = reproducibility_ci,
    method = "mls"
  )
}


#' MLS CI for a single variance component
#'
#' @noRd
#' @keywords internal
.ci_mls_single <- function(variance, df, alpha) {
  
  if (df <= 0 || !is.finite(df) || !is.finite(variance)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  # When variance is 0 (constrained), return [0, 0]
  if (variance <= 0) {
    return(c(lower = 0, upper = 0))
  }
  
  G1 <- 1 - stats::qchisq(alpha / 2, df) / df
  G2 <- stats::qchisq(1 - alpha / 2, df) / df - 1
  
  H1 <- (G1 * variance)^2
  H2 <- (G2 * variance)^2
  
  ci_lower <- max(0, variance - sqrt(H1))
  ci_upper <- variance + sqrt(H2)
  
  c(lower = ci_lower, upper = ci_upper)
}


#' MLS CI for a sum of variance components
#'
#' @noRd
#' @keywords internal
.ci_mls_sum <- function(variances, dfs, alpha) {
  
  valid <- dfs > 0 & is.finite(variances) & is.finite(dfs)
  
  if (sum(valid) == 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  variances <- variances[valid]
  dfs <- dfs[valid]
  
  L <- sum(variances)
  
  # When sum of variances is 0 (constrained), return [0, 0]
  if (L <= 0) {
    return(c(lower = 0, upper = 0))
  }
  
  G1 <- 1 - stats::qchisq(alpha / 2, dfs) / dfs
  G2 <- stats::qchisq(1 - alpha / 2, dfs) / dfs - 1
  
  H1 <- sum((G1 * variances)^2)
  H2 <- sum((G2 * variances)^2)
  
  ci_lower <- max(0, L - sqrt(H1))
  ci_upper <- L + sqrt(H2)
  
  c(lower = ci_lower, upper = ci_upper)
}


# Bootstrap CI ----

#' Confidence intervals using BCa bootstrap
#'
#' Resamples the data preserving nested structure and re-estimates
#' variance components.
#'
#' @noRd
#' @keywords internal
.ci_bootstrap <- function(vc_result, conf_level, boot_n, data, factors,
                          value_col, method) {
  
  alpha <- 1 - conf_level
  
  has_site <- !is.null(factors$site)
  has_day <- !is.null(factors$day)
  
  if (has_site) {
    resample_col <- factors$site
  } else if (has_day) {
    resample_col <- factors$day
  } else {
    resample_col <- factors$day
  }
  
  units <- unique(data[[resample_col]])
  n_units <- length(units)
  
  boot_repeatability <- numeric(boot_n)
  boot_intermediate <- numeric(boot_n)
  boot_reproducibility <- numeric(boot_n)
  boot_between_day <- numeric(boot_n)
  boot_between_run <- numeric(boot_n)
  boot_between_site <- numeric(boot_n)
  
  for (b in seq_len(boot_n)) {
    boot_units <- sample(units, n_units, replace = TRUE)
    
    boot_data <- do.call(rbind, lapply(seq_along(boot_units), function(i) {
      unit_data <- data[data[[resample_col]] == boot_units[i], , drop = FALSE]
      unit_data[[resample_col]] <- paste0(unit_data[[resample_col]], "_", i)
      unit_data
    }))
    
    boot_data <- .prepare_factors(boot_data, factors)
    
    tryCatch({
      if (method == "anova") {
        vc_boot <- .estimate_vc_anova(boot_data, factors, value_col)
      } else {
        vc_boot <- .estimate_vc_reml(boot_data, factors, value_col)
      }
      
      vc_b <- vc_boot$variance_components
      
      var_error_b <- vc_b$variance[vc_b$component == "error"]
      var_day_b <- if ("between_day" %in% vc_b$component) {
        vc_b$variance[vc_b$component == "between_day"]
      } else 0
      var_run_b <- if ("between_run" %in% vc_b$component) {
        vc_b$variance[vc_b$component == "between_run"]
      } else 0
      var_site_b <- if ("between_site" %in% vc_b$component) {
        vc_b$variance[vc_b$component == "between_site"]
      } else 0
      
      boot_repeatability[b] <- var_error_b
      boot_between_day[b] <- var_day_b
      boot_between_run[b] <- var_run_b
      boot_between_site[b] <- var_site_b
      boot_intermediate[b] <- var_day_b + var_run_b + var_error_b
      boot_reproducibility[b] <- var_site_b + var_day_b + var_run_b + var_error_b
      
    }, error = function(e) {
      boot_repeatability[b] <- NA_real_
      boot_intermediate[b] <- NA_real_
      boot_reproducibility[b] <- NA_real_
      boot_between_day[b] <- NA_real_
      boot_between_run[b] <- NA_real_
      boot_between_site[b] <- NA_real_
    })
  }
  
  vc <- vc_result$variance_components
  
  var_error_orig <- vc$variance[vc$component == "error"]
  var_day_orig <- if ("between_day" %in% vc$component) {
    vc$variance[vc$component == "between_day"]
  } else 0
  var_run_orig <- if ("between_run" %in% vc$component) {
    vc$variance[vc$component == "between_run"]
  } else 0
  var_site_orig <- if ("between_site" %in% vc$component) {
    vc$variance[vc$component == "between_site"]
  } else 0
  
  repeatability_ci <- .bca_precision_ci(boot_repeatability, var_error_orig, alpha)
  between_day_ci <- .bca_precision_ci(boot_between_day, var_day_orig, alpha)
  between_run_ci <- .bca_precision_ci(boot_between_run, var_run_orig, alpha)
  between_site_ci <- .bca_precision_ci(boot_between_site, var_site_orig, alpha)
  
  intermediate_orig <- var_day_orig + var_run_orig + var_error_orig
  intermediate_ci <- .bca_precision_ci(boot_intermediate, intermediate_orig, alpha)
  
  reproducibility_orig <- var_site_orig + var_day_orig + var_run_orig + var_error_orig
  reproducibility_ci <- .bca_precision_ci(boot_reproducibility, reproducibility_orig, alpha)
  
  list(
    repeatability = repeatability_ci,
    between_day = between_day_ci,
    between_run = between_run_ci,
    between_site = between_site_ci,
    intermediate = intermediate_ci,
    reproducibility = reproducibility_ci,
    method = "bootstrap",
    boot_n = boot_n
  )
}


#' BCa confidence interval for precision
#'
#' @param boot_stat Vector of bootstrap statistics
#' @param orig_stat Original point estimate
#' @param alpha Significance level (1 - conf_level)
#'
#' @noRd
#' @keywords internal
.bca_precision_ci <- function(boot_stat, orig_stat, alpha) {
  
  boot_stat <- boot_stat[is.finite(boot_stat)]
  
  if (length(boot_stat) < 100) {
    warning("Too few valid bootstrap samples for reliable CI.", call. = FALSE)
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  prop_less <- mean(boot_stat < orig_stat)
  prop_less <- max(0.001, min(0.999, prop_less))
  z0 <- stats::qnorm(prop_less)
  
  a <- 0  # Simplified (no jackknife acceleration)
  
  z_alpha_lower <- stats::qnorm(alpha / 2)
  z_alpha_upper <- stats::qnorm(1 - alpha / 2)
  
  alpha1 <- stats::pnorm(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)))
  alpha2 <- stats::pnorm(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)))
  
  alpha1 <- max(0.001, min(0.999, alpha1))
  alpha2 <- max(0.001, min(0.999, alpha2))
  
  ci <- stats::quantile(boot_stat, probs = c(alpha1, alpha2), na.rm = TRUE)
  
  ci_lower <- max(0, ci[[1]])
  ci_upper <- ci[[2]]
  
  c(lower = ci_lower, upper = ci_upper)
}


# Precision Summary ----

#' Build precision summary data frame with confidence intervals
#' @noRd
#' @keywords internal
.build_precision_summary <- function(vc_result, ci_result, grand_mean, factors) {
  
  vc <- vc_result$variance_components
  has_site <- "between_site" %in% vc$component
  has_run <- "between_run" %in% vc$component
  
  # Extract variances
  var_error <- vc$variance[vc$component == "error"]
  var_day <- vc$variance[vc$component == "between_day"]
  var_run <- if (has_run) vc$variance[vc$component == "between_run"] else 0
  var_site <- if (has_site) vc$variance[vc$component == "between_site"] else 0
  
  # Calculate composite precision measures (as variances)
  sd_repeatability <- sqrt(var_error)
  sd_intermediate <- sqrt(var_day + var_run + var_error)
  sd_reproducibility <- sqrt(var_site + var_day + var_run + var_error)
  
  # Build output vectors
  measures <- c()
  sds <- c()
  ci_lowers <- c()
  ci_uppers <- c()
  
  # Helper to safely extract CI and convert variance CI to SD CI
  get_sd_ci <- function(ci_name) {
    if (!is.null(ci_result[[ci_name]])) {
      ci_var <- ci_result[[ci_name]]
      lower <- if (is.finite(ci_var["lower"])) sqrt(max(0, ci_var["lower"])) else NA_real_
      upper <- if (is.finite(ci_var["upper"])) sqrt(ci_var["upper"]) else NA_real_
      names(lower) <- names(upper) <- NULL 
      return(c(lower = lower, upper = upper))
    }
    c(lower = NA_real_, upper = NA_real_)
  }
  
  # Repeatability (always present)
  measures <- c(measures, "Repeatability")
  sds <- c(sds, sd_repeatability)
  rep_ci <- get_sd_ci("repeatability")
  ci_lowers <- c(ci_lowers, rep_ci["lower"])
  ci_uppers <- c(ci_uppers, rep_ci["upper"])
  
  # Between-run (if present)
  if (has_run) {
    measures <- c(measures, "Between-run")
    sds <- c(sds, sqrt(var_run))
    run_ci <- get_sd_ci("between_run")
    ci_lowers <- c(ci_lowers, run_ci["lower"])
    ci_uppers <- c(ci_uppers, run_ci["upper"])
  }
  
  # Between-day
  measures <- c(measures, "Between-day")
  sds <- c(sds, sqrt(var_day))
  day_ci <- get_sd_ci("between_day")
  ci_lowers <- c(ci_lowers, day_ci["lower"])
  ci_uppers <- c(ci_uppers, day_ci["upper"])
  
  # Intermediate precision (within-lab)
  measures <- c(measures, "Intermediate precision")
  sds <- c(sds, sd_intermediate)
  int_ci <- get_sd_ci("intermediate")
  ci_lowers <- c(ci_lowers, int_ci["lower"])
  ci_uppers <- c(ci_uppers, int_ci["upper"])
  
  # Between-site and reproducibility (if multi-site)
  if (has_site) {
    measures <- c(measures, "Between-site")
    sds <- c(sds, sqrt(var_site))
    site_ci <- get_sd_ci("between_site")
    ci_lowers <- c(ci_lowers, site_ci["lower"])
    ci_uppers <- c(ci_uppers, site_ci["upper"])
    
    measures <- c(measures, "Reproducibility")
    sds <- c(sds, sd_reproducibility)
    repro_ci <- get_sd_ci("reproducibility")
    ci_lowers <- c(ci_lowers, repro_ci["lower"])
    ci_uppers <- c(ci_uppers, repro_ci["upper"])
  }
  
  # Calculate CVs (as percentage)
  cv_values <- 100 * sds / grand_mean
  cv_ci_lower <- 100 * ci_lowers / grand_mean
  cv_ci_upper <- 100 * ci_uppers / grand_mean
  
  data.frame(
    measure = measures,
    sd = sds,
    cv_pct = cv_values,
    ci_lower = ci_lowers,
    ci_upper = ci_uppers,
    cv_ci_lower = cv_ci_lower,
    cv_ci_upper = cv_ci_upper,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# S3 Methods ----

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