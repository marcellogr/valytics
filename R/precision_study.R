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


# Core Computation (Placeholder - to be implemented in Phase 1b) ----

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
    anova_table = vc_result$anova_table
  )
}


#' Estimate variance components using ANOVA method
#' @noRd
#' @keywords internal
.estimate_vc_anova <- function(data, factors, value_col) {

  # This will be fully implemented in Phase 1b
  # For now, return placeholder structure

  # Build formula for nested ANOVA
  # Hierarchy: site > day > run > replicate (error)

  # Placeholder - actual implementation coming in Phase 1b
  warning("ANOVA estimation not yet implemented. Returning placeholder.",
          call. = FALSE)

  # Return structure that matches expected output
  list(
    variance_components = data.frame(
      component = c("day", "run", "error", "total"),
      variance = c(NA_real_, NA_real_, NA_real_, NA_real_),
      sd = c(NA_real_, NA_real_, NA_real_, NA_real_),
      pct_total = c(NA_real_, NA_real_, NA_real_, 100),
      df = c(NA_real_, NA_real_, NA_real_, NA_real_),
      stringsAsFactors = FALSE
    ),
    anova_table = data.frame(
      source = character(0),
      df = numeric(0),
      ss = numeric(0),
      ms = numeric(0),
      stringsAsFactors = FALSE
    ),
    grand_mean = mean(data[[value_col]], na.rm = TRUE)
  )
}


#' Estimate variance components using REML
#' @noRd
#' @keywords internal
.estimate_vc_reml <- function(data, factors, value_col) {

  # This will be fully implemented in Phase 1d
  # Requires lme4 package

  warning("REML estimation not yet implemented. Returning placeholder.",
          call. = FALSE)

  # Return same structure as ANOVA
  list(
    variance_components = data.frame(
      component = c("day", "run", "error", "total"),
      variance = c(NA_real_, NA_real_, NA_real_, NA_real_),
      sd = c(NA_real_, NA_real_, NA_real_, NA_real_),
      pct_total = c(NA_real_, NA_real_, NA_real_, 100),
      df = c(NA_real_, NA_real_, NA_real_, NA_real_),
      stringsAsFactors = FALSE
    ),
    anova_table = NULL,
    grand_mean = mean(data[[value_col]], na.rm = TRUE)
  )
}


#' Calculate confidence intervals for precision estimates
#' @noRd
#' @keywords internal
.calculate_precision_ci <- function(vc_result, conf_level, ci_method, boot_n,
                                     data, factors, value_col, method) {

  # This will be fully implemented in Phase 1c
  # For now, return placeholder

  # Placeholder structure
  list(
    repeatability_ci = c(lower = NA_real_, upper = NA_real_),
    intermediate_ci = c(lower = NA_real_, upper = NA_real_),
    reproducibility_ci = c(lower = NA_real_, upper = NA_real_)
  )
}


#' Build precision summary data frame
#' @noRd
#' @keywords internal
.build_precision_summary <- function(vc_result, ci_result, grand_mean, factors) {

  # This will be refined in Phase 1b/1c
  # For now, return placeholder structure

  has_site <- !is.null(factors$site)
  has_run <- !is.null(factors$run)

  measures <- c("Repeatability")
  if (has_run) measures <- c(measures, "Between-run")
  measures <- c(measures, "Between-day", "Intermediate precision")
  if (has_site) measures <- c(measures, "Between-site", "Reproducibility")

  data.frame(
    measure = measures,
    sd = rep(NA_real_, length(measures)),
    cv_pct = rep(NA_real_, length(measures)),
    ci_lower = rep(NA_real_, length(measures)),
    ci_upper = rep(NA_real_, length(measures)),
    stringsAsFactors = FALSE
  )
}
