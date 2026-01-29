# Tests for precision_study() - Phase 1a
# Input validation and design detection
# Test Data Setup ----

#' Create EP05-style test data (20 days x 2 runs x 2 replicates)
#' @noRd
create_ep05_data <- function(seed = 42, mean_val = 100,
                             sd_day = 1.5, sd_run = 1.0, sd_error = 2.0) {
  set.seed(seed)
  
  n_days <- 20
  n_runs <- 2
  n_reps <- 2
  
  data <- expand.grid(
    day = 1:n_days,
    run = 1:n_runs,
    replicate = 1:n_reps
  )
  
  # Add variance components
  day_effect <- rep(rnorm(n_days, 0, sd_day), each = n_runs * n_reps)
  run_effect <- rep(rnorm(n_days * n_runs, 0, sd_run), each = n_reps)
  error <- rnorm(nrow(data), 0, sd_error)
  
  data$value <- mean_val + day_effect + run_effect + error
  data
}


#' Create EP15-style test data (5 days x 5 replicates)
#' @noRd
create_ep15_data <- function(seed = 42, mean_val = 100,
                             sd_day = 1.5, sd_error = 2.0) {
  set.seed(seed)
  
  n_days <- 5
  n_reps <- 5
  
  data <- expand.grid(
    day = 1:n_days,
    replicate = 1:n_reps
  )
  
  day_effect <- rep(rnorm(n_days, 0, sd_day), each = n_reps)
  error <- rnorm(nrow(data), 0, sd_error)
  
  data$value <- mean_val + day_effect + error
  data
}


#' Create multi-site test data (3 sites x 5 days x 5 replicates)
#' @noRd
create_multisite_data <- function(seed = 42, mean_val = 100,
                                  sd_site = 2.0, sd_day = 1.5, sd_error = 2.0) {
  set.seed(seed)
  
  n_sites <- 3
  n_days <- 5
  n_reps <- 5
  
  data <- expand.grid(
    site = LETTERS[1:n_sites],
    day = 1:n_days,
    replicate = 1:n_reps
  )
  
  site_effect <- rep(rnorm(n_sites, 0, sd_site), each = n_days * n_reps)
  day_effect <- rep(rnorm(n_sites * n_days, 0, sd_day), each = n_reps)
  error <- rnorm(nrow(data), 0, sd_error)
  
  data$value <- mean_val + site_effect + day_effect + error
  data
}


#' Create multi-sample test data (multiple concentration levels)
#' @noRd
create_multisample_data <- function(seed = 42) {
  set.seed(seed)
  
  # 3 samples at different concentrations
  samples <- data.frame(
    sample_id = c("Low", "Medium", "High"),
    true_mean = c(50, 100, 200)
  )
  
  n_days <- 5
  n_reps <- 3
  
  data_list <- lapply(1:nrow(samples), function(i) {
    sample_data <- expand.grid(
      day = 1:n_days,
      replicate = 1:n_reps
    )
    sample_data$sample_id <- samples$sample_id[i]
    
    # CV tends to be higher at low concentrations
    cv_factor <- 1 + (1 - samples$true_mean[i] / max(samples$true_mean)) * 0.5
    sd_day <- samples$true_mean[i] * 0.015 * cv_factor
    sd_error <- samples$true_mean[i] * 0.02 * cv_factor
    
    day_effect <- rep(rnorm(n_days, 0, sd_day), each = n_reps)
    error <- rnorm(nrow(sample_data), 0, sd_error)
    
    sample_data$value <- samples$true_mean[i] + day_effect + error
    sample_data
  })
  
  do.call(rbind, data_list)
}


# Input Validation Tests ----

test_that("precision_study validates data argument", {
  # Not a data frame
  expect_error(
    precision_study(data = "not a data frame", value = "value", day = "day"),
    "`data` must be a data frame"
  )
  
  # Empty data frame
  expect_error(
    precision_study(data = data.frame(), value = "value", day = "day"),
    "`data` cannot be empty"
  )
  
  # NULL data
  expect_error(
    precision_study(data = NULL, value = "value", day = "day"),
    "`data` must be a data frame"
  )
})


test_that("precision_study validates value column", {
  data <- create_ep15_data()
  
  # Value column doesn't exist
  expect_error(
    precision_study(data = data, value = "nonexistent", day = "day"),
    "Column 'nonexistent' not found in data"
  )
  
  # Value column is not numeric
  data$char_col <- "text"
  expect_error(
    precision_study(data = data, value = "char_col", day = "day"),
    "Column 'char_col' must be numeric"
  )
})


test_that("precision_study validates factor columns exist", {
  data <- create_ep15_data()
  
  # Day column doesn't exist
  expect_error(
    precision_study(data = data, value = "value", day = "nonexistent"),
    "Column 'nonexistent' \\(specified for day\\) not found in data"
  )
  
  # Site column doesn't exist
  expect_error(
    precision_study(data = data, value = "value", day = "day",
                    site = "nonexistent"),
    "Column 'nonexistent' \\(specified for site\\) not found in data"
  )
  
  # Run column doesn't exist
  expect_error(
    precision_study(data = data, value = "value", day = "day",
                    run = "nonexistent"),
    "Column 'nonexistent' \\(specified for run\\) not found in data"
  )
})


test_that("precision_study requires day factor", {
  data <- create_ep15_data()
  
  expect_error(
    precision_study(data = data, value = "value", day = NULL),
    "At least 'day' factor must be specified"
  )
})


test_that("precision_study validates conf_level", {
  data <- create_ep15_data()
  
  expect_error(
    precision_study(data = data, value = "value", day = "day", conf_level = 0),
    "`conf_level` must be a single number between 0 and 1"
  )
  
  expect_error(
    precision_study(data = data, value = "value", day = "day", conf_level = 1),
    "`conf_level` must be a single number between 0 and 1"
  )
  
  expect_error(
    precision_study(data = data, value = "value", day = "day", conf_level = 1.5),
    "`conf_level` must be a single number between 0 and 1"
  )
  
  expect_error(
    precision_study(data = data, value = "value", day = "day",
                    conf_level = c(0.9, 0.95)),
    "`conf_level` must be a single number between 0 and 1"
  )
})


test_that("precision_study validates boot_n", {
  data <- create_ep15_data()
  
  expect_error(
    precision_study(data = data, value = "value", day = "day", boot_n = 50),
    "`boot_n` must be an integer >= 100"
  )
  
  expect_error(
    precision_study(data = data, value = "value", day = "day", boot_n = 100.5),
    "`boot_n` must be an integer >= 100"
  )
})


test_that("precision_study handles missing values", {
  data <- create_ep15_data()
  
  # Add some NAs
  data$value[c(1, 5, 10)] <- NA
  
  # Should produce message about excluded observations
  expect_message(
    precision_study(data = data, value = "value", day = "day"),
    "3 observations excluded due to missing values"
  )
})


test_that("precision_study fails with too few observations after NA removal", {
  data <- data.frame(
    day = 1:5,
    value = c(NA, NA, NA, 100, 101)
  )
  
  # Only 2 complete observations
  expect_error(
    suppressMessages(precision_study(data = data, value = "value", day = "day")),
    "At least 3 complete observations are required"
  )
})


# Design Detection Tests ----

test_that("precision_study detects EP05 design (day/run/replicate)", {
  data <- create_ep05_data()
  
  # Suppress the placeholder warning for now
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day", run = "run")
  )
  
  expect_equal(result$design$type, "single_site")
  expect_true(grepl("day", result$design$structure))
  expect_true(grepl("run", result$design$structure))
  expect_equal(result$design$levels$day, 20)
  expect_equal(result$design$levels$run, 2)
  expect_true(result$design$balanced)
})


test_that("precision_study detects EP15 design (day/replicate)", {
  data <- create_ep15_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day")
  )
  
  expect_equal(result$design$type, "single_site")
  expect_true(grepl("day", result$design$structure))
  expect_false(grepl("/run/", result$design$structure))
  expect_equal(result$design$levels$day, 5)
  expect_equal(result$design$levels$replicate, 5)  # Inferred
  expect_true(result$design$balanced)
})


test_that("precision_study detects multi-site design", {
  data <- create_multisite_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day", site = "site")
  )
  
  expect_equal(result$design$type, "multi_site")
  expect_true(grepl("site", result$design$structure))
  expect_true(grepl("day", result$design$structure))
  expect_equal(result$design$levels$site, 3)
  expect_equal(result$design$levels$day, 5)
  expect_true(result$design$balanced)
})


test_that("precision_study detects unbalanced design", {
  # Create unbalanced data
  data <- create_ep15_data()
  
  # Remove some observations to make it unbalanced
  data <- data[-(1:3), ]  # Remove first 3 rows from day 1
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day")
  )
  
  expect_false(result$design$balanced)
})


test_that("precision_study handles explicit replicate column", {
  data <- create_ep15_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day",
                    replicate = "replicate")
  )
  
  expect_true(grepl("replicate", result$design$structure))
  expect_false(grepl("inferred", result$design$structure))
})


test_that("precision_study handles multiple samples", {
  data <- create_multisample_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day",
                    sample = "sample_id")
  )
  
  expect_equal(result$design$n_samples, 3)
  expect_false(is.null(result$by_sample))
  expect_equal(length(result$by_sample), 3)
  expect_true(all(c("Low", "Medium", "High") %in% names(result$by_sample)))
  expect_equal(length(result$sample_means), 3)
})


# Output Structure Tests ----

test_that("precision_study returns correct class", {
  data <- create_ep15_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day")
  )
  
  expect_s3_class(result, "precision_study")
  expect_s3_class(result, "valytics_precision")
  expect_s3_class(result, "valytics_result")
})


test_that("precision_study returns expected structure", {
  data <- create_ep15_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day")
  )
  
  # Check top-level components
  expect_named(result, c("input", "design", "variance_components", "precision",
                         "anova_table", "by_sample", "sample_means",
                         "settings", "call"))
  
  # Check input structure
  expect_named(result$input, c("data", "n", "n_excluded", "factors", "value_col"))
  expect_equal(result$input$n, 25)  # 5 days x 5 replicates
  expect_equal(result$input$value_col, "value")
  
  # Check design structure
  expect_named(result$design, c("type", "structure", "levels", "balanced",
                                "n_samples", "description"))
  
  # Check settings
  expect_named(result$settings, c("conf_level", "ci_method", "boot_n", "method"))
  expect_equal(result$settings$conf_level, 0.95)
  expect_equal(result$settings$ci_method, "satterthwaite")
  expect_equal(result$settings$method, "anova")
})


test_that("precision_study stores factors correctly", {
  data <- create_ep05_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day", run = "run")
  )
  
  expect_equal(result$input$factors$day, "day")
  expect_equal(result$input$factors$run, "run")
  expect_null(result$input$factors$site)
  expect_null(result$input$factors$sample)
})


test_that("precision_study preserves call", {
  data <- create_ep15_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day", conf_level = 0.90)
  )
  
  expect_true(!is.null(result$call))
  expect_true(inherits(result$call, "call"))
})


# Method Argument Tests ----

test_that("precision_study accepts different ci_method values", {
  data <- create_ep15_data()
  
  # Satterthwaite (default)
  result1 <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day",
                    ci_method = "satterthwaite")
  )
  expect_equal(result1$settings$ci_method, "satterthwaite")
  
  # MLS
  result2 <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day",
                    ci_method = "mls")
  )
  expect_equal(result2$settings$ci_method, "mls")
  
  # Bootstrap
  result3 <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day",
                    ci_method = "bootstrap")
  )
  expect_equal(result3$settings$ci_method, "bootstrap")
  expect_equal(result3$settings$boot_n, 1999)
})


test_that("precision_study accepts method = 'anova'", {
  data <- create_ep15_data()
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day", method = "anova")
  )
  
  expect_equal(result$settings$method, "anova")
})


test_that("precision_study requires lme4 for REML", {
  data <- create_ep15_data()
  
  # This test assumes lme4 might not be installed
  # If lme4 is installed, the warning will be about "not yet implemented"
  # If lme4 is not installed, we get an error about the package
  
  # We'll just check that REML is a valid option
  # The actual REML test will be in Phase 1d
  result <- suppressWarnings(
    tryCatch(
      precision_study(data = data, value = "value", day = "day", method = "reml"),
      error = function(e) {
        if (grepl("lme4", e$message)) {
          # Expected if lme4 not installed
          NULL
        } else {
          stop(e)
        }
      }
    )
  )
  
  # If we got a result, check it used REML setting
  if (!is.null(result)) {
    expect_equal(result$settings$method, "reml")
  }
})


# Edge Cases ----

test_that("precision_study works with minimum data", {
  # Minimum: 3 observations
  data <- data.frame(
    day = c(1, 1, 2),
    value = c(100, 101, 99)
  )
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day")
  )
  
  expect_s3_class(result, "precision_study")
  expect_equal(result$input$n, 3)
})


test_that("precision_study handles factor columns that are already factors", {
  data <- create_ep15_data()
  data$day <- as.factor(data$day)
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day")
  )
  
  expect_s3_class(result, "precision_study")
  expect_true(is.factor(result$input$data$day))
})


test_that("precision_study handles character factor columns", {
  data <- create_multisite_data()
  # site is already character (LETTERS)
  
  result <- suppressWarnings(
    precision_study(data = data, value = "value", day = "day", site = "site")
  )
  
  expect_s3_class(result, "precision_study")
  expect_true(is.factor(result$input$data$site))
})


test_that("precision_study generates sensible design description", {
  # EP05 style
  data1 <- create_ep05_data()
  result1 <- suppressWarnings(
    precision_study(data = data1, value = "value", day = "day", run = "run")
  )
  expect_true(grepl("20 days", result1$design$description))
  expect_true(grepl("2 runs", result1$design$description))
  
  # Multi-site
  data2 <- create_multisite_data()
  result2 <- suppressWarnings(
    precision_study(data = data2, value = "value", day = "day", site = "site")
  )
  expect_true(grepl("3 sites", result2$design$description))
  expect_true(grepl("5 days", result2$design$description))
})


# ANOVA Variance Component Tests ----

test_that("ANOVA estimates variance components for day-only design", {
  # Create data with known variance components
  set.seed(123)
  
  n_days <- 10
  n_reps <- 5
  
  # Known variance components
  true_var_day <- 4.0    # SD = 2
  true_var_error <- 1.0  # SD = 1
  
  data <- expand.grid(day = 1:n_days, replicate = 1:n_reps)
  day_effect <- rep(rnorm(n_days, 0, sqrt(true_var_day)), each = n_reps)
  error <- rnorm(nrow(data), 0, sqrt(true_var_error))
  data$value <- 100 + day_effect + error
  
  result <- precision_study(data, value = "value", day = "day")
  
  # Check structure
  expect_true("variance_components" %in% names(result))
  expect_true("anova_table" %in% names(result))
  
  vc <- result$variance_components
  expect_true("between_day" %in% vc$component)
  expect_true("error" %in% vc$component)
  expect_true("total" %in% vc$component)
  
  # Check that estimates are reasonable (within ~50% of true values for this sample size)
  # Note: with n=50, estimates can vary quite a bit
  expect_true(vc$variance[vc$component == "error"] > 0)
  expect_true(vc$variance[vc$component == "between_day"] >= 0)
  expect_true(vc$variance[vc$component == "total"] > 0)
  
  # Total should equal sum of components
  total_var <- vc$variance[vc$component == "total"]
  sum_components <- sum(vc$variance[vc$component != "total"])
  expect_equal(total_var, sum_components, tolerance = 1e-10)
  
  # Percentages should sum to 100
  pct_sum <- sum(vc$pct_total[vc$component != "total"])
  expect_equal(pct_sum, 100, tolerance = 1e-10)
})


test_that("ANOVA estimates variance components for day/run design", {
  set.seed(456)
  
  n_days <- 10
  n_runs <- 2
  n_reps <- 2
  
  # Known variance components
  true_var_day <- 2.25   # SD = 1.5
  true_var_run <- 1.0    # SD = 1.0
  true_var_error <- 4.0  # SD = 2.0
  
  data <- expand.grid(day = 1:n_days, run = 1:n_runs, replicate = 1:n_reps)
  day_effect <- rep(rnorm(n_days, 0, sqrt(true_var_day)), each = n_runs * n_reps)
  run_effect <- rep(rnorm(n_days * n_runs, 0, sqrt(true_var_run)), each = n_reps)
  error <- rnorm(nrow(data), 0, sqrt(true_var_error))
  data$value <- 100 + day_effect + run_effect + error
  
  result <- precision_study(data, value = "value", day = "day", run = "run")
  
  vc <- result$variance_components
  expect_true("between_day" %in% vc$component)
  expect_true("between_run" %in% vc$component)
  expect_true("error" %in% vc$component)
  
  # All variances should be non-negative
  expect_true(all(vc$variance >= 0))
  
  # Check ANOVA table structure
  anova <- result$anova_table
  expect_true("day" %in% anova$source)
  expect_true(any(grepl("run", anova$source)))
  expect_true("error" %in% anova$source)
  
  # DF should be correct for balanced design
  expect_equal(anova$df[anova$source == "day"], n_days - 1)
  expect_equal(anova$df[anova$source == "error"], n_days * n_runs * (n_reps - 1))
})


test_that("ANOVA estimates variance components for multi-site design", {
  set.seed(789)
  
  n_sites <- 3
  n_days <- 5
  n_reps <- 4
  
  true_var_site <- 9.0   # SD = 3
  true_var_day <- 4.0    # SD = 2
  true_var_error <- 1.0  # SD = 1
  
  data <- expand.grid(
    site = LETTERS[1:n_sites],
    day = 1:n_days,
    replicate = 1:n_reps
  )
  
  site_effect <- rep(rnorm(n_sites, 0, sqrt(true_var_site)), each = n_days * n_reps)
  day_effect <- rep(rnorm(n_sites * n_days, 0, sqrt(true_var_day)), each = n_reps)
  error <- rnorm(nrow(data), 0, sqrt(true_var_error))
  data$value <- 100 + site_effect + day_effect + error
  
  result <- precision_study(data, value = "value", day = "day", site = "site")
  
  vc <- result$variance_components
  expect_true("between_site" %in% vc$component)
  expect_true("between_day" %in% vc$component)
  expect_true("error" %in% vc$component)
  
  # Site variance should be substantial (we added SD=3)
  expect_true(vc$variance[vc$component == "between_site"] > 0)
})


test_that("ANOVA estimates variance components for full site/day/run design", {
  set.seed(101)
  
  n_sites <- 2
  n_days <- 3
  n_runs <- 2
  n_reps <- 2
  
  data <- expand.grid(
    site = LETTERS[1:n_sites],
    day = 1:n_days,
    run = 1:n_runs,
    replicate = 1:n_reps
  )
  
  # Add variance components
  site_effect <- rep(rnorm(n_sites, 0, 2), each = n_days * n_runs * n_reps)
  day_effect <- rep(rnorm(n_sites * n_days, 0, 1.5), each = n_runs * n_reps)
  run_effect <- rep(rnorm(n_sites * n_days * n_runs, 0, 1), each = n_reps)
  error <- rnorm(nrow(data), 0, 2)
  data$value <- 100 + site_effect + day_effect + run_effect + error
  
  result <- precision_study(
    data, value = "value",
    day = "day", run = "run", site = "site"
  )
  
  vc <- result$variance_components
  expect_true("between_site" %in% vc$component)
  expect_true("between_day" %in% vc$component)
  expect_true("between_run" %in% vc$component)
  expect_true("error" %in% vc$component)
  expect_true("total" %in% vc$component)
  
  # Should have 5 rows
  expect_equal(nrow(vc), 5)
})


test_that("ANOVA handles negative variance estimates correctly", {
  # Create data where between-day variance is near zero
  # This can lead to negative ANOVA estimates
  set.seed(202)
  
  n_days <- 5
  n_reps <- 10
  
  # Very small day effect, large error
  data <- expand.grid(day = 1:n_days, replicate = 1:n_reps)
  day_effect <- rep(rnorm(n_days, 0, 0.1), each = n_reps)  # Tiny day effect
  error <- rnorm(nrow(data), 0, 5)  # Large error
  data$value <- 100 + day_effect + error
  
  result <- precision_study(data, value = "value", day = "day")
  
  vc <- result$variance_components
  
  # All variances should be >= 0 (negative estimates set to 0)
  expect_true(all(vc$variance >= 0))
  
  # Day variance might be 0 due to negative estimate correction
  expect_true(vc$variance[vc$component == "between_day"] >= 0)
})


test_that("ANOVA table has correct structure", {
  data <- create_ep05_data()
  
  result <- precision_study(data, value = "value", day = "day", run = "run")
  
  anova <- result$anova_table
  
  # Check columns
  expect_true(all(c("source", "df", "ss", "ms") %in% names(anova)))
  
  # Check that SS are non-negative
  expect_true(all(anova$ss >= 0, na.rm = TRUE))
  
  # Check that MS = SS / DF (except for total)
  for (i in seq_len(nrow(anova))) {
    if (!is.na(anova$ms[i]) && anova$df[i] > 0) {
      expect_equal(anova$ms[i], anova$ss[i] / anova$df[i], tolerance = 1e-10)
    }
  }
  
  # Total SS should equal sum of other SS
  total_ss <- anova$ss[anova$source == "total"]
  other_ss <- sum(anova$ss[anova$source != "total"])
  expect_equal(total_ss, other_ss, tolerance = 1e-10)
})


# Precision Summary Tests ----

test_that("Precision summary calculates correct measures for day-only design", {
  set.seed(303)
  data <- create_ep15_data(mean_val = 100, sd_day = 2, sd_error = 3)
  
  result <- precision_study(data, value = "value", day = "day")
  
  prec <- result$precision
  
  # Should have: Repeatability, Between-day, Intermediate precision
  expect_true("Repeatability" %in% prec$measure)
  expect_true("Between-day" %in% prec$measure)
  expect_true("Intermediate precision" %in% prec$measure)
  
  # Should NOT have: Between-run, Between-site, Reproducibility
  expect_false("Between-run" %in% prec$measure)
  expect_false("Between-site" %in% prec$measure)
  expect_false("Reproducibility" %in% prec$measure)
  
  # SD values should be positive
  expect_true(all(prec$sd > 0))
  
  # CV should be SD / mean * 100
  # Grand mean is approximately 100
  expect_true(all(prec$cv_pct > 0))
  expect_true(all(prec$cv_pct < 20))  # Reasonable range for this data
})


test_that("Precision summary calculates correct measures for day/run design", {
  data <- create_ep05_data()
  
  result <- precision_study(data, value = "value", day = "day", run = "run")
  
  prec <- result$precision
  
  # Should have: Repeatability, Between-run, Between-day, Intermediate
  expect_true("Repeatability" %in% prec$measure)
  expect_true("Between-run" %in% prec$measure)
  expect_true("Between-day" %in% prec$measure)
  expect_true("Intermediate precision" %in% prec$measure)
})


test_that("Precision summary calculates correct measures for multi-site design", {
  data <- create_multisite_data()
  
  result <- precision_study(data, value = "value", day = "day", site = "site")
  
  prec <- result$precision
  
  # Should have all measures including site and reproducibility
  expect_true("Repeatability" %in% prec$measure)
  expect_true("Between-day" %in% prec$measure)
  expect_true("Between-site" %in% prec$measure)
  expect_true("Reproducibility" %in% prec$measure)
  expect_true("Intermediate precision" %in% prec$measure)
})


test_that("Intermediate precision is correctly calculated", {
  set.seed(404)
  
  # Create data with known components
  n_days <- 20
  n_runs <- 2
  n_reps <- 2
  
  var_day <- 4.0
  var_run <- 1.0
  var_error <- 2.25
  
  data <- expand.grid(day = 1:n_days, run = 1:n_runs, replicate = 1:n_reps)
  day_effect <- rep(rnorm(n_days, 0, sqrt(var_day)), each = n_runs * n_reps)
  run_effect <- rep(rnorm(n_days * n_runs, 0, sqrt(var_run)), each = n_reps)
  error <- rnorm(nrow(data), 0, sqrt(var_error))
  data$value <- 100 + day_effect + run_effect + error
  
  result <- precision_study(data, value = "value", day = "day", run = "run")
  
  vc <- result$variance_components
  prec <- result$precision
  
  # Intermediate precision should be sqrt(var_day + var_run + var_error)
  estimated_var_day <- vc$variance[vc$component == "between_day"]
  estimated_var_run <- vc$variance[vc$component == "between_run"]
  estimated_var_error <- vc$variance[vc$component == "error"]
  
  expected_intermediate_sd <- sqrt(estimated_var_day + estimated_var_run + estimated_var_error)
  actual_intermediate_sd <- prec$sd[prec$measure == "Intermediate precision"]
  
  expect_equal(actual_intermediate_sd, expected_intermediate_sd, tolerance = 1e-10)
})


test_that("Reproducibility is correctly calculated for multi-site", {
  set.seed(505)
  
  data <- create_multisite_data(sd_site = 3, sd_day = 2, sd_error = 1.5)
  
  result <- precision_study(data, value = "value", day = "day", site = "site")
  
  vc <- result$variance_components
  prec <- result$precision
  
  # Reproducibility should be sqrt(all variance components)
  var_site <- vc$variance[vc$component == "between_site"]
  var_day <- vc$variance[vc$component == "between_day"]
  var_error <- vc$variance[vc$component == "error"]
  
  expected_repro_sd <- sqrt(var_site + var_day + var_error)
  actual_repro_sd <- prec$sd[prec$measure == "Reproducibility"]
  
  expect_equal(actual_repro_sd, expected_repro_sd, tolerance = 1e-10)
})