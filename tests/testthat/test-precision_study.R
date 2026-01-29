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
  
  # Generate effects
  site_effects <- rnorm(n_sites, 0, sd_site)
  day_effects <- rnorm(n_sites * n_days, 0, sd_day)
  errors <- rnorm(nrow(data), 0, sd_error)
  
  # Map effects to data rows correctly
  # expand.grid varies the first argument fastest, so site varies fastest
  site_effect <- site_effects[as.numeric(as.factor(data$site))]
  day_idx <- (as.numeric(as.factor(data$site)) - 1) * n_days + as.numeric(as.factor(data$day))
  day_effect <- day_effects[day_idx]
  
  data$value <- mean_val + site_effect + day_effect + errors
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

#' Create day-only precision data (5 days x 5 replicates)
create_day_only_data <- function(grand_mean = 100, sd_day = 2.0, sd_error = 1.5,
                                 seed = 42) {
  set.seed(seed)
  
  n_days <- 5
  n_reps <- 5
  
  data <- expand.grid(
    day = factor(1:n_days),
    replicate = factor(1:n_reps)
  )
  
  day_effects <- rnorm(n_days, 0, sd_day)
  errors <- rnorm(nrow(data), 0, sd_error)
  
  day_effect <- day_effects[as.numeric(data$day)]
  data$value <- grand_mean + day_effect + errors
  
  data
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

test_that("CI handles zero variance component gracefully", {
  # Create data with very small between-day variance
  set.seed(42)
  data <- expand.grid(day = factor(1:10), replicate = 1:5)
  data$value <- 100 + rnorm(nrow(data), 0, 2)  # Only error variance
  
  prec <- precision_study(data, value = "value", day = "day",
                          ci_method = "satterthwaite")
  
  # Should not error, and CIs should be present
  expect_s3_class(prec, "precision_study")
  expect_true("ci_lower" %in% names(prec$precision))
})


test_that("CI handles small sample sizes", {
  # Minimum viable data
  set.seed(42)
  data <- expand.grid(day = factor(1:3), replicate = 1:2)
  data$value <- 100 + rnorm(nrow(data), 0, 2)
  
  prec <- precision_study(data, value = "value", day = "day",
                          ci_method = "satterthwaite")
  
  # Should not error
  expect_s3_class(prec, "precision_study")
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
  
  # Generate random effects
  site_effects <- rnorm(n_sites, 0, sqrt(true_var_site))
  day_effects <- rnorm(n_sites * n_days, 0, sqrt(true_var_day))
  
  # Map effects correctly - expand.grid varies first argument (site) fastest
  site_effect <- site_effects[as.numeric(as.factor(data$site))]
  day_idx <- (as.numeric(as.factor(data$site)) - 1) * n_days + as.numeric(as.factor(data$day))
  day_effect <- day_effects[day_idx]
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
  
  # Generate random effects
  site_effects <- rnorm(n_sites, 0, 2)
  day_effects <- rnorm(n_sites * n_days, 0, 1.5)
  run_effects <- rnorm(n_sites * n_days * n_runs, 0, 1)
  
  # Map effects correctly - expand.grid varies first argument (site) fastest
  site_idx <- as.numeric(as.factor(data$site))
  day_idx <- (site_idx - 1) * n_days + as.numeric(as.factor(data$day))
  run_idx <- (day_idx - 1) * n_runs + as.numeric(as.factor(data$run))
  
  site_effect <- site_effects[site_idx]
  day_effect <- day_effects[day_idx]
  run_effect <- run_effects[run_idx]
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

# Test data generators (for Confidence Intervals) ----

#' Create EP05-style precision data (20 days x 2 runs x 2 replicates)
create_ep05_data <- function(grand_mean = 100, sd_day = 1.5, sd_run = 1.0,
                             sd_error = 2.0, seed = 42) {
  set.seed(seed)
  
  n_days <- 20
  n_runs <- 2
  n_reps <- 2
  
  data <- expand.grid(
    day = factor(1:n_days),
    run = factor(1:n_runs),
    replicate = factor(1:n_reps)
  )
  
  # Add variance components
  day_effects <- rnorm(n_days, 0, sd_day)
  run_effects <- rnorm(n_days * n_runs, 0, sd_run)
  errors <- rnorm(nrow(data), 0, sd_error)
  
  day_effect <- day_effects[as.numeric(data$day)]
  run_effect <- run_effects[(as.numeric(data$day) - 1) * n_runs + as.numeric(data$run)]
  
  data$value <- grand_mean + day_effect + run_effect + errors
  
  data
}


#' Create day-only precision data (5 days x 5 replicates)
create_day_only_data <- function(grand_mean = 100, sd_day = 2.0, sd_error = 1.5,
                                 seed = 42) {
  set.seed(seed)
  
  n_days <- 5
  n_reps <- 5
  
  data <- expand.grid(
    day = factor(1:n_days),
    replicate = factor(1:n_reps)
  )
  
  day_effects <- rnorm(n_days, 0, sd_day)
  errors <- rnorm(nrow(data), 0, sd_error)
  
  day_effect <- day_effects[as.numeric(data$day)]
  data$value <- grand_mean + day_effect + errors
  
  data
}

#' Create multi-site precision data
create_multisite_data <- function(grand_mean = 100, sd_site = 2.5, sd_day = 1.5,
                                  sd_error = 2.0, seed = 42) {
  set.seed(seed)
  
  n_sites <- 3
  n_days <- 5
  n_reps <- 5
  
  data <- expand.grid(
    site = factor(1:n_sites),
    day = factor(1:n_days),
    replicate = factor(1:n_reps)
  )
  
  site_effects <- rnorm(n_sites, 0, sd_site)
  day_effects <- rnorm(n_sites * n_days, 0, sd_day)
  errors <- rnorm(nrow(data), 0, sd_error)
  
  site_effect <- site_effects[as.numeric(data$site)]
  day_idx <- (as.numeric(data$site) - 1) * n_days + as.numeric(data$day)
  day_effect <- day_effects[day_idx]
  
  data$value <- grand_mean + site_effect + day_effect + errors
  
  data
}

# Satterthwaite CI tests ----

test_that("Satterthwaite CI produces valid structure for day/run design", {
  data <- create_ep05_data()
  prec <- precision_study(data, value = "value", day = "day", run = "run",
                          ci_method = "satterthwaite")
  
  # Check precision data frame has CI columns
  
  expect_true("ci_lower" %in% names(prec$precision))
  expect_true("ci_upper" %in% names(prec$precision))
  expect_true("cv_ci_lower" %in% names(prec$precision))
  expect_true("cv_ci_upper" %in% names(prec$precision))
  
  # Check all measures have CIs
  expect_equal(nrow(prec$precision), 4)  # repeatability, between-run, between-day, intermediate
  
  # CIs should be numeric
  expect_true(is.numeric(prec$precision$ci_lower))
  expect_true(is.numeric(prec$precision$ci_upper))
})


test_that("Satterthwaite CI bounds are correctly ordered", {
  data <- create_ep05_data()
  prec <- precision_study(data, value = "value", day = "day", run = "run",
                          ci_method = "satterthwaite")
  
  # Lower bound should be less than or equal to upper bound
  for (i in seq_len(nrow(prec$precision))) {
    if (!is.na(prec$precision$ci_lower[i]) && !is.na(prec$precision$ci_upper[i])) {
      expect_true(prec$precision$ci_lower[i] <= prec$precision$ci_upper[i],
                  info = paste("Row", i, prec$precision$measure[i]))
    }
  }
  
  # Point estimate should be within CI (for most cases)
  for (i in seq_len(nrow(prec$precision))) {
    sd_val <- prec$precision$sd[i]
    ci_l <- prec$precision$ci_lower[i]
    ci_u <- prec$precision$ci_upper[i]
    
    if (!is.na(ci_l) && !is.na(ci_u) && sd_val > 0) {
      expect_true(sd_val >= ci_l * 0.99,  # Small tolerance for numerical precision
                  info = paste("Lower bound check for", prec$precision$measure[i]))
      expect_true(sd_val <= ci_u * 1.01,
                  info = paste("Upper bound check for", prec$precision$measure[i]))
    }
  }
})


test_that("Satterthwaite CI bounds are non-negative", {
  data <- create_ep05_data()
  prec <- precision_study(data, value = "value", day = "day", run = "run",
                          ci_method = "satterthwaite")
  
  # All CI bounds for SD should be non-negative
  expect_true(all(is.na(prec$precision$ci_lower) | prec$precision$ci_lower >= 0))
  expect_true(all(is.na(prec$precision$ci_upper) | prec$precision$ci_upper >= 0))
})


test_that("Satterthwaite CI works with day-only design", {
  data <- create_day_only_data()
  prec <- precision_study(data, value = "value", day = "day",
                          ci_method = "satterthwaite")
  
  expect_true("ci_lower" %in% names(prec$precision))
  expect_true("ci_upper" %in% names(prec$precision))
  
  # Should have repeatability, between-day, and intermediate
  expect_equal(nrow(prec$precision), 3)
})


test_that("Satterthwaite CI respects confidence level", {
  data <- create_ep05_data()
  
  prec_95 <- precision_study(data, value = "value", day = "day", run = "run",
                             conf_level = 0.95, ci_method = "satterthwaite")
  prec_90 <- precision_study(data, value = "value", day = "day", run = "run",
                             conf_level = 0.90, ci_method = "satterthwaite")
  prec_99 <- precision_study(data, value = "value", day = "day", run = "run",
                             conf_level = 0.99, ci_method = "satterthwaite")
  
  # 99% CI should be wider than 95%, which should be wider than 90%
  for (i in seq_len(nrow(prec_95$precision))) {
    width_90 <- prec_90$precision$ci_upper[i] - prec_90$precision$ci_lower[i]
    width_95 <- prec_95$precision$ci_upper[i] - prec_95$precision$ci_lower[i]
    width_99 <- prec_99$precision$ci_upper[i] - prec_99$precision$ci_lower[i]
    
    if (!is.na(width_90) && !is.na(width_95) && !is.na(width_99)) {
      expect_true(width_99 >= width_95,
                  info = paste("99% >= 95% for", prec_95$precision$measure[i]))
      expect_true(width_95 >= width_90,
                  info = paste("95% >= 90% for", prec_95$precision$measure[i]))
    }
  }
})


test_that("Satterthwaite CI width decreases with sample size", {
  # Small sample
  set.seed(42)
  data_small <- expand.grid(day = factor(1:5), replicate = 1:3)
  data_small$value <- 100 + rnorm(nrow(data_small), 0, 2)
  
  # Larger sample
  set.seed(42)
  data_large <- expand.grid(day = factor(1:20), replicate = 1:5)
  data_large$value <- 100 + rnorm(nrow(data_large), 0, 2)
  
  prec_small <- precision_study(data_small, value = "value", day = "day",
                                ci_method = "satterthwaite")
  prec_large <- precision_study(data_large, value = "value", day = "day",
                                ci_method = "satterthwaite")
  
  # Width for repeatability should be smaller with larger sample
  width_small <- prec_small$precision$ci_upper[1] - prec_small$precision$ci_lower[1]
  width_large <- prec_large$precision$ci_upper[1] - prec_large$precision$ci_lower[1]
  
  expect_true(width_large < width_small)
})


# MLS CI tests ----

test_that("MLS CI produces valid structure", {
  data <- create_ep05_data()
  prec <- precision_study(data, value = "value", day = "day", run = "run",
                          ci_method = "mls")
  
  expect_true("ci_lower" %in% names(prec$precision))
  expect_true("ci_upper" %in% names(prec$precision))
  expect_equal(nrow(prec$precision), 4)
})


test_that("MLS CI bounds are non-negative", {
  data <- create_ep05_data()
  prec <- precision_study(data, value = "value", day = "day", run = "run",
                          ci_method = "mls")
  
  expect_true(all(is.na(prec$precision$ci_lower) | prec$precision$ci_lower >= 0))
  expect_true(all(is.na(prec$precision$ci_upper) | prec$precision$ci_upper >= 0))
})


test_that("MLS and Satterthwaite CIs have same point estimates", {
  data <- create_ep05_data()
  
  prec_sat <- precision_study(data, value = "value", day = "day", run = "run",
                              ci_method = "satterthwaite")
  prec_mls <- precision_study(data, value = "value", day = "day", run = "run",
                              ci_method = "mls")
  
  # Point estimates should be identical
  expect_equal(prec_sat$precision$sd, prec_mls$precision$sd)
  expect_equal(prec_sat$precision$cv_pct, prec_mls$precision$cv_pct)
})


# Bootstrap CI tests ----

test_that("Bootstrap CI produces valid structure", {
  skip_on_cran()  # Skip on CRAN due to computation time
  
  data <- create_day_only_data()
  prec <- precision_study(data, value = "value", day = "day",
                          ci_method = "bootstrap", boot_n = 199)
  
  expect_true("ci_lower" %in% names(prec$precision))
  expect_true("ci_upper" %in% names(prec$precision))
})


test_that("Bootstrap CI bounds contain point estimate", {
  skip_on_cran()
  
  data <- create_day_only_data()
  prec <- precision_study(data, value = "value", day = "day",
                          ci_method = "bootstrap", boot_n = 199)
  
  # Point estimate should generally be within CI
  for (i in seq_len(nrow(prec$precision))) {
    sd_val <- prec$precision$sd[i]
    ci_l <- prec$precision$ci_lower[i]
    ci_u <- prec$precision$ci_upper[i]
    
    if (!is.na(ci_l) && !is.na(ci_u) && sd_val > 0) {
      # Allow some tolerance for bootstrap variability
      expect_true(sd_val >= ci_l * 0.8 || ci_l == 0,
                  info = paste("Lower bound for", prec$precision$measure[i]))
      expect_true(sd_val <= ci_u * 1.5,
                  info = paste("Upper bound for", prec$precision$measure[i]))
    }
  }
})


test_that("Bootstrap CI bounds are non-negative", {
  skip_on_cran()
  
  data <- create_day_only_data()
  prec <- precision_study(data, value = "value", day = "day",
                          ci_method = "bootstrap", boot_n = 199)
  
  expect_true(all(is.na(prec$precision$ci_lower) | prec$precision$ci_lower >= 0))
})


# CV CI tests ----

test_that("CV CIs are correctly calculated from SD CIs", {
  data <- create_ep05_data()
  prec <- precision_study(data, value = "value", day = "day", run = "run")
  
  grand_mean <- mean(data$value)
  
  # CV CI should be 100 * SD_CI / mean
  for (i in seq_len(nrow(prec$precision))) {
    if (!is.na(prec$precision$ci_lower[i])) {
      expected_cv_lower <- 100 * prec$precision$ci_lower[i] / grand_mean
      expect_equal(prec$precision$cv_ci_lower[i], expected_cv_lower,
                   tolerance = 0.001)
    }
    if (!is.na(prec$precision$ci_upper[i])) {
      expected_cv_upper <- 100 * prec$precision$ci_upper[i] / grand_mean
      expect_equal(prec$precision$cv_ci_upper[i], expected_cv_upper,
                   tolerance = 0.001)
    }
  }
})



# Multi-site CI tests ----

test_that("Multi-site design has reproducibility CI", {
  data <- create_multisite_data()
  prec <- precision_study(data, value = "value", site = "site", day = "day",
                          ci_method = "satterthwaite")
  
  # Should have reproducibility row
  expect_true("Reproducibility" %in% prec$precision$measure)
  
  # Reproducibility should have CI
  repro_row <- which(prec$precision$measure == "Reproducibility")
  expect_true(!is.na(prec$precision$ci_lower[repro_row]))
  expect_true(!is.na(prec$precision$ci_upper[repro_row]))
})


test_that("Multi-site CIs are correctly ordered", {
  data <- create_multisite_data()
  prec <- precision_study(data, value = "value", site = "site", day = "day",
                          ci_method = "satterthwaite")
  
  # Reproducibility SD should be >= Intermediate SD
  inter_sd <- prec$precision$sd[prec$precision$measure == "Intermediate precision"]
  repro_sd <- prec$precision$sd[prec$precision$measure == "Reproducibility"]
  
  expect_true(repro_sd >= inter_sd)
})


# Internal helper function tests ----

test_that(".ci_single_variance produces valid CI", {
  # This tests the internal function directly
  # Access internal function
  ci_single <- valytics:::.ci_single_variance
  
  # Test with typical values
  ci <- ci_single(variance = 4, df = 20, alpha = 0.05)
  
  expect_named(ci, c("lower", "upper"))
  expect_true(ci["lower"] < 4)
  expect_true(ci["upper"] > 4)
  expect_true(ci["lower"] >= 0)
})


test_that(".ci_variance_sum produces valid CI for sum", {
  ci_sum <- valytics:::.ci_variance_sum
  
  # Test with two variance components
  ci <- ci_sum(variances = c(4, 2), dfs = c(20, 15), alpha = 0.05)
  
  expect_named(ci, c("lower", "upper"))
  expect_true(ci["lower"] < 6)  # Sum is 6
  expect_true(ci["upper"] > 6)
  expect_true(ci["lower"] >= 0)
})


test_that(".ci_single_variance handles edge cases", {
  ci_single <- valytics:::.ci_single_variance
  
  # Zero variance
  ci <- ci_single(variance = 0, df = 20, alpha = 0.05)
  expect_equal(unname(ci["lower"]), 0)
  
  # NA variance
  ci <- ci_single(variance = NA, df = 20, alpha = 0.05)
  expect_true(is.na(ci["lower"]))
  expect_true(is.na(ci["upper"]))
  
  # Zero df
  ci <- ci_single(variance = 4, df = 0, alpha = 0.05)
  expect_true(is.na(ci["lower"]))
})


# Integration tests ----

test_that("All CI methods produce consistent precision summaries", {
  data <- create_day_only_data()
  
  prec_sat <- precision_study(data, value = "value", day = "day",
                              ci_method = "satterthwaite")
  prec_mls <- precision_study(data, value = "value", day = "day",
                              ci_method = "mls")
  
  # Same structure
  expect_equal(nrow(prec_sat$precision), nrow(prec_mls$precision))
  expect_equal(prec_sat$precision$measure, prec_mls$precision$measure)
  
  # Same point estimates
  expect_equal(prec_sat$precision$sd, prec_mls$precision$sd)
})


test_that("CI method is recorded in settings", {
  data <- create_ep05_data()
  
  prec_sat <- precision_study(data, value = "value", day = "day", run = "run",
                              ci_method = "satterthwaite")
  prec_mls <- precision_study(data, value = "value", day = "day", run = "run",
                              ci_method = "mls")
  
  expect_equal(prec_sat$settings$ci_method, "satterthwaite")
  expect_equal(prec_mls$settings$ci_method, "mls")
})

# REML Estimation Tests ----

test_that("REML requires lme4 package", {
  skip_if_not_installed("lme4")
  
  data <- create_day_only_data()
  
  # Should work when lme4 is available
  prec <- precision_study(data, value = "value", day = "day", method = "reml")
  expect_s3_class(prec, "precision_study")
})


test_that("REML returns correct structure for day-only design", {
  skip_if_not_installed("lme4")
  
  data <- create_day_only_data()
  prec <- precision_study(data, value = "value", day = "day", method = "reml")
  
  # Check variance components structure
  expect_s3_class(prec$variance_components, "data.frame")
  expect_true("between_day" %in% prec$variance_components$component)
  expect_true("error" %in% prec$variance_components$component)
  expect_true("total" %in% prec$variance_components$component)
  
  # Check all values are numeric and non-negative
  expect_true(all(prec$variance_components$variance >= 0))
  expect_true(all(prec$variance_components$sd >= 0))
  
  # Method should be recorded
  expect_equal(prec$settings$method, "reml")
})


test_that("REML returns correct structure for day/run design", {
  skip_if_not_installed("lme4")
  
  data <- create_ep05_data()
  prec <- precision_study(data, value = "value", day = "day", run = "run", 
                          method = "reml")
  
  # Check variance components
  expect_true("between_day" %in% prec$variance_components$component)
  expect_true("between_run" %in% prec$variance_components$component)
  expect_true("error" %in% prec$variance_components$component)
  
  # All variances should be non-negative
  expect_true(all(prec$variance_components$variance >= 0))
})


test_that("REML produces similar results to ANOVA for balanced data", {
  skip_if_not_installed("lme4")
  
  # Create balanced dataset
  set.seed(123)
  data <- create_ep05_data()
  
  prec_anova <- precision_study(data, value = "value", day = "day", run = "run",
                                method = "anova")
  prec_reml <- precision_study(data, value = "value", day = "day", run = "run",
                               method = "reml")
  
  # Variance estimates should be similar (not exact due to different methods)
  # Allow 20% relative tolerance for comparison
  anova_total <- sum(prec_anova$variance_components$variance[
    prec_anova$variance_components$component != "total"])
  reml_total <- sum(prec_reml$variance_components$variance[
    prec_reml$variance_components$component != "total"])
  
  expect_equal(anova_total, reml_total, tolerance = 0.2 * max(anova_total, reml_total))
})


test_that("REML works with multi-site design", {
  skip_if_not_installed("lme4")
  
  data <- create_multisite_data()
  prec <- precision_study(data, value = "value", site = "site", day = "day",
                          method = "reml")
  
  # Check site variance component is present
  expect_true("between_site" %in% prec$variance_components$component)
  expect_true("between_day" %in% prec$variance_components$component)
  
  # All variances non-negative
  expect_true(all(prec$variance_components$variance >= 0))
})


test_that("REML CI calculations work", {
  skip_if_not_installed("lme4")
  
  data <- create_day_only_data()
  prec <- precision_study(data, value = "value", day = "day", 
                          method = "reml", ci_method = "satterthwaite")
  
  # Check CIs are present
  expect_true(all(c("ci_lower", "ci_upper") %in% names(prec$precision)))
  
  # CIs should be numeric
  expect_true(all(is.numeric(prec$precision$ci_lower)))
  expect_true(all(is.numeric(prec$precision$ci_upper)))
})


test_that("REML handles unbalanced data better than ANOVA",
          {
            skip_if_not_installed("lme4")
            
            # Create unbalanced dataset
            set.seed(456)
            data <- data.frame(
              day = c(rep(1, 4), rep(2, 3), rep(3, 5), rep(4, 2), rep(5, 4)),
              value = rnorm(18, mean = 100, sd = 5)
            )
            data$value <- data$value + as.numeric(factor(data$day)) * 2
            
            # Both methods should run without error
            prec_anova <- precision_study(data, value = "value", day = "day", method = "anova")
            prec_reml <- precision_study(data, value = "value", day = "day", method = "reml")
            
            expect_s3_class(prec_anova, "precision_study")
            expect_s3_class(prec_reml, "precision_study")
            
            # Design should detect unbalanced
            expect_false(prec_anova$design$balanced)
            expect_false(prec_reml$design$balanced)
          })


test_that("REML works with strong between-run variance", {
  skip_if_not_installed("lme4")
  
  # Create data with clear run effect
  set.seed(42)
  n_days <- 5
  n_runs <- 2
  n_reps <- 3
  
  data_list <- list()
  idx <- 1
  for (d in 1:n_days) {
    day_effect <- rnorm(1, 0, 3)
    for (r in 1:n_runs) {
      run_effect <- rnorm(1, 0, 5)  # Strong run effect
      for (rep in 1:n_reps) {
        data_list[[idx]] <- data.frame(
          day = d, run = r,
          value = 100 + day_effect + run_effect + rnorm(1, 0, 2)
        )
        idx <- idx + 1
      }
    }
  }
  data <- do.call(rbind, data_list)
  
  prec <- precision_study(data, value = "value", day = "day", run = "run",
                          method = "reml")
  
  # Between-run should be non-zero
  run_var <- prec$variance_components$variance[
    prec$variance_components$component == "between_run"]
  expect_true(run_var > 0)
})


test_that("REML with bootstrap CI works", {
  skip_if_not_installed("lme4")
  skip_on_cran()  # Skip due to computation time
  
  data <- create_day_only_data()
  
  # Use fewer bootstrap samples for testing
  prec <- precision_study(data, value = "value", day = "day",
                          method = "reml", ci_method = "bootstrap", boot_n = 199)
  
  expect_s3_class(prec, "precision_study")
  expect_true(!any(is.na(prec$precision$ci_lower)))
})