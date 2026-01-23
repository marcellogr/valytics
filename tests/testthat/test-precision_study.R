# Tests for precision_study() - Phase 1a
# Input validation and design detection
# =============================================================================

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
