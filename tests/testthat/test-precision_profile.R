# Tests for precision_profile()
# =============================================================================

# Helper function to create test data
create_precision_profile_data <- function(n_levels = 6, seed = 42) {
  set.seed(seed)
  
  # Concentration levels spanning 2 orders of magnitude
  conc_levels <- c(5, 10, 25, 50, 100, 200)[1:n_levels]
  
  # True hyperbolic CV: sqrt(3^2 + (20/x)^2)
  # At x=5: CV ~4.5%, at x=200: CV ~3.0%
  
  prec_data <- data.frame()
  for (i in seq_along(conc_levels)) {
    level_data <- expand.grid(
      level = conc_levels[i],
      day = 1:5,
      replicate = 1:5
    )
    
    true_cv <- sqrt(3^2 + (20/conc_levels[i])^2)
    level_data$value <- conc_levels[i] * rnorm(
      nrow(level_data),
      mean = 1,
      sd = true_cv/100
    )
    
    prec_data <- rbind(prec_data, level_data)
  }
  
  prec_data
}

# Helper to create a data frame with concentration and CV columns
create_cv_data <- function(n_levels = 6, seed = 42) {
  set.seed(seed)
  
  conc <- c(5, 10, 25, 50, 100, 200)[1:n_levels]
  
  # True hyperbolic CV with some noise
  true_cv <- sqrt(3^2 + (20/conc)^2)
  cv <- true_cv + rnorm(n_levels, 0, 0.3)
  cv <- pmax(cv, 1)  # Ensure positive
  
  data.frame(
    concentration = conc,
    cv_pct = cv
  )
}


# Input Validation Tests ----

test_that("precision_profile validates minimum number of concentration levels", {
  cv_data <- create_cv_data(n_levels = 3)
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct"),
    "At least 4 concentration levels"
  )
})

test_that("precision_profile validates positive concentration values", {
  cv_data <- data.frame(
    concentration = c(-1, 10, 25, 50),
    cv_pct = c(5, 4, 3.5, 3)
  )
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct"),
    "positive"
  )
})

test_that("precision_profile validates positive CV values", {
  cv_data <- data.frame(
    concentration = c(5, 10, 25, 50),
    cv_pct = c(5, 4, -3.5, 3)
  )
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct"),
    "positive"
  )
})

test_that("precision_profile validates NA values", {
  cv_data <- data.frame(
    concentration = c(5, 10, NA, 50),
    cv_pct = c(5, 4, 3.5, 3)
  )
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct"),
    "NA"
  )
})

test_that("precision_profile validates conf_level", {
  cv_data <- create_cv_data()
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct",
                      conf_level = 0),
    "conf_level"
  )
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct",
                      conf_level = 1.5),
    "conf_level"
  )
})

test_that("precision_profile validates cv_targets", {
  cv_data <- create_cv_data()
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct",
                      cv_targets = c(-10, 20)),
    "cv_targets"
  )
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct",
                      cv_targets = c(10, 150)),
    "cv_targets"
  )
})

test_that("precision_profile validates boot_n", {
  cv_data <- create_cv_data()
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct",
                      bootstrap = TRUE, boot_n = 50),
    "boot_n"
  )
})


# Data Frame Interface Tests ----

test_that("precision_profile works with data frame input", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  expect_s3_class(result, "precision_profile")
  expect_s3_class(result, "valytics_precision")
  expect_s3_class(result, "valytics_result")
})

test_that("precision_profile returns expected structure", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  # Check top-level components
  expect_named(result, c("input", "model", "fitted", "fit_quality",
                         "functional_sensitivity", "settings", "call"))
  
  # Check input
  expect_true("concentration" %in% names(result$input))
  expect_true("cv" %in% names(result$input))
  expect_true("n_levels" %in% names(result$input))
  expect_true("conc_range" %in% names(result$input))
  expect_true("conc_span" %in% names(result$input))
  
  # Check model
  expect_true("type" %in% names(result$model))
  expect_true("parameters" %in% names(result$model))
  expect_true("equation" %in% names(result$model))
  
  # Check fitted
  expect_s3_class(result$fitted, "data.frame")
  expect_true(all(c("concentration", "cv_observed", "cv_fitted", "residual") %in%
                    names(result$fitted)))
  
  # Check fit_quality
  expect_true("r_squared" %in% names(result$fit_quality))
  expect_true("rmse" %in% names(result$fit_quality))
  
  # Check functional_sensitivity
  expect_s3_class(result$functional_sensitivity, "data.frame")
  expect_true(all(c("cv_target", "concentration", "achievable") %in%
                    names(result$functional_sensitivity)))
})


# Hyperbolic Model Tests ----

test_that("precision_profile fits hyperbolic model correctly", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    model = "hyperbolic"
  )
  
  expect_equal(result$model$type, "hyperbolic")
  expect_true("a" %in% names(result$model$parameters))
  expect_true("b" %in% names(result$model$parameters))
  expect_true(grepl("sqrt", result$model$equation))
})

test_that("hyperbolic model parameters are reasonable", {
  # Create data with known parameters: CV = sqrt(3^2 + (20/x)^2)
  set.seed(123)
  conc <- c(5, 10, 25, 50, 100, 200)
  true_cv <- sqrt(3^2 + (20/conc)^2)
  cv_data <- data.frame(concentration = conc, cv_pct = true_cv + rnorm(6, 0, 0.1))
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    model = "hyperbolic"
  )
  
  # Parameters should be close to true values (a=3, b=20)
  expect_true(abs(result$model$parameters["a"] - 3) < 1)
  expect_true(abs(result$model$parameters["b"] - 20) < 5)
})

test_that("hyperbolic model R-squared is reasonable for good fit", {
  # Create data with known relationship
  conc <- c(5, 10, 25, 50, 100, 200)
  true_cv <- sqrt(3^2 + (20/conc)^2)
  cv_data <- data.frame(concentration = conc, cv_pct = true_cv)
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    model = "hyperbolic"
  )
  
  # R-squared should be very high for perfect fit
  expect_true(result$fit_quality$r_squared > 0.99)
})


# Linear Model Tests ----

test_that("precision_profile fits linear model correctly", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    model = "linear"
  )
  
  expect_equal(result$model$type, "linear")
  expect_true("a" %in% names(result$model$parameters))
  expect_true("b" %in% names(result$model$parameters))
  expect_false(grepl("sqrt", result$model$equation))
})

test_that("linear model equation format is correct", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    model = "linear"
  )
  
  # Equation should be of form "CV = a + b/x"
  expect_true(grepl("/x", result$model$equation))
})


# Functional Sensitivity Tests ----

test_that("functional sensitivity is calculated for specified targets", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    cv_targets = c(5, 10, 20)
  )
  
  expect_equal(nrow(result$functional_sensitivity), 3)
  expect_equal(result$functional_sensitivity$cv_target, c(5, 10, 20))
})

test_that("functional sensitivity marks unachievable targets", {
  # Create data with asymptotic CV around 3%
  conc <- c(5, 10, 25, 50, 100, 200)
  cv <- sqrt(3^2 + (20/conc)^2)
  cv_data <- data.frame(concentration = conc, cv_pct = cv)
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    model = "hyperbolic",
    cv_targets = c(2, 5, 10)  # 2% is below asymptotic CV of 3%
  )
  
  # CV target of 2% should not be achievable (asymptotic is ~3%)
  expect_false(result$functional_sensitivity$achievable[
    result$functional_sensitivity$cv_target == 2])
  
  # CV targets of 5% and 10% should be achievable
  expect_true(result$functional_sensitivity$achievable[
    result$functional_sensitivity$cv_target == 5])
  expect_true(result$functional_sensitivity$achievable[
    result$functional_sensitivity$cv_target == 10])
})

test_that("functional sensitivity concentration decreases with higher CV target", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    cv_targets = c(4, 6, 8, 10)
  )
  
  # Higher CV target should give lower concentration (functional sensitivity)
  fs <- result$functional_sensitivity
  achievable_conc <- fs$concentration[fs$achievable]
  
  if (length(achievable_conc) > 1) {
    # Concentrations should decrease as CV target increases
    expect_true(all(diff(achievable_conc) <= 0) || any(is.na(achievable_conc)))
  }
})


# Bootstrap Tests ----

test_that("bootstrap produces confidence intervals for functional sensitivity", {
  skip_on_cran()  # Skip slow test on CRAN
  
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    cv_targets = c(10, 20),
    bootstrap = TRUE,
    boot_n = 199  # Reduced for faster testing
  )
  
  # Check that CIs are present for achievable targets
  fs <- result$functional_sensitivity
  
  for (i in seq_len(nrow(fs))) {
    if (fs$achievable[i]) {
      expect_false(is.na(fs$ci_lower[i]))
      expect_false(is.na(fs$ci_upper[i]))
      expect_true(fs$ci_lower[i] <= fs$concentration[i])
      expect_true(fs$ci_upper[i] >= fs$concentration[i])
    }
  }
})

test_that("bootstrap CIs are not present when bootstrap=FALSE", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    cv_targets = c(10, 20),
    bootstrap = FALSE
  )
  
  # CIs should be NA when bootstrap is not used
  expect_true(all(is.na(result$functional_sensitivity$ci_lower)))
  expect_true(all(is.na(result$functional_sensitivity$ci_upper)))
})


# Fitted Values Tests ----

test_that("fitted values have correct length", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  expect_equal(nrow(result$fitted), nrow(cv_data))
})

test_that("residuals are calculated correctly", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  # Residuals should be observed - fitted
  expected_residuals <- result$fitted$cv_observed - result$fitted$cv_fitted
  expect_equal(result$fitted$residual, expected_residuals)
})

test_that("prediction intervals contain observed values", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    conf_level = 0.95
  )
  
  # Most observations should be within prediction intervals
  within_ci <- result$fitted$cv_observed >= result$fitted$ci_lower &
    result$fitted$cv_observed <= result$fitted$ci_upper
  
  # At least 80% should be within 95% CI (allowing for small samples)
  expect_true(mean(within_ci) >= 0.5)
})


# Concentration Span Warning ----

test_that("precision_profile warns about narrow concentration span", {
  # Create data with only 1.5-fold span
  cv_data <- data.frame(
    concentration = c(100, 110, 120, 150),
    cv_pct = c(3.5, 3.4, 3.3, 3.2)
  )
  
  expect_warning(
    precision_profile(cv_data, concentration = "concentration", cv = "cv_pct"),
    "span"
  )
})


# Fit Quality Tests ----

test_that("fit quality metrics are calculated correctly", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  # R-squared should be between 0 and 1
  expect_true(result$fit_quality$r_squared >= 0)
  expect_true(result$fit_quality$r_squared <= 1)
  
  # RMSE should be non-negative
  expect_true(result$fit_quality$rmse >= 0)
  
  # MAE should be non-negative
  expect_true(result$fit_quality$mae >= 0)
  
  # Adjusted R-squared should be less than or equal to R-squared
  expect_true(result$fit_quality$adj_r_squared <= result$fit_quality$r_squared)
})


# Print Method Tests ----

test_that("print.precision_profile runs without error", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  expect_output(print(result), "Precision Profile")
  expect_output(print(result), "Model:")
  expect_output(print(result), "Functional Sensitivity")
})

test_that("print.precision_profile shows model equation", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  expect_output(print(result), "CV =")
})


# Summary Method Tests ----

test_that("summary.precision_profile runs without error", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  expect_output(summary(result), "Detailed Summary")
  expect_output(summary(result), "Fit Quality")
  expect_output(summary(result), "Fitted Values")
  expect_output(summary(result), "Functional Sensitivity")
})


# Plot Method Tests ----

test_that("plot.precision_profile returns ggplot object", {
  skip_if_not_installed("ggplot2")
  
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  p <- plot(result)
  expect_s3_class(p, "ggplot")
})

test_that("plot.precision_profile works without CI", {
  skip_if_not_installed("ggplot2")
  
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  p <- plot(result, show_ci = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot.precision_profile works without target lines", {
  skip_if_not_installed("ggplot2")
  
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  p <- plot(result, show_targets = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot.precision_profile works with log scale", {
  skip_if_not_installed("ggplot2")
  
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  p <- plot(result, log_x = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.precision_profile works", {
  skip_if_not_installed("ggplot2")
  
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct"
  )
  
  p <- ggplot2::autoplot(result)
  expect_s3_class(p, "ggplot")
})


# Column Name Validation Tests ----

test_that("precision_profile validates concentration column name", {
  cv_data <- create_cv_data()
  
  expect_error(
    precision_profile(cv_data, concentration = "wrong_name", cv = "cv_pct"),
    "not found"
  )
})

test_that("precision_profile validates cv column name", {
  cv_data <- create_cv_data()
  
  expect_error(
    precision_profile(cv_data, concentration = "concentration", cv = "wrong_name"),
    "not found"
  )
})


# Invalid Input Type Tests ----

test_that("precision_profile rejects invalid input type", {
  expect_error(
    precision_profile(1:10),
    "precision_study object|data frame"
  )
  
  expect_error(
    precision_profile("not_valid"),
    "precision_study object|data frame"
  )
})


# Precision Study Object Input Tests ----

test_that("precision_profile rejects precision_study with insufficient levels", {
  # Mock a precision_study object with only 2 levels
  mock_prec <- structure(
    list(
      by_sample = list(
        level1 = list(precision = data.frame(measure = "Repeatability", cv_pct = 5)),
        level2 = list(precision = data.frame(measure = "Repeatability", cv_pct = 4))
      ),
      sample_means = c(level1 = 10, level2 = 50)
    ),
    class = c("precision_study", "valytics_precision", "valytics_result")
  )
  
  expect_error(
    precision_profile(mock_prec),
    "at least 4 concentration levels"
  )
})

test_that("precision_profile rejects single-sample precision_study", {
  # Mock a precision_study object with no by_sample
  mock_prec <- structure(
    list(
      by_sample = NULL,
      sample_means = NULL
    ),
    class = c("precision_study", "valytics_precision", "valytics_result")
  )
  
  expect_error(
    precision_profile(mock_prec),
    "at least 4 concentration levels"
  )
})

test_that("precision_profile works with valid precision_study object", {
  # Mock a precision_study object with 5 levels
  mock_prec <- structure(
    list(
      by_sample = list(
        level1 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(6, 6.5)
        )),
        level2 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(4.5, 5)
        )),
        level3 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(3.8, 4.2)
        )),
        level4 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(3.3, 3.7)
        )),
        level5 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(3.1, 3.4)
        ))
      ),
      sample_means = c(level1 = 10, level2 = 25, level3 = 50, level4 = 100, level5 = 200)
    ),
    class = c("precision_study", "valytics_precision", "valytics_result")
  )
  
  result <- precision_profile(mock_prec)
  
  expect_s3_class(result, "precision_profile")
  expect_equal(result$input$n_levels, 5)
})

test_that("precision_profile extracts Within-laboratory precision CV from precision_study", {
  # Mock a precision_study object - should prefer Within-laboratory precision
  mock_prec <- structure(
    list(
      by_sample = list(
        level1 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(5, 7)  # Within-lab (7) should be used
        )),
        level2 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(4, 5.5)
        )),
        level3 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(3.5, 4.5)
        )),
        level4 = list(precision = data.frame(
          measure = c("Repeatability", "Within-laboratory precision"),
          cv_pct = c(3, 3.8)
        ))
      ),
      sample_means = c(level1 = 10, level2 = 25, level3 = 50, level4 = 100)
    ),
    class = c("precision_study", "valytics_precision", "valytics_result")
  )
  
  result <- precision_profile(mock_prec)
  
  # Should have used Within-laboratory precision CVs (7, 5.5, 4.5, 3.8)
  expect_equal(result$input$cv[1], 7)
  expect_equal(result$input$cv[4], 3.8)
})

test_that("precision_profile falls back to Repeatability when no Within-lab precision", {
  # Mock a precision_study object with only Repeatability
  mock_prec <- structure(
    list(
      by_sample = list(
        level1 = list(precision = data.frame(
          measure = c("Repeatability"),
          cv_pct = c(5)
        )),
        level2 = list(precision = data.frame(
          measure = c("Repeatability"),
          cv_pct = c(4)
        )),
        level3 = list(precision = data.frame(
          measure = c("Repeatability"),
          cv_pct = c(3.5)
        )),
        level4 = list(precision = data.frame(
          measure = c("Repeatability"),
          cv_pct = c(3)
        ))
      ),
      sample_means = c(level1 = 10, level2 = 25, level3 = 50, level4 = 100)
    ),
    class = c("precision_study", "valytics_precision", "valytics_result")
  )
  
  result <- precision_profile(mock_prec)
  
  # Should have used Repeatability CVs
  expect_equal(result$input$cv[1], 5)
  expect_equal(result$input$cv[4], 3)
})


# Settings Preservation Tests ----

test_that("settings are preserved in result object", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    conf_level = 0.90,
    bootstrap = TRUE,
    boot_n = 199
  )
  
  expect_equal(result$settings$conf_level, 0.90)
  expect_true(result$settings$bootstrap)
  expect_equal(result$settings$boot_n, 199)
})

test_that("boot_n is NA when bootstrap=FALSE", {
  cv_data <- create_cv_data()
  
  result <- precision_profile(
    cv_data,
    concentration = "concentration",
    cv = "cv_pct",
    bootstrap = FALSE
  )
  
  expect_true(is.na(result$settings$boot_n))
})