# Tests for verify_precision() - Phase 4
# Precision verification against manufacturer claims

# Test Data Setup ----

#' Create simple test data for verification
#' @noRd
create_verification_data <- function(seed = 42, n = 25, mean_val = 100, sd = 3.5) {
  set.seed(seed)
  rnorm(n, mean = mean_val, sd = sd)
}


#' Create precision study data for verification testing
#' @noRd
create_prec_study_data <- function(seed = 42, mean_val = 100,
                                   sd_day = 1.5, sd_error = 2.5) {
  set.seed(seed)
  
  n_days <- 5
  n_reps <- 5
  
  data <- expand.grid(
    day = 1:n_days,
    replicate = 1:n_reps
  )
  
  day_effects <- rnorm(n_days, 0, sd_day)
  errors <- rnorm(nrow(data), 0, sd_error)
  
  day_effect <- day_effects[as.numeric(as.factor(data$day))]
  data$value <- mean_val + day_effect + errors
  
  data
}


# Input Validation Tests ----

test_that("verify_precision validates claimed values", {
  x <- create_verification_data()
  
  # Must provide at least one of claimed_cv or claimed_sd
  expect_error(
    verify_precision(x, mean_value = 100),
    "Either `claimed_cv` or `claimed_sd` must be provided"
  )
  
  # Invalid claimed_cv
  expect_error(
    verify_precision(x, claimed_cv = -5, mean_value = 100),
    "`claimed_cv` must be a single positive number"
  )
  
  expect_error(
    verify_precision(x, claimed_cv = 0, mean_value = 100),
    "`claimed_cv` must be a single positive number"
  )
  
  expect_error(
    verify_precision(x, claimed_cv = c(3, 4), mean_value = 100),
    "`claimed_cv` must be a single positive number"
  )
  
  # Invalid claimed_sd
  expect_error(
    verify_precision(x, claimed_sd = -3, mean_value = 100),
    "`claimed_sd` must be a single positive number"
  )
  
  expect_error(
    verify_precision(x, claimed_sd = 0, mean_value = 100),
    "`claimed_sd` must be a single positive number"
  )
})


test_that("verify_precision validates alpha", {
  x <- create_verification_data()
  
  expect_error(
    verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = 0),
    "`alpha` must be a single number between 0 and 1"
  )
  
  expect_error(
    verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = 1),
    "`alpha` must be a single number between 0 and 1"
  )
  
  expect_error(
    verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = -0.1),
    "`alpha` must be a single number between 0 and 1"
  )
  
  expect_error(
    verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = c(0.05, 0.1)),
    "`alpha` must be a single number between 0 and 1"
  )
})


test_that("verify_precision validates conf_level", {
  x <- create_verification_data()
  
  expect_error(
    verify_precision(x, claimed_cv = 5, mean_value = 100, conf_level = 0),
    "`conf_level` must be a single number between 0 and 1"
  )
  
  expect_error(
    verify_precision(x, claimed_cv = 5, mean_value = 100, conf_level = 1),
    "`conf_level` must be a single number between 0 and 1"
  )
  
  expect_error(
    verify_precision(x, claimed_cv = 5, mean_value = 100, conf_level = 1.5),
    "`conf_level` must be a single number between 0 and 1"
  )
})


test_that("verify_precision validates numeric vector input", {
  # Non-numeric input
  expect_error(
    verify_precision("not numeric", claimed_cv = 5, mean_value = 100),
    "`x` must be a numeric vector"
  )
  
  # Too few observations
  expect_error(
    verify_precision(c(1, 2), claimed_cv = 5, mean_value = 100),
    "At least 3 observations are required"
  )
  
  # Empty vector after NA removal
  expect_error(
    verify_precision(rep(NA_real_, 5), claimed_cv = 5, mean_value = 100),
    "At least 3 observations are required"
  )
})


# Basic Functionality Tests ----

test_that("verify_precision returns correct class", {
  x <- create_verification_data()
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  expect_s3_class(result, "verify_precision")
  expect_s3_class(result, "valytics_precision")
  expect_s3_class(result, "valytics_result")
})


test_that("verify_precision returns expected structure", {
  x <- create_verification_data()
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  # Check main components
  expect_named(result, c("input", "observed", "claimed", "test",
                         "verification", "ci", "settings", "call"))
  
  # Check input
  expect_named(result$input, c("n", "df", "mean_value", "source"))
  expect_equal(result$input$n, 25)
  expect_equal(result$input$df, 24)
  
  # Check observed
  expect_named(result$observed, c("sd", "cv_pct", "variance"))
  expect_true(is.numeric(result$observed$sd))
  expect_true(is.numeric(result$observed$cv_pct))
  
  # Check claimed
  expect_named(result$claimed, c("sd", "cv_pct", "variance"))
  expect_equal(result$claimed$cv_pct, 5)
  
  # Check test
  expect_named(result$test, c("statistic", "df", "p_value",
                              "alternative", "method"))
  expect_true(is.numeric(result$test$statistic))
  expect_true(is.numeric(result$test$p_value))
  
  # Check verification
  expect_named(result$verification, c("verified", "ratio", "cv_ratio",
                                      "upper_verification_limit"))
  expect_true(is.logical(result$verification$verified))
  
  # Check CI
  expect_named(result$ci, c("variance_ci", "sd_ci", "cv_ci"))
  expect_named(result$ci$sd_ci, c("lower", "upper"))
})


test_that("verify_precision works with claimed_cv", {
  set.seed(42)
  x <- rnorm(25, mean = 100, sd = 3.5)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  # Claimed CV should be as specified
  expect_equal(result$claimed$cv_pct, 5)
  
  # Claimed SD should be derived from CV and mean
  expect_equal(result$claimed$sd, 5)  # 5% of 100
  
  # Observed CV should be reasonable
  expect_true(result$observed$cv_pct > 0)
  expect_true(result$observed$cv_pct < 20)  # Reasonable bound
})


test_that("verify_precision works with claimed_sd", {
  set.seed(42)
  x <- rnorm(25, mean = 100, sd = 3.5)
  
  result <- verify_precision(x, claimed_sd = 5, mean_value = 100)
  
  # Claimed SD should be as specified
  expect_equal(result$claimed$sd, 5)
  
  # Claimed CV should be derived from SD and mean
  expect_equal(result$claimed$cv_pct, 5)  # 5/100 * 100%
})


test_that("verify_precision works with both claimed_cv and claimed_sd", {
  x <- create_verification_data()
  
  result <- verify_precision(x, claimed_cv = 5, claimed_sd = 4, mean_value = 100)
  
  # Both should be used as given
  expect_equal(result$claimed$cv_pct, 5)
  expect_equal(result$claimed$sd, 4)
})


# Verification Logic Tests ----

test_that("verify_precision correctly verifies when observed < claimed", {
  # Generate data with SD = 3, claim CV = 5%
  set.seed(123)
  x <- rnorm(30, mean = 100, sd = 3)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  # Should be verified because observed CV is lower than claimed
  expect_true(result$verification$verified)
  expect_lt(result$observed$cv_pct, result$verification$upper_verification_limit)
})


test_that("verify_precision correctly fails when observed >> claimed", {
  # Generate data with SD = 8, claim CV = 3%
  set.seed(123)
  x <- rnorm(30, mean = 100, sd = 8)
  
  result <- verify_precision(x, claimed_cv = 3, mean_value = 100)
  
  # Should NOT be verified because observed CV is much higher than claimed
  expect_false(result$verification$verified)
  expect_gt(result$observed$cv_pct, result$verification$upper_verification_limit)
})


test_that("verify_precision upper verification limit is correct", {
  set.seed(42)
  x <- rnorm(25, mean = 100, sd = 4)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = 0.05)
  
  # Calculate expected UVL manually
  df <- 24
  chi_sq_crit <- qchisq(0.95, df = df)
  claimed_var <- 25  # (5/100 * 100)^2
  uvl_var <- claimed_var * chi_sq_crit / df
  uvl_cv <- 100 * sqrt(uvl_var) / 100
  
  expect_equal(result$verification$upper_verification_limit, uvl_cv,
               tolerance = 0.001)
})


test_that("verify_precision chi-square statistic is correct", {
  set.seed(42)
  x <- rnorm(25, mean = 100, sd = 4)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  # Calculate expected chi-square manually
  observed_var <- var(x)
  claimed_var <- 25  # (5/100 * 100)^2
  expected_chisq <- 24 * observed_var / claimed_var
  
  expect_equal(result$test$statistic, expected_chisq, tolerance = 0.001)
})


# Alternative Hypothesis Tests ----

test_that("verify_precision respects alternative argument", {
  set.seed(42)
  x <- rnorm(30, mean = 100, sd = 4)
  
  result_less <- verify_precision(x, claimed_cv = 5, mean_value = 100,
                                  alternative = "less")
  result_greater <- verify_precision(x, claimed_cv = 5, mean_value = 100,
                                     alternative = "greater")
  result_two <- verify_precision(x, claimed_cv = 5, mean_value = 100,
                                 alternative = "two.sided")
  
  # All should have different p-values (generally)
  expect_equal(result_less$test$alternative, "less")
  expect_equal(result_greater$test$alternative, "greater")
  expect_equal(result_two$test$alternative, "two.sided")
  
  # p-values should differ
  # For two-sided: p = 2 * min(p_less, p_greater)
  expect_true(result_two$test$p_value <= 1)
})


# Confidence Interval Tests ----

test_that("verify_precision confidence intervals are valid", {
  x <- create_verification_data()
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100, conf_level = 0.95)
  
  # Lower < observed < upper
  expect_lt(result$ci$sd_ci["lower"], result$observed$sd)
  expect_gt(result$ci$sd_ci["upper"], result$observed$sd)
  
  expect_lt(result$ci$cv_ci["lower"], result$observed$cv_pct)
  expect_gt(result$ci$cv_ci["upper"], result$observed$cv_pct)
  
  expect_lt(result$ci$variance_ci["lower"], result$observed$variance)
  expect_gt(result$ci$variance_ci["upper"], result$observed$variance)
})


test_that("verify_precision CI width changes with conf_level", {
  x <- create_verification_data()
  
  result_90 <- verify_precision(x, claimed_cv = 5, mean_value = 100,
                                conf_level = 0.90)
  result_95 <- verify_precision(x, claimed_cv = 5, mean_value = 100,
                                conf_level = 0.95)
  result_99 <- verify_precision(x, claimed_cv = 5, mean_value = 100,
                                conf_level = 0.99)
  
  # CI width should increase with conf_level
  width_90 <- diff(result_90$ci$sd_ci)
  width_95 <- diff(result_95$ci$sd_ci)
  width_99 <- diff(result_99$ci$sd_ci)
  
  expect_lt(width_90, width_95)
  expect_lt(width_95, width_99)
})


# Input Type Tests ----

test_that("verify_precision works with precision_study object", {
  skip_if_not_installed("valytics")
  
  data <- create_prec_study_data()
  prec <- precision_study(data, value = "value", day = "day")
  
  result <- verify_precision(prec, claimed_cv = 5)
  
  expect_s3_class(result, "verify_precision")
  expect_equal(result$input$source, "precision_study object")
  expect_true(result$input$n > 0)
})


test_that("verify_precision works with data frame input", {
  skip_if_not_installed("valytics")
  
  data <- create_prec_study_data()
  
  result <- verify_precision(
    data,
    claimed_cv = 5,
    value = "value",
    day = "day"
  )
  
  expect_s3_class(result, "verify_precision")
  expect_equal(result$input$source, "data frame (via precision_study)")
})


test_that("verify_precision handles NA values in numeric vector", {
  x <- c(rnorm(20, 100, 3.5), NA, NA, NA)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  expect_equal(result$input$n, 20)
  expect_equal(result$input$df, 19)
})


test_that("verify_precision uses provided mean_value for CV calculation", {
  set.seed(42)
  x <- rnorm(25, mean = 100, sd = 4)
  
  # Override mean with a different value
  result <- verify_precision(x, claimed_cv = 5, mean_value = 200)
  
  # CV should be calculated using 200 as the mean
  expected_cv <- 100 * sd(x) / 200
  expect_equal(result$observed$cv_pct, expected_cv, tolerance = 0.001)
})


test_that("verify_precision calculates mean from data when not provided", {
  set.seed(42)
  x <- rnorm(25, mean = 100, sd = 4)
  
  # When using claimed_sd, mean can be inferred from data
  # But when using claimed_cv, we need mean to convert
  # Test with claimed_sd:
  result <- verify_precision(x, claimed_sd = 4)
  
  expect_equal(result$input$mean_value, mean(x), tolerance = 0.001)
})


# Print and Summary Tests ----

test_that("print.verify_precision runs without error", {
  x <- create_verification_data()
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  expect_output(print(result), "Precision Verification")
  expect_output(print(result), "VERIFIED|NOT VERIFIED")
  expect_output(print(result), "Chi-square")
})


test_that("summary.verify_precision returns correct class", {
  x <- create_verification_data()
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  summ <- summary(result)
  
  expect_s3_class(summ, "summary.verify_precision")
  expect_named(summ, c("call", "input", "observed", "claimed", "test",
                       "verification", "ci", "settings"))
})


test_that("print.summary.verify_precision runs without error", {
  x <- create_verification_data()
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  summ <- summary(result)
  
  expect_output(print(summ), "Precision Verification - Detailed Summary")
  expect_output(print(summ), "Confidence Intervals")
  expect_output(print(summ), "Verification Outcome")
})


# Edge Cases ----

test_that("verify_precision handles minimum sample size (n=3)", {
  x <- c(100, 102, 98)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  expect_s3_class(result, "verify_precision")
  expect_equal(result$input$n, 3)
  expect_equal(result$input$df, 2)
})


test_that("verify_precision handles perfect precision (zero variance)", {
  x <- rep(100, 10)  # No variation
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  expect_equal(result$observed$sd, 0)
  expect_equal(result$observed$cv_pct, 0)
  expect_true(result$verification$verified)  # Zero variance is always verified
})


test_that("verify_precision handles large sample sizes", {
  set.seed(42)
  x <- rnorm(500, mean = 100, sd = 4)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  expect_s3_class(result, "verify_precision")
  expect_equal(result$input$n, 500)
  expect_equal(result$input$df, 499)
})


test_that("verify_precision handles small mean values", {
  set.seed(42)
  x <- rnorm(25, mean = 1, sd = 0.04)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 1)
  
  expect_s3_class(result, "verify_precision")
  expect_true(result$observed$cv_pct > 0)
})


test_that("verify_precision handles large mean values", {
  set.seed(42)
  x <- rnorm(25, mean = 10000, sd = 400)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 10000)
  
  expect_s3_class(result, "verify_precision")
  expect_true(result$observed$cv_pct > 0)
})


# Numerical Consistency Tests ----

test_that("verify_precision ratios are consistent", {
  x <- create_verification_data()
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  # CV ratio should equal SD ratio (when mean is same)
  expect_equal(result$verification$cv_ratio,
               result$observed$cv_pct / result$claimed$cv_pct,
               tolerance = 0.001)
  
  # Variance ratio should be SD ratio squared
  expect_equal(result$verification$ratio,
               (result$observed$sd / result$claimed$sd)^2,
               tolerance = 0.001)
})


test_that("verify_precision variance calculations are consistent", {
  x <- create_verification_data()
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  # Observed variance should equal sd^2
  expect_equal(result$observed$variance, result$observed$sd^2, tolerance = 0.001)
  
  # Claimed variance should equal sd^2
  expect_equal(result$claimed$variance, result$claimed$sd^2, tolerance = 0.001)
})


# Degrees of Freedom Tests ----

test_that("verify_precision uses correct df for numeric vector", {
  x <- rnorm(30, 100, 4)
  
  result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
  
  expect_equal(result$input$df, 29)
  expect_equal(result$test$df, 29)
})


test_that("verify_precision df equals input n-1 for vector input", {
  for (n in c(5, 10, 25, 50)) {
    x <- rnorm(n, 100, 4)
    result <- verify_precision(x, claimed_cv = 5, mean_value = 100)
    
    expect_equal(result$input$df, n - 1,
                 info = sprintf("Failed for n = %d", n))
  }
})


# Significance Level Tests ----

test_that("verify_precision alpha affects verification limit", {
  x <- create_verification_data()
  
  result_01 <- verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = 0.01)
  result_05 <- verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = 0.05)
  result_10 <- verify_precision(x, claimed_cv = 5, mean_value = 100, alpha = 0.10)
  
  # More stringent alpha -> higher UVL (easier to verify)
  expect_gt(result_01$verification$upper_verification_limit,
            result_05$verification$upper_verification_limit)
  expect_gt(result_05$verification$upper_verification_limit,
            result_10$verification$upper_verification_limit)
})


# Integration Tests ----

test_that("verify_precision results are reproducible with same seed", {
  set.seed(123)
  x1 <- rnorm(25, 100, 4)
  result1 <- verify_precision(x1, claimed_cv = 5, mean_value = 100)
  
  set.seed(123)
  x2 <- rnorm(25, 100, 4)
  result2 <- verify_precision(x2, claimed_cv = 5, mean_value = 100)
  
  expect_equal(result1$observed$cv_pct, result2$observed$cv_pct)
  expect_equal(result1$test$statistic, result2$test$statistic)
  expect_equal(result1$verification$verified, result2$verification$verified)
})


test_that("verify_precision handles realistic clinical lab scenario", {
  # Simulate verification study: 5 days x 5 replicates
  # Manufacturer claims CV = 4%
  set.seed(42)
  
  n_days <- 5
  n_reps <- 5
  true_sd <- 3.2  # Actual CV ~ 3.2%
  mean_val <- 100
  
  data <- expand.grid(day = 1:n_days, replicate = 1:n_reps)
  day_effects <- rnorm(n_days, 0, 1)
  errors <- rnorm(nrow(data), 0, 2)
  
  data$value <- mean_val + day_effects[data$day] + errors
  
  # Test with data frame input
  result <- verify_precision(
    data,
    claimed_cv = 4,
    value = "value",
    day = "day"
  )
  
  expect_s3_class(result, "verify_precision")
  # With these parameters, should likely be verified
  expect_true(is.logical(result$verification$verified))
})