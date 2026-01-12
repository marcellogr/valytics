# Tests for sigma_metric()
# =============================================================================

# Basic functionality ----

test_that("sigma_metric returns correct class", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  expect_s3_class(sm, "sigma_metric")
  expect_s3_class(sm, "valytics_ate")
  expect_s3_class(sm, "valytics_result")
})


test_that("sigma_metric returns expected structure", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  # Check top-level names
  expect_named(sm, c("sigma", "input", "interpretation"))

  # Check input
  expect_named(sm$input, c("bias", "cv", "tea"))

  # Check interpretation
  expect_named(sm$interpretation, c("category", "defect_rate"))

  # Sigma should be numeric
  expect_true(is.numeric(sm$sigma))
})


test_that("sigma_metric calculates correctly", {
  # Sigma = (TEa - |Bias|) / CV
  # Sigma = (10 - 1.5) / 2.0 = 8.5 / 2.0 = 4.25
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  expect_equal(sm$sigma, 4.25)
})


test_that("sigma_metric uses absolute value of bias", {
  # Negative bias should give same result as positive
  sm_pos <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)
  sm_neg <- sigma_metric(bias = -1.5, cv = 2.0, tea = 10)

  expect_equal(sm_pos$sigma, sm_neg$sigma)
})


test_that("sigma_metric handles zero bias", {
  # Sigma = (10 - 0) / 2.0 = 5.0
  sm <- sigma_metric(bias = 0, cv = 2.0, tea = 10)

  expect_equal(sm$sigma, 5.0)
})


test_that("sigma_metric stores input values correctly", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  expect_equal(sm$input$bias, 1.5)
  expect_equal(sm$input$cv, 2.0)
  expect_equal(sm$input$tea, 10)
})


# Interpretation categories ----

test_that("sigma_metric interprets >= 6 as World Class", {
  # Sigma = (20 - 2) / 3 = 6.0
  sm <- sigma_metric(bias = 2, cv = 3, tea = 20)

  expect_equal(sm$interpretation$category, "World Class")
  expect_equal(sm$interpretation$defect_rate, 3.4)
})


test_that("sigma_metric interprets >= 5 as Excellent", {
  # Sigma = (15 - 2.5) / 2.5 = 5.0
  sm <- sigma_metric(bias = 2.5, cv = 2.5, tea = 15)

  expect_equal(sm$sigma, 5.0)
  expect_equal(sm$interpretation$category, "Excellent")
  expect_equal(sm$interpretation$defect_rate, 230)
})


test_that("sigma_metric interprets >= 4 as Good", {
  # Sigma = (10 - 2) / 2 = 4.0
  sm <- sigma_metric(bias = 2, cv = 2, tea = 10)

  expect_equal(sm$sigma, 4.0)
  expect_equal(sm$interpretation$category, "Good")
  expect_equal(sm$interpretation$defect_rate, 6210)
})


test_that("sigma_metric interprets >= 3 as Marginal", {
  # Sigma = (10 - 4) / 2 = 3.0
  sm <- sigma_metric(bias = 4, cv = 2, tea = 10)

  expect_equal(sm$sigma, 3.0)
  expect_equal(sm$interpretation$category, "Marginal")
  expect_equal(sm$interpretation$defect_rate, 66800)
})


test_that("sigma_metric interprets >= 2 as Poor", {
  # Sigma = (10 - 6) / 2 = 2.0
  sm <- sigma_metric(bias = 6, cv = 2, tea = 10)

  expect_equal(sm$sigma, 2.0)
  expect_equal(sm$interpretation$category, "Poor")
  expect_equal(sm$interpretation$defect_rate, 308500)
})


test_that("sigma_metric interprets >= 1 as Unacceptable", {
  # Sigma = (10 - 8) / 2 = 1.0
  sm <- sigma_metric(bias = 8, cv = 2, tea = 10)

  expect_equal(sm$sigma, 1.0)
  expect_equal(sm$interpretation$category, "Unacceptable")
  expect_equal(sm$interpretation$defect_rate, 690000)
})


test_that("sigma_metric handles sigma < 1", {
  # Sigma = (10 - 9) / 2 = 0.5
  sm <- sigma_metric(bias = 9, cv = 2, tea = 10)

  expect_equal(sm$sigma, 0.5)
  expect_equal(sm$interpretation$category, "Unacceptable")
  expect_true(is.na(sm$interpretation$defect_rate))
})


test_that("sigma_metric handles negative sigma", {
  # When bias > TEa, sigma becomes negative
  # Sigma = (10 - 15) / 2 = -2.5
  sm <- sigma_metric(bias = 15, cv = 2, tea = 10)

  expect_equal(sm$sigma, -2.5)
  expect_equal(sm$interpretation$category, "Unacceptable")
})


# Input validation ----

test_that("sigma_metric validates bias", {
  # Non-numeric
  expect_error(sigma_metric(bias = "1.5", cv = 2, tea = 10),
               "`bias` must be a single numeric value")

  # Vector
  expect_error(sigma_metric(bias = c(1.5, 2.0), cv = 2, tea = 10),
               "`bias` must be a single numeric value")

  # NA
  expect_error(sigma_metric(bias = NA, cv = 2, tea = 10),
               "`bias` must be a single numeric value")
})


test_that("sigma_metric validates cv", {
  # Non-numeric
  expect_error(sigma_metric(bias = 1.5, cv = "2", tea = 10),
               "`cv` must be a single numeric value")

  # Vector
  expect_error(sigma_metric(bias = 1.5, cv = c(2, 3), tea = 10),
               "`cv` must be a single numeric value")

  # NA
  expect_error(sigma_metric(bias = 1.5, cv = NA, tea = 10),
               "`cv` must be a single numeric value")

  # Zero
  expect_error(sigma_metric(bias = 1.5, cv = 0, tea = 10),
               "`cv` must be a positive number")

  # Negative
  expect_error(sigma_metric(bias = 1.5, cv = -2, tea = 10),
               "`cv` must be a positive number")
})


test_that("sigma_metric validates tea", {
  # Non-numeric
  expect_error(sigma_metric(bias = 1.5, cv = 2, tea = "10"),
               "`tea` must be a single numeric value")

  # Vector
  expect_error(sigma_metric(bias = 1.5, cv = 2, tea = c(10, 15)),
               "`tea` must be a single numeric value")

  # NA
  expect_error(sigma_metric(bias = 1.5, cv = 2, tea = NA),
               "`tea` must be a single numeric value")

  # Zero
  expect_error(sigma_metric(bias = 1.5, cv = 2, tea = 0),
               "`tea` must be a positive number")

  # Negative
  expect_error(sigma_metric(bias = 1.5, cv = 2, tea = -10),
               "`tea` must be a positive number")
})


# Edge cases ----

test_that("sigma_metric handles very small values", {
  sm <- sigma_metric(bias = 0.01, cv = 0.1, tea = 1)

  expect_equal(sm$sigma, (1 - 0.01) / 0.1)
})


test_that("sigma_metric handles very large values", {
  sm <- sigma_metric(bias = 10, cv = 50, tea = 100)

  expect_equal(sm$sigma, (100 - 10) / 50)
})


test_that("sigma_metric handles borderline category values", {
  # Just below 6
  sm_5_99 <- sigma_metric(bias = 0.02, cv = 1, tea = 6)
  expect_equal(sm_5_99$interpretation$category, "Excellent")

  # Just at 6
  sm_6 <- sigma_metric(bias = 0, cv = 1, tea = 6)
  expect_equal(sm_6$interpretation$category, "World Class")
})


# Integration with ate_from_bv ----

test_that("sigma_metric works with ate_from_bv output", {
  # Calculate TEa from biological variation
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  # Use TEa in sigma calculation
  sm <- sigma_metric(bias = 1.5, cv = 2.5, tea = ate$specifications$tea)

  expect_s3_class(sm, "sigma_metric")
  expect_true(is.numeric(sm$sigma))
})


# Print and summary methods ----

test_that("print.sigma_metric runs without error", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  expect_output(print(sm), "Six Sigma Metric")
  expect_output(print(sm), "Observed bias")
  expect_output(print(sm), "Observed CV")
  expect_output(print(sm), "Sigma:")
  expect_output(print(sm), "Performance:")
})


test_that("print.sigma_metric respects digits argument", {
  sm <- sigma_metric(bias = 1.555, cv = 2.333, tea = 10.999)

  # Default digits = 2
  output_default <- capture.output(print(sm))

  # Custom digits = 4
  output_custom <- capture.output(print(sm, digits = 4))

  expect_false(identical(output_default, output_custom))
})


test_that("print.sigma_metric returns object invisibly", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  result <- capture.output(returned <- print(sm))
  expect_identical(returned, sm)
})


test_that("summary.sigma_metric runs without error", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  expect_output(summary(sm), "Detailed Summary")
  expect_output(summary(sm), "Formula")
  expect_output(summary(sm), "Sigma Scale Reference")
  expect_output(summary(sm), "World Class")
  expect_output(summary(sm), "Current performance level")
})


test_that("summary.sigma_metric shows calculation steps", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  output <- capture.output(summary(sm))
  output_text <- paste(output, collapse = "\n")

  # Should show the formula breakdown
  expect_true(grepl("TEa - \\|Bias\\|", output_text))
})


test_that("summary.sigma_metric returns object invisibly", {
  sm <- sigma_metric(bias = 1.5, cv = 2.0, tea = 10)

  result <- capture.output(returned <- summary(sm))
  expect_identical(returned, sm)
})


test_that("summary.sigma_metric handles negative sigma", {
  sm <- sigma_metric(bias = 15, cv = 2, tea = 10)

  # Should not error
  expect_output(summary(sm), "Unacceptable")
})


# Real-world examples ----

test_that("sigma_metric produces reasonable values for typical lab data", {
  # Good performance: low bias, low CV, reasonable TEa
  sm_good <- sigma_metric(bias = 1, cv = 2, tea = 12)
  expect_true(sm_good$sigma > 4)

  # Marginal performance: higher bias and CV
  sm_marginal <- sigma_metric(bias = 3, cv = 3, tea = 12)
  expect_true(sm_marginal$sigma >= 3 && sm_marginal$sigma < 4)
})
