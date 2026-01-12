# Tests for ate_from_bv()
# =============================================================================

# Basic functionality ----

test_that("ate_from_bv returns correct class", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  expect_s3_class(ate, "ate_specs")
  expect_s3_class(ate, "valytics_ate")
  expect_s3_class(ate, "valytics_result")
})


test_that("ate_from_bv returns expected structure", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  # Check top-level names
 expect_named(ate, c("specifications", "input", "multipliers"))

  # Check specifications
  expect_named(ate$specifications, c("allowable_cv", "allowable_bias", "tea"))

  # Check input
  expect_named(ate$input, c("cvi", "cvg", "level", "k"))

  # Check multipliers
  expect_named(ate$multipliers, c("imprecision", "bias"))
})


test_that("ate_from_bv calculates desirable level correctly", {
  # Using known values for verification
  # CV_I = 5.6, CV_G = 7.5
  # Desirable CV_A = 0.50 * 5.6 = 2.8
  # Total BV = sqrt(5.6^2 + 7.5^2) = sqrt(31.36 + 56.25) = sqrt(87.61) = 9.36
 # Desirable Bias = 0.25 * 9.36 = 2.34
  # TEa = 1.65 * 2.8 + 2.34 = 4.62 + 2.34 = 6.96

  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "desirable")

  expect_equal(ate$specifications$allowable_cv, 0.50 * 5.6)
  expect_equal(ate$specifications$allowable_bias,
               0.25 * sqrt(5.6^2 + 7.5^2), tolerance = 1e-10)
  expect_equal(ate$specifications$tea,
               1.65 * (0.50 * 5.6) + 0.25 * sqrt(5.6^2 + 7.5^2),
               tolerance = 1e-10)
})


test_that("ate_from_bv calculates optimal level correctly", {
  # Optimal = 0.5 * desirable multipliers
  # Optimal CV_A = 0.25 * CV_I
  # Optimal Bias = 0.125 * sqrt(CV_I^2 + CV_G^2)

  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "optimal")

  expect_equal(ate$specifications$allowable_cv, 0.25 * 5.6)
  expect_equal(ate$specifications$allowable_bias,
               0.125 * sqrt(5.6^2 + 7.5^2), tolerance = 1e-10)

  # Check multipliers
  expect_equal(ate$multipliers$imprecision, 0.25)
  expect_equal(ate$multipliers$bias, 0.125)
})


test_that("ate_from_bv calculates minimum level correctly", {
  # Minimum = 1.5 * desirable multipliers
  # Minimum CV_A = 0.75 * CV_I
  # Minimum Bias = 0.375 * sqrt(CV_I^2 + CV_G^2)

  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "minimum")

  expect_equal(ate$specifications$allowable_cv, 0.75 * 5.6)
  expect_equal(ate$specifications$allowable_bias,
               0.375 * sqrt(5.6^2 + 7.5^2), tolerance = 1e-10)

  # Check multipliers
  expect_equal(ate$multipliers$imprecision, 0.75)
  expect_equal(ate$multipliers$bias, 0.375)
})


test_that("ate_from_bv handles cvg = NULL correctly", {
  ate <- ate_from_bv(cvi = 5.6, cvg = NULL)

  # CV_A should be calculated
  expect_equal(ate$specifications$allowable_cv, 0.50 * 5.6)

  # Bias and TEa should be NULL
  expect_null(ate$specifications$allowable_bias)
  expect_null(ate$specifications$tea)

  # Input should reflect NULL cvg
  expect_null(ate$input$cvg)
})


test_that("ate_from_bv respects custom k value", {
  # Default k = 1.65
  ate_default <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  # Custom k = 2.0
  ate_custom <- ate_from_bv(cvi = 5.6, cvg = 7.5, k = 2.0)

  # CV_A and Bias should be the same
  expect_equal(ate_default$specifications$allowable_cv,
               ate_custom$specifications$allowable_cv)
  expect_equal(ate_default$specifications$allowable_bias,
               ate_custom$specifications$allowable_bias)

  # TEa should differ
  expect_true(ate_custom$specifications$tea > ate_default$specifications$tea)

  # Verify formula: TEa = k * CV_A + Bias
  expect_equal(ate_custom$specifications$tea,
               2.0 * ate_custom$specifications$allowable_cv +
                 ate_custom$specifications$allowable_bias)
})


test_that("ate_from_bv stores input values correctly", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "minimum", k = 2.0)

  expect_equal(ate$input$cvi, 5.6)
  expect_equal(ate$input$cvg, 7.5)
  expect_equal(ate$input$level, "minimum")
  expect_equal(ate$input$k, 2.0)
})


# Input validation ----

test_that("ate_from_bv validates cvi", {
  # Non-numeric
  expect_error(ate_from_bv(cvi = "5.6"), "`cvi` must be a single numeric value")

  # Vector
  expect_error(ate_from_bv(cvi = c(5.6, 6.0)), "`cvi` must be a single numeric value")

  # NA
 expect_error(ate_from_bv(cvi = NA), "`cvi` must be a single numeric value")

  # Zero
  expect_error(ate_from_bv(cvi = 0), "`cvi` must be a positive number")

  # Negative
  expect_error(ate_from_bv(cvi = -5.6), "`cvi` must be a positive number")
})


test_that("ate_from_bv validates cvg", {
  # Non-numeric
  expect_error(ate_from_bv(cvi = 5.6, cvg = "7.5"),
               "`cvg` must be a single numeric value or NULL")

  # Vector
  expect_error(ate_from_bv(cvi = 5.6, cvg = c(7.5, 8.0)),
               "`cvg` must be a single numeric value or NULL")

  # NA
  expect_error(ate_from_bv(cvi = 5.6, cvg = NA),
               "`cvg` must be a single numeric value or NULL")

  # Zero
  expect_error(ate_from_bv(cvi = 5.6, cvg = 0),
               "`cvg` must be a positive number")

  # Negative
  expect_error(ate_from_bv(cvi = 5.6, cvg = -7.5),
               "`cvg` must be a positive number")
})


test_that("ate_from_bv validates k", {
  # Non-numeric
  expect_error(ate_from_bv(cvi = 5.6, k = "1.65"),
               "`k` must be a single numeric value")

  # Vector
  expect_error(ate_from_bv(cvi = 5.6, k = c(1.65, 2.0)),
               "`k` must be a single numeric value")

  # NA
  expect_error(ate_from_bv(cvi = 5.6, k = NA),
               "`k` must be a single numeric value")

  # Zero
  expect_error(ate_from_bv(cvi = 5.6, k = 0),
               "`k` must be a positive number")

  # Negative
  expect_error(ate_from_bv(cvi = 5.6, k = -1.65),
               "`k` must be a positive number")
})


test_that("ate_from_bv validates level", {
  # Invalid level
  expect_error(ate_from_bv(cvi = 5.6, level = "extreme"))
  expect_error(ate_from_bv(cvi = 5.6, level = "Desirable"))  # case-sensitive
})


# Hierarchy validation ----

test_that("optimal < desirable < minimum for all specifications", {
  ate_opt <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "optimal")
  ate_des <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "desirable")
  ate_min <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "minimum")

  # CV_A: optimal < desirable < minimum
  expect_lt(ate_opt$specifications$allowable_cv,
            ate_des$specifications$allowable_cv)
  expect_lt(ate_des$specifications$allowable_cv,
            ate_min$specifications$allowable_cv)

  # Bias: optimal < desirable < minimum
  expect_lt(ate_opt$specifications$allowable_bias,
            ate_des$specifications$allowable_bias)
  expect_lt(ate_des$specifications$allowable_bias,
            ate_min$specifications$allowable_bias)

  # TEa: optimal < desirable < minimum
  expect_lt(ate_opt$specifications$tea, ate_des$specifications$tea)
  expect_lt(ate_des$specifications$tea, ate_min$specifications$tea)
})


# Edge cases ----

test_that("ate_from_bv handles very small CV values", {
  ate <- ate_from_bv(cvi = 0.1, cvg = 0.2)

  expect_true(ate$specifications$allowable_cv > 0)
  expect_true(ate$specifications$allowable_bias > 0)
  expect_true(ate$specifications$tea > 0)
})


test_that("ate_from_bv handles very large CV values", {
  ate <- ate_from_bv(cvi = 50, cvg = 100)

  expect_true(ate$specifications$allowable_cv > 0)
  expect_true(ate$specifications$allowable_bias > 0)
  expect_true(ate$specifications$tea > 0)
})


test_that("ate_from_bv handles equal cvi and cvg", {
  ate <- ate_from_bv(cvi = 5.0, cvg = 5.0)

  # Total BV = sqrt(5^2 + 5^2) = sqrt(50) = 7.07
  total_bv <- sqrt(5^2 + 5^2)
  expect_equal(ate$specifications$allowable_bias, 0.25 * total_bv, tolerance = 1e-10)
})


test_that("ate_from_bv handles cvi >> cvg", {
  ate <- ate_from_bv(cvi = 20, cvg = 2)

  # Total BV dominated by cvi
  total_bv <- sqrt(20^2 + 2^2)
  expect_equal(ate$specifications$allowable_bias, 0.25 * total_bv, tolerance = 1e-10)
})


test_that("ate_from_bv handles cvg >> cvi", {
  ate <- ate_from_bv(cvi = 2, cvg = 20)

  # Total BV dominated by cvg
  total_bv <- sqrt(2^2 + 20^2)
  expect_equal(ate$specifications$allowable_bias, 0.25 * total_bv, tolerance = 1e-10)
})


# Print and summary methods ----

test_that("print.ate_specs runs without error", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  expect_output(print(ate), "Analytical Performance Specifications")
  expect_output(print(ate), "Within-subject CV")
  expect_output(print(ate), "Allowable imprecision")
  expect_output(print(ate), "Total allowable error")
})


test_that("print.ate_specs handles cvg = NULL", {
  ate <- ate_from_bv(cvi = 5.6)

  expect_output(print(ate), "Between-subject CV.*not provided")
  expect_output(print(ate), "Allowable bias.*requires CV_G")
})


test_that("print.ate_specs respects digits argument", {
  ate <- ate_from_bv(cvi = 5.666, cvg = 7.555)

  # Default digits = 2
  output_default <- capture.output(print(ate))

  # Custom digits = 4
  output_custom <- capture.output(print(ate, digits = 4))

  # Output should differ in precision
  expect_false(identical(output_default, output_custom))
})


test_that("print.ate_specs returns object invisibly", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  result <- capture.output(returned <- print(ate))
  expect_identical(returned, ate)
})


test_that("summary.ate_specs runs without error", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  expect_output(summary(ate), "Detailed Summary")
  expect_output(summary(ate), "Biological Variation Data")
  expect_output(summary(ate), "Fraser & Petersen")
  expect_output(summary(ate), "Comparison Across Performance Levels")
  expect_output(summary(ate), "biologicalvariation.eu")
})


test_that("summary.ate_specs shows all three levels in comparison", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "desirable")

  output <- capture.output(summary(ate))
  output_text <- paste(output, collapse = "\n")

  expect_true(grepl("optimal", output_text))
  expect_true(grepl("desirable", output_text))
  expect_true(grepl("minimum", output_text))
})


test_that("summary.ate_specs marks selected level", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5, level = "optimal")

  output <- capture.output(summary(ate))
  output_text <- paste(output, collapse = "\n")

  # Selected level should have asterisk
  expect_true(grepl("optimal \\*", output_text))
})


test_that("summary.ate_specs returns object invisibly", {
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  result <- capture.output(returned <- summary(ate))

  expect_s3_class(returned, "summary.ate_specs")
})


# Real-world example values ----

test_that("ate_from_bv produces reasonable values for glucose",
{
  # Glucose BV from EFLM database (approximate values for illustration)
  # CV_I ~ 5.6%, CV_G ~ 7.5%
  ate <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  # Desirable specifications should be reasonable for glucose
  # CV_A around 2-3%, Bias around 2-3%, TEa around 6-8%
  expect_true(ate$specifications$allowable_cv > 2 &&
                ate$specifications$allowable_cv < 4)
  expect_true(ate$specifications$allowable_bias > 1 &&
                ate$specifications$allowable_bias < 4)
  expect_true(ate$specifications$tea > 5 && ate$specifications$tea < 10)
})


test_that("ate_from_bv produces reasonable values for creatinine", {
  # Creatinine BV (approximate values for illustration)
  # CV_I ~ 5.95%, CV_G ~ 14.7%
  ate <- ate_from_bv(cvi = 5.95, cvg = 14.7)

  # TEa should reflect higher total BV
  expect_true(ate$specifications$tea > 8)
})
