# Tests for ate_assessment()
# =============================================================================

# Basic functionality ----

test_that("ate_assessment returns correct class", {
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  expect_s3_class(assess, "ate_assessment")
  expect_s3_class(assess, "valytics_ate")
  expect_s3_class(assess, "valytics_result")
})


test_that("ate_assessment returns expected structure", {
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  # Check top-level names
  expect_named(assess, c("assessment", "observed", "specifications",
                         "sigma", "settings"))

  # Check assessment
  expect_named(assess$assessment, c("bias_acceptable", "cv_acceptable",
                                    "tea_acceptable", "overall"))

  # Check observed
  expect_named(assess$observed, c("bias", "cv", "te"))

  # Check specifications
  expect_named(assess$specifications, c("allowable_bias", "allowable_cv", "tea"))

  # Check sigma
  expect_named(assess$sigma, c("value", "category"))

  # Check settings
  expect_named(assess$settings, c("k"))
})


test_that("ate_assessment calculates observed TE correctly", {
  # TE = k * CV + |Bias| = 1.65 * 2.5 + 1.5 = 4.125 + 1.5 = 5.625
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  expect_equal(assess$observed$te, 1.65 * 2.5 + 1.5)
})


test_that("ate_assessment uses custom k value", {
  # TE = 2.0 * 2.5 + 1.5 = 5.0 + 1.5 = 6.5
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10, k = 2.0)

  expect_equal(assess$observed$te, 2.0 * 2.5 + 1.5)
  expect_equal(assess$settings$k, 2.0)
})


test_that("ate_assessment uses absolute value of bias for TE", {
  assess_pos <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)
  assess_neg <- ate_assessment(bias = -1.5, cv = 2.5, tea = 10)

  expect_equal(assess_pos$observed$te, assess_neg$observed$te)
})


test_that("ate_assessment stores input values correctly", {
  assess <- ate_assessment(
    bias = 1.5,
    cv = 2.5,
    tea = 10,
    allowable_bias = 3.0,
    allowable_cv = 4.0,
    k = 1.96
  )

  expect_equal(assess$observed$bias, 1.5)
  expect_equal(assess$observed$cv, 2.5)
  expect_equal(assess$specifications$tea, 10)
  expect_equal(assess$specifications$allowable_bias, 3.0)
  expect_equal(assess$specifications$allowable_cv, 4.0)
  expect_equal(assess$settings$k, 1.96)
})


# TEa assessment ----

test_that("ate_assessment TEa passes when TE <= TEa", {
  # TE = 1.65 * 2 + 1 = 4.3, TEa = 10
  assess <- ate_assessment(bias = 1, cv = 2, tea = 10)

  expect_true(assess$assessment$tea_acceptable)
  expect_true(assess$assessment$overall)
})


test_that("ate_assessment TEa fails when TE > TEa", {
  # TE = 1.65 * 5 + 3 = 11.25, TEa = 10
  assess <- ate_assessment(bias = 3, cv = 5, tea = 10)

  expect_false(assess$assessment$tea_acceptable)
  expect_false(assess$assessment$overall)
})


test_that("ate_assessment TEa passes at boundary", {
  # Set up so TE exactly equals TEa
  # TE = 1.65 * 2 + 1 = 4.3
  assess <- ate_assessment(bias = 1, cv = 2, tea = 4.3)

  expect_true(assess$assessment$tea_acceptable)
})


# Bias assessment ----

test_that("ate_assessment bias passes when |bias| <= allowable_bias", {
  assess <- ate_assessment(
    bias = 2.0, cv = 2.5, tea = 10,
    allowable_bias = 3.0
  )

  expect_true(assess$assessment$bias_acceptable)
})


test_that("ate_assessment bias fails when |bias| > allowable_bias", {
  assess <- ate_assessment(
    bias = 4.0, cv = 2.5, tea = 10,
    allowable_bias = 3.0
  )

  expect_false(assess$assessment$bias_acceptable)
})


test_that("ate_assessment bias uses absolute value", {
  assess_pos <- ate_assessment(
    bias = 2.0, cv = 2.5, tea = 10,
    allowable_bias = 3.0
  )
  assess_neg <- ate_assessment(
    bias = -2.0, cv = 2.5, tea = 10,
    allowable_bias = 3.0
  )

  expect_equal(assess_pos$assessment$bias_acceptable,
               assess_neg$assessment$bias_acceptable)
})


test_that("ate_assessment bias is NA when not specified", {
  assess <- ate_assessment(bias = 2.0, cv = 2.5, tea = 10)

  expect_true(is.na(assess$assessment$bias_acceptable))
})


# CV assessment ----

test_that("ate_assessment CV passes when cv <= allowable_cv", {
  assess <- ate_assessment(
    bias = 1.5, cv = 3.0, tea = 10,
    allowable_cv = 4.0
  )

  expect_true(assess$assessment$cv_acceptable)
})


test_that("ate_assessment CV fails when cv > allowable_cv", {
  assess <- ate_assessment(
    bias = 1.5, cv = 5.0, tea = 10,
    allowable_cv = 4.0
  )

  expect_false(assess$assessment$cv_acceptable)
})


test_that("ate_assessment CV is NA when not specified", {
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  expect_true(is.na(assess$assessment$cv_acceptable))
})


# Overall assessment ----

test_that("ate_assessment overall passes when all components pass", {
  assess <- ate_assessment(
    bias = 1.0, cv = 2.0, tea = 10,
    allowable_bias = 3.0,
    allowable_cv = 4.0
  )

  expect_true(assess$assessment$bias_acceptable)
  expect_true(assess$assessment$cv_acceptable)
  expect_true(assess$assessment$tea_acceptable)
  expect_true(assess$assessment$overall)
})


test_that("ate_assessment overall fails if bias fails", {
  assess <- ate_assessment(
    bias = 5.0, cv = 2.0, tea = 10,
    allowable_bias = 3.0,
    allowable_cv = 4.0
  )

  expect_false(assess$assessment$bias_acceptable)
  expect_true(assess$assessment$cv_acceptable)
  expect_false(assess$assessment$overall)
})


test_that("ate_assessment overall fails if cv fails", {
  assess <- ate_assessment(
    bias = 1.0, cv = 5.0, tea = 15,
    allowable_bias = 3.0,
    allowable_cv = 4.0
  )

  expect_true(assess$assessment$bias_acceptable)
  expect_false(assess$assessment$cv_acceptable)
  expect_false(assess$assessment$overall)
})


test_that("ate_assessment overall fails if tea fails", {
  assess <- ate_assessment(
    bias = 1.0, cv = 3.0, tea = 5,
    allowable_bias = 3.0,
    allowable_cv = 4.0
  )

  expect_true(assess$assessment$bias_acceptable)
  expect_true(assess$assessment$cv_acceptable)
  expect_false(assess$assessment$tea_acceptable)
  expect_false(assess$assessment$overall)
})


test_that("ate_assessment overall based on TEa when no component specs", {
  # Passing case
  assess_pass <- ate_assessment(bias = 1.0, cv = 2.0, tea = 10)
  expect_true(assess_pass$assessment$overall)

  # Failing case
  assess_fail <- ate_assessment(bias = 5.0, cv = 5.0, tea = 10)
  expect_false(assess_fail$assessment$overall)
})


# Sigma calculation ----

test_that("ate_assessment calculates sigma correctly", {
  # Sigma = (TEa - |Bias|) / CV = (10 - 1.5) / 2.5 = 3.4
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  expect_equal(assess$sigma$value, (10 - 1.5) / 2.5)
})


test_that("ate_assessment sigma category is correct", {
  # High sigma (>= 6)
  assess_high <- ate_assessment(bias = 1, cv = 1, tea = 10)
  expect_equal(assess_high$sigma$category, "World Class")

  # Low sigma (>= 2 but < 3) -> Poor
  # Sigma = (10 - 3) / 3 = 7/3 = 2.33
  assess_low <- ate_assessment(bias = 3, cv = 3, tea = 10)
  expect_equal(assess_low$sigma$category, "Poor")
})


# Input validation ----

test_that("ate_assessment validates bias", {
  expect_error(ate_assessment(bias = "1.5", cv = 2.5, tea = 10),
               "`bias` must be a single numeric value")
  expect_error(ate_assessment(bias = NA, cv = 2.5, tea = 10),
               "`bias` must be a single numeric value")
  expect_error(ate_assessment(bias = c(1, 2), cv = 2.5, tea = 10),
               "`bias` must be a single numeric value")
})


test_that("ate_assessment validates cv", {
  expect_error(ate_assessment(bias = 1.5, cv = "2.5", tea = 10),
               "`cv` must be a single numeric value")
  expect_error(ate_assessment(bias = 1.5, cv = NA, tea = 10),
               "`cv` must be a single numeric value")
  expect_error(ate_assessment(bias = 1.5, cv = 0, tea = 10),
               "`cv` must be a positive number")
  expect_error(ate_assessment(bias = 1.5, cv = -2.5, tea = 10),
               "`cv` must be a positive number")
})


test_that("ate_assessment validates tea", {
  expect_error(ate_assessment(bias = 1.5, cv = 2.5, tea = "10"),
               "`tea` must be a single numeric value")
  expect_error(ate_assessment(bias = 1.5, cv = 2.5, tea = NA),
               "`tea` must be a single numeric value")
  expect_error(ate_assessment(bias = 1.5, cv = 2.5, tea = 0),
               "`tea` must be a positive number")
  expect_error(ate_assessment(bias = 1.5, cv = 2.5, tea = -10),
               "`tea` must be a positive number")
})


test_that("ate_assessment validates allowable_bias", {
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, allowable_bias = "3"),
    "`allowable_bias` must be a single numeric value or NULL"
  )
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, allowable_bias = NA),
    "`allowable_bias` must be a single numeric value or NULL"
  )
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, allowable_bias = 0),
    "`allowable_bias` must be a positive number"
  )
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, allowable_bias = -3),
    "`allowable_bias` must be a positive number"
  )
})


test_that("ate_assessment validates allowable_cv", {
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, allowable_cv = "4"),
    "`allowable_cv` must be a single numeric value or NULL"
  )
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, allowable_cv = NA),
    "`allowable_cv` must be a single numeric value or NULL"
  )
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, allowable_cv = 0),
    "`allowable_cv` must be a positive number"
  )
})


test_that("ate_assessment validates k", {
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, k = "1.65"),
    "`k` must be a single numeric value"
  )
  expect_error(
    ate_assessment(bias = 1.5, cv = 2.5, tea = 10, k = 0),
    "`k` must be a positive number"
  )
})


# Integration with ate_from_bv ----

test_that("ate_assessment works with ate_from_bv output", {
  specs <- ate_from_bv(cvi = 5.6, cvg = 7.5)

  assess <- ate_assessment(
    bias = 1.5,
    cv = 2.5,
    tea = specs$specifications$tea,
    allowable_bias = specs$specifications$allowable_bias,
    allowable_cv = specs$specifications$allowable_cv
  )

  expect_s3_class(assess, "ate_assessment")
  expect_true(is.logical(assess$assessment$overall))
})


# Print and summary methods ----

test_that("print.ate_assessment runs without error", {
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  expect_output(print(assess), "Analytical Performance Assessment")
  expect_output(print(assess), "METHOD")
  expect_output(print(assess), "Sigma Metric")
})


test_that("print.ate_assessment shows ACCEPTABLE for passing method", {
  assess <- ate_assessment(bias = 1, cv = 2, tea = 10)

  expect_output(print(assess), "METHOD ACCEPTABLE")
})


test_that("print.ate_assessment shows NOT ACCEPTABLE for failing method", {
  assess <- ate_assessment(bias = 5, cv = 5, tea = 8)

  expect_output(print(assess), "NOT ACCEPTABLE")
})


test_that("print.ate_assessment respects digits argument", {
  assess <- ate_assessment(bias = 1.555, cv = 2.333, tea = 10.999)

  output_default <- capture.output(print(assess))
  output_custom <- capture.output(print(assess, digits = 4))

  expect_false(identical(output_default, output_custom))
})


test_that("print.ate_assessment returns object invisibly", {
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  result <- capture.output(returned <- print(assess))
  expect_identical(returned, assess)
})


test_that("summary.ate_assessment runs without error", {
  assess <- ate_assessment(
    bias = 1.5, cv = 2.5, tea = 10,
    allowable_bias = 3.0, allowable_cv = 4.0
  )

  expect_output(summary(assess), "Detailed Summary")
  expect_output(summary(assess), "Observed Performance")
  expect_output(summary(assess), "Allowable Specifications")
  expect_output(summary(assess), "Component Assessment")
  expect_output(summary(assess), "Sigma Metric")
  expect_output(summary(assess), "Interpretation")
})


test_that("summary.ate_assessment shows calculation breakdown", {
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  output <- capture.output(summary(assess))
  output_text <- paste(output, collapse = "\n")

  expect_true(grepl("TE =", output_text))
  expect_true(grepl("Sigma =", output_text))
})


test_that("summary.ate_assessment returns object invisibly", {
  assess <- ate_assessment(bias = 1.5, cv = 2.5, tea = 10)

  result <- capture.output(returned <- summary(assess))
  expect_identical(returned, assess)
})


test_that("summary.ate_assessment handles failing method with recommendations", {
  assess <- ate_assessment(
    bias = 5.0, cv = 6.0, tea = 10,
    allowable_bias = 3.0, allowable_cv = 4.0
  )

  output <- capture.output(summary(assess))
  output_text <- paste(output, collapse = "\n")

  expect_true(grepl("does not meet specifications", output_text))
})


# Edge cases ----

test_that("ate_assessment handles zero bias", {
  assess <- ate_assessment(bias = 0, cv = 2, tea = 10)

  expect_equal(assess$observed$te, 1.65 * 2 + 0)
  expect_true(assess$assessment$tea_acceptable)
})


test_that("ate_assessment handles borderline sigma values", {
  # Sigma exactly at 4
  assess <- ate_assessment(bias = 2, cv = 2, tea = 10)
  # Sigma = (10 - 2) / 2 = 4

  expect_equal(assess$sigma$value, 4)
  expect_equal(assess$sigma$category, "Good")
})
