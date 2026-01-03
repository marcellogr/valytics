# Tests for ba_analysis() ----

test_that("ba_analysis works with basic vector input", {
  set.seed(123)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- x + rnorm(50, mean = 2, sd = 5)  # +2 bias

  result <- ba_analysis(x, y)

  # Check class structure

  expect_s3_class(result, "ba_analysis")
  expect_s3_class(result, "valytics_comparison")
  expect_s3_class(result, "valytics_result")

  # Check structure components
  expect_named(result, c("input", "results", "settings", "call"))
  expect_named(result$input, c("x", "y", "n", "n_excluded", "var_names"))
  expect_named(result$results, c("differences", "averages", "bias", "bias_se",
                                  "bias_ci", "sd_diff", "loa_lower", "loa_upper",
                                  "loa_lower_ci", "loa_upper_ci"))

  # Check basic values
  expect_equal(result$input$n, 50)
  expect_equal(result$input$n_excluded, 0)
  expect_equal(length(result$results$differences), 50)
  expect_equal(length(result$results$averages), 50)

  # Bias should be close to 2 (the true bias we added)
  expect_true(abs(result$results$bias - 2) < 2)  # Within 2 units
})


test_that("ba_analysis works with formula interface", {
  set.seed(456)
  df <- data.frame(
    method_a = rnorm(30, 50, 10),
    method_b = rnorm(30, 52, 10)
  )

  result <- ba_analysis(method_a ~ method_b, data = df)

  expect_s3_class(result, "ba_analysis")
  expect_equal(result$input$n, 30)
  expect_equal(result$input$var_names[["x"]], "method_a")
  expect_equal(result$input$var_names[["y"]], "method_b")
})


test_that("ba_analysis handles percent differences correctly", {
  x <- c(100, 200, 300, 400, 500)
  y <- c(105, 210, 315, 420, 525)  # 5% higher

  result <- ba_analysis(x, y, type = "percent")

  expect_equal(result$settings$type, "percent")
  # Differences should be around 5%
  expect_true(all(abs(result$results$differences - 5) < 1))
})


test_that("ba_analysis handles NA values correctly with na_action = 'omit'", {
  x <- c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10)
  y <- c(1.1, 2.1, 3.1, NA, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1)

  result <- ba_analysis(x, y, na_action = "omit")

  expect_equal(result$input$n, 8)
  expect_equal(result$input$n_excluded, 2)
})


test_that("ba_analysis fails with na_action = 'fail' when NAs present", {
  x <- c(1, 2, NA, 4, 5)
  y <- c(1.1, 2.1, 3.1, 4.1, 5.1)

  expect_error(
    ba_analysis(x, y, na_action = "fail"),
    "Missing values detected"
  )
})


test_that("ba_analysis validates input correctly", {
  # Non-numeric input
  expect_error(ba_analysis("a", "b"), "`x` must be a numeric vector")

  # Unequal lengths
  expect_error(
    ba_analysis(1:5, 1:3),
    "`x` and `y` must have the same length"
  )

  # Invalid conf_level
  expect_error(
    ba_analysis(1:10, 1:10, conf_level = 1.5),
    "`conf_level` must be a single number between 0 and 1"
  )
  expect_error(
    ba_analysis(1:10, 1:10, conf_level = 0),
    "`conf_level` must be a single number between 0 and 1"
  )

  # Too few observations
  expect_error(
    ba_analysis(1:2, 1:2),
    "At least 3 complete paired observations are required"
  )

  # Missing y when not using formula
  expect_error(
    ba_analysis(1:10),
    "Either provide a formula or both `x` and `y` vectors"
  )
})


test_that("ba_analysis confidence intervals contain point estimates", {
  set.seed(789)
  x <- rnorm(100, 50, 10)
  y <- x + rnorm(100, 0, 3)

  result <- ba_analysis(x, y)

  # Bias CI should contain bias
  expect_true(result$results$bias >= result$results$bias_ci["lower"])
  expect_true(result$results$bias <= result$results$bias_ci["upper"])

  # LoA CIs should contain LoA values
  expect_true(result$results$loa_lower >= result$results$loa_lower_ci["lower"])
  expect_true(result$results$loa_lower <= result$results$loa_lower_ci["upper"])

  expect_true(result$results$loa_upper >= result$results$loa_upper_ci["lower"])
  expect_true(result$results$loa_upper <= result$results$loa_upper_ci["upper"])
})


test_that("ba_analysis produces symmetric LoA when bias is near zero", {
  set.seed(101)
  x <- rnorm(100, 100, 20)
  y <- x + rnorm(100, 0, 5)  # No systematic bias

  result <- ba_analysis(x, y)

  # LoA should be roughly symmetric around bias
  distance_lower <- result$results$bias - result$results$loa_lower
  distance_upper <- result$results$loa_upper - result$results$bias

  expect_equal(distance_lower, distance_upper, tolerance = 0.001)
})


test_that("ba_analysis handles different confidence levels", {
  set.seed(202)
  x <- rnorm(50, 100, 10)
  y <- x + rnorm(50, 5, 3)

  result_95 <- ba_analysis(x, y, conf_level = 0.95)
  result_99 <- ba_analysis(x, y, conf_level = 0.99)
  result_90 <- ba_analysis(x, y, conf_level = 0.90)

  # Point estimates should be identical
  expect_equal(result_95$results$bias, result_99$results$bias)
  expect_equal(result_95$results$bias, result_90$results$bias)

  # 99% CI should be wider than 95%, which should be wider than 90%
  width_90 <- diff(result_90$results$bias_ci)
  width_95 <- diff(result_95$results$bias_ci)
  width_99 <- diff(result_99$results$bias_ci)

  expect_true(width_99 > width_95)
  expect_true(width_95 > width_90)
})


# Tests for print and summary methods
# =============================================================================

test_that("print.ba_analysis runs without error", {
  set.seed(303)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 2)

  result <- ba_analysis(x, y)

  expect_output(print(result), "Bland-Altman Analysis")
  expect_output(print(result), "Bias")
  expect_output(print(result), "Limits of Agreement")
})


test_that("summary.ba_analysis produces correct structure", {
  set.seed(404)
  x <- rnorm(40, 100, 15)
  y <- x + rnorm(40, 2, 4)

  result <- ba_analysis(x, y)
  summ <- summary(result)

  expect_s3_class(summ, "summary.ba_analysis")
  expect_named(summ, c("call", "n", "n_excluded", "var_names", "type",
                       "conf_level", "descriptives", "agreement",
                       "normality_test", "sd_diff"))

  # Check descriptives data frame
  expect_s3_class(summ$descriptives, "data.frame")
  expect_equal(nrow(summ$descriptives), 2)

  # Check agreement data frame
  expect_s3_class(summ$agreement, "data.frame")
  expect_equal(nrow(summ$agreement), 3)

  # Normality test should be present
  expect_true(!is.null(summ$normality_test))
})


test_that("print.summary.ba_analysis runs without error", {
  set.seed(505)
  x <- rnorm(25, 80, 12)
  y <- x + rnorm(25, 0, 3)

  result <- ba_analysis(x, y)
  summ <- summary(result)

  expect_output(print(summ), "Bland-Altman Analysis - Detailed Summary")
  expect_output(print(summ), "Descriptive Statistics")
  expect_output(print(summ), "Agreement Statistics")
  expect_output(print(summ), "Shapiro-Wilk")
})


# Edge cases
# =============================================================================

test_that("ba_analysis handles perfect agreement",
{
  x <- 1:20
  y <- 1:20  # Perfect agreement

  result <- ba_analysis(x, y)

  expect_equal(result$results$bias, 0)
  expect_equal(result$results$sd_diff, 0)
  expect_equal(result$results$loa_lower, 0)
  expect_equal(result$results$loa_upper, 0)
})


test_that("ba_analysis handles constant bias", {
  x <- 1:30
  y <- x + 5  # Constant +5 bias

  result <- ba_analysis(x, y)

  expect_equal(result$results$bias, 5)
  expect_equal(result$results$sd_diff, 0)
})


test_that("ba_analysis handles minimum sample size (n=3)", {
  x <- c(10, 20, 30)
  y <- c(11, 21, 31)

  result <- ba_analysis(x, y)

  expect_s3_class(result, "ba_analysis")
  expect_equal(result$input$n, 3)
})
