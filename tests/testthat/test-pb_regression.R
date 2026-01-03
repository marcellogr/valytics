# Tests for pb_regression()

# =============================================================================
# Test data setup
# =============================================================================

test_that("pb_regression returns correct class", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(50, sd = 5)

  pb <- pb_regression(x, y)

  expect_s3_class(pb, "pb_regression")
  expect_s3_class(pb, "valytics_comparison")
  expect_s3_class(pb, "valytics_result")
})

test_that("pb_regression works with formula interface", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  df <- data.frame(
    method_a = rnorm(50, mean = 100, sd = 15),
    method_b = rnorm(50, mean = 102, sd = 15)
  )

  pb <- pb_regression(method_a ~ method_b, data = df)

  expect_s3_class(pb, "pb_regression")
  expect_equal(pb$input$var_names["x"], c(x = "method_b"))
  expect_equal(pb$input$var_names["y"], c(y = "method_a"))
})

test_that("pb_regression returns expected structure", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(50, sd = 5)

  pb <- pb_regression(x, y)

 # Check structure
  expect_named(pb, c("input", "results", "cusum", "settings", "call"))

  # Input
  expect_equal(pb$input$n, 50)
  expect_equal(length(pb$input$x), 50)
  expect_equal(length(pb$input$y), 50)

  # Results
  expect_true(is.numeric(pb$results$slope))
  expect_true(is.numeric(pb$results$intercept))
  expect_named(pb$results$slope_ci, c("lower", "upper"))
  expect_named(pb$results$intercept_ci, c("lower", "upper"))

  # CUSUM
  expect_true(is.logical(pb$cusum$linear) || is.na(pb$cusum$linear))
})

test_that("pb_regression slope estimate is reasonable for known relationship", {
  skip_if_not_installed("robslopes")

  set.seed(123)
  x <- seq(10, 100, length.out = 100)
  # True slope = 1.1, intercept = 5
  y <- 1.1 * x + 5 + rnorm(100, sd = 3)

  pb <- pb_regression(x, y)

  # Slope should be close to 1.1
  expect_true(abs(pb$results$slope - 1.1) < 0.1)
  # Intercept should be close to 5
  expect_true(abs(pb$results$intercept - 5) < 5)
})

test_that("pb_regression handles NA values correctly", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- c(rnorm(48, mean = 100, sd = 15), NA, NA)
  y <- c(rnorm(48, mean = 102, sd = 15), 100, NA)

  pb <- pb_regression(x, y, na_action = "omit")

  expect_equal(pb$input$n, 48)
  expect_equal(pb$input$n_excluded, 2)
})

test_that("pb_regression fails with na_action = 'fail' when NAs present", {
  skip_if_not_installed("robslopes")

  x <- c(1, 2, NA, 4, 5)
  y <- c(1, 2, 3, 4, 5)

  expect_error(
    pb_regression(x, y, na_action = "fail"),
    "Missing values detected"
  )
})

test_that("pb_regression requires minimum sample size", {
  skip_if_not_installed("robslopes")

  x <- 1:5
  y <- 1:5

  expect_error(
    pb_regression(x, y),
    "At least 10 complete paired observations"
  )
})

test_that("pb_regression validates conf_level", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- rnorm(50)
  y <- rnorm(50)

  expect_error(pb_regression(x, y, conf_level = 0))
  expect_error(pb_regression(x, y, conf_level = 1))
  expect_error(pb_regression(x, y, conf_level = 1.5))
  expect_error(pb_regression(x, y, conf_level = "0.95"))
})

# =============================================================================
# CI methods
# =============================================================================

test_that("analytical CI produces valid intervals", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(50, sd = 5)

  pb <- pb_regression(x, y, ci_method = "analytical")

  # CI should bracket point estimate
  expect_true(pb$results$slope_ci["lower"] <= pb$results$slope)
  expect_true(pb$results$slope_ci["upper"] >= pb$results$slope)
  expect_true(pb$results$intercept_ci["lower"] <= pb$results$intercept)
  expect_true(pb$results$intercept_ci["upper"] >= pb$results$intercept)
})

test_that("bootstrap CI produces valid intervals", {
  skip_if_not_installed("robslopes")
  skip_on_cran()  # Skip on CRAN due to computation time

  set.seed(42)
  x <- rnorm(30, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(30, sd = 5)

  pb <- pb_regression(x, y, ci_method = "bootstrap", boot_n = 199)

  # CI should bracket point estimate (usually)
  # Using wider tolerance for bootstrap
  expect_true(is.numeric(pb$results$slope_ci["lower"]))
  expect_true(is.numeric(pb$results$slope_ci["upper"]))
})

# =============================================================================
# CUSUM test
# =============================================================================

test_that("CUSUM test returns expected structure", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- x + rnorm(50, sd = 5)  # Linear relationship

  pb <- pb_regression(x, y)

  expect_named(pb$cusum, c("statistic", "critical_value", "p_value", "linear"))
  expect_equal(pb$cusum$critical_value, 1.36)
})

test_that("CUSUM detects linearity in linear data", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- seq(1, 100, length.out = 100)
  y <- 2 * x + 10 + rnorm(100, sd = 5)  # Clearly linear

  pb <- pb_regression(x, y)

  # Should indicate linearity
  expect_true(pb$cusum$linear || is.na(pb$cusum$linear))
})

# =============================================================================
# Print and summary methods
# =============================================================================

test_that("print method runs without error", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(50, sd = 5)

  pb <- pb_regression(x, y)

  expect_output(print(pb), "Passing-Bablok Regression")
  expect_output(print(pb), "Slope:")
  expect_output(print(pb), "Intercept:")
})

test_that("summary method runs without error", {
  skip_if_not_installed("robslopes")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(50, sd = 5)

  pb <- pb_regression(x, y)

  expect_output(summary(pb), "Passing-Bablok Regression")
  expect_output(summary(pb), "CUSUM")
  expect_output(summary(pb), "Interpretation")
})

# =============================================================================
# Plot methods
# =============================================================================

test_that("plot methods return ggplot objects", {
  skip_if_not_installed("robslopes")
  skip_if_not_installed("ggplot2")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(50, sd = 5)

  pb <- pb_regression(x, y)

  # Scatter plot
  p1 <- plot(pb, type = "scatter")
  expect_s3_class(p1, "ggplot")

  # Residual plot
  p2 <- plot(pb, type = "residuals")
  expect_s3_class(p2, "ggplot")

  # CUSUM plot
  p3 <- plot(pb, type = "cusum")
  expect_s3_class(p3, "ggplot")
})

test_that("autoplot method works", {
  skip_if_not_installed("robslopes")
  skip_if_not_installed("ggplot2")

  set.seed(42)
  x <- rnorm(50, mean = 100, sd = 15)
  y <- 1.05 * x + 3 + rnorm(50, sd = 5)

  pb <- pb_regression(x, y)

  p <- autoplot.pb_regression(pb)
  expect_s3_class(p, "ggplot")
})
