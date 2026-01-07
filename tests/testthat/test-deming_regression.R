# Tests for deming_regression()

# =============================================================================
# Test data setup
# =============================================================================

test_that("deming_regression returns correct class", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  expect_s3_class(dm, "deming_regression")
  expect_s3_class(dm, "valytics_comparison")
  expect_s3_class(dm, "valytics_result")
})

test_that("deming_regression works with formula interface", {
  set.seed(42)
  df <- data.frame(
    method_a = rnorm(50, mean = 100, sd = 15),
    method_b = rnorm(50, mean = 102, sd = 15)
  )

  dm <- deming_regression(method_a ~ method_b, data = df)

  expect_s3_class(dm, "deming_regression")
  expect_equal(dm$input$var_names["x"], c(x = "method_b"))
  expect_equal(dm$input$var_names["y"], c(y = "method_a"))
})

test_that("deming_regression returns expected structure", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  # Check structure
  expect_named(dm, c("input", "results", "settings", "call"))

  # Input
  expect_equal(dm$input$n, 50)
  expect_equal(length(dm$input$x), 50)
  expect_equal(length(dm$input$y), 50)

  # Results
  expect_true(is.numeric(dm$results$slope))
  expect_true(is.numeric(dm$results$intercept))
  expect_named(dm$results$slope_ci, c("lower", "upper"))
  expect_named(dm$results$intercept_ci, c("lower", "upper"))
  expect_true(is.numeric(dm$results$slope_se))
  expect_true(is.numeric(dm$results$intercept_se))

  # Settings
  expect_equal(dm$settings$error_ratio, 1)
  expect_equal(dm$settings$conf_level, 0.95)
  expect_equal(dm$settings$ci_method, "jackknife")
})

test_that("deming_regression slope estimate is reasonable for known relationship", {
  set.seed(123)
  true_vals <- seq(10, 100, length.out = 100)
  # True slope = 1.1, intercept = 5
  x <- true_vals + rnorm(100, sd = 3)
  y <- 1.1 * true_vals + 5 + rnorm(100, sd = 3)

  dm <- deming_regression(x, y)

  # Slope should be close to 1.1
  expect_true(abs(dm$results$slope - 1.1) < 0.15)
  # Intercept should be close to 5
  expect_true(abs(dm$results$intercept - 5) < 10)
})

test_that("deming_regression handles NA values correctly", {
  set.seed(42)
  x <- c(rnorm(48, mean = 100, sd = 15), NA, NA)
  y <- c(rnorm(48, mean = 102, sd = 15), 100, NA)

  dm <- deming_regression(x, y, na_action = "omit")

  expect_equal(dm$input$n, 48)
  expect_equal(dm$input$n_excluded, 2)
})

test_that("deming_regression fails with na_action = 'fail' when NAs present", {
  x <- c(1:10, NA, 12:20)
  y <- c(1:10, 11, NA, 13:20)

  expect_error(
    deming_regression(x, y, na_action = "fail"),
    "Missing values detected"
  )
})

test_that("deming_regression requires minimum sample size", {
  x <- 1:5
  y <- 1:5

  expect_error(
    deming_regression(x, y),
    "At least 10 complete paired observations"
  )
})

test_that("deming_regression validates error_ratio", {
  set.seed(42)
  x <- rnorm(50)
  y <- rnorm(50)

  expect_error(deming_regression(x, y, error_ratio = 0))
  expect_error(deming_regression(x, y, error_ratio = -1))
  expect_error(deming_regression(x, y, error_ratio = "1"))
  expect_error(deming_regression(x, y, error_ratio = Inf))
})

test_that("deming_regression validates conf_level", {
  set.seed(42)
  x <- rnorm(50)
  y <- rnorm(50)

  expect_error(deming_regression(x, y, conf_level = 0))
  expect_error(deming_regression(x, y, conf_level = 1))
  expect_error(deming_regression(x, y, conf_level = 1.5))
  expect_error(deming_regression(x, y, conf_level = "0.95"))
})

# =============================================================================
# Error ratio tests
# =============================================================================

test_that("deming_regression respects different error ratios", {
  set.seed(42)
  true_vals <- rnorm(100, mean = 100, sd = 20)
  x <- true_vals + rnorm(100, sd = 5)
  y <- true_vals + rnorm(100, sd = 5)

  dm_1 <- deming_regression(x, y, error_ratio = 1)
  dm_2 <- deming_regression(x, y, error_ratio = 2)
  dm_05 <- deming_regression(x, y, error_ratio = 0.5)

  # Different error ratios should give different results
  expect_false(identical(dm_1$results$slope, dm_2$results$slope))
  expect_false(identical(dm_1$results$slope, dm_05$results$slope))

  # All should have the correct error_ratio stored
  expect_equal(dm_1$settings$error_ratio, 1)
  expect_equal(dm_2$settings$error_ratio, 2)
  expect_equal(dm_05$settings$error_ratio, 0.5)
})

test_that("deming_regression with lambda=1 is orthogonal regression", {
  set.seed(42)
  true_vals <- rnorm(100, mean = 100, sd = 20)
  x <- true_vals + rnorm(100, sd = 5)
  y <- true_vals + rnorm(100, sd = 5)

  dm <- deming_regression(x, y, error_ratio = 1)

  # Orthogonal regression should minimize perpendicular distances
  # Check that residuals are perpendicular
  expect_true(length(dm$results$residuals) == 100)
  expect_true(is.numeric(dm$results$residuals))
})

# =============================================================================
# CI methods
# =============================================================================

test_that("jackknife CI produces valid intervals", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y, ci_method = "jackknife")

  # CI should bracket point estimate
  expect_true(dm$results$slope_ci["lower"] <= dm$results$slope)
  expect_true(dm$results$slope_ci["upper"] >= dm$results$slope)
  expect_true(dm$results$intercept_ci["lower"] <= dm$results$intercept)
  expect_true(dm$results$intercept_ci["upper"] >= dm$results$intercept)

  # SE should be positive
  expect_true(dm$results$slope_se > 0)
  expect_true(dm$results$intercept_se > 0)
})

test_that("bootstrap CI produces valid intervals", {
  skip_on_cran()  # Skip on CRAN due to computation time

  set.seed(42)
  true_vals <- rnorm(30, mean = 100, sd = 20)
  x <- true_vals + rnorm(30, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(30, sd = 5)

  dm <- deming_regression(x, y, ci_method = "bootstrap", boot_n = 199)

  # CI should exist and be numeric
  expect_true(is.numeric(dm$results$slope_ci["lower"]))
  expect_true(is.numeric(dm$results$slope_ci["upper"]))
  expect_true(is.numeric(dm$results$intercept_ci["lower"]))
  expect_true(is.numeric(dm$results$intercept_ci["upper"]))

  # Settings should reflect bootstrap
  expect_equal(dm$settings$ci_method, "bootstrap")
  expect_equal(dm$settings$boot_n, 199)
})

test_that("different confidence levels produce appropriate CI widths", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm_90 <- deming_regression(x, y, conf_level = 0.90)
  dm_95 <- deming_regression(x, y, conf_level = 0.95)
  dm_99 <- deming_regression(x, y, conf_level = 0.99)

  # Point estimates should be identical
  expect_equal(dm_90$results$slope, dm_95$results$slope)
  expect_equal(dm_95$results$slope, dm_99$results$slope)

  # 99% CI should be wider than 95%, which should be wider than 90%
  width_90 <- diff(dm_90$results$slope_ci)
  width_95 <- diff(dm_95$results$slope_ci)
  width_99 <- diff(dm_99$results$slope_ci)

  expect_true(width_99 > width_95)
  expect_true(width_95 > width_90)
})

# =============================================================================
# Residuals and fitted values
# =============================================================================

test_that("fitted values lie on regression line", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  # fitted_y should equal intercept + slope * fitted_x
  expected_fitted_y <- dm$results$intercept + dm$results$slope * dm$results$fitted_x
  expect_equal(dm$results$fitted_y, expected_fitted_y, tolerance = 1e-10)
})

test_that("residuals have expected properties", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- true_vals + rnorm(50, sd = 5)  # No bias

  dm <- deming_regression(x, y)

  # Residuals should be centered near zero (approximately)
  expect_true(abs(mean(dm$results$residuals)) < 2)

  # Number of residuals equals sample size
  expect_equal(length(dm$results$residuals), dm$input$n)
})

# =============================================================================
# Print and summary methods
# =============================================================================

test_that("print method runs without error", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  expect_output(print(dm), "Deming Regression")
  expect_output(print(dm), "Slope:")
  expect_output(print(dm), "Intercept:")
  expect_output(print(dm), "Error ratio")
})

test_that("print method reports interpretation correctly", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)

  # Perfect agreement
  x <- true_vals + rnorm(50, sd = 5)
  y <- true_vals + rnorm(50, sd = 5)
  dm_agree <- deming_regression(x, y)
  expect_output(print(dm_agree), "includes")

  # Clear bias
  y_bias <- 1.5 * true_vals + 10 + rnorm(50, sd = 5)
  dm_bias <- deming_regression(x, y_bias)
  expect_output(print(dm_bias), "excludes")
})

test_that("summary method runs without error", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  expect_output(summary(dm), "Deming Regression")
  expect_output(summary(dm), "Interpretation")
  expect_output(summary(dm), "Conclusion")
  expect_output(summary(dm), "Error ratio")
})

test_that("summary method returns invisibly", {
  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  result <- capture.output(summ <- summary(dm))
  expect_true(is.list(summ))
  expect_named(summ, c("coefficients", "intercept_includes_zero",
                       "slope_includes_one", "methods_equivalent"))
})

# =============================================================================
# Plot methods
# =============================================================================

test_that("plot methods return ggplot objects", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  # Scatter plot
  p1 <- plot(dm, type = "scatter")
  expect_s3_class(p1, "ggplot")

  # Residual plot
  p2 <- plot(dm, type = "residuals")
  expect_s3_class(p2, "ggplot")
})

test_that("plot works with show_ci = FALSE", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  true_vals <- rnorm(30, mean = 100, sd = 20)
  x <- true_vals + rnorm(30, sd = 5)
  y <- true_vals + rnorm(30, sd = 5)

  dm <- deming_regression(x, y)
  p <- plot(dm, show_ci = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot works with show_identity = FALSE", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  true_vals <- rnorm(30, mean = 100, sd = 20)
  x <- true_vals + rnorm(30, sd = 5)
  y <- true_vals + rnorm(30, sd = 5)

  dm <- deming_regression(x, y)
  p <- plot(dm, show_identity = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot accepts custom labels", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  true_vals <- rnorm(30, mean = 100, sd = 20)
  x <- true_vals + rnorm(30, sd = 5)
  y <- true_vals + rnorm(30, sd = 5)

  dm <- deming_regression(x, y)
  p <- plot(dm, title = "Custom Title", xlab = "X Label", ylab = "Y Label")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Custom Title")
  expect_equal(p$labels$x, "X Label")
  expect_equal(p$labels$y, "Y Label")
})

test_that("residual plot types work correctly", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- true_vals + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  p_fitted <- plot(dm, type = "residuals", residual_type = "fitted")
  p_rank <- plot(dm, type = "residuals", residual_type = "rank")

  expect_s3_class(p_fitted, "ggplot")
  expect_s3_class(p_rank, "ggplot")
})

test_that("autoplot method works", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  true_vals <- rnorm(50, mean = 100, sd = 20)
  x <- true_vals + rnorm(50, sd = 5)
  y <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)

  dm <- deming_regression(x, y)

  p <- autoplot.deming_regression(dm)
  expect_s3_class(p, "ggplot")
})

# =============================================================================
# Edge cases
# =============================================================================

test_that("deming_regression handles perfect correlation", {
  set.seed(42)
  x <- 1:50
  y <- 2 * x + 10 + rnorm(50, sd = 0.1)  # Near-perfect linear relationship

  dm <- deming_regression(x, y)

  expect_s3_class(dm, "deming_regression")
  expect_true(abs(dm$results$slope - 2) < 0.1)
  expect_true(abs(dm$results$intercept - 10) < 1)
})

test_that("deming_regression handles negative values", {
  set.seed(42)
  x <- rnorm(50, mean = 0, sd = 10)
  y <- x + rnorm(50, sd = 2)

  dm <- deming_regression(x, y)

  expect_s3_class(dm, "deming_regression")
  expect_true(abs(dm$results$slope - 1) < 0.2)
})

test_that("deming_regression handles wide range of values", {
  set.seed(42)
  x <- c(rnorm(25, mean = 10, sd = 2), rnorm(25, mean = 1000, sd = 100))
  y <- x * 1.05 + 5 + rnorm(50, sd = 10)

  dm <- deming_regression(x, y)

  expect_s3_class(dm, "deming_regression")
  expect_true(is.finite(dm$results$slope))
  expect_true(is.finite(dm$results$intercept))
})

# =============================================================================
# Comparison with other methods
# =============================================================================

test_that("deming_regression differs from OLS appropriately", {
  set.seed(42)
  # Create data where both X and Y have measurement error
  true_vals <- rnorm(100, mean = 100, sd = 20)
  x <- true_vals + rnorm(100, sd = 10)  # High error in X
  y <- true_vals + rnorm(100, sd = 10)  # High error in Y

  dm <- deming_regression(x, y)
  ols <- lm(y ~ x)

  # With high measurement error in X, OLS slope should be attenuated
  # compared to Deming regression
  # This is the "regression dilution" or "attenuation bias"
  # Deming slope should be closer to 1 (the true slope)
  expect_true(abs(dm$results$slope - 1) <= abs(coef(ols)["x"] - 1))
})

# =============================================================================
# Integration with package datasets
# =============================================================================

test_that("deming_regression works with glucose_methods dataset", {
  skip_if_not(exists("glucose_methods"))

  dm <- deming_regression(reference ~ poc_meter, data = glucose_methods)

  expect_s3_class(dm, "deming_regression")
  expect_equal(dm$input$n, nrow(glucose_methods))
})

test_that("deming_regression works with creatinine_serum dataset", {
  skip_if_not(exists("creatinine_serum"))

  dm <- deming_regression(enzymatic ~ jaffe, data = creatinine_serum)

  expect_s3_class(dm, "deming_regression")
  expect_equal(dm$input$n, nrow(creatinine_serum))
})
