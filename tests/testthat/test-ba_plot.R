# Tests for plot.ba_analysis()
# =============================================================================

test_that("plot.ba_analysis returns a ggplot object", {
  set.seed(123)
  x <- rnorm(50, 100, 15)
  y <- x + rnorm(50, 2, 5)

  ba <- ba_analysis(x, y)
  p <- plot(ba)

  expect_s3_class(p, "ggplot")
})


test_that("plot.ba_analysis works with show_ci = FALSE", {
  set.seed(234)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 3)

  ba <- ba_analysis(x, y)
  p <- plot(ba, show_ci = FALSE)

  expect_s3_class(p, "ggplot")
})


test_that("plot.ba_analysis works with show_points = FALSE", {
  set.seed(345)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 3)

  ba <- ba_analysis(x, y)
  p <- plot(ba, show_points = FALSE)

  expect_s3_class(p, "ggplot")
})


test_that("plot.ba_analysis accepts custom colors", {
  set.seed(456)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 3)

  ba <- ba_analysis(x, y)

  custom_colors <- c(bias = "darkgreen", loa = "purple", ci = "lightblue")
  p <- plot(ba, line_colors = custom_colors)

  expect_s3_class(p, "ggplot")
})


test_that("plot.ba_analysis accepts partial custom colors", {
  set.seed(567)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 3)

  ba <- ba_analysis(x, y)

  # Only override bias color
  p <- plot(ba, line_colors = c(bias = "orange"))

  expect_s3_class(p, "ggplot")
})


test_that("plot.ba_analysis accepts custom labels", {
  set.seed(678)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 3)

  ba <- ba_analysis(x, y)

  p <- plot(ba,
            title = "Custom Title",
            xlab = "Custom X Label",
            ylab = "Custom Y Label")

  expect_s3_class(p, "ggplot")

  # Check that labels are set (access ggplot internals)
  expect_equal(p$labels$title, "Custom Title")
  expect_equal(p$labels$x, "Custom X Label")
  expect_equal(p$labels$y, "Custom Y Label")
})


test_that("plot.ba_analysis works with percent differences", {
  x <- c(100, 200, 300, 400, 500)
  y <- c(105, 210, 315, 420, 525)

  ba <- ba_analysis(x, y, type = "percent")
  p <- plot(ba)

  expect_s3_class(p, "ggplot")

  # Y-axis label should mention "Percent"
  expect_true(grepl("Percent|percent", p$labels$y, ignore.case = TRUE))
})


test_that("plot.ba_analysis respects point_alpha and point_size", {
  set.seed(789)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 3)

  ba <- ba_analysis(x, y)

  # These should not error
  p1 <- plot(ba, point_alpha = 0.3, point_size = 4)
  p2 <- plot(ba, point_alpha = 1, point_size = 1)


  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("plot.ba_analysis auto-generates appropriate axis labels", {
  set.seed(890)
  df <- data.frame(
    glucose_ref = rnorm(30, 100, 20),
    glucose_test = rnorm(30, 102, 20)
  )

  ba <- ba_analysis(glucose_ref ~ glucose_test, data = df)
  p <- plot(ba)

  # X label should mention both variable names
  expect_true(grepl("glucose_ref", p$labels$x))
  expect_true(grepl("glucose_test", p$labels$x))

  # Y label should show difference direction
  expect_true(grepl("glucose_test", p$labels$y))
  expect_true(grepl("glucose_ref", p$labels$y))
})


test_that("plot.ba_analysis can be extended with ggplot2 functions", {
  skip_if_not_installed("ggplot2")

  set.seed(901)
  x <- rnorm(30, 50, 10)
  y <- x + rnorm(30, 1, 3)

  ba <- ba_analysis(x, y)
  p <- plot(ba)

  # Add ggplot2 layers
  p2 <- p + ggplot2::theme_minimal()
  p3 <- p + ggplot2::ggtitle("Modified Title")

  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})


test_that("plot.ba_analysis handles edge case with small sample", {
  x <- c(10, 20, 30)
  y <- c(11, 21, 31)

  ba <- ba_analysis(x, y)
  p <- plot(ba)

  expect_s3_class(p, "ggplot")
})


test_that("plot.ba_analysis handles perfect agreement", {
  x <- 1:20
  y <- 1:20

  ba <- ba_analysis(x, y)
  p <- plot(ba)

  expect_s3_class(p, "ggplot")
})
