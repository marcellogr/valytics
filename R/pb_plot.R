#' Plot method for pb_regression objects
#'
#' @description
#' Creates publication-ready plots for Passing-Bablok regression results.
#' Multiple plot types are available: scatter plot with regression line,
#' residual plot, and CUSUM plot for linearity assessment.
#'
#' @param x An object of class `pb_regression`.
#' @param type Character; type of plot to create:
#'   \itemize{
#'     \item `"scatter"` (default): Scatter plot with regression line, CI band,
#'       and identity line
#'     \item `"residuals"`: Residuals vs. fitted values or rank
#'     \item `"cusum"`: CUSUM plot for linearity assessment
#'   }
#' @param show_ci Logical; if `TRUE` (default), displays confidence band for
#'   the regression line (only for `type = "scatter"`).
#' @param show_identity Logical; if `TRUE` (default), displays the identity
#'   line (y = x) for reference.
#' @param residual_type Character; for `type = "residuals"`, plot residuals
#'   against `"fitted"` (default) or `"rank"` (ordered by x).
#' @param point_alpha Numeric; transparency of points (0-1, default: 0.6).
#' @param point_size Numeric; size of points (default: 2).
#' @param line_colors Named character vector with colors for `"regression"`,
#'   `"identity"`, and `"ci"`. Defaults to a clean color scheme.
#' @param title Character; plot title. If `NULL` (default), generates an
#'   automatic title.
#' @param xlab,ylab Character; axis labels. If `NULL`, auto-generates based
#'   on variable names.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A `ggplot` object that can be further customized.
#'
#' @details
#' **Scatter plot** (`type = "scatter"`):
#' Displays the raw data with the fitted Passing-Bablok regression line and
#' optional confidence band. The identity line (y = x) is shown for reference.
#' If the regression line overlaps substantially with the identity line, the
#' methods are in good agreement.
#'
#' **Residual plot** (`type = "residuals"`):
#' Displays perpendicular residuals. Look for:
#' \itemize{
#'   \item Random scatter around zero (good)
#'   \item Patterns or trends (suggests non-linearity)
#'   \item Funnel shape (suggests heteroscedasticity)
#' }
#'
#' **CUSUM plot** (`type = "cusum"`):
#' Displays the cumulative sum of residual signs, used to assess linearity.
#' The CUSUM should stay within the critical bounds if the linearity assumption
#' holds.
#'
#' @examples
#' set.seed(42)
#' method_a <- rnorm(50, mean = 100, sd = 15)
#' method_b <- 1.05 * method_a + 3 + rnorm(50, sd = 5)
#' pb <- pb_regression(method_a, method_b)
#'
#' # Scatter plot with regression line
#' plot(pb)
#'
#' # Without identity line
#' plot(pb, show_identity = FALSE)
#'
#' # Residual plot
#' plot(pb, type = "residuals")
#'
#' # Residuals by rank
#' plot(pb, type = "residuals", residual_type = "rank")
#'
#' # CUSUM plot
#' plot(pb, type = "cusum")
#'
#' # Customized appearance
#' plot(pb, point_size = 3, title = "Glucose: POC vs Reference")
#'
#' @seealso [pb_regression()] for performing the analysis,
#'   [summary.pb_regression()] for detailed results
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_ribbon geom_hline
#'   geom_line geom_segment labs theme_bw theme element_text coord_cartesian
#'   scale_x_continuous scale_y_continuous
#' @export
plot.pb_regression <- function(x,
                               type = c("scatter", "residuals", "cusum"),
                               show_ci = TRUE,
                               show_identity = TRUE,
                               residual_type = c("fitted", "rank"),
                               point_alpha = 0.6,
                               point_size = 2,
                               line_colors = NULL,
                               title = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               ...) {

  # Check ggplot2 availability
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. ",
         "Please install it with install.packages('ggplot2').",
         call. = FALSE)
  }

  type <- match.arg(type)
  residual_type <- match.arg(residual_type)

  # Default colors
  default_colors <- c(
    regression = "#2166AC",
    identity = "#666666",
    ci = "#B2182B",
    zero = "#999999"
  )

  if (is.null(line_colors)) {
    line_colors <- default_colors
  } else {
    line_colors <- modifyList(as.list(default_colors), as.list(line_colors))
    line_colors <- unlist(line_colors)
  }

  # Dispatch to appropriate plot function
  switch(type,
         scatter = .plot_pb_scatter(x, show_ci, show_identity, point_alpha,
                                    point_size, line_colors, title, xlab, ylab),
         residuals = .plot_pb_residuals(x, residual_type, point_alpha,
                                        point_size, line_colors, title,
                                        xlab, ylab),
         cusum = .plot_pb_cusum(x, line_colors, title, xlab, ylab)
  )
}


# Plot Helper Functions ----

#' Scatter plot for Passing-Bablok regression
#' @noRd
.plot_pb_scatter <- function(x, show_ci, show_identity, point_alpha,
                             point_size, line_colors, title, xlab, ylab) {

  # Extract data
  res <- x$results
  input <- x$input
  settings <- x$settings

  # Prepare plot data
  plot_data <- data.frame(
    x = input$x,
    y = input$y
  )

  # Axis labels
  if (is.null(xlab)) xlab <- input$var_names["x"]
  if (is.null(ylab)) ylab <- input$var_names["y"]

  if (is.null(title)) {
    title <- "Passing-Bablok Regression"
  }

  ci_pct <- paste0(settings$conf_level * 100, "%")

  # Build plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y))

  # Add CI band for regression line
  if (show_ci && !any(is.na(res$slope_ci)) && !any(is.na(res$intercept_ci))) {
    # Create prediction band
    x_range <- range(input$x)
    x_seq <- seq(x_range[1], x_range[2], length.out = 100)

    # Upper and lower bounds using CI extremes
    # Conservative approach: use combinations that give widest band
    y_upper <- pmax(
      res$intercept_ci["upper"] + res$slope_ci["upper"] * x_seq,
      res$intercept_ci["upper"] + res$slope_ci["lower"] * x_seq,
      res$intercept_ci["lower"] + res$slope_ci["upper"] * x_seq
    )
    y_lower <- pmin(
      res$intercept_ci["lower"] + res$slope_ci["lower"] * x_seq,
      res$intercept_ci["lower"] + res$slope_ci["upper"] * x_seq,
      res$intercept_ci["upper"] + res$slope_ci["lower"] * x_seq
    )

    ci_data <- data.frame(x = x_seq, y_lower = y_lower, y_upper = y_upper)

    p <- p +
      ggplot2::geom_ribbon(
        data = ci_data,
        ggplot2::aes(x = .data$x, ymin = .data$y_lower, ymax = .data$y_upper),
        fill = line_colors["ci"],
        alpha = 0.2,
        inherit.aes = FALSE
      )
  }

  # Add identity line (y = x)
  if (show_identity) {
    p <- p +
      ggplot2::geom_abline(
        intercept = 0, slope = 1,
        color = line_colors["identity"],
        linetype = "dashed",
        linewidth = 0.7
      )
  }

  # Add regression line
  p <- p +
    ggplot2::geom_abline(
      intercept = res$intercept, slope = res$slope,
      color = line_colors["regression"],
      linewidth = 1
    )

  # Add points
  p <- p +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size,
      color = "black"
    )

  # Create equation label
  eq_label <- sprintf("y = %.3f + %.3f x", res$intercept, res$slope)

  # Add labels and theme
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d, %s CI (%s)",
                         input$n, ci_pct,
                         if (settings$ci_method == "analytical") "analytical" else "bootstrap"),
      x = xlab,
      y = ylab,
      caption = eq_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0),
      axis.title = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      aspect.ratio = 1
    )

  p
}


#' Residual plot for Passing-Bablok regression
#' @noRd
.plot_pb_residuals <- function(x, residual_type, point_alpha,
                               point_size, line_colors, title, xlab, ylab) {

  # Extract data
  res <- x$results
  input <- x$input

  # Prepare plot data based on residual type
  if (residual_type == "fitted") {
    plot_data <- data.frame(
      x_val = res$fitted_y,
      residual = res$residuals
    )
    if (is.null(xlab)) xlab <- "Fitted values"
  } else {
    # Rank order by x
    ord <- order(input$x)
    plot_data <- data.frame(
      x_val = seq_along(input$x),
      residual = res$residuals[ord]
    )
    if (is.null(xlab)) xlab <- "Rank (ordered by X)"
  }

  if (is.null(ylab)) ylab <- "Perpendicular residual"
  if (is.null(title)) title <- "Passing-Bablok Residuals"

  # Build plot
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = .data$x_val, y = .data$residual))

  # Add zero reference line
  p <- p +
    ggplot2::geom_hline(
      yintercept = 0,
      color = line_colors["zero"],
      linetype = "dashed",
      linewidth = 0.7
    )

  # Add points
  p <- p +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size,
      color = "black"
    )

  # Add smooth line to detect trends (optional, using loess)
  if (nrow(plot_data) >= 10) {
    p <- p +
      ggplot2::geom_smooth(
        method = "loess",
        formula = y ~ x,
        se = FALSE,
        color = line_colors["regression"],
        linewidth = 0.8,
        alpha = 0.7
      )
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d", input$n),
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      axis.title = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}


#' CUSUM plot for Passing-Bablok regression
#' @noRd
.plot_pb_cusum <- function(x, line_colors, title, xlab, ylab) {

  # Extract data
  input <- x$input
  res <- x$results
  cusum_res <- x$cusum

  n <- input$n

  # Recalculate CUSUM values for plotting
  fitted_y <- res$intercept + res$slope * input$x
  residuals <- input$y - fitted_y

  n_pos <- sum(residuals > 0)
  n_neg <- sum(residuals < 0)

  if (n_pos == 0 || n_neg == 0) {
    stop("Cannot create CUSUM plot: all residuals have the same sign.",
         call. = FALSE)
  }

  score_pos <- sqrt(n_neg / n_pos)
  score_neg <- -sqrt(n_pos / n_neg)

  scores <- ifelse(residuals > 0, score_pos,
                   ifelse(residuals < 0, score_neg, 0))

  # Sort by x values
  ord <- order(input$x)
  scores_sorted <- scores[ord]

  # Calculate cumulative sums
  cusum <- cumsum(scores_sorted)

  # Critical bounds (approximate, based on K-S distribution)
  # At alpha = 0.05, critical value is 1.36
  critical <- 1.36 * sqrt(n)

  # Prepare plot data
  plot_data <- data.frame(
    rank = seq_along(cusum),
    cusum = cusum
  )

  if (is.null(title)) title <- "CUSUM Linearity Test"
  if (is.null(xlab)) xlab <- "Rank (ordered by X)"
  if (is.null(ylab)) ylab <- "Cumulative sum"

  # Build plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$rank, y = .data$cusum))

  # Add critical bounds
  p <- p +
    ggplot2::geom_hline(
      yintercept = c(-critical, critical),
      color = line_colors["ci"],
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = line_colors["zero"],
      linetype = "dotted",
      linewidth = 0.5
    )

  # Add CUSUM line and points
  p <- p +
    ggplot2::geom_line(
      color = line_colors["regression"],
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      color = line_colors["regression"],
      size = 1.5,
      alpha = 0.7
    )

  # Add annotation for test result
  result_text <- if (!is.na(cusum_res$linear) && cusum_res$linear) {
    sprintf("Linearity OK (H = %.3f, p = %.3f)", cusum_res$statistic, cusum_res$p_value)
  } else if (!is.na(cusum_res$linear)) {
    sprintf("Linearity violated (H = %.3f, p = %.3f)", cusum_res$statistic, cusum_res$p_value)
  } else {
    "Test not available"
  }

  # Labels and theme
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = result_text,
      x = xlab,
      y = ylab,
      caption = sprintf("Critical bounds at alpha = 0.05: +/-%.2f", critical)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0),
      axis.title = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}


#' @rdname plot.pb_regression
#' @param object An object of class `pb_regression`.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.pb_regression <- function(object,
                                   type = c("scatter", "residuals", "cusum"),
                                   show_ci = TRUE,
                                   show_identity = TRUE,
                                   residual_type = c("fitted", "rank"),
                                   point_alpha = 0.6,
                                   point_size = 2,
                                   line_colors = NULL,
                                   title = NULL,
                                   xlab = NULL,
                                   ylab = NULL,
                                   ...) {
  plot.pb_regression(
    x = object,
    type = type,
    show_ci = show_ci,
    show_identity = show_identity,
    residual_type = residual_type,
    point_alpha = point_alpha,
    point_size = point_size,
    line_colors = line_colors,
    title = title,
    xlab = xlab,
    ylab = ylab,
    ...
  )
}