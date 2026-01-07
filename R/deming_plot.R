#' Plot method for deming_regression objects
#'
#' @description
#' Creates publication-ready plots for Deming regression results.
#' Multiple plot types are available: scatter plot with regression line
#' and residual plot.
#'
#' @param x An object of class `deming_regression`.
#' @param type Character; type of plot to create:
#'   \itemize{
#'     \item `"scatter"` (default): Scatter plot with regression line, CI band,
#'       and identity line
#'     \item `"residuals"`: Residuals vs. fitted values or rank
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
#' Displays the raw data with the fitted Deming regression line and
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
#' @examples
#' set.seed(42)
#' true_vals <- rnorm(50, 100, 20)
#' method_a <- true_vals + rnorm(50, sd = 5)
#' method_b <- 1.05 * true_vals + 3 + rnorm(50, sd = 5)
#' dm <- deming_regression(method_a, method_b)
#'
#' # Scatter plot with regression line
#' plot(dm)
#'
#' # Without identity line
#' plot(dm, show_identity = FALSE)
#'
#' # Residual plot
#' plot(dm, type = "residuals")
#'
#' # Residuals by rank
#' plot(dm, type = "residuals", residual_type = "rank")
#'
#' # Customized appearance
#' plot(dm, point_size = 3, title = "Glucose: POC vs Reference")
#'
#' @seealso [deming_regression()] for performing the analysis,
#'   [summary.deming_regression()] for detailed results
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_ribbon geom_hline
#'   geom_line labs theme_bw theme element_text coord_cartesian
#'   scale_x_continuous scale_y_continuous geom_smooth
#' @export
plot.deming_regression <- function(x,
                                   type = c("scatter", "residuals"),
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
         scatter = .plot_deming_scatter(x, show_ci, show_identity, point_alpha,
                                        point_size, line_colors, title, xlab, ylab),
         residuals = .plot_deming_residuals(x, residual_type, point_alpha,
                                            point_size, line_colors, title,
                                            xlab, ylab)
  )
}


# Plot Helper Functions ----

#' Scatter plot for Deming regression
#' @noRd
.plot_deming_scatter <- function(x, show_ci, show_identity, point_alpha,
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
    title <- "Deming Regression"
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

  # Create equation and lambda label
  eq_label <- sprintf("y = %.3f + %.3f x  (lambda = %.2f)",
                      res$intercept, res$slope, settings$error_ratio)

  # Add labels and theme
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d, %s CI (%s)",
                         input$n, ci_pct, settings$ci_method),
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


#' Residual plot for Deming regression
#' @noRd
.plot_deming_residuals <- function(x, residual_type, point_alpha,
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
  if (is.null(title)) title <- "Deming Regression Residuals"

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

  # Add smooth line to detect trends (using loess)
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
      subtitle = sprintf("n = %d, lambda = %.2f",
                         input$n, x$settings$error_ratio),
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


#' @rdname plot.deming_regression
#' @param object An object of class `deming_regression`.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.deming_regression <- function(object,
                                       type = c("scatter", "residuals"),
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
  plot.deming_regression(
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
