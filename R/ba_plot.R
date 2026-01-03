#' Plot method for ba_analysis objects
#'
#' @description
#' Creates a Bland-Altman plot (difference vs. average) for visualizing
#' agreement between two measurement methods. The plot displays the bias
#' (mean difference) and limits of agreement with optional confidence intervals.
#'
#' @param x An object of class `ba_analysis`.
#' @param show_ci Logical; if `TRUE` (default), displays confidence interval
#'   bands for bias and limits of agreement.
#' @param show_points Logical; if `TRUE` (default), displays individual
#'   data points.
#' @param point_alpha Numeric; transparency of points (0-1, default: 0.6).
#' @param point_size Numeric; size of points (default: 2).
#' @param line_colors Named character vector with colors for `"bias"`,
#'   `"loa"`, and `"ci"`. Defaults to a clean color scheme.
#' @param title Character; plot title. If `NULL` (default), generates an
#'   automatic title.
#' @param xlab Character; x-axis label. If `NULL`, uses "Mean of methods".
#' @param ylab Character; y-axis label. If `NULL`, auto-generates based on
#'   difference type.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A `ggplot` object that can be further customized.
#'
#' @details
#' The Bland-Altman plot displays:
#' \itemize{
#'   \item **Points**: Each point represents a paired observation, plotted as
#'     the difference (y - x) against the average ((x + y) / 2).
#'   \item **Bias line**: Solid horizontal line at the mean difference.
#'   \item **Limits of agreement**: Dashed horizontal lines at bias +/- 1.96 x SD.
#'   \item **Confidence intervals**: Shaded bands showing the uncertainty in
#'     the bias and LoA estimates.
#' }
#'
#' Patterns to look for:
#' \itemize{
#'   \item **Funnel shape**: Suggests proportional bias (variance increases
#'     with magnitude).
#'   \item **Trend**: Suggests systematic relationship between difference and
#'     magnitude.
#'   \item **Outliers**: Points outside the LoA may warrant investigation.
#' }
#'
#' @examples
#' # Basic Bland-Altman plot
#' set.seed(42)
#' method_a <- rnorm(50, mean = 100, sd = 15)
#' method_b <- method_a + rnorm(50, mean = 2, sd = 5)
#'
#' ba <- ba_analysis(method_a, method_b)
#' plot(ba)
#'
#' # Without confidence intervals
#' plot(ba, show_ci = FALSE)
#'
#' # Customized appearance
#' plot(ba,
#'      point_alpha = 0.8,
#'      point_size = 3,
#'      title = "Method Comparison: A vs B")
#'
#' # Further customization with ggplot2
#' library(ggplot2)
#' plot(ba) +
#'   theme_minimal() +
#'   scale_color_brewer(palette = "Set1")
#'
#' # Using autoplot (ggplot2-style)
#' autoplot(ba)
#'
#' @seealso
#' [ba_analysis()] for performing the analysis,
#' [summary.ba_analysis()] for detailed results
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_rect annotate
#'   labs theme_bw theme element_text scale_y_continuous coord_cartesian
#' @export
plot.ba_analysis <- function(x,
                             show_ci = TRUE,
                             show_points = TRUE,
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

  # Setup ----

  # Default colors
  default_colors <- c(
    bias = "#2166AC",
    loa = "#B2182B",
    ci = "#DDDDDD"
  )

  if (is.null(line_colors)) {
    line_colors <- default_colors
  } else {
    # Merge with defaults for any missing colors
    line_colors <- modifyList(as.list(default_colors), as.list(line_colors))
    line_colors <- unlist(line_colors)
  }

  # Extract results
  res <- x$results
  settings <- x$settings

  # Prepare plot data
  plot_data <- data.frame(
    average = res$averages,
    difference = res$differences
  )

  # Labels
  if (is.null(xlab)) {
    xlab <- sprintf("Mean of %s and %s",
                    x$input$var_names["x"],
                    x$input$var_names["y"])
  }

  if (is.null(ylab)) {
    if (settings$type == "absolute") {
      ylab <- sprintf("Difference (%s - %s)",
                      x$input$var_names["y"],
                      x$input$var_names["x"])
    } else {
      ylab <- sprintf("Percent difference (%s - %s)",
                      x$input$var_names["y"],
                      x$input$var_names["x"])
    }
  }

  if (is.null(title)) {
    title <- "Bland-Altman Plot"
  }

  # CI level for annotation
  ci_pct <- paste0(settings$conf_level * 100, "%")

  # Build plot ----

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$average,
                                               y = .data$difference))

  # Add CI bands first (so they're behind everything)
  # Use -Inf/Inf for full-width bands that extend across the entire plot
  if (show_ci) {
    # CI band for bias
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = res$bias_ci["lower"], ymax = res$bias_ci["upper"],
        fill = line_colors["ci"], alpha = 0.5
      )

    # CI band for lower LoA
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = res$loa_lower_ci["lower"], ymax = res$loa_lower_ci["upper"],
        fill = line_colors["ci"], alpha = 0.5
      )

    # CI band for upper LoA
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = res$loa_upper_ci["lower"], ymax = res$loa_upper_ci["upper"],
        fill = line_colors["ci"], alpha = 0.5
      )
  }

  # Add horizontal lines for bias and LoA
  p <- p +
    # Bias line (solid)
    ggplot2::geom_hline(
      yintercept = res$bias,
      color = line_colors["bias"],
      linewidth = 0.8
    ) +
    # Lower LoA (dashed)
    ggplot2::geom_hline(
      yintercept = res$loa_lower,
      color = line_colors["loa"],
      linetype = "dashed",
      linewidth = 0.7
    ) +
    # Upper LoA (dashed)
    ggplot2::geom_hline(
      yintercept = res$loa_upper,
      color = line_colors["loa"],
      linetype = "dashed",
      linewidth = 0.7
    ) +
    # Zero reference line (subtle)
    ggplot2::geom_hline(
      yintercept = 0,
      color = "gray50",
      linetype = "dotted",
      linewidth = 0.5
    )

  # Add points
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        alpha = point_alpha,
        size = point_size,
        color = "black"
      )
  }

  # Annotations ----

  # Calculate x position for annotations (right edge of data + margin)
  x_range <- range(plot_data$average)
  x_margin <- diff(x_range) * 0.02
  x_annot <- x_range[2] + x_margin

  p <- p +
    # Bias annotation
    ggplot2::annotate(
      "text",
      x = x_annot,
      y = res$bias,
      label = sprintf("Bias: %.2f", res$bias),
      hjust = 0,
      vjust = -0.5,
      size = 3,
      color = line_colors["bias"]
    ) +
    # Upper LoA annotation
    ggplot2::annotate(
      "text",
      x = x_annot,
      y = res$loa_upper,
      label = sprintf("+%.2f SD: %.2f", settings$multiplier, res$loa_upper),
      hjust = 0,
      vjust = -0.5,
      size = 3,
      color = line_colors["loa"]
    ) +
    # Lower LoA annotation
    ggplot2::annotate(
      "text",
      x = x_annot,
      y = res$loa_lower,
      label = sprintf("-%.2f SD: %.2f", settings$multiplier, res$loa_lower),
      hjust = 0,
      vjust = 1.5,
      size = 3,
      color = line_colors["loa"]
    )

  # Labels and theme
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d, %s CI", x$input$n, ci_pct),
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      axis.title = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      # Add right margin to accommodate annotations
      plot.margin = ggplot2::margin(5.5, 40, 5.5, 5.5, "pt")
    ) +
    # Set axis limits with clipping disabled for annotations
    ggplot2::coord_cartesian(
      xlim = c(x_range[1] - diff(x_range) * 0.05,
               x_range[2] + diff(x_range) * 0.15),
      clip = "off"
    )

  p
}


#' @rdname plot.ba_analysis
#' @param object An object of class `ba_analysis`.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.ba_analysis <- function(object,
                                 show_ci = TRUE,
                                 show_points = TRUE,
                                 point_alpha = 0.6,
                                 point_size = 2,
                                 line_colors = NULL,
                                 title = NULL,
                                 xlab = NULL,
                                 ylab = NULL,
                                 ...) {
  plot.ba_analysis(
    x = object,
    show_ci = show_ci,
    show_points = show_points,
    point_alpha = point_alpha,
    point_size = point_size,
    line_colors = line_colors,
    title = title,
    xlab = xlab,
    ylab = ylab,
    ...
  )
}
