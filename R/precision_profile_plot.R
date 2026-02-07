#' Plot method for precision_profile objects
#'
#' @description
#' Creates publication-ready visualization of precision profile results,
#' showing CV vs concentration with the fitted model curve.
#'
#' @param x An object of class `precision_profile`.
#' @param show_ci Logical; if `TRUE` (default), displays prediction interval
#'   bands for the fitted curve.
#' @param show_targets Logical; if `TRUE` (default), displays horizontal lines
#'   at functional sensitivity target CV values.
#' @param show_points Logical; if `TRUE` (default), displays the observed
#'   data points.
#' @param point_alpha Numeric; transparency of points (0-1, default: 0.8).
#' @param point_size Numeric; size of points (default: 3).
#' @param line_colors Named character vector with colors for `"fitted"`,
#'   `"ci"`, and `"target"`. Defaults to a clean color scheme.
#' @param title Character; plot title. If `NULL` (default), generates an
#'   automatic title.
#' @param xlab Character; x-axis label. If `NULL`, uses "Concentration".
#' @param ylab Character; y-axis label. If `NULL`, uses "CV (%)".
#' @param log_x Logical; if `TRUE`, uses logarithmic scale for x-axis
#'   (default: `FALSE`).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A `ggplot` object that can be further customized.
#'
#' @details
#' The precision profile plot displays:
#' \itemize{
#'   \item **Observed points**: CV values at each tested concentration
#'   \item **Fitted curve**: Model-predicted CV across the concentration range
#'   \item **Prediction intervals**: Confidence bands showing uncertainty
#'   \item **Target lines**: Horizontal lines at functional sensitivity thresholds
#' }
#'
#' The plot helps visualize:
#' - How measurement precision changes with concentration
#' - Model fit quality (points should follow the curve)
#' - Functional sensitivity estimates (intersection of curve with target lines)
#'
#' @examples
#' # See ?precision_profile for complete examples
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_hline
#'   labs theme_bw theme element_text scale_x_continuous scale_y_continuous
#'   scale_x_log10 annotation_logticks
#' @export
plot.precision_profile <- function(x,
                                   show_ci = TRUE,
                                   show_targets = TRUE,
                                   show_points = TRUE,
                                   point_alpha = 0.8,
                                   point_size = 3,
                                   line_colors = NULL,
                                   title = NULL,
                                   xlab = NULL,
                                   ylab = NULL,
                                   log_x = FALSE,
                                   ...) {
  
  # Check ggplot2 availability
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. ",
         "Please install it with install.packages('ggplot2').",
         call. = FALSE)
  }
  
  # Default colors
  default_colors <- c(
    fitted = "#2166AC",
    ci = "#B2182B",
    target = "#666666",
    points = "#000000"
  )
  
  if (is.null(line_colors)) {
    line_colors <- default_colors
  } else {
    line_colors <- modifyList(as.list(default_colors), as.list(line_colors))
    line_colors <- unlist(line_colors)
  }
  
  # Labels
  if (is.null(xlab)) xlab <- "Concentration"
  if (is.null(ylab)) ylab <- "CV (%)"
  if (is.null(title)) {
    title <- sprintf("Precision Profile (%s model)",
                     tools::toTitleCase(x$model$type))
  }
  
  # Prepare plotting data ----
  
  # Generate smooth curve for fitted line
  conc_range <- range(x$input$concentration)
  conc_smooth <- seq(conc_range[1], conc_range[2], length.out = 200)
  
  # Calculate fitted CV for smooth curve
  if (x$model$type == "hyperbolic") {
    a <- x$model$parameters["a"]
    b <- x$model$parameters["b"]
    cv_smooth <- sqrt(a^2 + (b / conc_smooth)^2)
  } else {
    a <- x$model$parameters["a"]
    b <- x$model$parameters["b"]
    cv_smooth <- a + b / conc_smooth
  }
  
  smooth_data <- data.frame(
    concentration = conc_smooth,
    cv_fitted = cv_smooth
  )
  
  # Build plot ----
  p <- ggplot2::ggplot()
  
  # Add prediction interval band
  if (show_ci) {
    p <- p +
      ggplot2::geom_ribbon(
        data = x$fitted,
        ggplot2::aes(x = .data$concentration,
                     ymin = .data$ci_lower,
                     ymax = .data$ci_upper),
        fill = line_colors["ci"],
        alpha = 0.2
      )
  }
  
  # Add target lines for functional sensitivity
  if (show_targets) {
    for (i in seq_len(nrow(x$functional_sensitivity))) {
      fs <- x$functional_sensitivity[i, ]
      
      if (fs$achievable && !is.na(fs$concentration)) {
        # Horizontal line at target CV
        p <- p +
          ggplot2::geom_hline(
            yintercept = fs$cv_target,
            color = line_colors["target"],
            linetype = "dashed",
            linewidth = 0.5
          )
        
        # Vertical line at functional sensitivity concentration
        p <- p +
          ggplot2::geom_vline(
            xintercept = fs$concentration,
            color = line_colors["target"],
            linetype = "dotted",
            linewidth = 0.5
          )
      }
    }
  }
  
  # Add fitted curve
  p <- p +
    ggplot2::geom_line(
      data = smooth_data,
      ggplot2::aes(x = .data$concentration, y = .data$cv_fitted),
      color = line_colors["fitted"],
      linewidth = 1
    )
  
  # Add observed points
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        data = x$fitted,
        ggplot2::aes(x = .data$concentration, y = .data$cv_observed),
        color = line_colors["points"],
        size = point_size,
        alpha = point_alpha
      )
  }
  
  # Add labels and theme
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("R² = %.3f, %s",
                         x$fit_quality$r_squared,
                         x$model$equation),
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
  
  # Apply log scale if requested
  if (log_x) {
    p <- p +
      ggplot2::scale_x_log10() +
      ggplot2::annotation_logticks(sides = "b")
  }
  
  # Ensure y-axis starts at 0
  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0))
  
  p
}


#' @rdname plot.precision_profile
#' @param object An object of class `precision_profile`.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.precision_profile <- function(object,
                                       show_ci = TRUE,
                                       show_targets = TRUE,
                                       show_points = TRUE,
                                       point_alpha = 0.8,
                                       point_size = 3,
                                       line_colors = NULL,
                                       title = NULL,
                                       xlab = NULL,
                                       ylab = NULL,
                                       log_x = FALSE,
                                       ...) {
  plot.precision_profile(
    x = object,
    show_ci = show_ci,
    show_targets = show_targets,
    show_points = show_points,
    point_alpha = point_alpha,
    point_size = point_size,
    line_colors = line_colors,
    title = title,
    xlab = xlab,
    ylab = ylab,
    log_x = log_x,
    ...
  )
}