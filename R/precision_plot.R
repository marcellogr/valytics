#' Plot method for precision_study objects
#'
#' @description
#' Creates visualizations for precision study results. Multiple plot types
#' are available: variance component chart, CV profile across samples,
#' and precision estimates summary.
#'
#' @param x An object of class `precision_study`.
#' @param type Character; type of plot to create:
#'   \itemize{
#'     \item `"variance"` (default): Bar chart showing variance components
#'       as percentage of total variance
#'     \item `"cv"`: CV profile across samples (requires multi-sample data)
#'     \item `"precision"`: Forest plot of precision estimates with CIs
#'   }
#' @param show_ci Logical; if `TRUE` (default), displays confidence intervals
#'   where applicable.
#' @param colors Character vector of colors for the plot elements. If `NULL`,
#'   uses a default color palette.
#' @param title Character; plot title. If `NULL` (default), generates an
#'   automatic title.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A `ggplot` object that can be further customized.
#'
#' @details
#' **Variance component chart** (`type = "variance"`):
#' Displays the proportion of total variance attributable to each source
#' (between-day, between-run, error/repeatability). Helps identify which
#' factors contribute most to measurement variability.
#'
#' **CV profile** (`type = "cv"`):
#' For multi-sample studies, displays how CV varies across concentration
#' levels. Typically CV is higher at low concentrations. Requires data
#' from multiple samples/levels.
#'
#' **Precision summary** (`type = "precision"`):
#' Forest plot showing precision estimates (repeatability, intermediate
#' precision, reproducibility) with confidence intervals.
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' prec_data <- data.frame(
#'   day = rep(1:5, each = 6),
#'   run = rep(rep(1:2, each = 3), 5),
#'   value = rnorm(30, mean = 100, sd = 5)
#' )
#' prec_data$value <- prec_data$value + rep(rnorm(5, 0, 3), each = 6)
#'
#' prec <- precision_study(prec_data, value = "value", day = "day", run = "run")
#'
#' # Variance component chart (default)
#' plot(prec)
#' plot(prec, type = "variance")
#'
#' # Precision estimates with CIs
#' plot(prec, type = "precision")
#'
#' @seealso [precision_study()] for performing the analysis,
#'   [summary.precision_study()] for detailed results
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_errorbar
#'   geom_errorbarh geom_line geom_ribbon geom_text labs theme_bw theme
#'   element_text element_blank coord_flip scale_fill_manual scale_color_manual
#'   scale_y_continuous scale_x_continuous guides guide_legend
#' @export
plot.precision_study <- function(x,
                                 type = c("variance", "cv", "precision"),
                                 show_ci = TRUE,
                                 colors = NULL,
                                 title = NULL,
                                 ...) {
  
  # Check ggplot2 availability
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. ",
         "Please install it with install.packages('ggplot2').",
         call. = FALSE)
  }
  
  type <- match.arg(type)
  
  # Dispatch to appropriate plot function
  switch(type,
         variance = .plot_prec_variance(x, colors, title),
         cv = .plot_prec_cv(x, show_ci, colors, title),
         precision = .plot_prec_forest(x, show_ci, colors, title)
  )
}


#' Variance component bar chart
#' @noRd
.plot_prec_variance <- function(x, colors, title) {
  
  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "between_site" = "#1B9E77",
      "between_day" = "#D95F02",
      "between_run" = "#7570B3",
      "error" = "#E7298A"
    )
  }
  
  # Get variance components (exclude 'total')
  vc <- x$variance_components
  vc <- vc[vc$component != "total", ]
  
  # Order components logically (nested order)
  component_order <- c("between_site", "between_day", "between_run", "error")
  vc$component <- factor(vc$component, 
                         levels = intersect(component_order, vc$component))
  
  # Create display labels
  label_map <- c(
    "between_site" = "Between-site",
    "between_day" = "Between-day",
    "between_run" = "Between-run",
    "error" = "Repeatability (error)"
  )
  vc$label <- label_map[as.character(vc$component)]
  vc$label <- factor(vc$label, levels = label_map[levels(vc$component)])
  
  # Default title
  if (is.null(title)) {
    title <- "Variance Components"
  }
  
  # Build plot
  p <- ggplot2::ggplot(vc, ggplot2::aes(x = .data$label, 
                                        y = .data$pct_total,
                                        fill = .data$component)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", .data$pct_total)),
                       hjust = -0.1, size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(limits = c(0, max(vc$pct_total) * 1.15),
                                expand = c(0, 0)) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Design: %s, n = %d", x$design$structure, x$input$n),
      x = NULL,
      y = "Percentage of total variance (%)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      axis.title = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  p
}


#' CV profile across samples
#' @noRd
.plot_prec_cv <- function(x, show_ci, colors, title) {
  
  # Check for multi-sample data
  if (is.null(x$by_sample) || length(x$by_sample) < 2) {
    stop("CV profile plot requires multi-sample data (multiple concentration levels).\n",
         "Use type = 'variance' or type = 'precision' for single-sample studies.",
         call. = FALSE)
  }
  
  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "Repeatability" = "#2166AC",
      "Intermediate precision" = "#B2182B",
      "Reproducibility" = "#1B7837"
    )
  }
  
  # Extract CV data from each sample
  sample_names <- names(x$by_sample)
  sample_means <- x$sample_means
  
  cv_data <- do.call(rbind, lapply(seq_along(x$by_sample), function(i) {
    samp <- x$by_sample[[i]]
    prec <- samp$precision
    
    data.frame(
      sample = sample_names[i],
      mean_conc = if (!is.null(sample_means)) sample_means[sample_names[i]] else i,
      measure = prec$measure,
      cv_pct = prec$cv_pct,
      ci_lower = if ("cv_ci_lower" %in% names(prec)) prec$cv_ci_lower else NA_real_,
      ci_upper = if ("cv_ci_upper" %in% names(prec)) prec$cv_ci_upper else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  
  # Filter to key measures (ones that make sense to compare across samples)
  key_measures <- c("Repeatability", "Intermediate precision", "Reproducibility")
  cv_data <- cv_data[cv_data$measure %in% key_measures, ]
  cv_data$measure <- factor(cv_data$measure, levels = key_measures)
  
  # Default title
  if (is.null(title)) {
    title <- "CV Profile Across Concentration Levels"
  }
  
  # Build plot
  p <- ggplot2::ggplot(cv_data, ggplot2::aes(x = .data$mean_conc, 
                                             y = .data$cv_pct,
                                             color = .data$measure,
                                             group = .data$measure)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3)
  
  # Add CI ribbons if requested and available
  if (show_ci && !all(is.na(cv_data$ci_lower))) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, 
                                        ymax = .data$ci_upper,
                                        fill = .data$measure),
                           alpha = 0.2, color = NA)
  }
  
  p <- p +
    ggplot2::scale_color_manual(values = colors, name = "Precision measure") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("%d concentration levels", length(x$by_sample)),
      x = "Mean concentration",
      y = "CV (%)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      axis.title = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))
  
  p
}


#' Forest plot of precision estimates
#' @noRd
.plot_prec_forest <- function(x, show_ci, colors, title) {
  
  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "Repeatability" = "#2166AC",
      "Between-run" = "#7570B3",
      "Between-day" = "#D95F02",
      "Between-site" = "#1B9E77",
      "Intermediate precision" = "#B2182B",
      "Reproducibility" = "#1B7837"
    )
  }
  
  # Get precision data
  prec <- x$precision
  
  # Order measures logically (bottom to top in forest plot after coord_flip)
  measure_order <- c("Reproducibility", "Intermediate precision", 
                     "Between-site", "Between-day", "Between-run",
                     "Repeatability")
  prec$measure <- factor(prec$measure, 
                         levels = rev(intersect(measure_order, prec$measure)))
  
  # Default title
  if (is.null(title)) {
    title <- "Precision Estimates"
  }
  
  ci_pct <- sprintf("%g%%", x$settings$conf_level * 100)
  
  # Calculate x limit
  x_max <- if ("cv_ci_upper" %in% names(prec) && !all(is.na(prec$cv_ci_upper))) {
    max(prec$cv_ci_upper, prec$cv_pct, na.rm = TRUE)
  } else {
    max(prec$cv_pct, na.rm = TRUE)
  }
  
  # Build plot
  p <- ggplot2::ggplot(prec, ggplot2::aes(x = .data$cv_pct, 
                                          y = .data$measure,
                                          color = .data$measure))
  
  # Add CI error bars if requested and available
  if (show_ci && "cv_ci_lower" %in% names(prec) && !all(is.na(prec$cv_ci_lower))) {
    p <- p +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$cv_ci_lower, 
                                           xmax = .data$cv_ci_upper),
                              height = 0.2, linewidth = 0.8)
  }
  
  p <- p +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%%", .data$cv_pct)),
                       hjust = -0.3, size = 3.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max * 1.25),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Design: %s, n = %d, %s CI", 
                         x$design$structure, x$input$n, ci_pct),
      x = "CV (%)",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      axis.title = ggplot2::element_text(size = 10),
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  p
}


#' @rdname plot.precision_study
#' @param object An object of class `precision_study`.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.precision_study <- function(object,
                                     type = c("variance", "cv", "precision"),
                                     show_ci = TRUE,
                                     colors = NULL,
                                     title = NULL,
                                     ...) {
  plot.precision_study(
    x = object,
    type = type,
    show_ci = show_ci,
    colors = colors,
    title = title,
    ...
  )
}