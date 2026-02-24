#' Plot method for linearity_study objects
#'
#' @description
#' Creates publication-ready plots for linearity evaluation results.
#' Multiple plot types are available: linearity plot with fitted models,
#' deviation plot, and residual plot.
#'
#' @param x An object of class `linearity_study`.
#' @param type Character; type of plot to create:
#'   \itemize{
#'     \item `"linearity"` (default): Measured vs target with linear and
#'       polynomial fits
#'     \item `"deviation"`: Deviation from linearity at each level, with
#'       optional ADL limits
#'     \item `"residuals"`: Residuals from the linear model vs target
#'       concentration
#'     \item `"dilution"`: Dilution recovery plot showing observed means vs
#'       dilution proportion with WLS regression line. Only available for
#'       dilution-based studies.
#'   }
#' @param deviation_scale Character; for `type = "deviation"`, controls the
#'   y-axis scale: `"auto"` (default) uses absolute deviations, `"absolute"`
#'   forces absolute deviations, `"percent"` forces percent deviations. When
#'   ADL is specified as percent and plotted on the absolute scale, the ADL
#'   limits form a fan shape that widens with concentration.
#' @param show_points Logical; if `TRUE` (default), shows individual
#'   replicate measurements (for `type = "linearity"` and `"residuals"`).
#' @param show_means Logical; if `TRUE` (default), shows level means with
#'   error bars (for `type = "linearity"`).
#' @param show_identity Logical; if `TRUE` (default), shows the identity
#'   line y = x (for `type = "linearity"`).
#' @param point_alpha Numeric; transparency of points (0-1, default: 0.4).
#' @param point_size Numeric; size of points (default: 2).
#' @param line_colors Named character vector with colors for `"linear"`,
#'   `"polynomial"`, `"identity"`, `"adl"`, and `"zero"`.
#' @param title Character; plot title. If `NULL` (default), auto-generates.
#' @param xlab,ylab Character; axis labels. If `NULL`, auto-generates.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A `ggplot` object that can be further customized.
#'
#' @details
#' **Linearity plot** (`type = "linearity"`):
#' Displays measured values against target concentrations with the linear fit
#' and (if non-linear) the best polynomial fit. Level means are shown with
#' error bars representing +/- 1 SD. The identity line provides a visual
#' reference.
#'
#' **Deviation plot** (`type = "deviation"`):
#' Shows the deviation from linearity at each concentration level, calculated
#' as the difference between the polynomial best-fit and the linear fit.
#' Individual replicate deviations are shown as scatter points, with level
#' means as connected square markers. If ADL is provided, the acceptance
#' region is shown as a shaded ribbon with dashed boundary lines. For
#' percentage-based ADL plotted on the absolute scale, the ribbon forms a
#' fan shape that widens proportionally with concentration.
#'
#' **Residual plot** (`type = "residuals"`):
#' Displays residuals from the linear model against target concentration.
#' Systematic patterns suggest non-linearity.
#'
#' **Dilution plot** (`type = "dilution"`):
#' For dilution-based studies only. Displays observed values against the
#' proportion of the high pool, with individual replicate points, level means,
#' and a weighted least squares (WLS) regression line with equation annotation.
#' This provides a visual assessment of dilution recovery.
#'
#' @examples
#' set.seed(42)
#' target <- rep(c(10, 50, 100, 200, 300, 400, 500), each = 4)
#' measured <- target * rnorm(28, 1.0, 0.02) + rnorm(28, 0, 2)
#'
#' lin <- linearity_study(target, measured, adl = 5, adl_type = "percent")
#'
#' # Linearity plot
#' plot(lin)
#'
#' # Deviation plot (key diagnostic)
#' plot(lin, type = "deviation")
#'
#' # Residual plot
#' plot(lin, type = "residuals")
#'
#' @seealso [linearity_study()] for performing the analysis,
#'   [summary.linearity_study()] for detailed results
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_abline geom_hline
#'   geom_errorbar geom_bar geom_col geom_ribbon labs theme_bw theme
#'   element_text scale_color_manual scale_fill_manual coord_cartesian
#' @export
plot.linearity_study <- function(x,
                                 type = c("linearity", "deviation",
                                          "residuals", "dilution"),
                                 deviation_scale = c("auto", "absolute",
                                                     "percent"),
                                 show_points = TRUE,
                                 show_means = TRUE,
                                 show_identity = TRUE,
                                 point_alpha = 0.4,
                                 point_size = 2,
                                 line_colors = NULL,
                                 title = NULL,
                                 xlab = NULL,
                                 ylab = NULL,
                                 ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. ",
         "Please install it with install.packages('ggplot2').",
         call. = FALSE)
  }
  
  type <- match.arg(type)
  deviation_scale <- match.arg(deviation_scale)
  
  # Default colors
  default_colors <- c(
    linear = "#2166AC",
    polynomial = "#B2182B",
    identity = "#999999",
    adl = "#D6604D",
    zero = "#999999",
    mean = "#333333",
    pass = "#4DAF4A",
    fail = "#E41A1C"
  )
  
  if (is.null(line_colors)) {
    line_colors <- default_colors
  } else {
    line_colors <- modifyList(as.list(default_colors), as.list(line_colors))
    line_colors <- unlist(line_colors)
  }
  
  switch(type,
         linearity = .plot_linearity(x, show_points, show_means, show_identity,
                                     point_alpha, point_size, line_colors,
                                     title, xlab, ylab),
         deviation = .plot_deviation(x, deviation_scale, point_size, line_colors,
                                     title, xlab, ylab),
         residuals = .plot_linearity_residuals(x, show_points, point_alpha,
                                               point_size, line_colors,
                                               title, xlab, ylab),
         dilution = .plot_dilution_recovery(x, show_points, point_alpha,
                                            point_size, line_colors,
                                            title, xlab, ylab)
  )
}


# Plot Helper Functions ----

#' Linearity plot (measured vs target)
#' @noRd
.plot_linearity <- function(x, show_points, show_means, show_identity,
                            point_alpha, point_size, line_colors,
                            title, xlab, ylab) {
  
  input <- x$input
  results <- x$results
  
  # Raw data
  raw_data <- data.frame(target = input$x, measured = input$y)
  
  # Level means for error bars
  level_data <- input$level_stats
  level_data$se <- level_data$sd / sqrt(level_data$n)
  
  # Prediction lines
  x_seq <- seq(min(input$x), max(input$x), length.out = 200)
  pred_linear <- stats::predict(results$models$linear,
                                newdata = data.frame(x = x_seq))
  
  line_data <- data.frame(x = x_seq, linear = pred_linear)
  
  # Labels
  if (is.null(xlab)) xlab <- input$var_names["x"]
  if (is.null(ylab)) ylab <- input$var_names["y"]
  if (is.null(title)) title <- "Linearity Evaluation"
  
  # Build plot
  p <- ggplot2::ggplot()
  
  # Identity line
  if (show_identity) {
    p <- p +
      ggplot2::geom_abline(
        intercept = 0, slope = 1,
        color = line_colors["identity"],
        linetype = "dotted",
        linewidth = 0.5
      )
  }
  
  # Individual points
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        data = raw_data,
        ggplot2::aes(x = .data$target, y = .data$measured),
        alpha = point_alpha,
        size = point_size,
        color = "gray40"
      )
  }
  
  # Level means with error bars
  if (show_means) {
    p <- p +
      ggplot2::geom_errorbar(
        data = level_data,
        ggplot2::aes(x = .data$target,
                     ymin = .data$mean - .data$sd,
                     ymax = .data$mean + .data$sd),
        width = diff(range(input$x)) * 0.015,
        color = line_colors["mean"],
        linewidth = 0.5
      ) +
      ggplot2::geom_point(
        data = level_data,
        ggplot2::aes(x = .data$target, y = .data$mean),
        size = point_size + 1,
        color = line_colors["mean"],
        shape = 16
      )
  }
  
  # Linear fit
  p <- p +
    ggplot2::geom_line(
      data = line_data,
      ggplot2::aes(x = .data$x, y = .data$linear),
      color = line_colors["linear"],
      linewidth = 0.9,
      linetype = "solid"
    )
  
  # Polynomial fit (only if non-linear)
  if (results$best_model != "linear") {
    pred_poly <- stats::predict(results$models[[results$best_model]],
                                newdata = data.frame(x = x_seq))
    poly_data <- data.frame(x = x_seq, poly = pred_poly)
    
    p <- p +
      ggplot2::geom_line(
        data = poly_data,
        ggplot2::aes(x = .data$x, y = .data$poly),
        color = line_colors["polynomial"],
        linewidth = 0.9,
        linetype = "dashed"
      )
  }
  
  # Build subtitle
  subtitle_parts <- sprintf("n = %d, %d levels", input$n, input$n_levels)
  if (results$best_model != "linear") {
    subtitle_parts <- paste0(subtitle_parts,
                             sprintf(" | Best fit: %s", results$best_model))
  }
  
  # Labels and theme
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle_parts,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      axis.title = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      aspect.ratio = 1
    )
  
  p
}


#' Deviation from linearity plot
#' @noRd
.plot_deviation <- function(x, deviation_scale, point_size, line_colors,
                            title, xlab, ylab) {
  
  results <- x$results
  settings <- x$settings
  
  # Determine scale for y-axis ----
  # "auto": use percent if adl_type is percent, else absolute
  # "absolute": always plot absolute deviations
  # "percent": always plot percent deviations
  if (deviation_scale == "auto") {
    if (!is.null(settings$adl) && settings$adl_type == "percent") {
      use_pct <- FALSE  # Plot absolute deviations with fan-shaped ADL
    } else {
      use_pct <- FALSE  # Default to absolute
    }
  } else {
    use_pct <- (deviation_scale == "percent")
  }
  
  if (use_pct) {
    dev_col <- "deviation_pct"
    dev_label <- "Deviation from linearity (%)"
  } else {
    dev_col <- "deviation_abs"
    dev_label <- "Deviation from linearity"
  }
  
  adl_val <- settings$adl
  
  dev_data <- results$deviations
  dev_data$deviation <- dev_data[[dev_col]]
  
  # Determine pass/fail for coloring
  if (!is.null(adl_val)) {
    dev_data$status <- ifelse(abs(dev_data$deviation) <= adl_val,
                              "Within ADL", "Exceeds ADL")
  } else {
    dev_data$status <- "No ADL"
  }
  
  if (is.null(title)) title <- "Deviation from Linearity"
  if (is.null(xlab)) xlab <- x$input$var_names["x"]
  if (is.null(ylab)) ylab <- dev_label
  
  # Build plot
  p <- ggplot2::ggplot(dev_data,
                       ggplot2::aes(x = .data$target, y = .data$deviation))
  
  # ADL shaded region ----
  if (!is.null(adl_val)) {
    
    if (settings$adl_type == "percent" && !use_pct) {
      # Percent ADL plotted on absolute scale: fan-shaped region
      # ADL limits in absolute units = +/- (adl_pct / 100) * linear_predicted
      x_seq <- seq(min(dev_data$target), max(dev_data$target),
                   length.out = 100)
      linear_pred <- stats::predict(results$models$linear,
                                    newdata = data.frame(x = x_seq))
      adl_abs <- adl_val / 100 * abs(linear_pred)
      adl_ribbon <- data.frame(
        x = x_seq,
        ymin = -adl_abs,
        ymax = adl_abs
      )
    } else if (settings$adl_type == "percent" && use_pct) {
      # Percent ADL plotted on percent scale: flat horizontal
      adl_ribbon <- data.frame(
        x = range(dev_data$target),
        ymin = rep(-adl_val, 2),
        ymax = rep(adl_val, 2)
      )
    } else if (settings$adl_type == "absolute" && use_pct) {
      # Absolute ADL plotted on percent scale: inverse fan
      x_seq <- seq(min(dev_data$target), max(dev_data$target),
                   length.out = 100)
      linear_pred <- stats::predict(results$models$linear,
                                    newdata = data.frame(x = x_seq))
      adl_pct <- ifelse(abs(linear_pred) > .Machine$double.eps,
                        100 * adl_val / abs(linear_pred), NA_real_)
      adl_ribbon <- data.frame(
        x = x_seq,
        ymin = -adl_pct,
        ymax = adl_pct
      )
      # Remove NAs
      adl_ribbon <- adl_ribbon[complete.cases(adl_ribbon), ]
    } else {
      # Absolute ADL on absolute scale: flat horizontal
      adl_ribbon <- data.frame(
        x = range(dev_data$target),
        ymin = rep(-adl_val, 2),
        ymax = rep(adl_val, 2)
      )
    }
    
    p <- p +
      ggplot2::geom_ribbon(
        data = adl_ribbon,
        ggplot2::aes(x = .data$x, ymin = .data$ymin, ymax = .data$ymax),
        fill = line_colors["adl"],
        alpha = 0.15,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = adl_ribbon,
        ggplot2::aes(x = .data$x, y = .data$ymax),
        color = line_colors["adl"],
        linetype = "dashed",
        linewidth = 0.7,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = adl_ribbon,
        ggplot2::aes(x = .data$x, y = .data$ymin),
        color = line_colors["adl"],
        linetype = "dashed",
        linewidth = 0.7,
        inherit.aes = FALSE
      )
  }
  
  # Zero line
  p <- p +
    ggplot2::geom_hline(
      yintercept = 0,
      color = line_colors["zero"],
      linetype = "solid",
      linewidth = 0.5
    )
  
  # Individual replicate deviations (scatter) ----
  # Compute per-observation deviation from linear fit
  raw_x <- x$input$x
  raw_y <- x$input$y
  linear_pred_all <- stats::predict(results$models$linear,
                                    newdata = data.frame(x = raw_x))
  poly_pred_all <- stats::predict(results$models[[results$best_model]],
                                  newdata = data.frame(x = raw_x))
  
  if (use_pct) {
    # Individual replicate: deviation of observed from linear, in percent
    raw_dev <- ifelse(abs(linear_pred_all) > .Machine$double.eps,
                      100 * (raw_y - linear_pred_all) / linear_pred_all,
                      NA_real_)
  } else {
    raw_dev <- raw_y - linear_pred_all
  }
  
  raw_data <- data.frame(
    target = raw_x,
    deviation = raw_dev
  )
  
  p <- p +
    ggplot2::geom_point(
      data = raw_data,
      ggplot2::aes(x = .data$target, y = .data$deviation),
      size = point_size,
      alpha = 0.5,
      color = "grey40",
      shape = 16,
      inherit.aes = FALSE
    )
  
  # Level mean deviations (connected points) ----
  if (!is.null(adl_val)) {
    dev_data$color_group <- dev_data$status
    p <- p +
      ggplot2::geom_point(
        data = dev_data,
        ggplot2::aes(x = .data$target, y = .data$deviation,
                     color = .data$color_group),
        size = point_size + 1.5,
        shape = 15
      ) +
      ggplot2::geom_line(
        data = dev_data,
        ggplot2::aes(x = .data$target, y = .data$deviation),
        color = line_colors["mean"],
        linewidth = 0.6,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = c("Within ADL" = line_colors["pass"],
                   "Exceeds ADL" = line_colors["fail"]),
        name = NULL
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        data = dev_data,
        ggplot2::aes(x = .data$target, y = .data$deviation),
        size = point_size + 1.5,
        color = line_colors["mean"],
        shape = 15
      ) +
      ggplot2::geom_line(
        data = dev_data,
        ggplot2::aes(x = .data$target, y = .data$deviation),
        color = line_colors["mean"],
        linewidth = 0.6,
        inherit.aes = FALSE
      )
  }
  
  # Subtitle
  md <- results$max_deviation
  subtitle <- sprintf("Best model: %s | Max deviation: %s",
                      results$best_model,
                      if (!is.null(settings$adl) &&
                          settings$adl_type == "percent") {
                        sprintf("%.2f%%", md$max_pct)
                      } else {
                        sprintf("%.3f", md$max_abs)
                      })
  
  # Caption for ADL
  caption <- NULL
  if (!is.null(adl_val)) {
    adl_str <- if (settings$adl_type == "percent") {
      sprintf("+/-%.1f%%", adl_val)
    } else {
      sprintf("+/-%.3f", adl_val)
    }
    verdict <- if (isTRUE(results$linear)) "LINEAR" else "NON-LINEAR"
    caption <- sprintf("ADL: %s | Assessment: %s", adl_str, verdict)
  }
  
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0),
      axis.title = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (!is.null(adl_val)) "bottom" else "none"
    )
  
  p
}


#' Residual plot for linearity
#' @noRd
.plot_linearity_residuals <- function(x, show_points, point_alpha, point_size,
                                      line_colors, title, xlab, ylab) {
  
  input <- x$input
  results <- x$results
  
  # Residuals from linear model
  resid_vals <- stats::residuals(results$models$linear)
  
  raw_data <- data.frame(
    target = input$x,
    residual = resid_vals
  )
  
  # Level means of residuals
  levels <- sort(unique(round(input$x, digits = 10)))
  level_resid <- data.frame(
    target = numeric(0),
    mean_resid = numeric(0)
  )
  for (lev in levels) {
    idx <- which(abs(input$x - lev) < .Machine$double.eps * 100)
    level_resid <- rbind(level_resid, data.frame(
      target = lev,
      mean_resid = mean(resid_vals[idx])
    ))
  }
  
  if (is.null(title)) title <- "Residuals from Linear Fit"
  if (is.null(xlab)) xlab <- x$input$var_names["x"]
  if (is.null(ylab)) ylab <- "Residual"
  
  p <- ggplot2::ggplot()
  
  # Zero line
  p <- p +
    ggplot2::geom_hline(
      yintercept = 0,
      color = line_colors["zero"],
      linetype = "dashed",
      linewidth = 0.6
    )
  
  # Individual points
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        data = raw_data,
        ggplot2::aes(x = .data$target, y = .data$residual),
        alpha = point_alpha,
        size = point_size,
        color = "gray40"
      )
  }
  
  # Level mean residuals
  p <- p +
    ggplot2::geom_point(
      data = level_resid,
      ggplot2::aes(x = .data$target, y = .data$mean_resid),
      size = point_size + 1,
      color = line_colors["polynomial"],
      shape = 18
    ) +
    ggplot2::geom_line(
      data = level_resid,
      ggplot2::aes(x = .data$target, y = .data$mean_resid),
      color = line_colors["polynomial"],
      linewidth = 0.6,
      alpha = 0.7
    )
  
  # Smooth to detect trend
  if (nrow(raw_data) >= 10) {
    p <- p +
      ggplot2::geom_smooth(
        data = raw_data,
        ggplot2::aes(x = .data$target, y = .data$residual),
        method = "loess",
        formula = y ~ x,
        se = FALSE,
        color = line_colors["linear"],
        linewidth = 0.7,
        alpha = 0.5
      )
  }
  
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d | Best model: %s",
                         input$n, results$best_model),
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


# Dilution Recovery Plot ----

#' Dilution recovery plot for linearity studies
#' @noRd
.plot_dilution_recovery <- function(x, show_points, point_alpha, point_size,
                                    line_colors, title, xlab, ylab) {
  
  input <- x$input
  results <- x$results
  
  # Check that this is a dilution-based study
  if (is.null(input$dilution)) {
    stop("Dilution plot requires a dilution-based linearity study. ",
         "Use `dilution` parameter in linearity_study().", call. = FALSE)
  }
  
  dilution_vec <- input$dilution
  raw_y <- input$y
  
  # Use the stored dilution-scale WLS fit
  dil_fit <- results$dilution_fit
  intercept_wls <- dil_fit$intercept
  slope_wls <- dil_fit$slope
  
  # Compute level means from raw data for plotting
  dil_rounded <- round(dilution_vec, 8)
  unique_dils <- sort(unique(dil_rounded))
  
  level_data <- data.frame(
    dilution = numeric(length(unique_dils)),
    mean = numeric(length(unique_dils))
  )
  for (i in seq_along(unique_dils)) {
    idx <- which(dil_rounded == unique_dils[i])
    level_data$dilution[i] <- unique_dils[i]
    level_data$mean[i] <- mean(raw_y[idx], na.rm = TRUE)
  }
  
  # Individual replicate data
  raw_data <- data.frame(
    dilution = dilution_vec,
    measured = raw_y
  )
  
  # Regression line data
  dil_seq <- seq(0, 1, length.out = 100)
  line_data <- data.frame(
    dilution = dil_seq,
    predicted = intercept_wls + slope_wls * dil_seq
  )
  
  # Labels
  if (is.null(title)) title <- "Dilution Recovery"
  if (is.null(xlab)) xlab <- "Proportion of HIGH pool"
  if (is.null(ylab)) ylab <- "Mean"
  
  # Equation label
  eq_label <- sprintf("y = %.2f + %.2f * x",
                      intercept_wls, slope_wls)
  
  # Build plot
  p <- ggplot2::ggplot()
  
  # Individual replicate points
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        data = raw_data,
        ggplot2::aes(x = .data$dilution, y = .data$measured),
        alpha = point_alpha,
        size = point_size,
        color = "grey40",
        shape = 16
      )
  }
  
  # Level means
  p <- p +
    ggplot2::geom_point(
      data = level_data,
      ggplot2::aes(x = .data$dilution, y = .data$mean),
      size = point_size + 1.5,
      color = line_colors["mean"],
      shape = 16
    )
  
  # WLS regression line
  p <- p +
    ggplot2::geom_line(
      data = line_data,
      ggplot2::aes(x = .data$dilution, y = .data$predicted),
      color = line_colors["linear"],
      linewidth = 0.8,
      linetype = "dashed"
    )
  
  # Equation annotation
  p <- p +
    ggplot2::annotate(
      "text",
      x = 0.02,
      y = max(raw_data$measured, na.rm = TRUE) * 0.98,
      label = eq_label,
      hjust = 0,
      vjust = 1,
      size = 3.5,
      color = line_colors["linear"],
      fontface = "italic"
    )
  
  # Theme and labels
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(round(level_data$dilution, 4))),
      limits = c(-0.02, 1.02)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d observations, %d levels | WLS regression line",
                         input$n, input$n_levels),
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


#' @rdname plot.linearity_study
#' @param object An object of class `linearity_study`.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.linearity_study <- function(object,
                                     type = c("linearity", "deviation",
                                              "residuals", "dilution"),
                                     deviation_scale = c("auto", "absolute",
                                                         "percent"),
                                     show_points = TRUE,
                                     show_means = TRUE,
                                     show_identity = TRUE,
                                     point_alpha = 0.4,
                                     point_size = 2,
                                     line_colors = NULL,
                                     title = NULL,
                                     xlab = NULL,
                                     ylab = NULL,
                                     ...) {
  plot.linearity_study(
    x = object,
    type = type,
    deviation_scale = deviation_scale,
    show_points = show_points,
    show_means = show_means,
    show_identity = show_identity,
    point_alpha = point_alpha,
    point_size = point_size,
    line_colors = line_colors,
    title = title,
    xlab = xlab,
    ylab = ylab,
    ...
  )
}