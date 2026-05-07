# ===============================================================
# BSA-Seq Plotting Functions
# ===============================================================
#
# Description:
# This script provides modular functions for visualization of
# bulked-segregant analysis sequencing (BSA-Seq) data. It supports
# SNP-based and window-based analyses and produces publication-quality
# chromosome-scale plots.
#
# Features:
# - Line plots with optional smoothing (rolling median or locfit)
# - Point plots for SNP-level statistics (e.g., -log10 p-values)
# - Histogram plots for SNP density
# - Automatic chromosome ordering (natural sort)
# - Automatic detection of pre-smoothed window metrics (_wmd, _rmd, _lft)
# - High-resolution export for publication figures
#
# Dependencies:
# ggplot2, dplyr, rlang, zoo, scales, stringr, grid
# Optional: locfit (required for smooth_type = "locfit")
#
# Input requirements:
# - data must contain columns: CHROM, POS
# - POS is assumed to be in base pairs
#
# Notes:
# - POS is converted to megabases (Mb) for plotting
# - Window-based metrics are not smoothed again to avoid over-smoothing
#'
#'
make_plot_theme <- function(base_size = 60, legendit = "none") {
  require(grid)
  line_scale <- base_size / 20
  theme_linedraw(base_size = base_size) + 
    theme(
      plot.title = element_text(size = base_size, face = "bold", 
                                hjust = 0.5, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.7, 
                                 face = "bold", color = "black"),
      axis.text.y = element_text(size = base_size, face = "bold", 
                                 colour = "black", hjust = 1),
      axis.title.x = element_text(size = base_size, face = "bold", color = "black"),
      axis.title.y = element_text(size = base_size, face = "bold", 
                                  angle = 90, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 1.5 * line_scale),
      axis.ticks = element_line(colour = "black", linewidth = 2 * line_scale),
      axis.ticks.length = unit(0.5, "cm"),
      panel.border = element_rect(linewidth = 3 * line_scale, fill = NA),
      panel.spacing.x = unit(2, "lines"),
      panel.spacing.y = unit(2, "lines"),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 2),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.0),
      strip.background = element_rect(fill = "grey80", colour = "black", 
                                      linewidth = 3 * line_scale),
      strip.text = element_text(colour = "black", face = "bold", size = base_size),
      strip.text.x = element_text(colour = "black", face = "bold", size = base_size, 
                                  margin = margin(0.6, 0.0, 0.6, 0.0, "cm")),
      strip.text.y = element_text(colour = "black", face = "bold", size = base_size, 
                                  margin = margin(0.0, 0.6, 0.0, 0.6, "cm")), 
      legend.position = legendit
    )
}

#' Create base BSA-Seq plot
#'
#' Builds faceted ggplot with chromosome panels, titles, and axes.
#' @param data Data containing CHROM and POS.
#' @param prefix Character. Reference name.
#' @param y_title Character. Y-axis label.
#' @param plot_title Character. Plot title.
#' @param facet_column Integer. Number of columns for facet_wrap.
#' @param plot_style Character. "wrap" or "grid".
#' @param remove_x_text Logical. Remove x-axis labels.
#' @return List with plot, plot_style, is_wrap, and inbred name.
#' @export
#' 
make_plot_base <- function(data, prefix, y_title, plot_title, facet_column = 5, 
                           plot_style = c("wrap", "grid"), remove_x_text = FALSE) {
  
  # Helper mini-function : Capitalize the first letter 
  capitalize_first <- function(text) {
    substr(text, 1, 1) <- toupper(substr(text, 1, 1))
    text
  }
  
  # Capitalize inbred name from prefix
  inbred <- capitalize_first(prefix)
  
  # plot_style
  plot_style <- match.arg(plot_style)
  is_wrap <- plot_style == "wrap"
  
  # Select faceting style
  facet_layer <- if (is_wrap) {
    facet_wrap(~ CHROM, ncol = facet_column, scales = "free_x")
  } else {
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x")
  }
  
  # Shared labels + theme bits
  title_txt <- paste0("Aligned to ", inbred, " :  ", plot_title)
  xlab_txt  <- "Chromosome Position (Mb)\n"
  ylab_txt  <- paste0("\n", y_title, "\n")
  
  axis_x_theme <- theme(axis.text.x  = if (remove_x_text) element_blank() else 
    element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = -1)),
    axis.ticks.x = if (remove_x_text) element_blank() else element_line(),
    axis.title.x = element_text(color = "black"))
  
  p <- ggplot(data) + facet_layer + 
    labs(title = title_txt, x = xlab_txt, y = ylab_txt) + make_plot_theme() + 
    axis_x_theme
  return(list(plot = p, plot_style = plot_style, is_wrap = is_wrap, inbred = inbred))
}

#' Add line layer to BSA-Seq plot
#' Adds a chromosome-colored line layer with optional smoothing.
#' @param plot ggplot object from make_plot_base().
#' @param column Character. Metric column to plot.
#' @param color_panel Character vector of colors.
#' @param smooth_type Character. "rollmedian", "locfit", or "none".
#' @param nn_prop Numeric. Locfit neighborhood proportion.
#' @param rollmedian Integer. Rolling median window size.
#' @param line_size Numeric. Line width.
#' @param threshold Optional numeric threshold.
#' @details
#' locfit smoothing requires the locfit package.
#' @return ggplot object with line layer.
#' @export

make_line_plot <- function(plot, column, color_panel, threshold = NULL,
                           smooth_type = c("rollmedian", "locfit", "none"), 
                           nn_prop = 0.1, rollmedian = 501, line_size = 4) {
  
  smooth_type <- match.arg(smooth_type)
  p <- plot + aes(x = PositionMb, y = !!rlang::sym(column), color = CHROM) +
    scale_color_manual(values = color_panel) + guides(color = "none")
  
  if (smooth_type == "locfit") {
      if (!requireNamespace("locfit", quietly = TRUE)) {
        stop("locfit package required.")
        }
    p <- p + stat_smooth(method = "locfit", formula = y ~ lp(x, nn = nn_prop), 
                         linewidth = line_size)
  } else if (smooth_type == "rollmedian") {
    p <- p + geom_line(aes(y = zoo::rollmedian(!!rlang::sym(column), 
                                               k = rollmedian, na.pad = TRUE)), 
                       linewidth = line_size)
  } else {
    p <- p + geom_line(linewidth = line_size)
  }
  
  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed",
                        color = "black", linewidth = line_size * 0.6)
  }
  return(p)
}

#' Add point layer to BSA-Seq plot
#' @param plot ggplot object from make_plot_base().
#' @param column Character. Metric column.
#' @param color_panel Character vector of colors.
#' @param point_size Numeric. Point size.
#' @param line_size Numeric. Line width for threshold.
#' @param threshold Optional numeric threshold.
#' @return ggplot object with point layer.
#' @export
make_point_plot <- function(plot, column, color_panel, point_size = 4, 
                            line_size = 4, threshold = NULL) {
  
  p <- plot + aes(x = PositionMb, y = !!rlang::sym(column), color = CHROM) +
    geom_point(size = point_size) + scale_color_manual(values = color_panel) +
    guides(color = "none")
  
  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed", 
                        color = "black", linewidth = line_size * 0.6)
  }
  return(p)
}

#' Add histogram layer to BSA-Seq plot
#' @param plot ggplot object from make_plot_base().
#' @param color_panel Character vector of colors.
#' @param bwidth Numeric. Bin width in base pairs.
#' @param alpha_size Numeric. Transparency.
#' @param line_size Numeric. Border thickness.
#' @return ggplot object with histogram layer.
#' @export
make_histogram_plot <- function(
    plot, color_panel, bwidth = 1000000, alpha_size = 1, line_size = 4) {
  
  bscale <- bwidth / 1e6
  p <- plot + aes(x = PositionMb, fill = CHROM) +
    geom_histogram(binwidth = bscale, alpha = alpha_size, 
                   linewidth = line_size * 0.15) + 
    scale_fill_manual(values = color_panel) + guides(fill = "none") + 
    scale_y_continuous(
      labels = scales::label_number(scale = 1e-3, accuracy = 0.01, trim = TRUE))
  return(p)
}

#' Save BSA-Seq plot to file
#'
#' @param plot ggplot object.
#' @param plots_dir Output directory.
#' @param inbred Character. Prefix used in filename.
#' @param file_suffix Character. File suffix.
#' @param plot_style Character. "wrap" or "grid".
#' @param device File format.
#' @param hwidth,hheight Dimensions for wrap.
#' @param width,height Dimensions for grid.
#' @param dpi Resolution.
#' @param is_wrap Logical.
#' @return NULL (invisible). Saves file to disk.
#' @export
save_bsa_plot <- function(plot, plots_dir, inbred, file_suffix, plot_style,
                          device = "png", hwidth = 30, hheight = 18,
                          width = 45, height = 15, dpi = 300, is_wrap = TRUE) {
  
  # Ensure plots directory exists
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Construct file path and save plot
  date_tag <- format(Sys.Date(), "%Y%m%d")
  
  file_path   <- file.path(
    plots_dir, paste0(date_tag, "_", inbred, "_", 
                      file_suffix, "_", plot_style, ".", device))
  plot_width  <- if (is_wrap) hwidth else width
  plot_height <- if (is_wrap) hheight else height
  
  ggsave(filename = file_path, plot = plot, device = device,
         width = plot_width, height = plot_height, dpi = dpi)
  
  message(paste0("Plot saved to: ", file_path))
}

#' Plot BSA-Seq signals across chromosomes
#'
#' High-level wrapper for generating BSA-Seq visualizations from SNP-based
#' or window-based data. Supports line, point, and histogram plots with
#' optional smoothing for line plots.
#'
#' @param data A data.frame or data.table containing at least CHROM and POS.
#' @param prefix Character. Reference name used in titles and filenames.
#' @param column Character. Name of metric column to plot (required for
#'   plot_type = "line" or "points").
#' @param y_title Character. Y-axis label.
#' @param plot_title Character. Plot title.
#' @param file_suffix Character. Output file suffix.
#' @param plot_type Character. One of "line", "points", or "histogram".
#' @param ylim Optional numeric vector specifying y-axis limits.
#' @param smooth_type Character. One of "rollmedian", "locfit", or "none".
#' @param threshold Optional numeric threshold displayed as a dashed line.
#' @param bwidth Numeric. Histogram bin width in base pairs.
#' @param rollmedian Integer. Window size for rolling median smoothing.
#' @param nn_prop Numeric. Neighborhood proportion for locfit smoothing.
#' @param point_size Numeric. Size of points for point plots.
#' @param line_size Numeric. Line width.
#' @param alpha_size Numeric. Histogram transparency.
#' @param facet_column Integer. Number of facet columns.
#' @param plot_style Character. One of "wrap" or "grid".
#' @param remove_x_text Logical. Whether to remove x-axis text.
#' @param color_panel Character vector of colors for chromosomes.
#' @param plots_dir Character. Output directory.
#' @param device Character. File format (e.g. "png", "pdf").
#' @param hwidth,hheight Numeric. Dimensions for wrapped plots.
#' @param width,height Numeric. Dimensions for grid plots.
#' @param dpi Numeric. Resolution.
#' @details
#' Position (POS) is converted to megabases (Mb) for plotting. Chromosomes
#' are ordered using natural sorting. Window-based metrics (ending in
#' "_wmd", "_rmd", "_lft") are treated as pre-smoothed and are not
#' smoothed again.
#' @return Invisibly returns a ggplot object. Plot is saved to disk.
#'
#' @export
#' 
plot_bsa <- function(
    data, prefix, column = NULL, y_title, plot_title, file_suffix,
    plot_type = c("line", "points", "histogram"), ylim = NULL, 
    smooth_type = c("rollmedian", "locfit", "none"), threshold = NULL, 
    bwidth = 1000000, rollmedian = 501, nn_prop = 0.1, point_size = 4, line_size = 4, 
    alpha_size = 1, facet_column = 5, plot_style = c("wrap", "grid"),
    remove_x_text = FALSE, color_panel = c("blue", "red"),
    plots_dir, device = "png", hwidth = 30, hheight = 18,
    width = 49, height = 13, dpi = 300) {
  
  plot_type <- match.arg(plot_type)
  smooth_type <- match.arg(smooth_type)
  
  # check required columns
  if (!all(c("CHROM", "POS") %in% names(data))) {
    stop("data must contain 'CHROM' and 'POS' columns.")
  }
  
  if (plot_type != "histogram") {
    if (is.null(column)) {
      stop("'column' must be provided for plot_type 'line' or 'points'.")
    }
    if (!column %in% names(data)) {
      stop(paste("column", column, "not found in data."))
    }
  }
  
  # Add Mb position for plotting
  data <- dplyr::mutate(data, PositionMb = POS / 1e6)
  data$CHROM <- factor(data$CHROM,
    levels = stringr::str_sort(unique(data$CHROM), numeric = TRUE)
  )
  
  # Prevent double smoothing for window-based line plots
  if (plot_type == "line" && grepl("_(wmd|rmd|lft)$", column)) {
    if (smooth_type != "none") {
      message("window-based column provide: ", column)
    }
    smooth_type <- "none"
  }
  
  # chromosome colors to match number of chromosomes
  n_chr <- length(unique(data$CHROM))
  color_panel <- rep(color_panel, length.out = n_chr)
  
  # Build base plot
  base <- make_plot_base(data = data, prefix = prefix, y_title = y_title,
                         plot_title = plot_title, facet_column = facet_column,
                         plot_style = plot_style, remove_x_text = remove_x_text
                         )
  
  p <- base$plot
  
  # Shared arguments
  shared_args <- list(plot = p, color_panel = color_panel)
  
  # Select plot function and add only type-specific arguments
  if (plot_type == "histogram") {
    plot_fun <- make_histogram_plot
    plot_args <- c(shared_args, 
                   list(bwidth = bwidth, alpha_size = alpha_size,
                        line_size = line_size)
    )
    
  } else if (plot_type == "points") {
    plot_fun <- make_point_plot
    plot_args <- c(shared_args, 
                   list(column = column, point_size = point_size,
                        line_size = line_size, threshold = threshold)
    )
    
  } else {
    plot_fun <- make_line_plot
    plot_args <- c(shared_args, 
                   list(column = column, smooth_type = smooth_type, 
                        nn_prop = nn_prop, rollmedian = rollmedian, 
                        line_size = line_size, threshold = threshold)
    )
  }
  
  # Add plot layer
  p <- do.call(plot_fun, plot_args)
  
  # Apply y limits if provided
  if (!is.null(ylim)) {p <- p + coord_cartesian(ylim = ylim)}
  
  # Save plot
  save_bsa_plot(
    plot = p, plots_dir = plots_dir, inbred = base$inbred, file_suffix = file_suffix,
    plot_style = base$plot_style, device = device, hwidth = hwidth, hheight = hheight, 
    width = width, height = height, dpi = dpi, is_wrap = base$is_wrap
  )
  
  invisible(p)
}