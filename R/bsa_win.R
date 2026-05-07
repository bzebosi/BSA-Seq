#' Sliding-window engine for genomic data
#' Applies a user-defined function to each genomic window.
#' @param data data frame or data.table.
#' @param window_size Window size in base pairs.
#' @param step_size Step size between adjacent windows in base pairs.
#' @param fun function applied to each window. Must return a data.table or NULL.
#' @return A data.table with one row per window, including standard
#'   window coordinates and any columns returned by fun.
#' @export
slide_windows <- function(data, window_size = 2e6, step_size = 1e5, fun) {
  
  if (window_size <= 0 || step_size <= 0) {
    stop("window_size and step_size must be > 0")
  }
  
  datax <- data.table::as.data.table(data)
  
  if (!all(c("CHROM", "POS") %in% names(datax))) {
    stop("Missing required columns: CHROM and POS")
  }
  
  datax[, POS := as.integer(POS)]
  datax <- datax[is.finite(POS)]
  
  if (nrow(datax) == 0L) {
    return(data.table::data.table())
  }
  
  chr_levels <- naturalsort(unique(as.character(datax$CHROM)))
  results <- list()
  
  for (chr in chr_levels) {
    chr_snps <- datax[CHROM == chr]
    data.table::setorder(chr_snps, POS)
    
    chr_len <- max(chr_snps$POS, na.rm = TRUE)
    half <- floor(window_size / 2)
    
    centers <- if (chr_len <= window_size) {
      as.integer(chr_len / 2)
    } else {
      seq(from = half, to = chr_len - half, by = step_size)
    }
    
    out_rows <- lapply(centers, function(center) {
      center <- as.integer(center)
      start <- as.integer(max(1L, center - half))
      end   <- as.integer(min(chr_len, center + half))
      
      win   <- chr_snps[POS >= start & POS <= end]
      
      extra <- fun(win, chr = chr, start = start, end = end, center = center)
      
      if (is.null(extra)) {return(NULL)}
      
      if (!data.table::is.data.table(extra)) {
        stop("fun must return a data.table or NULL")
      }
      
      data.table::data.table(CHROM = chr, start = start, end = end, 
                             POS = center,N = as.integer(nrow(win)), extra)
    })
    
    out_rows <- out_rows[lengths(out_rows) > 0]
    
    if (length(out_rows) > 0L) {
      results[[length(results) + 1L]] <- data.table::rbindlist(
        out_rows, use.names = TRUE, fill = TRUE)}
  }
  
  if (!length(results)) {
    out <- data.table::data.table()
  } else {
    out <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
    out[, CHROM := factor(as.character(CHROM), levels = chr_levels)]
    data.table::setorder(out, CHROM, POS)
  }
  out
}

#' Interval detection
#' finds intervals around the highest peak in each score column within each chromosome.
#' @param data A data frame or data.table containing CHROM, start, end, and POS.
#' @param offhold Fraction of peak score used as cutoff for interval expansion.
#' @param min_vsize Minimum interval width in base pairs.
#' @param use_cols Character vector of score columns to scan.
#' @return data.table of peak-centered intervals and summary statistics.
#' @export
peak_interval <- function(data, offhold=0.80, min_vsize = 0L, use_cols = c()) {
  
  datax <- data.table::as.data.table(data)
  need <- c("CHROM","start","end","POS")
  
  if (!all(need %in% names(datax))) stop("Input needs CHROM, start, end, POS")
  cols_use <- intersect(use_cols, names(datax))
  if (!length(cols_use)) stop("None of use_cols found in data")
  
  # empty result (data.table)
  out <- data.table::data.table(
    CHROM = character(), score_col = character(), cutoff = numeric(), 
    start = integer(), end = integer(), width_bp = integer(), peak_pos = integer(), 
    peak_score = numeric(), mean_score = numeric(), max_score = numeric(), 
    min_score = numeric(), area = numeric()
  )
  
  results <- list()
  for (chr in unique(datax$CHROM)) {
    dt <- datax[CHROM == chr][order(POS)]
    if (!nrow(dt)) next
    n <- nrow(dt)
    
    for (col in cols_use) {
      s <- dt[[col]]
      if (all(is.na(s))) next
      
      # peak and cutoff
      s2 <- ifelse(is.na(s), -Inf, s)
      peak_idx <- which.max(s2)
      peak_val <- s2[peak_idx]
      
      if (!is.finite(peak_val)) next
      cutoff <- peak_val * offhold
      
      # expand from peak while score >= cutoff
      L <- R <- peak_idx
      while (L > 1L && s2[L - 1L] >= cutoff) L <- L - 1L
      while (R < n  && s2[R + 1L] >= cutoff) R <- R + 1L
      
      # build segment + size filter
      seg <- dt[L:R]
      seg_start <- as.integer(min(seg$start, na.rm = TRUE))
      seg_end   <- as.integer(max(seg$end,   na.rm = TRUE))
      width_bp  <- as.integer(seg_end - seg_start + 1L)
      if (width_bp < as.integer(min_vsize)) next
      
      # collect
      mean_sc <- as.numeric(mean(seg[[col]], na.rm = TRUE))
      max_sc  <- as.numeric(max(seg[[col]],  na.rm = TRUE))
      min_sc  <- as.numeric(min(seg[[col]],  na.rm = TRUE))
      results[[length(results) + 1L]] <- data.table::data.table(
        CHROM = chr, score_col = col, cutoff = cutoff, start = seg_start, 
        end = seg_end, width_bp = width_bp, peak_pos = as.integer(dt$POS[peak_idx]), 
        peak_score = as.numeric(s[peak_idx]), mean_score = mean_sc,
        max_score  = max_sc, min_score  = min_sc, area = mean_sc * width_bp)
    }
  }
  
  if (!length(results)) return(out)
  
  data.table::rbindlist(results, use.names = TRUE)
}

#' Sliding-window summary statistics for BSA-Seq data
#'
#' Computes sliding-window summary tracks from metric or homozygosity data.
#'
#' @param data A data frame or data.table containing CHROM, POS,
#'   and required columns.
#' @param metric Name of numeric metric column.
#' @param af_col Name of allele-frequency column.
#' @param window_size Window size in base pairs.
#' @param step_size Step size between adjacent windows in base pairs.
#' @param rollmedian Rolling median window size.
#' @param nn_prop Smoothing parameter for locfit.
#' @param af_min Threshold for defining homozygous SNPs.
#' @param find_intervals Logical; if TRUE, return peak-based intervals.
#' @param offhold Fraction of peak score used for interval expansion.
#' @param min_vsize Minimum interval size in base pairs.
#' @param use_col Which tracks to return: "wmd", "rmd", "lft", or "all".
#' @return A list with $windows, and optionally $intervals.
#' @export
window_compute <- function(
    data, metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    only_mutant = FALSE, window_size = 2e6, step_size = 1e5,
    rollmedian = 101L, nn_prop = 0.1, af_min = 0.99, find_intervals = FALSE, 
    offhold = 0.90, min_vsize = 1e6, use_col = "all") {
  
  metric <- match.arg(metric)
  valid_tracks <- c("wmd", "lft", "rmd")
  tracks <- tolower(use_col)
  tracks <- if ("all" %in% tracks) valid_tracks else intersect(tracks, valid_tracks)
  if (!length(tracks)) tracks <- "wmd"
  
  if ("rmd" %in% tracks && !requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required when use_col includes 'rmd'.")
  }
  
  if ("lft" %in% tracks && !requireNamespace("locfit", quietly = TRUE)) {
    stop("Package 'locfit' is required when use_col includes 'lft'.")
  }
  
  datax <- data.table::as.data.table(data)
  datax[, POS := as.integer(POS)]
  
  if (!all(c("CHROM", "POS") %in% names(datax))) {
    stop("Missing required columns: CHROM and POS.")
  }
  
  all_chr <- naturalsort(unique(as.character(datax$CHROM)))
  
  run_one_metric <- function(datax, colname, base_name, is_homo = FALSE) {
    
    if (!(colname %in% names(datax))) return(NULL)
    
    if (!is_homo && !is.numeric(datax[[colname]])) {
      stop(paste(colname, "must be numeric."))
    }
    
    if (is_homo) {
      af <- datax[[colname]]
      datax[, Hom := as.integer(af >= af_min | af <= (1 - af_min))]
      
      stat_fun <- function(win, chr, start, end, center) {
        data.table::data.table(
          Stat = sum(win$Hom, na.rm = TRUE) / nrow(win))
      }
      
    } else {
      stat_fun <- function(win, chr, start, end, center) {
        data.table::data.table(
          Stat = stats::median(win[[colname]], na.rm = TRUE))
      }
    }
    
    out <- slide_windows(data = datax, window_size = window_size, 
                         step_size = step_size, fun = stat_fun)
    
    if (nrow(out) == 0L) return(NULL)
    
    out[, CHROM := factor(as.character(CHROM), levels = all_chr)]
    data.table::setorder(out, CHROM, POS)
    
    stat_col <- paste0(base_name, "_wmd")
    rmd_col  <- paste0(base_name, "_rmd")
    lft_col  <- paste0(base_name, "_lft")
    
    data.table::setnames(out, "Stat", stat_col)
    
    if ("rmd" %in% tracks) {
      k <- max(1L, as.integer(rollmedian))
      out[, (rmd_col) := zoo::rollapply(
        get(stat_col), width = k, FUN = stats::median,
        na.rm = TRUE, fill = NA_real_, align = "center"), by = CHROM]
    }
    
    if ("lft" %in% tracks) {
      out[, (lft_col) := {
        x <- POS
        y <- get(stat_col)
        keep <- is.finite(x) & is.finite(y)
        pred <- rep(NA_real_, .N)
        
        if (sum(keep) >= 5L) {
          df_fit <- data.frame(x = x[keep], y = y[keep])
          df_fit <- df_fit[order(df_fit$x), , drop = FALSE]
          df_fit <- stats::aggregate(y ~ x, data = df_fit, FUN = stats::median)
          
          if (nrow(df_fit) >= 5L) {
            fit <- tryCatch(
              locfit::locfit.raw(df_fit$x, df_fit$y, alpha = c(nn = nn_prop)), 
              error = function(e) NULL)
            
            if (!is.null(fit)) {
              pred <- tryCatch(as.numeric(stats::predict(fit, newdata = x)),
                               error = function(e) rep(NA_real_, .N))
            }
          }
        }
        pred <- if (is_homo) pmin(1, pmax(0, pred)) else pmax(0, pred)
        pred
      }, by = CHROM]
    }
    
    if (!isTRUE(find_intervals)) {
      return(list(windows = out))
    }
    
    if (!exists("peak_interval", mode = "function")) {
      stop("peak_interval() is not available.")
    }
    
    candidate_cols <- paste0(base_name, "_", tracks)
    candidate_cols <- candidate_cols[candidate_cols %in% names(out)]
    if (!length(candidate_cols)) candidate_cols <- stat_col
    
    intervals <- peak_interval(data = out, offhold = offhold, min_vsize = min_vsize,
                               use_cols = candidate_cols)
    
    list(windows = out, intervals = intervals)
  }
  
  if (metric %in% c("AF", "homozygosity")) {
    
    af_cols <- if (only_mutant) "mt_AF" else c("wt_AF", "mt_AF")
    out <- list()
    
    for (col in af_cols) {
      base_name <- if (metric == "homozygosity") paste0("homozygosity_", col) else col
      
      res <- run_one_metric(datax = data.table::copy(datax), colname = col,
                            base_name = base_name, is_homo = metric == "homozygosity")
      
      if (!is.null(res)) out[[col]] <- res
    }
    return(out)
  }
  
  if (only_mutant && metric %in% c("AFD", "ED", "ED4", "G")) {
    message(metric, " requires both WT and mutant data. Skipping.")
    return(NULL)
  }
  
  run_one_metric(datax = datax, colname = metric, base_name = metric, is_homo = FALSE)
}
#' Sliding-window Fisher test for BSA-Seq data
#' Performs sliding-window enrichment testing using Fisher's exact test on
#' metric-based or homozygosity-based signals.
#' @param data data frame or data.table containing CHROM, POS, and required metric columns.
#' @param metric_col name of numeric metric column for stat_type = "metric".
#' @param af_col name of allele-frequency column for stat_type = "homozygosity".
#' @param window_size Window size in base pairs.
#' @param step_size Step size between adjacent windows in base pairs.
#' @param doorstep Threshold for defining high versus low metric values.
#' @param af_min Threshold for defining homozygous SNPs from allele frequency.
#' @param threshold significance threshold on the -10 * log10(adjusted p-value) scale.
#' @param find_intervals Logical; if TRUE, return peak-based intervals.
#' @param offhold Fraction of peak score used for interval expansion.
#' @param min_vsize Minimum interval size in base pairs.
#' @param use_cols Score column(s) passed to peak_interval().
#' @return A list with $windows, and optionally $intervals.
#' @export
window_pval <- function(
    data, metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    only_mutant = FALSE, window_size = 2e6, step_size = 1e5,
    doorstep = NULL, af_min = 0.99, threshold = -10 * log10(0.05),
    find_intervals = FALSE, offhold = 0.80, min_vsize = 0L,
    use_cols = "log.pval") {
  
  metric <- match.arg(metric)
  
  datax <- data.table::as.data.table(data)
  
  if (!all(c("CHROM", "POS") %in% names(datax))) {
    stop("Missing columns: CHROM and POS")
  }
  
  run_one_pval <- function(datax, colname, is_homo = FALSE) {
    
    if (!(colname %in% names(datax))) return(NULL)
    
    if (!is_homo && !is.numeric(datax[[colname]])) {
      stop(paste(colname, "must be numeric."))
    }
    
    if (is_homo) {
      
      af <- datax[[colname]]
      datax[, Hom := as.integer(af >= af_min | af <= (1 - af_min))]
      
      bg_high <- sum(datax$Hom == 1L, na.rm = TRUE)
      bg_low  <- sum(datax$Hom == 0L, na.rm = TRUE)
      
      count_fun <- function(win) {
        list(high = sum(win$Hom == 1L, na.rm = TRUE),
             low  = sum(win$Hom == 0L, na.rm = TRUE))
      }
    } else {
      
      if (is.null(doorstep)) {
        stop("doorstep is required for metric p-value analysis.")
      }
      
      bg_high <- sum(datax[[colname]] > doorstep, na.rm = TRUE)
      bg_low  <- sum(datax[[colname]] <= doorstep, na.rm = TRUE)
      
      count_fun <- function(win) {
        list(high = sum(win[[colname]] > doorstep, na.rm = TRUE),
             low  = sum(win[[colname]] <= doorstep, na.rm = TRUE))
      }
    }
    
    fisher_window <- function(win, chr, start, end, center) {
      counts <- count_fun(win)
      
      bg_out_high <- bg_high - counts$high
      bg_out_low  <- bg_low  - counts$low
      
      mat <- matrix(c(counts$high, counts$low, bg_out_high, bg_out_low), 
                    nrow = 2, byrow = TRUE)
      
      p_val <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
      
      data.table::data.table(high = counts$high, low = counts$low,
                             bg_high = bg_out_high, bg_low = bg_out_low,
                             pval = p_val)
    }
    
    out <- slide_windows(data = datax, window_size = window_size,
                         step_size = step_size, fun = fisher_window)
    
    if (nrow(out) > 0L) {
      out[, adj.pval := p.adjust(pval, method = "bonferroni")]
      out[, log.pval := -10 * log10(ifelse(adj.pval == 0, 1e-10, adj.pval))]
      out[, sig := log.pval > threshold]
    }
    
    if (isTRUE(find_intervals)) {
      if (!exists("peak_interval", mode = "function")) {
        stop("peak_interval() is not available")
      }
      
      intervals <- peak_interval(data = out, offhold = offhold,
                                 min_vsize = min_vsize, use_cols = use_cols)
      
      return(list(windows = out, intervals = intervals))
    }
    
    list(windows = out)
  }
  
  if (metric %in% c("AF", "homozygosity")) {
    
    af_cols <- if (only_mutant) "mt_AF" else c("wt_AF", "mt_AF")
    out <- list()
    
    for (col in af_cols) {
      
      res <- run_one_pval(datax = data.table::copy(datax), colname = col,
                          is_homo = metric == "homozygosity")
      
      if (!is.null(res)) out[[col]] <- res
    }
    
    return(out)
  }
  
  if (only_mutant && metric %in% c("AFD", "ED", "ED4", "G")) {
    message(metric, " requires both WT and mutant data. Skipping.")
    return(NULL)
  }
  
  run_one_pval(datax = datax, colname = metric, is_homo = FALSE)
}
#' Sliding-window BSA analysis
#' Wrapper for running sliding-window summaries and optional enrichment
#' p-value analysis.
#' @param data SNP-level data frame.
#' @param stat_type Analysis type: `"metric"`, `"homozygosity"`, or `"all"`.
#' @param metric_col Metric column name (e.g. `"ED4"`). Required for metric analysis.
#' @param af_col Allele frequency column name. Required for homozygosity analysis.
#' @param window_size Sliding-window size.
#' @param step_size Sliding-window step size.
#' @param nn_prop Neighbor proportion used for smoothing.
#' @param af_min Homozygosity cutoff. Default is `0.99`.
#' @param rollmedian Rolling median window size.
#' @param doorstep Metric threshold for enrichment analysis.
#' @param threshold Significance threshold for p-value analysis.
#' @param run_pval Logical; if `TRUE`, run `window_pval()`.
#' @param find_intervals Logical; if `TRUE`, identify significant intervals.
#' @param offhold Interval boundary threshold.
#' @param min_vsize Minimum interval size.
#' @param use_col Output column mode passed to `window_compute()`.
#' @return A list containing `compute` and optional `pval` results.
#' @export
window_analysis <- function(
    data, metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    only_mutant = FALSE, window_size = 2e6, step_size = 1e5, nn_prop = 0.1, 
    af_min = 0.99, rollmedian = 101L, doorstep = NULL, run_pval = TRUE,
    threshold = -10 * log10(0.05), find_intervals = FALSE, offhold = 0.90, 
    min_vsize = 1e6, use_col = "all") {
  
  metric <- match.arg(metric)
  
  if (run_pval && metric %in% c("AF","AFD", "ED", "ED4", "G") &&
      is.null(doorstep)) {
    stop("doorstep is required for metric p-value analysis.")
  }
  
  common_args <- list(
    data = data, metric = metric, only_mutant = only_mutant, 
    window_size = window_size, step_size = step_size, af_min = af_min,
    find_intervals = find_intervals, offhold = offhold, min_vsize = min_vsize)
  
  out_compute <- do.call(
    window_compute, c(common_args, list(rollmedian = rollmedian,
                                        nn_prop = nn_prop, use_col = use_col)))
  
  out_pval <- NULL
  
  if (run_pval) {
    
    pval_args <- c(common_args, list(threshold = threshold, 
                                     use_cols = "log.pval"))
    
    if (metric %in% c("AF", "AFD", "ED", "ED4", "G")) {
      pval_args$doorstep <- doorstep
    }
    
    out_pval <- do.call(window_pval, pval_args)
  }
  
  list(compute = out_compute, pval = out_pval)
}

run_window_analysis <- function(
    data, metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    only_mutant = FALSE, window_size = 2e6, step_size = 1e5, nn_prop = 0.1,
    af_min = 0.99, rollmedian = 101L, doorstep = NULL, threshold = -10 * log10(0.05),
    run_pval = TRUE, find_intervals = FALSE, offhold = 0.90, min_vsize = 1e6,
    use_col = "all") {
  
  metric <- match.arg(metric, several.ok = TRUE)
  out <- list()
  
  for (m in metric) {
    
    # allow one doorstep or metric-specific doorsteps
    ds <- if (is.list(doorstep)) doorstep[[m]] else doorstep
    
    if (run_pval && m %in% c("AFD", "ED", "ED4", "G") && is.null(ds)) {
      message(m, " skipped: doorstep is required for p-value analysis.")
      next
    }
    
    out[[m]] <- window_analysis(
      data = data, metric = m, only_mutant = only_mutant,
      window_size = window_size, step_size = step_size, nn_prop = nn_prop,
      af_min = af_min, rollmedian = rollmedian, doorstep = ds, threshold = threshold, 
      run_pval = run_pval, find_intervals = find_intervals, offhold = offhold,
      min_vsize = min_vsize, use_col = use_col)
  }
  out
}

make_window_parameters <- function(
    data,
    wt = "wildtype", mt = "mutant",
    metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    only_mutant = FALSE, use_ems = FALSE,
    use_pval = TRUE,
    use_cols = c("wmd", "rmd", "lft", "log.pval")) {
  
  metric <- match.arg(metric, several.ok = TRUE)
  
  snp <- if (use_ems) "ems snps only" else "all snps"
  v <- if (use_ems) "ems" else "all"
  
  outp <- list()
  
  add_params <- function(tbl, cols, y_title, plot_title, plotid_base) {
    
    cols <- cols[cols %in% names(tbl)]
    
    for (col in cols) {
      outp[[paste(plotid_base, col, sep = "_")]] <<- list(
        column = col,
        data = tbl,
        y_title = y_title,
        plot_title = plot_title,
        plotid = paste(plotid_base, col, sep = "_")
      )
    }
  }
  
  for (m in metric) {
    
    if (m %in% c("AF", "homozygosity")) {
      
      gens <- if (only_mutant) "mt_AF" else c("wt_AF", "mt_AF")
      
      for (g in gens) {
        
        who <- if (g == "mt_AF") mt else wt
        base <- if (m == "homozygosity") paste0("homozygosity_", g) else g
        
        if (!is.null(data[[m]]$compute[[g]]$windows)) {
          
          tbl <- data[[m]]$compute[[g]]$windows
          cols <- paste0(base, "_", use_cols[use_cols != "log.pval"])
          
          add_params(
            tbl = tbl,
            cols = cols,
            y_title = m,
            plot_title = sprintf("%s : %s (%s)", m, who, snp),
            plotid_base = sprintf("%s_%s_%s_window", who, m, v)
          )
        }
        
        if (isTRUE(use_pval) && !is.null(data[[m]]$pval[[g]]$windows)) {
          
          tbl <- data[[m]]$pval[[g]]$windows
          
          add_params(
            tbl = tbl,
            cols = "log.pval",
            y_title = "-log10(p-value)",
            plot_title = sprintf("%s enrichment : %s (%s)", m, who, snp),
            plotid_base = sprintf("%s_%s_%s_pval", who, m, v)
          )
        }
      }
      
      next
    }
    
    if (!is.null(data[[m]]$compute$windows)) {
      
      tbl <- data[[m]]$compute$windows
      cols <- paste0(m, "_", use_cols[use_cols != "log.pval"])
      
      add_params(
        tbl = tbl,
        cols = cols,
        y_title = m,
        plot_title = sprintf("%s : %s - %s (%s)", m, wt, mt, snp),
        plotid_base = sprintf("%s_%s_%s_%s_window", wt, mt, m, v)
      )
    }
    
    if (isTRUE(use_pval) && !is.null(data[[m]]$pval$windows)) {
      
      tbl <- data[[m]]$pval$windows
      
      add_params(
        tbl = tbl,
        cols = "log.pval",
        y_title = "-log10(p-value)",
        plot_title = sprintf("%s enrichment : %s - %s (%s)", m, wt, mt, snp),
        plotid_base = sprintf("%s_%s_%s_%s_pval", wt, mt, m, v)
      )
    }
  }
  
  outp
}

run_window_metrics <- function(
    data, wt = "wildtype", mt = "mutant",
    metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    prefix = "sample", plots_dir = "plots", only_mutant = FALSE,
    plot_data = TRUE, plot_style = c("wrap", "grid"), ylim = NULL,
    threshold = NULL, line_size = 4, bwidth = 1e6,
    color_panel = c("blue", "red"), device = "png",
    width = 45, height = 13, hwidth = 30, hheight = 18, dpi = 300,
    use_ems = TRUE, use_pval = TRUE,
    use_cols = c("wmd", "rmd", "lft", "log.pval")) {
  
  metric <- match.arg(metric, several.ok = TRUE)
  plot_style <- match.arg(plot_style)
  
  subfolder <- if (only_mutant || is.null(wt)) {
    paste(prefix, mt, sep = "_")
  } else {
    paste(prefix, wt, mt, sep = "_")
  }
  
  plots_dir <- file.path(plots_dir, "window_plots", subfolder)
  ensure_dir(plots_dir)
  
  parameters <- list()
  
  if (!is.null(data$all)) {
    parameters <- c(
      parameters,
      make_window_parameters(
        data = data$all, wt = wt, mt = mt, metric = metric,
        only_mutant = only_mutant, use_ems = FALSE,
        use_pval = use_pval, use_cols = use_cols
      )
    )
  }
  
  if (isTRUE(use_ems) && !is.null(data$ems)) {
    parameters <- c(
      parameters,
      make_window_parameters(
        data = data$ems, wt = wt, mt = mt, metric = metric,
        only_mutant = only_mutant, use_ems = TRUE,
        use_pval = use_pval, use_cols = use_cols
      )
    )
  }
  
  if (length(parameters) == 0) {
    message("No window plotting data found.")
    return(invisible(NULL))
  }
  
  if (plot_data) {
    
    line_parameters <- parameters[!grepl("_pval", names(parameters))]
    pval_parameters <- parameters[grepl("_pval", names(parameters))]
    
    if (length(line_parameters) > 0) {
      plot_parameter_list(
        parameters = line_parameters, prefix = prefix, plots_dir = plots_dir,
        plot_type = "line", plot_mode = "none", plot_style = plot_style,
        ylim = ylim, threshold = NULL, line_size = line_size,
        bwidth = bwidth, color_panel = color_panel, device = device,
        width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi
      )
    }
    
    if (length(pval_parameters) > 0) {
      plot_parameter_list(
        parameters = pval_parameters, prefix = prefix, plots_dir = plots_dir,
        plot_type = "points", plot_mode = "none", plot_style = plot_style,
        ylim = ylim, threshold = threshold, line_size = line_size,
        bwidth = bwidth, color_panel = color_panel, device = device,
        width = width, height = height, hwidth = hwidth, hheight = hheight, dpi = dpi
      )
    }
  }
  
  message("Window plotting completed.")
  invisible(parameters)
}



bsa_window <- function(
    vcf_dir, prefix = "b73", pattern = "snps\\.tsv$",
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE, plot_data = FALSE, 
    Genotypes = list(wt = "wildtype", mt = "mutant"),
    save_results = FALSE, save_interval = FALSE, output_dir = NULL,
    metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    use_ems = TRUE, window_size = 2e6, step_size = 1e5,
    nn_prop = 0.1, af_min = 0.99, rollmedian = 101L,
    doorstep = NULL, threshold = -10 * log10(0.05),
    run_pval = TRUE, find_intervals = FALSE,
    offhold = 0.90, min_vsize = 1e6, use_col = "all",
    plot_style = c("grid", "wrap"), ylim = NULL,
    color_panel = c("blue", "red"), device = "png",
    width = 48, height = 13, hwidth = 30, hheight = 18,
    dpi = 300, line_size = 4, bwidth = 1e6,
    use_pval = TRUE,
    use_cols = c("wmd", "rmd", "lft", "log.pval")) {
  
  metric <- match.arg(metric, several.ok = TRUE)
  plot_style <- match.arg(plot_style)
  
  if (is.null(output_dir)) {
    project_dir <- dirname(dirname(vcf_dir))
    output_dir <- file.path(project_dir, "post_analysis")
  }
  
  wt_label <- Genotypes$wt
  mt_label <- Genotypes$mt
  
  subfolder <- if (only_mutant || is.null(wt_label)) {
    paste(prefix, mt_label, sep = "_")
  } else {
    paste(prefix, wt_label, mt_label, sep = "_")
  }
  
  message("Step 1: Reading SNP tables")
  
  geno_data <- read_vcf(
    vcf_dir = vcf_dir, prefix = prefix, pattern = pattern,
    Genotypes = Genotypes, min_DP = min_DP,
    min_QUAL = min_QUAL, only_mutant = only_mutant)
  
  message("Step 2: Computing BSA statistics")
  
  tables_dir <- file.path(output_dir, "tables", subfolder)
  ensure_dir(tables_dir)
  
  res <- compute_bsa(
    geno_data = geno_data, prefix = prefix,
    save_results = save_results, output_dir = tables_dir,
    only_mutant = only_mutant
  )
  
  message("Step 3: Running window-based analysis")
  
  data_all <- if (only_mutant) res$mt_all else res$wt_mt_all
  data_ems <- if (only_mutant) res$mt_ems else res$wt_mt_ems
  
  common_args <- list(
    metric = metric, only_mutant = only_mutant, window_size = window_size,
    step_size = step_size, nn_prop = nn_prop, af_min = af_min,
    rollmedian = rollmedian, doorstep = doorstep, threshold = threshold,
    run_pval = run_pval, find_intervals = find_intervals, offhold = offhold,
    min_vsize = min_vsize, use_col = use_col)
  
  win <- list()
  
  win$all <- do.call(
    run_window_analysis, c(list(data = data_all), common_args))
  
  if (use_ems && !is.null(data_ems)) {
    win$ems <- do.call(
      run_window_analysis,
      c(list(data = data_ems), common_args)
    )
  }
  
  if (plot_data) {
    
    message("Step 4: Plotting window-based metrics")
    
    plot_out <- file.path(output_dir, "plots")
    ensure_dir(plot_out)
    
    win$plot_results <- run_window_metrics(
      data = win, wt = wt_label, mt = mt_label, metric = metric, 
      prefix = prefix, plots_dir = plot_out, only_mutant = only_mutant,
      plot_data = TRUE, plot_style = plot_style, ylim = ylim, threshold = threshold,
      line_size = line_size, bwidth = bwidth, color_panel = color_panel,
      device = device, width = width, height = height, hwidth = hwidth,
      hheight = hheight, dpi = dpi, use_ems = use_ems, use_pval = use_pval,
      use_cols = use_cols)
  }
  
  if (isTRUE(save_interval) && isTRUE(find_intervals)) {
    
    message("Step 5: Saving interval Excel")
    
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      
      intervals_dir <- file.path(output_dir, "intervals", subfolder)
      ensure_dir(intervals_dir)
      
      file_name <- if (only_mutant) {
        sprintf("%s_%s_window_bsa_intervals.xlsx", prefix, mt_label)
      } else {
        sprintf("%s_%s_vs_%s_window_bsa_intervals.xlsx", prefix, wt_label, mt_label)
      }
      
      wb <- openxlsx::createWorkbook()
      
      add_interval_sheets <- function(x, label = "all") {
        
        for (m in names(x)) {
          
          obj <- x[[m]]
          
          if (m %in% c("AF", "homozygosity")) {
            
            for (g in names(obj$compute)) {
              
              iv <- obj$compute[[g]]$intervals
              
              if (!is.null(iv) && is.data.frame(iv) && nrow(iv) > 0L) {
                sheet <- substr(paste(label, m, g, sep = "_"), 1, 31)
                openxlsx::addWorksheet(wb, sheet)
                openxlsx::writeData(wb, sheet, iv)
              }
            }
            
          } else {
            
            iv <- obj$compute$intervals
            
            if (!is.null(iv) && is.data.frame(iv) && nrow(iv) > 0L) {
              sheet <- substr(paste(label, m, sep = "_"), 1, 31)
              openxlsx::addWorksheet(wb, sheet)
              openxlsx::writeData(wb, sheet, iv)
            }
          }
        }
      }
      
      add_interval_sheets(win$all, "all")
      
      if (!is.null(win$ems)) {
        add_interval_sheets(win$ems, "ems")
      }
      
      xlsx_path <- file.path(intervals_dir, file_name)
      
      if (length(wb$worksheets) > 0L) {
        openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
        message("Saved intervals to: ", xlsx_path)
      } else {
        message("No interval rows to save.")
      }
      
    } else {
      warning("save_interval = TRUE but 'openxlsx' is not installed. Skipping Excel write.")
    }
  }
  
  win$bsa_results <- res
  
  invisible(win)
}

bsa_window_all <- function(
    vcf_dir, prefix_list, wt_list = NULL, mt_list, pattern = "snps\\.tsv$", 
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE, save_results = FALSE, 
    save_interval = FALSE, output_dir = NULL,
    metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    use_ems = TRUE, window_size = 2e6, step_size = 1e5,
    nn_prop = 0.1, af_min = 0.99, rollmedian = 101L,
    doorstep = NULL, threshold = -10 * log10(0.05),
    run_pval = TRUE, find_intervals = FALSE,
    offhold = 0.90, min_vsize = 1e6, use_col = "all",
    plot_data = FALSE, plot_style = c("grid", "wrap"), ylim = NULL,
    color_panel = c("blue", "red"), device = "png",
    width = 48, height = 13, hwidth = 30, hheight = 18,
    dpi = 300, line_size = 4, bwidth = 1e6,
    use_pval = TRUE,
    use_cols = c("wmd", "rmd", "lft", "log.pval")) {
  
  metric <- match.arg(metric, several.ok = TRUE)
  plot_style <- match.arg(plot_style)
  
  if (!only_mutant && is.null(wt_list)) {
    stop("wt_list is required when only_mutant = FALSE.")
  }
  
  results <- list()
  seen_combos <- character()
  
  if (is.null(output_dir)) {
    project_dir <- dirname(dirname(vcf_dir))
    output_dir <- file.path(project_dir, "post_analysis")
  }
  
  ensure_dir(output_dir)
  
  common_args <- list(
    vcf_dir = vcf_dir,
    pattern = pattern,
    min_DP = min_DP,
    min_QUAL = min_QUAL,
    save_results = save_results,
    save_interval = save_interval,
    output_dir = output_dir,
    metric = metric,
    use_ems = use_ems,
    window_size = window_size,
    step_size = step_size,
    nn_prop = nn_prop,
    af_min = af_min,
    rollmedian = rollmedian,
    doorstep = doorstep,
    threshold = threshold,
    run_pval = run_pval,
    find_intervals = find_intervals,
    offhold = offhold,
    min_vsize = min_vsize,
    use_col = use_col,
    plot_data = plot_data,
    plot_style = plot_style,
    ylim = ylim,
    color_panel = color_panel,
    device = device,
    width = width,
    height = height,
    hwidth = hwidth,
    hheight = hheight,
    dpi = dpi,
    line_size = line_size,
    bwidth = bwidth,
    use_pval = use_pval,
    use_cols = use_cols
  )
  
  for (pfx in prefix_list) {
    
    message("\n==== Running window-based BSA pipeline for prefix: ", pfx, " ====")
    
    if (only_mutant) {
      
      for (mt in mt_list) {
        
        combo_key <- paste(pfx, mt, sep = "_")
        if (combo_key %in% seen_combos) next
        seen_combos <- c(seen_combos, combo_key)
        
        message("Running mutant-only: ", combo_key)
        
        rex <- do.call(
          bsa_window,
          c(common_args, list(
            prefix = pfx,
            only_mutant = TRUE,
            Genotypes = list(wt = NULL, mt = mt)
          ))
        )
        
        results[[combo_key]] <- rex
      }
      
    } else {
      
      for (wt in wt_list) {
        for (mt in mt_list) {
          
          if (wt == mt) next
          
          pair_name <- paste(sort(c(wt, mt)), collapse = "_")
          combo_key <- paste(pfx, pair_name, sep = "_")
          
          if (combo_key %in% seen_combos) next
          seen_combos <- c(seen_combos, combo_key)
          
          message("Running: ", combo_key,
                  "  (wildtype = ", wt, ", mutant = ", mt, ")")
          
          rex <- do.call(
            bsa_window,
            c(common_args, list(
              prefix = pfx,
              only_mutant = FALSE,
              Genotypes = list(wt = wt, mt = mt)
            ))
          )
          
          results[[combo_key]] <- rex
        }
      }
    }
  }
  
  message("\n==== bsa_window_all finished ====")
  invisible(results)
}