make_metric_parameters <- function(
    wt_mt = NULL, wt_mt_ems = NULL, ant_wt = NULL, ant_mt = NULL,
    ant_wt_ems = NULL, ant_mt_ems = NULL, wt = "wildtype", mt = "mutant",
    metric = c("AFD", "ED", "ED4", "G", "AF"), only_mutant = FALSE,
    use_ems = TRUE) {
  
  # metrics to generate
  metric <- match.arg(metric, several.ok = TRUE)
  vars <- if (use_ems) c("all", "ems") else "all"
  
  # input datasets
  tables <- list(wt_mt_all = wt_mt, wt_mt_ems = wt_mt_ems,
                 wt_all = ant_wt, wt_ems = ant_wt_ems,
                 mt_all = ant_mt, mt_ems = ant_mt_ems)
  
  outp <- list()
  for (m in metric) {
    if (m == "AF") {
      gens <- if (only_mutant) "mt" else c("wt", "mt")
      for (g in gens) {
        for (v in vars) {
          id <- paste(g, v, sep = "_")
          col <- if (g == "mt") "mt_AF" else "wt_AF"
          who <- if (g == "mt") mt else wt
          tbl <- check_cols(tables[[id]], col)
          if (is.null(tbl)) next
          snp <- if (v == "ems") "ems snps only" else "all snps"
          outp[[paste("AF", g, v, sep = "_")]] <- list(
            column = col, data = tbl, y_title = "AF",
            plot_title = sprintf("%s AF (%s)", who, snp),
            plotid = sprintf("%s_AF_%s", who, v)
          )
        }
      }
      next
    }
    
    # comparative metrics require WT + MT
    if (only_mutant) next
    for (v in vars) {
      id <- paste("wt_mt", v, sep = "_")
      tbl <- check_cols(tables[[id]], m)
      if (is.null(tbl)) next
      snp <- if (v == "ems") "ems snps only" else "all snps"
      outp[[paste(m, v, sep = "_")]] <- list(
        column = m, data = tbl, y_title = m,
        plot_title = sprintf("%s : %s - %s (%s)", m, wt, mt, snp),
        plotid = sprintf("%s_%s_%s_%s", wt, mt, m, v)
      )
    }
  }
  outp
}

make_hist_parameters <- function(
    ant_wt = NULL, ant_mt = NULL, ant_wt_ems = NULL, ant_mt_ems = NULL,
    wt = "wildtype", mt = "mutant", only_mutant = FALSE,
    use_ems = TRUE, af_min = 0.99) {
  
  gens <- if (only_mutant) "mt" else c("mt", "wt")
  vars <- if (use_ems) c("all", "ems") else "all"

  tables <- list(mt_all = ant_mt, mt_ems = ant_mt_ems,
                 wt_all = ant_wt, wt_ems = ant_wt_ems)
  
  outp <- list()
  for (v in vars) {
    for (g in gens) {
      id <- paste(g, v, sep = "_")
      col <- if (g == "mt") "mt_AF" else "wt_AF"
      tbl <- check_cols(tables[[id]], col, min_value = af_min)
      if (is.null(tbl)) next
      who <- if (g == "mt") mt else wt
      snp <- if (v == "ems") "ems snps only" else "all snps"
      outp[[id]] <- list(
        column = col, data = tbl, y_title = "SNPs / Mb (×10³)",
        plot_title = sprintf("%s unique %s", who, snp),
        plotid = sprintf("%s_AF_unique_%s", who, v)
      )
    }
  }
  outp
}


plot_parameter_list <- function(
    parameters, prefix = "sample", plots_dir = "plots",
    plot_type = c("line", "points", "histogram"),
    plot_mode = c("both", "rollmedian", "locfit", "none"),
    plot_style = c("wrap", "grid"), ylim = NULL, threshold = NULL, 
    rollmedian = 501, nn_prop = 0.1, line_size = 4, bwidth = 1e6, 
    color_panel = c("blue", "red"), device = "png", width = 45, height = 13, 
    hwidth = 30, hheight = 18, dpi = 300) {
  
  plot_type <- match.arg(plot_type)
  plot_mode <- match.arg(plot_mode)
  plot_style <- match.arg(plot_style)
  
  # Histograms do not use smoothing
  smooth_types <- if (plot_type == "histogram") {
    "none"
  } else {
    switch(plot_mode,
           both = c("rollmedian", "locfit"),
           rollmedian = "rollmedian",
           locfit = "locfit", none = "none")
  }
  
  # Shared plotting arguments
  common_args <- list(
    prefix = prefix, plot_type = plot_type, ylim = ylim, threshold = threshold, 
    rollmedian = rollmedian, nn_prop = nn_prop, line_size = line_size,
    bwidth = bwidth, plot_style = plot_style, color_panel = color_panel,
    plots_dir = plots_dir, device = device, width = width, height = height,
    hwidth = hwidth, hheight = hheight, dpi = dpi)
  
  for (pmt in parameters) {
    for (sm in smooth_types) {
      
      message(sprintf("Generating plot: %s", pmt$plot_title))
      
      # build output filename
      file_suffix <- if (plot_type == "histogram") {
        sprintf("histogram_%s", pmt$plotid)
      } else {
        sprintf("%s_%s", pmt$plotid, sm)
      }
      
      args <- c(
        pmt[names(pmt) != "plotid"], common_args,
        list(file_suffix = file_suffix, smooth_type = sm)
      )
      do.call(plot_bsa, args)
    }
  }
  invisible(parameters)
}

run_snp_metrics <- function(
    data, wt = "wildtype", mt = "mutant",
    metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    prefix = "sample", plots_dir = "plots", only_mutant = FALSE, 
    plot_data = TRUE, plot_mode = c("both", "rollmedian", "locfit", "none"),
    plot_style = c("wrap", "grid"), ylim = NULL, threshold = NULL,
    rollmedian = 501, nn_prop = 0.1, line_size = 4, bwidth = 1e6,
    color_panel = c("blue", "red"), device = "png", width = 45, 
    height = 13, hwidth = 30, hheight = 18, dpi = 300,
    use_ems = TRUE, af_min = 0.99) {
  
  result <- data 
  metric <- match.arg(metric, several.ok = TRUE)
  
  subfolder <- if (only_mutant || is.null(wt)) {
    paste(prefix, mt, sep = "_")
  } else {
    paste(prefix, wt, mt, sep = "_")
  }
  
  plots_dir <- file.path(plots_dir, "snp_plots", subfolder)
  
  ensure_dir(plots_dir)
  
  common_args <- if (only_mutant) {
    list(ant_mt = result$mt_all, ant_mt_ems = result$mt_ems, wt = wt, mt = mt, 
         only_mutant = only_mutant, use_ems = use_ems)
  } else {
    list(ant_wt = result$ant_wt_all, ant_mt = result$ant_mt_all,
         ant_wt_ems = result$ant_wt_ems, ant_mt_ems = result$ant_mt_ems,
         wt = wt, mt = mt, only_mutant = only_mutant, use_ems = use_ems)
  }
  
  plot_args <- list(
    prefix = prefix, plots_dir = plots_dir, plot_style = plot_style, ylim = ylim,
    bwidth = bwidth, color_panel = color_panel, device = device, width = width,
    height = height, hwidth = hwidth, hheight = hheight, dpi = dpi
  )
  
  out <- list()
  line_metrics <- metric[metric != "homozygosity"]
  
  if (length(line_metrics) > 0) {
    
    out$metric <- do.call(
      make_metric_parameters,
      c(common_args,
        if (only_mutant) {
          list(metric = line_metrics)
        } else {
          list(wt_mt = result$wt_mt_all, wt_mt_ems = result$wt_mt_ems,
               metric = line_metrics)
        }
      )
    )
    
    if (length(out$metric) > 0 && plot_data) {
      do.call(plot_parameter_list, 
              c(list(parameters = out$metric, plot_type = "line", 
                     plot_mode = plot_mode, threshold = threshold, 
                     rollmedian = rollmedian, nn_prop = nn_prop,
                     line_size = line_size), plot_args))
    }
  }
  
  if ("homozygosity" %in% metric) {
    
    out$homozygosity <- do.call(make_hist_parameters,
                                c(common_args, list(af_min = af_min)))
    
    if (length(out$homozygosity) > 0 && plot_data) {
      do.call(plot_parameter_list, 
              c(list(parameters = out$homozygosity, plot_type = "histogram"), 
                plot_args))
    }
  }
  message("SNP plotting completed.")
  invisible(out)
}


bsa_snp <- function(
    vcf_dir, prefix = "b73", pattern = "snps\\.tsv$", min_DP = 5, min_QUAL = 5, 
    only_mutant = FALSE, plot_data = FALSE, save_results = FALSE, output_dir = NULL,
    Genotypes = list(wt = "wildtype", mt = "mutant"),
    metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    plot_mode = c("both", "rollmedian", "locfit", "none"),
    plot_style = c("grid", "wrap"), use_ems = TRUE, af_min = 0.99, nn_prop = 0.1,
    rollmedian = 501, ylim = NULL, color_panel = c("blue", "red"), device = "png", 
    width = 48, height = 13, hwidth = 30, hheight = 18, dpi = 300, bwidth = 1e6) {
  
  metric <- match.arg(metric, several.ok = TRUE)
  plot_mode <- match.arg(plot_mode)
  plot_style <- match.arg(plot_style)
  
  if (is.null(output_dir)) {
    project_dir <- dirname(dirname(vcf_dir))
    output_dir <- file.path(project_dir, "post_analysis")
  }
  
  wt_label <- Genotypes$wt
  mt_label <- Genotypes$mt
  
  message("Step 1: Reading SNP tables")
  
  geno_data <- read_vcf(vcf_dir = vcf_dir, prefix = prefix, pattern = pattern,
                        Genotypes = Genotypes, min_DP = min_DP, 
                        min_QUAL = min_QUAL, only_mutant = only_mutant)
  
  message("Step 2: Computing BSA statistics")
  
  tables_dir <- file.path(output_dir, "tables")
  ensure_dir(tables_dir)
  
  res <- compute_bsa(geno_data = geno_data, prefix = prefix, 
                     save_results = save_results, output_dir = tables_dir,
                     only_mutant = only_mutant)
  
  if (plot_data) {
    
    message("Step 3: Plotting SNP-based metrics")
    
    plot_out <- file.path(output_dir, "plots")
    ensure_dir(plot_out)
    
    res$plot_results <- run_snp_metrics(
      data = res, wt = wt_label, mt = mt_label, metric = metric, prefix = prefix,
      plots_dir = plot_out, only_mutant = only_mutant, plot_data = TRUE,
      plot_mode = plot_mode, plot_style = plot_style, ylim = ylim, 
      rollmedian = rollmedian, nn_prop = nn_prop, bwidth = bwidth, 
      color_panel = color_panel, device = device, width = width, 
      height = height, hwidth = hwidth, hheight = hheight, dpi = dpi,
      use_ems = use_ems, af_min = af_min)
  }
  
  invisible(res)
}



bsa_snp_all <- function(
    vcf_dir, prefix_list, wt_list = NULL, mt_list, pattern = "snps\\.tsv$", 
    min_DP = 5, min_QUAL = 5, only_mutant = FALSE, save_results = FALSE, 
    output_dir = NULL, metric = c("AF", "AFD", "ED", "ED4", "G", "homozygosity"),
    plot_mode = c("both", "rollmedian", "locfit", "none"), nn_prop = 0.1, 
    af_min = 0.99, bwidth = 1e6, plot_style = c("grid", "wrap"), 
    use_ems = TRUE, plot_data = FALSE, rollmedian = 501, ylim = NULL,
    color_panel = c("blue", "red"), device = "png", width = 48, 
    height = 13, hwidth = 30, hheight = 18, dpi = 300) {
  
  metric <- match.arg(metric, several.ok = TRUE)
  plot_mode <- match.arg(plot_mode)
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
    vcf_dir = vcf_dir, pattern = pattern, min_DP = min_DP, min_QUAL = min_QUAL,
    save_results = save_results, output_dir = output_dir, metric = metric,
    plot_mode = plot_mode, plot_style = plot_style, use_ems = use_ems,
    plot_data = plot_data, rollmedian = rollmedian, ylim = ylim, 
    color_panel = color_panel, device = device, width = width,
    height = height, hwidth = hwidth, hheight = hheight, dpi = dpi,
    nn_prop = nn_prop, af_min = af_min, bwidth = bwidth
  )
  
  for (pfx in prefix_list) {
    
    message("\n==== Running SNP-based BSA pipeline for prefix: ", pfx, " ====")
    
    if (only_mutant) {
      
      for (mt in mt_list) {
        
        combo_key <- paste(pfx, mt, sep = "_")
        if (combo_key %in% seen_combos) next
        seen_combos <- c(seen_combos, combo_key)
        
        message("Running mutant-only: ", combo_key)
        
        rex <- do.call(
          bsa_snp, c(common_args, list(prefix = pfx, only_mutant = TRUE, 
                                       Genotypes = list(wt = NULL, mt = mt)))
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
            bsa_snp, c(common_args, list(prefix = pfx, only_mutant = FALSE,
                                         Genotypes = list(wt = wt, mt = mt)))
          )
          
          results[[combo_key]] <- rex
        }
      }
    }
  }
  
  message("\n==== bsa_snp_all finished ====")
  invisible(results)
}