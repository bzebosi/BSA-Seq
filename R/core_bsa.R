#' Read SNP tables and compute allele frequencies for BSA-Seq
#' This function reads SNP tables for wild-type and mutant bulks, filters SNPs 
#' based on depth and quality, and computes allele frequency (AF) and related 
#' summary statistics for downstream BSA-Seq analysis.
#' @param vcf_dir Directory containing SNP tables
#' @param prefix Prefix used to identify sample files
#' @param pattern File pattern to match SNP files (e.g., "snps\\.tsv$")
#' @param Genotypes Named list of genotype labels (wt, mt)
#' @param min_DP Minimum read depth threshold
#' @param min_QUAL Minimum quality threshold
#' @param only_mutant Logical; if TRUE, load only mutant data
#' @return A list containing filtered SNP tables for each genotype
#' @details
#' Input files must contain at least 10 columns:
#' CHROM, POS, REF, ALT, QUAL, DP, Fref, Rref, Falt, Ralt
#' @export 

read_vcf <- function(vcf_dir, prefix, pattern, 
                     Genotypes = list(wt = "wildtype", mt = "mutant"), 
                     min_DP = 5, min_QUAL = 5, only_mutant = FALSE) {
  
  # Find matching files
  vcf_list <- list.files(path = vcf_dir, pattern = pattern, full.names = TRUE)
  if (length(vcf_list) == 0) stop("No vcf files found.")
  
  # decide whether to load only mutant or both genotypes
  selected_genotypes <- if (only_mutant) {
    list(mt = Genotypes[["mt"]])
  } else {
    Genotypes[c("wt", "mt")]
  }
  
  geno_data <- list()
  
  # Loop through genotypes
  for (genotype in names(selected_genotypes)) {
    file_pattern <- paste0("^", prefix, "_", selected_genotypes[[genotype]])
    sample_name <- paste(prefix, selected_genotypes[[genotype]], sep = "_")
    matched_file <- vcf_list[grepl(file_pattern, basename(vcf_list))]
    
    if (length(matched_file) == 0) {
      stop(paste("File for", genotype, "not found."))
    }
    if (length(matched_file) > 1) {
      stop(paste0("Multiple files found for ", genotype, ":\n",
                  paste(basename(matched_file), collapse = "\n"))
      )
    }
    geno_file <- matched_file[1]
    message("Reading file for ", basename(geno_file), " (", genotype, ").")
    
    # Read file
    data <- data.table::fread(geno_file)
    if (is.null(data) || nrow(data) == 0) {stop(paste(genotype, " file is invalid"))}
    message(nrow(data), " rows loaded for ", sample_name, " (", genotype, ").")
    
    # check and rename columns
    if (ncol(data) >= 10) {
      exp_colnames <- c("CHROM", "POS", "REF", "ALT", "QUAL", "DP", "Fref", 
                        "Rref", "Falt", "Ralt")
      data.table::setnames(data, old = colnames(data)[1:10], new = exp_colnames)
    } else {
      stop("Input file does not contain the expected number of columns (10).")
    }
    
    # Ensure numeric columns are numeric
    num_cols <- c("POS", "QUAL", "DP", "Fref", "Rref", "Falt", "Ralt")
    data[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
    
    # Filter and Keep only SNPs (exclude indels)
    data <- data[nchar(REF) == 1 & nchar(ALT) == 1]
    
    # Filter by depth and quality
    data <- data[!is.na(DP) & !is.na(QUAL) & DP >= min_DP & QUAL >= min_QUAL]
    
    # Allele frequency (AF)
    data[, `:=`(Tref = Fref + Rref, Talt = Falt + Ralt)]
    data[, AF := data.table::fifelse((Tref + Talt) > 0, 
                                     Talt / (Tref + Talt), NA_real_)]
    data[, nonhom := as.integer(Tref > 0 & Talt > 0)]
    
    # Sort and prefix columns
    data[, POS := as.integer(POS)]
    data.table::setorder(data, CHROM, POS)
    cols_to_prefix <- setdiff(names(data), c("CHROM", "POS"))
    data.table::setnames(data, old = cols_to_prefix,
                         new = paste0(genotype, "_", cols_to_prefix))
    geno_data[[genotype]] <- data
    message(nrow(data), " rows kept for ", sample_name, ".\n")
  }
  
  # Ensure both genotypes present if not mutant-only mode
  if (!only_mutant && (is.null(geno_data$wt) || is.null(geno_data$mt))) {
    stop("Both wild-type and mutant data are required.")
  }
  
  message("Created list for:", paste(names(geno_data), collapse = ", "))
  return(geno_data)
}
#' 
#' Compute BSA-Seq statistics from filtered SNP tables
#' This function takes filtered wild-type and mutant SNP tables generated
#' by `read_vcf()`, identifies shared and unique SNPs, extracts EMS-type
#' variants, and computes BSA-Seq mapping statistics including allele
#' frequency difference (AFD), Euclidean distance (ED), ED4, and G statistic.
#' @param geno_data A list of filtered SNP tables returned by `read_vcf()`
#' @param prefix Prefix used for naming output files
#' @param save_results Logical; if TRUE, save output tables and RDS object
#' @param output_dir Directory for saving output files
#' @param only_mutant Logical; if TRUE, run mutant-only mode
#' @return A list of BSA-Seq result tables
#' @export
#' 
compute_bsa <- function(geno_data, prefix, save_results = FALSE, 
                        output_dir = "post_analysis", only_mutant = FALSE) {
  # Function to sort variants by chromosome and Position
  sort_variants <- function(dt) {
    dt <- data.table::as.data.table(dt)
    dt[, CHROM := factor(CHROM, levels = naturalsort(unique(CHROM)))]
    dt <- dt[order(CHROM, POS)]
    dt[, CHROM := as.character(CHROM)]
    dt
  }
  
  # Function to identify ems specific SNPs (G > A or C > T)
  is_ems <- function(ref, alt) {
    (ref == "G" & alt == "A") | (ref == "C" & alt == "T")
  }
  # filter and Extract EMS SNPs from mutant
  get_ems <- function(data, ref_col, alt_col) {
    sort_variants(data[ is_ems(get(ref_col), get(alt_col)) ])
  }
  
  # anti-join and identity unique SNPs
  id_unique_snps <- function(data1, data2) {
    sort_variants(data.table::as.data.table(
      dplyr::anti_join(data1, data2, by = c("CHROM", "POS"))))
  }
  
  # Calculate G-statistics for each SNP
  compute_G <- function(wt_Tref, wt_Talt, mt_Tref, mt_Talt) {
    obs <- c(mt_Tref, mt_Talt, wt_Tref, wt_Talt)
    total <- sum(obs)
    
    if (total == 0) return(NA_real_)
    
    # Compute marginal counts for 2x2 contingency tables
    row_ref <- mt_Tref + wt_Tref
    row_alt <- mt_Talt + wt_Talt
    col_mt  <- mt_Tref + mt_Talt
    col_wt  <- wt_Tref + wt_Talt
    
    # Expected counts under independence
    exp <- c(col_mt * row_ref / total, col_mt * row_alt / total,
             col_wt * row_ref / total, col_wt * row_alt / total)
    
    # Avoid division by zero or log(0)
    obs[obs == 0] <- 1e-10
    exp[exp == 0] <- 1e-10
    
    # compute G-statistic (livelihood ratio)
    return(2 * sum(obs * log(obs / exp), na.rm = TRUE))
  }
  
  # helper to save outputs
  save_output <- function(result, prefix, output_dir, rds_name) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(result, file.path(output_dir, rds_name))
    for (name in names(result)) {
      data.table::fwrite(
        result[[name]], file.path(output_dir, paste0(prefix, "_", name,".csv")))
    }
  }
  
  # Process mutant-only mode
  if (only_mutant) {
    if (is.null(geno_data$mt)) stop("Mutant data is required.")
    
    mt_data <- sort_variants(geno_data$mt)
    mt_ems <- get_ems(mt_data, "mt_REF", "mt_ALT")
    result <- list(mt_all = mt_data, mt_ems = mt_ems)
    
    # Save output
    if (save_results) {
      save_output(result, prefix, output_dir, 
                  paste0(prefix, "_mutant_only_results.rds"))
      message("Mutant-only results saved to ", output_dir)
    }
    
    return(result)
  } else {
    
    # WT vs MT analysis 
    if (is.null(geno_data$wt) || is.null(geno_data$mt)) {
      stop("Both wild-type and mutant datasets are required.")}
    
    wt_data <- sort_variants(geno_data$wt)
    mt_data <- sort_variants(geno_data$mt)
    
    # Merge wild-type and mutant by CHROM and POS
    message("Merging wild-type and mutant datasets...")
    wt_mt <- data.table::merge.data.table(wt_data, mt_data, by = c("CHROM", "POS"))
    wt_mt <- sort_variants(wt_mt)
    message("Number of shared SNPs (wt_mt): ", nrow(wt_mt))
    
    # Calculate Allele frequency difference (AFD)
    wt_mt[, AFD := abs(wt_AF - mt_AF)]
    
    # Calculate Euclidean Distance (ED) and its fourth power (ED4)
    wt_mt[, ED := sqrt(2) * AFD]
    wt_mt[, ED4 := ED^4]
    
    # Calculate G-statistics (G-stat)
    wt_mt[, G := mapply(compute_G, wt_Tref, wt_Talt, mt_Tref, mt_Talt)]
    
    # Unique SNPs and EMS filtering
    ant_wt <- id_unique_snps(wt_data, mt_data)
    ant_mt <- id_unique_snps(mt_data, wt_data)
    ant_wt_ems <- get_ems(ant_wt, "wt_REF", "wt_ALT")
    ant_mt_ems <- get_ems(ant_mt, "mt_REF", "mt_ALT")
    wt_mt_ems <- wt_mt[ is_ems(wt_REF, wt_ALT) & is_ems(mt_REF, mt_ALT) ]
    wt_mt_ems <- sort_variants(wt_mt_ems)
    
    result <- list(wt_mt_all = wt_mt, ant_mt_all = ant_mt, ant_wt_all = ant_wt, 
                   wt_mt_ems = wt_mt_ems, ant_mt_ems = ant_mt_ems, 
                   ant_wt_ems = ant_wt_ems)
    
    # Save output
    if (save_results) {
      save_output(result, prefix, output_dir, paste0(prefix, "_results.rds"))
      message("All results saved to ", output_dir)
    }
    return(result)
  }
}



# vcf_dir="/Users/zebosi/Documents/osu_postdoc/BSA/S004/data/snps"
# wt <- c("S004A")
# mt <- c("S004B")
# prefix = "b73"
# pattern = "snps\\.tsv$"
# min_DP = 5
# min_QUAL = 20
# 
# geno_data <- read_vcf(
#   vcf_dir, prefix, pattern,
#   Genotypes = list(wt = wt, mt = mt),
#   min_DP, min_QUAL, only_mutant = FALSE
# )


# result <- compute_bsa(
#   geno_data = geno_data, prefix = prefix, save_results = FALSE,
#   output_dir = output_dir, only_mutant = FALSE
# )