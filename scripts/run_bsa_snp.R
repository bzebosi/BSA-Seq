base <- "https://raw.githubusercontent.com/bzebosi/BSA-Seq/main/"
src_base <- paste0(base, "R/")

# Load pipeline functions
for (f in c("install_pkgs.R", "utils.R", 
            "plot_bsa.R", "core_bsa.R", "bsa_snp.R")) {
  source(paste0(src_base, f))
}

# Required packages
packages <- c(
  "reshape2", "readxl", "BiocManager", "zoo", "plyr", "GlobalOptions", 
  "shape", "scales", "tidyverse", "openxlsx", "stringr", "IRanges", "magrittr", 
  "data.table", "naturalsort", "locfit", "rlang"
)

# Install/load required packages
Install_pkgs(packages)


# -------------------------------------------------
# Load config
# -------------------------------------------------

# Option 1: use already loaded cfg/common_args
if (exists("cfg") && is.list(cfg) &&
    exists("common_args") && is.list(common_args)) {
  
  message("Using existing config in environment")
  
} else {
  
  # Option 2: choose config interactively
  message("Select config file")
  
  cfg_file <- file.choose()
  
  source(cfg_file)
}


message("BSA-Seq SNP environment ready")
message("Running analysis...")


# -------------------------------------------------
# Run BSA-SNP analysis
# -------------------------------------------------

results <- list()

for (s in names(cfg)) {
  
  message("\n=== Running ", s, " ===")
  
  args <- c(
    list(
      vcf_dir = cfg[[s]]$vcf_dir,
      output_dir = cfg[[s]]$output_dir,
      wt_list = cfg[[s]]$wt_list,
      mt_list = cfg[[s]]$mt_list,
      prefix_list = cfg[[s]]$prefix_list
    ),
    common_args
  )
  
  results[[s]] <- do.call(run_bsa_snp, args)
}


# Save results
saveRDS(results, "results_bsa_snp.rds")

message("\nAnalysis complete")