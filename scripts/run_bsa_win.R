src_base <- "https://raw.githubusercontent.com/bzebosi/BSA-Seq/main/R/"

for (f in c("install_pkgs.R", "utils.R", 
            "plot_bsa.R", "core_bsa.R", "bsa_win.R")) {
  source(paste0(src_base, f), local = TRUE)
}

message("BSA-Seq SNP environment ready. Call run_bsa_snp() when you’re set.")

# required packages
packages <- c(
  "reshape2", "readxl", "BiocManager", "zoo", "plyr", "GlobalOptions", 
  "shape", "scales", "tidyverse", "openxlsx", "stringr", "IRanges", "magrittr", 
  "data.table", "naturalsort", "locfit", "rlang"
)

# Install required packages
Install_multi_package_bz(packages)

message("Select a config file")

cfg_file <- file.choose()

source(cfg_file)


message("BSA-Seq SNP environment ready")
message("Running analysis...")

# --------------------------------
# Run BSA-SNP analysis
# --------------------------------
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
  
  results[[s]] <- do.call(run_bsa_win, args)
}