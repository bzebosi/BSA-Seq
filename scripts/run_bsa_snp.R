src_base <- "https://raw.githubusercontent.com/bzebosi/BSA-Seq/main/R/"

for (f in c("install_pkgs.R", "utils.R", 
            "plot_bsa.R", "core_bsa.R", "bsa_snp.R")) {
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