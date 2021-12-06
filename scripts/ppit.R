### Load required packages

# List required packages
cranPackages <- c("devtools", "BiocManager", "ggplot2", "dplyr", "tidyr", "ape")
biocPackages <- c("phyloseq", "Biostrings")

# Install missing CRAN packages
installedCRANPackages <- cranPackages %in% rownames(installed.packages())
if (any(installedCRANPackages == FALSE)) {
  install.packages(cranPackages[!installedCRANPackages],
                   repos='http://cran.us.r-project.org')
}

# Install missing Bioconductor packages
installedBioPackages <- biocPackages %in% rownames(installed.packages())
if (any(installedBioPackages == FALSE)) {
  BiocManager::install(biocPackages[!installedBioPackages])
}

# Install ppit if missing
ppitInstall <- "ppit" %in% rownames(installed.packages())
if (ppitInstall == FALSE) {
  devtools::install_github("BKapili/ppit", build_vignettes = TRUE)
}

