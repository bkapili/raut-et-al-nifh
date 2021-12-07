### Load required packages
# List required packages
cranPackages <- c("BiocManager", "ape", "devtools")
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
  devtools::install_github("BKapili/ppit", build_vignettes = FALSE)
}

# Load packages
lapply(c(cranPackages, biocPackages, "ppit"), library, character.only = TRUE)


### Write out SEPP input files
# Write out chimera-filtered ASV sequences
psRaw <- readRDS(file = "../robjects/psRaw.rds")
writeXStringSet(refseq(psRaw), filepath = "../robjects/psRaw_ASV_seqs.fasta", format = "fasta")

# Write out nifH reference alignment
writeXStringSet(DNAStringSet(ppit::nifH_reference_alignment_v2),
                filepath ="../robjects/nifH_reference_alignment_v2.fasta",
                format = "fasta")

# Write out nifH reference tree
write.tree(ppit::nifH_reference_tree_v2, file = "../robjects/nifH_reference_tree_v2.nwk")

# Write out RAxML info file
write.table(ppit::nifH_reference_RAxML_info_v2, file = "../robjects/nifH_reference_RAxML_info_v2.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)