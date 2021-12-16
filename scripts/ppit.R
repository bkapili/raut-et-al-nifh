# -------------------------------------------------------------
# Script purpose: Infer taxonomic identities for nifH ASVs using
#                 PPIT (v.1.2.0).
#
# Inputs:
#   * Phyloseq object before decontamination ("psRaw.rds")
#   * SEPP output tree ("sepp_placement.nwk")
#   * SEPP output alignment ("sepp_alignment.fasta")
#
# Output: Phyloseq object with taxonomic inferences (overwritten "psRaw.rds")
# -------------------------------------------------------------

### Load required packages
# List required packages
cranPackages <- c("devtools", "BiocManager", "dplyr", "tidyr", "ape")
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

# Load packages
lapply(c(cranPackages, biocPackages, "ppit"), library, character.only = TRUE)

### Run PPIT
# Load phyloseq object
psRaw <- readRDS(file = "../robjects/psRaw.rds")

# Set PPIT input
type <- "partial"
query <- names(refseq(psRaw))
tree <- read.tree(file = "../robjects/sepp_placement.nwk")
alignment <- readDNAMultipleAlignment("../robjects/sepp_alignment.fasta")
taxonomy <- ppit::nifH_reference_taxonomy_v2
cutoffs <- ppit::nifH_cutoffs_v2

# Optimize phylogenetic neighborhood
optThresh <- tree.partition(type, query, tree, alignment, taxonomy, cutoffs)
print(optThresh)
#optThresh = 0.42 for MIN_HAMMING=1
#optThresh <- cutoffs[1,5]

# Infer taxonomy
dfInferences <- ppit(type, query, tree, alignment,
                     taxonomy, optThresh, cutoffs)

# Requery those with pat distance fail
requery <- dfInferences %>%
  filter(Suspected_homolog == "X") %>%
  rownames

dfRequery <- ppit(type, requery, tree, alignment,
                  taxonomy, cutoffs[1,1], cutoffs)
dfInferences[match(rownames(dfRequery), rownames(dfInferences)),] <- dfRequery

# Rename homolog and unidentified ASV taxonomy
allHomologs <- dfInferences %>%
  filter(Suspected_homolog == "X" | Percent_id_fail == "X" | Pat_dist_fail == "X") %>%
  rownames

unidentified <- dfInferences %>%
  filter(Potential_HGT == "X") %>%
  rownames

tax_table <- dfInferences %>%
  select(Domain, Phylum, Class, Order, Family, Genus) %>%
  as.matrix

tax_table[match(allHomologs, rownames(tax_table)),] <- "suspected_homolog"
tax_table[match(unidentified, rownames(tax_table)),] <- "unidentified"

# Add to phyloseq tax_table slot
tax_table(psRaw) <- tax_table

# Save phyloseq object
saveRDS(psRaw, "../robjects/psRaw.rds")