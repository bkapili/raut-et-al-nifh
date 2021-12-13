# Portions of code adapted from decontam tutorial
# (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)

### Load required packages
# List required packages
cranPackages <- c("BiocManager", "dplyr", "tidyr", "ape", "ggplot2")
biocPackages <- c("phyloseq", "Biostrings", "decontam")

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

# Load packages
suppressPackageStartupMessages(lapply(c(cranPackages, biocPackages),
                                      library,
                                      character.only = TRUE))

### Prepare phyloseq object for decontam
# Load phyloseq object
psTaxa <- readRDS(file = "../robjects/psTaxa.rds")

# Add sample_data variable of sample type
blankMock <- c("Blank", "Mock")
sampleInd <- sample_data(psTaxa)$sample_name %>%
  grep(paste(blankMock, collapse = "|"), ., invert = TRUE)

blankInd <- sample_data(psTaxa)$sample_name %>%
  grep("Blank", .)

mockInd <- sample_data(psTaxa)$sample_name %>%
  grep("Mock", .)

sample_data(psTaxa)$sample_type[sampleInd] <- "Sample"
sample_data(psTaxa)$sample_type[blankInd] <- "Negative"
sample_data(psTaxa)$sample_type[mockInd] <- "Mock"

# Add sample_data variable of DNA quantitation data
quant <- read.csv(file = "../data/decontam_dna_quantitation.csv")
sample_data(psTaxa)$picogreen_ngul <- as.numeric(quant$picogreen_ngul)


### Plot library sizes
dfLibSize <- psTaxa %>%
  sample_data %>%
  data.frame %>%
  mutate(library_size = sample_sums(psTaxa)) %>%
  arrange(library_size)

pLibSizes <- ggplot(data = dfLibSize, aes(x = seq(nrow(dfLibSize)),
                              y = library_size,
                              color = sample_type)) +
  geom_point() +
  xlab("Rank") +
  ylab("Library size (reads)") +
  theme_bw()

ggsave(filename = "../supplemental/library_sizes.pdf", plot = pLibSizes,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)


### Run decontam with frequency method
# Remove negatives, mocks, and sample with no PicoGreen data
psDecontam <- subset_samples(psTaxa, sample_type != "Negative" &
                                   sample_type !="Mock" & is.na(picogreen_ngul) == FALSE)
dfDecontam <- isContaminant(psDecontam, method = "frequency",
                                conc = "picogreen_ngul", threshold = 0.1)

# Add total number of reads for each taxon
dfDecontam <- dfDecontam %>%
  mutate(reads = taxa_sums(psDecontam))

# Plot decontam score histograms
pHistASV <- ggplot(data = dfDecontam, aes(x = p)) + 
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  xlab("decontam score") +
  ylab("Number ASVs") +
  theme_bw()

pHistReads <- ggplot(data = dfDecontam, aes(x = p, weight = reads)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(trans = "log10",
                     minor_breaks = rep(1:9, 3)*(10^rep(0:5, each=9))) +
  xlab("decontam score") +
  ylab("Number reads") +
  theme_bw()

ggsave(filename = "../supplemental/decontam_score_asvs.pdf", plot = pHistASV,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)

ggsave(filename = "../supplemental/decontam_score_reads.pdf", plot = pHistReads,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)

# Add taxonomy to decontam data frame
dfDecontam <- dfDecontam %>%
  cbind(., psDecontam %>% tax_table %>% as.data.frame) %>%
  arrange(-contaminant, Domain, Phylum, Class, Order, Family, Genus)

# Remove putative contaminants
notContam <- dfDecontam %>%
  filter(contaminant == FALSE) %>%
  rownames

psNoContam <- prune_taxa(notContam, psTaxa)

# Write phyloseq object
saveRDS(psNoContam, file = "../robjects/psNoContam.rds")

### Add filtering summary to read and ASV retention logs
# Load files
readTrack <- read.csv(file = "../supplemental/read_retention.csv",
                      header = TRUE, row.names = 1)
asvTrack <- read.csv(file = "../supplemental/asv_retention.csv",
                     header = TRUE, row.names = 1)

# Add per-sample reads and total ASVs
readTrack <- mutate(readTrack, decontam = sample_sums(psNoContam))
asvTrack <- mutate(asvTrack, decontam = ncol(otu_table(psNoContam)))

# Write to csv
write.csv(readTrack, file = "../supplemental/read_retention.csv", quote = FALSE,
          row.names = TRUE, col.names = TRUE)
write.csv(asvTrack, file = "../supplemental/asv_retention.csv", quote = FALSE,
          row.names = TRUE, col.names = TRUE)