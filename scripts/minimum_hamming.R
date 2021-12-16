# -------------------------------------------------------------
# Script purpose: Plot mock community compositions with a specified
#                 MIN_HAMMING value during dada() ASV inferencing.
#
# Run as: Rscript minimum_hamming.R MIN_HAMMING
# Example to run DADA2 w/ minimum hamming distance=1: Rscript minimum_hamming.R 1
#
# Inputs:
#   * Primer-removed paired-end .fastq files
#   * FIGARO output of ranked trimming positions ("trimParameters.json")
#   * Sample metadata ("SampleID_ExperimentID.csv")
#   * Tab-delimited summary of cutadapt statistics ("cutadapt_log.txt")
#   * PicoGreen dsDNA quantitation data ("decontam_dna_quantitation.csv")
#   * BLAST database of expected mock nifH sequences ("db_nifh_mock")
#
# Outputs:
#   * Trimmed and quality-filtered reads in .fastq.gz format ("filtered" subdirectory)
#   * Plots of DADA2 error model ("learned_errors_F.pdf" and "learned_errors_R.pdf")
#   * Seqtabs of raw, length_filtered, and length_filtered+chimera_removed merged reads
#   * Phyloseq object before decontamination ("psRaw_min_ham_X.rds"; X = MIN_HAMMING)
#   * Decontaminated phyloseq object ("psNoContam_X.rds"; X = MIN_HAMMING)
#   * Plot of per-sample library sizes before decontamination ("library_sizes_X.pdf"; X = MIN_HAMMING)
#   * Histogram of decontam scores with num. of ASVs per bin ("decontam_score_asvs_X.pdf"; X = MIN_HAMMING)
#   * Histogram of decontam scores with num. of reads per bin ("decontam_score_reads_X.pdf"; X = MIN_HAMMING)
#   * Updated log of per-sample read tracking after contaminant removal ("read_retention_X.csv"; X = MIN_HAMMING)
#   * Log of total ASV tracking ("asv_retention_X.csv"; X = MIN_HAMMING)
#   * Plot of mock community compositions ("mock_plot_min_ham_X.pdf"; X = MIN_HAMMING)
#
# Notes: Portions of code adapted from DADA2 tutorial (v.1.16)
# (https://benjjneb.github.io/dada2/tutorial.html)
# -------------------------------------------------------------

### Load required packages
# List required packages
cranPackages <- c("BiocManager", "ggplot2", "dplyr", "tidyr",
                  "jsonlite", "ape", "ggplot2", "reshape2")
biocPackages <- c("phyloseq", "dada2", "Biostrings", "decontam")

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

# Install rblast if missing
rblastInstall <- "rBLAST" %in% rownames(installed.packages())
if (rblastInstall == FALSE) {
  devtools::install_github("mhahsler/rBLAST")
}

# Load packages
lapply(c(cranPackages, biocPackages, "rBLAST"), library, character.only = TRUE)


### Trim and quality-filter reads
# Set minimum hamming
minHam <- commandArgs(trailingOnly = TRUE) %>% as.integer

# Extract file and sample names
fnFs <- list.files("../data/cutadapt", pattern = "_L001_R1_001.fastq", full.names = TRUE) %>% sort
fnRs <- list.files("../data/cutadapt", pattern = "_L001_R2_001.fastq", full.names = TRUE) %>% sort
sampleNames <- sapply(strsplit(basename(fnFs), "-CUTADAPT"), `[`, 1)

# Create directory for filtered reads
filtPath <- paste0("../data/filtered/min_ham_", minHam)
dir.create(filtPath, recursive = TRUE)

# Set path for future filtered files
filtFs <- file.path(filtPath, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filtPath, paste0(sampleNames, "_R_filt.fastq.gz"))
names(filtFs) <- sampleNames
names(filtRs) <- sampleNames

# Read in optimal trim positions
trimParams <- fromJSON("../data/figaro/trimParameters.json", flatten = TRUE) %>%
  filter(row_number() == 1) %>%
  pull(trimPosition) %>%
  unlist

# Trim and quality-filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(trimParams[1], trimParams[2]),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = FALSE)


### Model error rates
# Dereplicate reads
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Learn error rates
errF <- learnErrors(derepFs, multithread = TRUE)
errR <- learnErrors(derepRs, multithread = TRUE)

# Export PDF of learned errors
dir.create("../supplemental")
pErrF <- plotErrors(errF, nominalQ = TRUE)
pErrR <- plotErrors(errR, nominalQ = TRUE)

ggsave(filename = "../supplemental/learned_errors_F.pdf", plot = pErrF,
       device = "pdf", units = "cm", width = 35, height = 20,
       dpi = 300, useDingbats=FALSE)

ggsave(filename = "../supplemental/learned_errors_R.pdf", plot = pErrR,
       device = "pdf", units = "cm", width = 35, height = 20,
       dpi = 300, useDingbats=FALSE)


### Infer ASVs
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE, MIN_HAMMING = minHam)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = TRUE, MIN_HAMMING = minHam)


### Create and filter ASV table
# Merge read pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Make ASV table
rawSeqtab <- makeSequenceTable(mergers)

# Remove sequences too short or too long
table(nchar(getSequences(rawSeqtab)))
lenFiltSeqtab <- rawSeqtab[,nchar(colnames(rawSeqtab)) %in% seq(330, 370)]

# Remove bimeras
chimFiltSeqtab <- removeBimeraDenovo(lenFiltSeqtab, method = "consensus",
                                     multithread = TRUE, verbose = TRUE)
dim(chimFiltSeqtab)
sum(chimFiltSeqtab)/sum(lenFiltSeqtab)

# Save seqtabs
dir.create("../robjects")
saveRDS(rawSeqtab, paste0("../robjects/rawSeqtab_", minHam, ".rds"))
saveRDS(lenFiltSeqtab, paste0("../robjects/lenFiltSeqtab_", minHam, ".rds"))
saveRDS(chimFiltSeqtab, paste0("../robjects/chimFiltSeqtab_", minHam, ".rds"))


### Track read and ASV retention
# Load cutadapt log
cutLog <- read.table(file = "../data/cutadapt/log/cutadapt_log.txt", header = TRUE) %>%
  arrange(sample)

# Define function to return the number of unique reads/sequences
getN <- function(x) sum(getUniques(x))

# Apply function to output from each step and store in matrix
readTrack <- cbind(pull(cutLog, in_reads), out, sapply(dadaFs, getN), sapply(dadaRs, getN),
                   sapply(mergers, getN), rowSums(lenFiltSeqtab), rowSums(chimFiltSeqtab))

# Rename column names to step names and row names to sample names
colnames(readTrack) <- c("raw", "cutadapt", "filtered", "denoisedF", "denoisedR",
                         "merged", "length_filtered", "chim_filtered")
rownames(readTrack) <- sampleNames

# Write read retention matrix to csv
write.csv(readTrack, file = paste0("../supplemental/read_retention_", minHam, ".csv"),
          quote = FALSE, row.names = TRUE, col.names = TRUE)

# Record number of ASVs in matrix
asvTrack <- cbind(ncol(rawSeqtab), ncol(lenFiltSeqtab), ncol(chimFiltSeqtab))

# Rename column names to step names
colnames(asvTrack) <- c("merged", "length_filtered", "chim_filtered")

# Write ASV retention matrix to csv
write.csv(asvTrack, file = paste0("../supplemental/asv_retention_", minHam, ".csv"),
          quote = FALSE, row.names = TRUE, col.names = TRUE)


### Create phyloseq object
# Import sample metadata
sampMetadata <- read.csv(file = "../data/SampleID_ExperimentID.csv", header = TRUE, row.names = 1)

# Format ASV table for otu_table class
asvSeqtab <- chimFiltSeqtab
colnames(asvSeqtab) <- paste0("ASV", 1:ncol(asvSeqtab))

# Format sequences for refseq class
seqs <- colnames(chimFiltSeqtab) %>% DNAStringSet()
names(seqs) <- paste0("ASV", 1:length(seqs))

# Format sample metadata for sample_data class
sampMetadata <- sampMetadata %>%
  separate(ExperimentID, sep = "-", c("project", "sample_name", "seaweed_species",
                                      "day_of_decay", "replicate", "nucleic_acid", "gene"))

# Build and save phyloseq object
psRaw <- phyloseq(otu_table(asvSeqtab, taxa_are_rows = FALSE),  
                  refseq(seqs),
                  sample_data(sampMetadata))

saveRDS(psRaw, paste0("../robjects/psRaw_min_ham_", minHam, ".rds"))


### decontaminate sequences
# Add sample_data variable of sample type
blankMock <- c("Blank", "Mock")
sampleInd <- sample_data(psRaw)$sample_name %>%
  grep(paste(blankMock, collapse = "|"), ., invert = TRUE)

blankInd <- sample_data(psRaw)$sample_name %>% grep("Blank", .)
mockInd <- sample_data(psRaw)$sample_name %>% grep("Mock", .)

sample_data(psRaw)$sample_type[sampleInd] <- "Sample"
sample_data(psRaw)$sample_type[blankInd] <- "Negative"
sample_data(psRaw)$sample_type[mockInd] <- "Mock"

# Add sample_data variable of DNA quantitation data
quant <- read.csv(file = "../data/decontam_dna_quantitation.csv")
sample_data(psRaw)$picogreen_ngul <- as.numeric(quant$picogreen_ngul)


# Plot library sizes
dfLibSize <- psRaw %>%
  sample_data %>%
  data.frame %>%
  mutate(library_size = sample_sums(psRaw)) %>%
  arrange(library_size)

pLibSizes <- ggplot(data = dfLibSize, aes(x = seq(nrow(dfLibSize)),
                                          y = library_size,
                                          color = sample_type)) +
  geom_point() +
  xlab("Rank") +
  ylab("Library size (reads)") +
  theme_bw()

ggsave(filename = paste0("../supplemental/library_sizes_", minHam, ".pdf"), plot = pLibSizes,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)

# Remove negatives, mocks, and sample with no PicoGreen data
psDecontam <- subset_samples(psRaw, sample_type != "Negative" &
                               sample_type !="Mock" & is.na(picogreen_ngul) == FALSE)

dfDecontam <- isContaminant(psDecontam, method = "frequency",
                            conc = "picogreen_ngul", threshold = 0.1)

# Add total number of reads for each taxon
dfDecontam <- dfDecontam %>% mutate(reads = taxa_sums(psDecontam))

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

# Save histograms
ggsave(filename = paste0("../supplemental/decontam_score_asvs_", minHam, ".pdf"), plot = pHistASV,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)

ggsave(filename = paste0("../supplemental/decontam_score_reads_", minHam, ".pdf"), plot = pHistReads,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)

# Remove putative contaminants
notContam <- dfDecontam %>%
  filter(contaminant == FALSE) %>%
  rownames

psNoContam <- prune_taxa(notContam, psRaw)

# Write phyloseq object
saveRDS(psNoContam, file = paste0("../robjects/psNoContam_", minHam, ".rds"))

### Add filtering summary to read and ASV retention logs
# Add per-sample reads and total ASVs
readTrack <- readTrack %>% as.data.frame %>% mutate(decontam = sample_sums(psNoContam))
asvTrack <- asvTrack %>% as.data.frame %>% mutate(decontam = ncol(otu_table(psNoContam)))

# Write to csv
write.csv(readTrack, file = paste0("../supplemental/read_retention_", minHam, ".csv"), quote = FALSE,
          row.names = TRUE, col.names = TRUE)
write.csv(asvTrack, file = paste0("../supplemental/asv_retention_", minHam, ".csv"), quote = FALSE,
          row.names = TRUE, col.names = TRUE)

### Plot mock community
# Remove singletons and subset to mocks
psMock <- psNoContam %>%
  prune_taxa(taxa_sums(.) > 1, .) %>%
  subset_samples(., sample_type == "Mock") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(., function(x) 100*x/sum(x))

# BLAST ASVs against known mock sequences
dbMock <- blast(db = "../data/rblast/db_nifh_mock")
dfBlast <- predict(dbMock, refseq(psMock), type = "blastn",
                   custom_format = "qseqid sseqid pident mismatch length qcovs")

# Identify exact, 1-mismatch, and 2-mismatch ASVs
exactMatch <- dfBlast %>% filter(pident == 100 & qcovs == 100) %>% pull(qseqid)
oneMismatch <- dfBlast %>% filter(mismatch == 1 & qcovs == 100) %>% pull(qseqid)
twoMismatch <- dfBlast %>% filter(mismatch == 2 & qcovs == 100) %>% pull(qseqid)
print(paste(length(exactMatch), "of 12 expected mock sequences recovered."))

# Calculate % of reads the exact matches comprise in each sample
psNoContam %>%
  subset_samples(., sample_type == "Mock") %>%
  transform_sample_counts(., function(x) 100*x/sum(x)) %>%
  prune_taxa(exactMatch, .) %>%
  sample_sums

### Create summary figure
# Prepare data frame for scatter plot
dfMock <- melt(otu_table(psMock)) %>%
  group_by(Var2) %>%
  mutate(prevalence = sum(value != 0),
         median = median(value),
         Prevalence_bin = case_when(prevalence == 6 ~ "6",
                                    prevalence %in% 2:5 ~ "2 - 5",
                                    prevalence == 1 ~ "1"),
         type = case_when(Var2 %in% exactMatch ~ "exact_match",
                          Var2 %in% oneMismatch ~ "one_mismatch",
                          Var2 %in% twoMismatch ~ "two_mismatch",
                          TRUE ~ "other")) %>%
  arrange(-median)

# Prepare mismatch data frame
numTaxa <- ncol(otu_table(psMock))
dfRect <- dfMock %>%
  ungroup %>%
  select(Var2, type) %>%
  unique %>%
  mutate(xmin = seq(0.5, numTaxa-0.5, 1), xmax = seq(1.5, numTaxa+0.5, 1))

# Export PDF of scatter plot overlaid with mismatch data
p <- ggplot() +
  geom_rect(data = dfRect, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 30,
                               fill = type), alpha = 0.3) +
  geom_point(data = dfMock %>% filter(value > 0),
             aes(x = reorder(Var2, -value), y = value, color = Prevalence_bin)) +
  scale_x_discrete() +
  scale_y_continuous(trans = "log10",
                     minor_breaks = rep(1:9, 5)*(10^rep(-3:1, each=9))) +
  theme_bw() +
  xlab("ASV") +
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_blank())

ggsave(filename = paste0("../supplemental/mock_plot_min_ham_", minHam, ".pdf"), plot = p,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)