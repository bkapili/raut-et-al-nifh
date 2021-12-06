### Load required packages

# List required packages
cranPackages <- c("BiocManager", "ggplot2", "dplyr", "tidyr", "jsonlite")
biocPackages <- c("phyloseq", "dada2", "Biostrings")

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
lapply(c(cranPackages, biocPackages), library, character.only = TRUE)


### Trim and quality-filter reads
# Extract file and sample names
fnFs <- list.files("../data/cutadapt", pattern = "_L001_R1_001.fastq", full.names = TRUE) %>% sort
fnRs <- list.files("../data/cutadapt", pattern = "_L001_R2_001.fastq", full.names = TRUE) %>% sort
sampleNames <- sapply(strsplit(basename(fnFs), "-CUTADAPT"), `[`, 1)

# Create directory for filtered reads
dir.create("../data/filtered")

# Set path for future filtered files
filtFs <- file.path("../data/filtered", paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path("../data/filtered", paste0(sampleNames, "_R_filt.fastq.gz"))
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
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = TRUE)


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
saveRDS(rawSeqtab, paste0("../robjects", "/rawSeqtab.rds"))
saveRDS(lenFiltSeqtab, paste0("../robjects", "/lenFiltSeqtab.rds"))
saveRDS(chimFiltSeqtab, paste0("../robjects", "/chimFiltSeqtab.rds"))


### Track read retention
# Load cutadapt log
cutLog <- read.table(file = "../data/cutadapt/log/cutadapt_log.txt", header = TRUE) %>%
  arrange(sample)

# Define function to return the number of unique reads/sequences
getN <- function(x) sum(getUniques(x))

# Apply function to output from each step and store in matrix
track <- cbind(pull(cutLog, in_reads), out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(lenFiltSeqtab), rowSums(chimFiltSeqtab))

# Rename column names to step names and row names to sample names
colnames(track) <- c("raw", "cutadapt", "filtered", "denoisedF", "denoisedR",
                     "merged", "length_filtered", "chim_filtered")
rownames(track) <- sampleNames

# Write read retention data frame to csv
write.csv(track, file = "../supplemental/read_retention.csv", quote = FALSE,
          row.names = TRUE, col.names = TRUE)


### Create phyloseq object
# Import sample metadata
sampMetadata <- read.csv(file = "../data/SampleID_ExperimentID.csv",
                         header = TRUE,
                         row.names = 1)

# Format ASV table for otu_table class
asvSeqtab <- chimFiltSeqtab
colnames(asvSeqtab) <- paste("ASV", 1:ncol(asvSeqtab), sep = ".")

# Format sequences for refseq class
seqs <- colnames(chimFiltSeqtab) %>%
  DNAStringSet()
names(seqs) <- paste("ASV", 1:length(seqs), sep = ".")

# Format sample metadata for sample_data class
sampMetadata <- sampMetadata %>%
  separate(ExperimentID, sep = "-", c("project", "site_depth",  "deployment", "core",
                           "sed_depth", "nucleic_acid", "replicate", "gene"))

# Build and save phyloseq object
psRaw <- phyloseq(otu_table(asvSeqtab, taxa_are_rows = FALSE),  
               refseq(seqs),
               sample_data(sampMetadata))
saveRDS(psRaw, paste0("../robjects", "/psRaw.rds"))