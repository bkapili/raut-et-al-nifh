### Load required packages
# List required packages
cranPackages <- c("BiocManager", "dplyr", "tidyr", "reshape2", "ggplot2")
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

# Install rblast if missing
rblastInstall <- "rBLAST" %in% rownames(installed.packages())
if (rblastInstall == FALSE) {
  devtools::install_github("mhahsler/rBLAST")
}

# Load packages
suppressPackageStartupMessages(lapply(c(cranPackages, biocPackages, "rBLAST"),
                                      library,
                                      character.only = TRUE))

### Evaluate mock
# Load in decontaminated phyloseq
psNoContam <- readRDS(file = "../robjects/psNoContam.rds")

# Remove singletons and subset to mocks
psMock <- psNoContam %>%
  prune_taxa(taxa_sums(.) > 1, .) %>%
  subset_samples(., sample_type == "Mock") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(., function(x) 100*x/sum(x))

# BLAST ASVs against known mock sequences
dbMock <- blast(db = "../data/db_nifh_mock")
dfBlast <- predict(dbMock, refseq(psMock), type = "blastn",
                   custom_format = "qseqid sseqid pident mismatch length qcovs")

# Identify exact, 1-mismatch, and 2-mismatch ASVs
exactMatch <- dfBlast %>% filter(pident == 100 & qcovs == 100) %>% pull(qseqid)
oneMismatch <- dfBlast %>% filter(mismatch == 1 & qcovs == 100) %>% pull(qseqid)
twoMismatch <- dfBlast %>% filter(mismatch == 2 & qcovs == 100) %>% pull(qseqid)
print(paste(nrow(exactMatch), "of 12 expected mock sequences recovered."))

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

# Create scatter plot overlaid with mismatch
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

ggsave(filename = "../supplemental/mock_plot.pdf", plot = p,
       device = "pdf", units = "cm", width = 15, height = 10,
       dpi = 300, useDingbats = FALSE)