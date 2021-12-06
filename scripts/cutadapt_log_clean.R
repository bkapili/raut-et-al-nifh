# Required packages
packages <- c("dplyr", "tidyr")

# Install missing packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
library(dplyr)

# Add percent reads written column
df <- read.table(file = "../data/cutadapt/log/cutadapt_log.txt", sep = "\t", header = TRUE) %>%
  filter(row_number() %% 2 != 0) %>%
  mutate(percent_written = round(as.integer(out_reads)/as.integer(in_reads)*100, 1)) %>%
  arrange(-percent_written)

# Overwrite existing log
write.table(df, file = "../data/cutadapt/log/cutadapt_log.txt", sep = "\t",quote = FALSE, col.names = TRUE, row.names = FALSE)