# -------------------------------------------------------------
# Script purpose: Summarize forward and reverse read quality scores 
#                 for each sample.
#
# Input:  Raw paired-end .fastq files
#
# Output: .html files containing summaries of forward and reverse
#         read quality scores per sample.
# -------------------------------------------------------------

# Change directory
cd ../data/raw

# Run fastqc
fastqc *.fastq

# Move to new directory
mkdir ./fastqc && mv *_fastqc* ./fastqc