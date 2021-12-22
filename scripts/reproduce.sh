# -------------------------------------------------------------
# Script purpose: Execute full pipeline for processing raw reads
#                 up through phyloseq object generation.
#
# Notes:          See individual scripts for more information
#                 about each step.
# -------------------------------------------------------------

# Run fastqc
bash fastqc.sh

# Trim primers
bash cutadapt.sh

# Infer ASVs
Rscript dada2.R

# Place ASVs of nifH reference tree
bash sepp.sh

# Infer taxonomy for nifH ASVs
Rscript ppit.R

# Remove contaminant ASVs
Rscript decontam.R
