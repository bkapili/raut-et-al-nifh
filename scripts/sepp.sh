# -------------------------------------------------------------
# Script purpose: Place ASVs on reference nifH tree from PPIT
#                 (v.1.2.0) using SEPP for subsequent taxonomic inferencing.
#
# Input:  Phyloseq object containing chimera-removed ASV sequences
#         ("psRaw.rds")
#
# Output: Tree in Newick format of ASVs placed on reference tree
#         ("sepp_placement.nwk").
# -------------------------------------------------------------

# Prepare necessary input files
Rscript sepp_input_prep.R

# Run SEPP
run_sepp.py \
  -t ../robjects/nifH_reference_tree_v2.nwk \
  -r ../robjects/nifH_reference_RAxML_info_v2.txt \
  -a ../robjects/nifH_reference_alignment_v2.fasta \
  -f ../robjects/psRaw_ASV_seqs.fasta \
  -o sepp \
  -d ../robjects \
  -seed 5935622

# Convert to Newick format
guppy tog --xml ../robjects/sepp_placement.json && mv sepp_placement.tog.xml ../robjects
gotree reformat newick \
  -i ../robjects/sepp_placement.tog.xml \
  -f phyloxml \
  -o ../robjects/sepp_placement.nwk