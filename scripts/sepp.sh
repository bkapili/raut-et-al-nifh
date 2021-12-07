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