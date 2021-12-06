# Change directories
cd ../data/cutadapt

# Extract max and min read lengths from all trimmed files
for f in `cat sample_names.txt`; do
  prinseq-lite.pl -fastq "$f"_R1_001_CUTADAPT.fastq -stats_len |\
  grep 'max\|min' >> ./log/prinseq_log.txt
  
  prinseq-lite.pl -fastq "$f"_R2_001_CUTADAPT.fastq -stats_len |\
  grep 'max\|min' >> ./log/prinseq_log.txt; done

# Check all read lengths are identical and match expected length
TRIMMED_LEN=$(cat ./log/prinseq_log.txt | cut -f3 | uniq)
EXPECTED_LEN=273

if [ "$TRIMMED_LEN" == "$EXPECTED_LEN" ]; then
  echo "Success! Length of trimmed reads are identical and equal to specified length. Proceeding to FIGARO."; else
  echo "Warning! Length of trimmed reads are not identical and/or not equal to specified length. Inspect cutadapt output before proceeding to FIGARO."; fi
  
# Rename files
rename _CUTADAPT.fastq .fastq  *.fastq
rename DL13_OC17_ DL13-OC17- *.fastq
rename DL13_MOCK_ DL13-MOCK- *.fastq
rename DL13_BLAN_ DL13-BLAN- *.fastq
rename _EP \\-EP-CUTADAPT *.fastq
for file in *; do mv "$file" "$(echo "$file" | sed 's/\\//g')"; done

# Make new directory
mkdir ../figaro

# Run FIGARO
figaro.py \
  -i . \
  -o ../figaro \
  -a 370 \
  -f 1 \
  -r 1 \
  -m 18