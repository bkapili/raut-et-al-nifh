# -------------------------------------------------------------
# Script purpose: Remove primer sequences using cutadapt.
#
# Input:  Raw paired-end .fastq files
#
# Outputs:
#   * Primer-removed .fastq files in a new directory named "cutadapt"
#   * A tab-delimited summary of run statistics ("cutadapt_log.txt")
# -------------------------------------------------------------

# Create new directory
mkdir -p ../data/cutadapt/log

# Set path variables
RAW_PATH=$(echo '../data/raw')
CUT_PATH=$(echo '../data/cutadapt')

# Create file of sample names to loop through
ls "$RAW_PATH"/*.fastq |\
  sed 's/_R1_001.fastq//' |\
  sed 's/_R2_001.fastq//' |\
  sed 's/.*\///' | uniq > "$CUT_PATH"/sample_names.txt

# Set read length filter
EXPECTED_LEN=273

# Loop through samples running cutadapt and print to log
for f in `cat "$CUT_PATH"/sample_names.txt`; do
  cutadapt \
    --report=minimal \
    -g GGHAARGGHGGHATHGGNAARTC \
    -G GGCATNGCRAANCCVCCRCANAC \
    --discard-untrimmed \
    --max-n=0 \
    --match-read-wildcards \
    -e 0.1 \
    -m $EXPECTED_LEN -M $EXPECTED_LEN \
    -o "$CUT_PATH"/"$f"_R1_001_CUTADAPT.fastq \
    -p "$CUT_PATH"/"$f"_R2_001_CUTADAPT.fastq \
    "$RAW_PATH"/"$f"_R1_001.fastq "$RAW_PATH"/"$f"_R2_001.fastq > tmp.txt
  
  sed -i '1 s/^/sample\t/' tmp.txt
  sed -i "2 s/^/$f\t/" tmp.txt
  
  cat tmp.txt >> "$CUT_PATH"/log/cutadapt_log.txt; done
  
# Remove tmp file
rm tmp.txt

# Run script to reformat cutadapt log
Rscript ./cutadapt_log_clean.R

# Rename files
rename _CUTADAPT.fastq .fastq  *.fastq
rename DL13_YR_ DL13-YR- *.fastq
rename _EP \\-EP-CUTADAPT *.fastq
for file in *; do mv "$file" "$(echo "$file" | sed 's/\\//g')"; done