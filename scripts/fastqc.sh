# Change directory
cd ../data/raw

# Run fastqc
fastqc *.fastq

# Move to new directory
mkdir ./fastqc && mv *_fastqc* ./fastqc