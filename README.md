
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

## Overview

This repo contains the code to reproduce the processing of raw PE300
Illumina MiSeq data as performed in Raut *et al.*, in prep. The sections
below provide instructions on how to:

  - Create a `conda` environment for data processing
  - Clone this repo
  - Download the raw data
  - Execute the scripts

### conda environment setup and installations

The code below assumes you already have
[conda](https://docs.conda.io/en/latest/) installed. It will create a
new `conda` environment called `raut_et_al_nifH`, activate it, and
install all the necessary packages (and the specific versions) used
during analysis. The environment with installations will require ~1.6 Gb
of disk space.

``` bash
# Set path for environment
CONDA_PATH=$(echo "SET_PATH_HERE")

# Create conda environment
conda create --prefix $CONDA_PATH/raut_et_al_nifh
conda activate $CONDA_PATH/raut_et_al_nifh

# Install preprocessing software
conda install -c bioconda fastqc=0.11.9 cutadapt=3.5 \
  prinseq=0.20.4 figaro=1.1.2 sepp=4.4.0 \
  gotree=0.4.2 blast=2.12.0 libgit2=1.3.0

# Install R
conda install -c conda-forge r-base=4.1.1
```

### Clone git repo

The code below will clone this repo into a subdirectory named
`raut-et-al-nifh` in the directory you specify in `REPO_PATH`.

``` bash
# Set path for repo
REPO_PATH=$(echo "SET_PATH_HERE")

# Clone repo
mkdir $REPO_PATH/raut-et-al-nifh
git clone https://github.com/bkapili/raut-et-al-nifh.git $REPO_PATH/raut-et-al-nifh
```

### Download raw data

The code below will download the raw PE300 Illumina MiSeq data into a
new subdirectory of `data` named `raw` and unzip the .fastq files
(required for `prinseq`).

``` bash
# Create raw data subdirectory
mkdir $REPO_PATH/raut-et-al-nifh/data/raw && cd "$_"

# Download data from SLIMS
wget -r -nH -nc -A "DL13_YR*" --cut-dirs 5    "http://slimsdata.genomecenter.ucdavis.edu/Data/by0xe2e70/210830_M00384_0012_MS3086631-600V3/Unaligned/Project_ADAS_Dekas13_MCRA_NIFH/"

# Unzip (req'd for prinseq)
gzip -d *.fastq.gz
```

### Run scripts

The code in this section will run the script `reproduce.sh`, which
executes all the individual steps in order. In order for it to run
properly, it should be executed directly from the scripts folder. After
it runs, the phyloseq object `psNoContam.rds` is saved in a new
subdirectory named `robjects` that contains the fully processed data.

``` bash
# Change directory
cd $REPO_PATH/raut-et-al-nifh/scripts

# Execute reproduce.sh
bash reproduce.sh
```

Additional information about each script is provided in the block of
comment code at the head of each script.
