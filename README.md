# bcltofrag
Pipeline TFM

#bwa-mem2 index ./data/reference/hg38.fa

#bsub -o out_bwa_index.txt -e err_bwa_index.txt -q bio -n 1 -W 1400 -M 100000 -hl -R 'rusage[mem=100000]' bash bwa-mem2 index ./data/reference/hg38.fa

* wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit -> reference for endmotifs


wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz


# About the repository

# Installation

- save Illumina output in root/data/ folder (root directory should have with >size available (x seq - total GB))

- clone repository in root directory

- install conda environment from yaml file and python kernel

# Usage

- generate reference genome file

```bash
# create containing folder within ../data
mkdir -p ../data/reference
# update path to reference of interest
wget -P ../data/reference https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# unzip file
gunzip ../data/reference/hg38.fa.gz
```

- snakefile1 - ejemplo codigo usado
- .sh - ejemplo codigo usado 
- index from referece - ejemplo codigo usado
- snakefile2 - ejemplo codigo usado

