#!/bin/bash
# Script to convert SRA files to FASTQ files

#SBATCH --job-name=sra_to_fastq
#SBATCH --time=04:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=../logs/%j_fastq.txt

module load StdEnv/2023 sra-toolkit/3.0.9

# Set directories
WORKING_DIR="$HOME/assignment3/data/raw"
OUT_DIR="$HOME/assignment3/data/fastq"
TEMP_DIR="/scratch/$USER/sra_tmp"

# Create directories
mkdir -p "$OUT_DIR"
mkdir -p "$TEMP_DIR"

# Loop through each SRR folder
for srr_dir in "$WORKING_DIR"/*/; do
    srr=$(basename "$srr_dir")

    echo "Processing $srr..."

    fasterq-dump "$srr_dir/${srr}.sra" \
        -O "$OUT_DIR" \
        --split-files \
        --threads 8 \
        --temp "$TEMP_DIR"
done

echo "FASTQ conversion complete!"
