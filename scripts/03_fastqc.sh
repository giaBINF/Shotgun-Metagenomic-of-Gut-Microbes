#!/bin/bash
# Script to check read quality for samples with FastQC

#SBATCH --job-name=fastqc_raw
#SBATCH --time=02:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=../logs/%j_fastqc.txt
#SBATCH --error=../logs/%j_fastqc_error.txt

# Load module
module load StdEnv/2023 fastqc/0.12.1

# Define variables
INPUT_DIR="$HOME/assignment3/data/fastq"
OUTPUT_DIR="$HOME/assignment3/results/qc_results/raw_reads"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# Run FastQC for all FASTQ files
fastqc "$INPUT_DIR"/*.fastq -o "$OUTPUT_DIR" -t 8
