#!/bin/bash
#SBATCH --job-name=trim_fastq
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=../logs/%x_%j_output.txt
#SBATCH --error=../logs/%x_%j_error.txt

module load StdEnv/2023 fastp/1.0.1

RAW_DIR="$HOME/assignment3/data/fastq"
TRIM_DIR="$HOME/assignment3/data/fastq_trimmed"
REPORT_DIR="$HOME/assignment3/results/fastp"

mkdir -p "$TRIM_DIR"
mkdir -p "$REPORT_DIR"

for i in 72 73 74 75 76 77
do
    echo "Trimming SRR81469${i}"

    fastp \
        -i "$RAW_DIR/SRR81469${i}_1.fastq" \
        -I "$RAW_DIR/SRR81469${i}_2.fastq" \
        -o "$TRIM_DIR/SRR81469${i}_1.trimmed.fastq" \
        -O "$TRIM_DIR/SRR81469${i}_2.trimmed.fastq" \
        -h "$REPORT_DIR/SRR81469${i}.html" \
        -j "$REPORT_DIR/SRR81469${i}.json" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread 8
done

echo "Trimming complete."
