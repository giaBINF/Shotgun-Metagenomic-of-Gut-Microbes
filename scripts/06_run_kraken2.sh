#!/bin/bash
#SBATCH --job-name=kraken2_run
#SBATCH --account=def-cottenie
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=../logs/%x_%j_output.txt
#SBATCH --error=../logs/%x_%j_error.txt

module load StdEnv/2023 kraken2/2.1.6

DB_DIR="/scratch/$USER/assignment3/kraken_db"
FASTQ_DIR="$HOME/assignment3/data/fastq"
OUT_DIR="$HOME/assignment3/results/kraken2"

mkdir -p "$OUT_DIR"

if [ ! -f "$DB_DIR/hash.k2d" ]; then
    echo "Kraken2 database not found in $DB_DIR"
    exit 1
fi

for i in 72 73 74 75 76 77
do
    echo "Running Kraken2 on SRR81469${i}"

    kraken2 \
        --db "$DB_DIR" \
        --confidence 0.15 \
        --memory-mapping \
        --paired \
        --threads 8 \
        --output "$OUT_DIR/SRR81469${i}.kraken" \
        --report "$OUT_DIR/SRR81469${i}.report" \
        "$FASTQ_DIR/SRR81469${i}_1.fastq" \
        "$FASTQ_DIR/SRR81469${i}_2.fastq"
done

echo "Kraken2 run complete."
