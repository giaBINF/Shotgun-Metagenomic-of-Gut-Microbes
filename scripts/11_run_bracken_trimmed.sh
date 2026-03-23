#!/bin/bash
#SBATCH --job-name=bracken_trimmed
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=../logs/%x_%j_output.txt
#SBATCH --error=../logs/%x_%j_error.txt

module load StdEnv/2023 bracken

DB_DIR="/scratch/$USER/assignment3/kraken_db"
KRAKEN_DIR="$HOME/assignment3/results/kraken2_trimmed"

for i in 72 73 74 75 76 77
do
    echo "Running Bracken on trimmed SRR81469${i}"

    bracken \
        -d "$DB_DIR" \
        -i "$KRAKEN_DIR/SRR81469${i}.report" \
        -o "$KRAKEN_DIR/SRR81469${i}.bracken" \
        -r 300 \
        -l S
done

echo "Bracken complete."
