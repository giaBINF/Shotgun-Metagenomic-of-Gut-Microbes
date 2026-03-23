#!/bin/bash
#SBATCH --job-name=download_srr
#SBATCH --time=04:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=../logs/%j_download.txt
#SBATCH --error=../logs/%j_download_error.txt

module load StdEnv/2023 sra-toolkit/3.0.9

OUT_DIR="$HOME/assignment3/data/raw"
mkdir -p "$OUT_DIR"

for i in {72..77}
do
    echo "Downloading SRR81469${i}"
    prefetch SRR81469${i} -O "$OUT_DIR"
done
