#!/bin/bash
#SBATCH --job-name=make_biom
#SBATCH --account=def-cottenie
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=../logs/%x_%j_output.txt
#SBATCH --error=../logs/%x_%j_error.txt

module load StdEnv/2023 python/3.11

BIOM_DIR="$HOME/assignment3/results/biom"
KRAKEN_DIR="$HOME/assignment3/results/kraken2_trimmed"

mkdir -p "$BIOM_DIR"
cd "$BIOM_DIR"

if [ ! -f kraken_biom.py ]; then
    wget https://raw.githubusercontent.com/smdabdoub/kraken-biom/master/kraken_biom.py
fi

python -m pip install --user biom-format h5py

python kraken_biom.py "$KRAKEN_DIR"/*.report --fmt json -o table.biom

echo "BIOM file created:"
ls -lh table.biom
