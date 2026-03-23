#!/bin/bash
#SBATCH --job-name=download_reference_database
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

mkdir -p /scratch/$USER/assignment3/kraken_db
cd /scratch/$USER/assignment3/kraken_db

wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz

# Extract
tar -xvf k2_standard_08_GB_20251015.tar.gz
