#!/bin/bash
#SBATCH --job-name=gene_fetch
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.log

# Activate env + credentials
source activate gene-fetch
source ~/.ncbi_credentials

# ---- Define taxa + genes ----
taxon=$1
taxid=$2
gene=$3

echo "Running: $taxon ($taxid) - $gene"

outdir="trees/${taxon}/${taxon}_${gene}_anchors/"
mkdir -p ${outdir}

gene-fetch \
    --gene "$gene" \
    -s "$taxid" \
    --nucleotide-size 100 \
    --type nucleotide \
    --out "$outdir" \
    --max-sequences 5000 \
    --email "$NCBI_EMAIL" \
    --api-key "$NCBI_API_KEY"
