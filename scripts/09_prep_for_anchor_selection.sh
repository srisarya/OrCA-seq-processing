#!/bin/bash
#SBATCH --job-name=phylo_anchor_filter
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log

set -euo pipefail

source activate phylo_r

aligned_fasta=$1
samples_fasta=$2
gene=$3

phylo_output=$(dirname "${aligned_fasta}")/${gene}
mkdir -p "${phylo_output}"

# Sanitise FASTA headers
clean_fasta="${phylo_output}/${gene}_cleaned.fa"
awk '/^>/ {
    header = substr($0, 2)
    gsub(/[^A-Za-z0-9._]/, "_", header)
    print ">" header
    next
} { print }' "${aligned_fasta}" > "${clean_fasta}"

# Extract and sanitise sample IDs from samples FASTA
sanitised_samples="${phylo_output}/sample_ids.tmp"
grep "^>" "${samples_fasta}" | sed 's/^>//' \
    | awk '{ gsub(/[^A-Za-z0-9._]/, "_", $0); print }' \
    | sort > "${sanitised_samples}"

# Build metadata: sample if in list, else anchor
metadata="${phylo_output}/${gene}_metadata.csv"
echo "label,type" > "${metadata}"
grep "^>" "${clean_fasta}" | sed 's/^>//' | awk '
    NR == FNR { sample_ids[$1] = 1; next }
    { print $0 "," ($0 in sample_ids ? "sample" : "anchor") }
' "${sanitised_samples}" - >> "${metadata}"

rm "${sanitised_samples}"

echo "  Anchors: $(grep -c ',anchor' "${metadata}" || true)"
echo "  Samples: $(grep -c ',sample' "${metadata}" || true)"
