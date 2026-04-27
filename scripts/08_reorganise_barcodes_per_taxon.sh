#!/bin/bash
#SBATCH --job-name=tree_filter
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log

# Expects input files:
# 1. trees/names_samples_for_treenames.csv from R analysis: this is a list of contigs which are correct based on manual blasting
# 2. concatenated  fasta files from each barcode (do this manually)
#   - all_COI_amplicons.fa
#   - all_18S_accurate_amplicons.fa
#   - all_28S_accurate_amplicons.fa

set -euo pipefail

mkdir -p trees/logs

echo "[$(date)] Starting tree_filter job $SLURM_JOB_ID"

CSV="trees/names_samples_for_treenames.csv"

FASTA_COI="all_COI_amplicons.fa"
FASTA_18S="all_18S_accurate_amplicons.fa"
FASTA_28S="all_28S_accurate_amplicons.fa"

# Check all inputs exist before starting
for f in "$CSV" "$FASTA_COI" "$FASTA_18S" "$FASTA_28S"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required input file not found: $f" >&2
        exit 1
    fi
done

echo "[$(date)] Step 1: Building wanted-header lists..."

tail -n +2 "$CSV" | tr -d '\r' | awk -F',' '
    { barcode=$3; header=$2 }
    barcode == "18S" { print "18S_rRNA::consensus_"header > "wanted_18S.txt" }
    barcode == "28S" { print "28S_rRNA::consensus_"header > "wanted_28S.txt" }
    barcode == "COI" { print "consensus_"header > "wanted_COI.txt" }
'

echo "Wanted-header counts:"
wc -l wanted_*.txt

echo "[$(date)] Step 2: Filtering FASTA files..."

for barcode in COI 18S 28S; do
    case "$barcode" in
        COI) fasta="$FASTA_COI" ;;
        18S) fasta="$FASTA_18S" ;;
        28S) fasta="$FASTA_28S" ;;
    esac

    awk -v hfile="wanted_${barcode}.txt" '
        BEGIN {
            while ((getline line < hfile) > 0)
                wanted[line] = 1
        }
        /^>/ {
            header = substr($0, 2)
            sub(/:[0-9].*$/, "", header)
            printing = (header in wanted)
        }
        printing { print }
    ' "$fasta" > "filtered_${barcode}.fasta"

    echo "$barcode: $(grep -c '^>' filtered_${barcode}.fasta) sequences extracted"
done

echo "[$(date)] Step 3: Renaming headers and splitting by taxon..."

for barcode in COI 18S 28S; do
    case "$barcode" in
        COI) fasta="filtered_COI.fasta" ;;
        18S) fasta="filtered_18S.fasta" ;;
        28S) fasta="filtered_28S.fasta" ;;
    esac

    awk -F',' -v b="$barcode" -v csvfile="$CSV" -v barcode="$barcode" '
        BEGIN {
            while ((getline line < csvfile) > 0) {
                n = split(line, fields, ",")
                if (fields[3] != barcode) continue
                fasta_header = fields[2]
                name = fields[5]
                expected_taxon = fields[4]
                sample = fields[1]
                ns = split(sample, parts, "_")
                dataset = parts[ns]
                adapter = ""
                for (i=1; i<ns; i++) adapter = (adapter == "" ? parts[i] : adapter"_"parts[i])
                lookup[fasta_header] = name"|"adapter"|"dataset
                taxon[fasta_header] = expected_taxon
            }
        }
        /^>/ {
            header = substr($0, 2)
            sub(/:[0-9].*$/, "", header)
            sub(/^[^:]*::/, "", header)
            sub(/^consensus_/, "", header)
            if (header in lookup) {
                new_header = lookup[header]
                out_taxon = taxon[header]
                printing = 1
                outfile = "trees/" out_taxon "/" barcode ".fasta"
                system("mkdir -p trees/" out_taxon)
                print ">" new_header > outfile
            } else {
                printing = 0
            }
        }
        !(/^>/) && printing {
            print >> outfile
        }
    ' "$fasta"

    echo "$barcode: done"
done

echo "[$(date)] Job $SLURM_JOB_ID complete."
