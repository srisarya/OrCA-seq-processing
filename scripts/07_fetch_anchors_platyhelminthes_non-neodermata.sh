# this was run in an interactive node, in the Platyhelminthes folder

while IFS=',' read -r taxon taxid; do
    gene-fetch \
        --gene 28S \
        -s "$taxid" \
        --nucleotide-size 100 \
        --type nucleotide \
        --out "${taxon}_28S" \
        --max-sequences 500 \
        --email "$NCBI_EMAIL" \
        --api-key "$NCBI_API_KEY"
done < non_neodermata.csv
