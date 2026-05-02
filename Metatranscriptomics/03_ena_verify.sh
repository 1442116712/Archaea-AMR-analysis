#!/bin/bash
# Step 5: Cross-validate via ENA to remove PCR amplicons mislabeled as metatranscriptomic
#
# Input:  runs_v3plus.txt (370 runs)
# Output: ena_verify_full.tsv, final_real_runs.tsv (325 runs)

set -u
cd ./mgnify_search

INPUT="runs_v3plus.txt"
ENA_FULL="ena_verify_full.tsv"
ENA_CKPT="ena_checkpoint.txt"
FINAL="final_real_runs.tsv"
FINAL_IDS="final_run_ids.txt"

[ ! -f "$ENA_CKPT" ] && touch "$ENA_CKPT"
[ ! -s "$ENA_FULL" ] && echo -e "run_accession\tlibrary_strategy\tlibrary_selection\tlibrary_source\tinstrument_model\tread_count\tbase_count\tstudy_accession\tstudy_title" > "$ENA_FULL"

while read RUN_ID; do
    grep -qxF "$RUN_ID" "$ENA_CKPT" 2>/dev/null && continue
    
    RESP=$(curl -s --max-time 30 --retry 3 \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN_ID}&result=read_run&fields=run_accession,library_strategy,library_selection,library_source,instrument_model,read_count,base_count,study_accession,study_title&format=tsv")
    
    LINE=$(echo "$RESP" | tail -n +2 | head -n 1)
    [ -z "$LINE" ] && LINE="${RUN_ID}\tNOT_FOUND\tNOT_FOUND\tNOT_FOUND\tNA\tNA\tNA\tNA\tNA"
    echo -e "$LINE" >> "$ENA_FULL"
    echo "$RUN_ID" >> "$ENA_CKPT"
    sleep 0.2
done < "$INPUT"

# Filter: remove PCR amplicons and not-found
(head -1 "$ENA_FULL"; \
 awk -F'\t' 'NR>1 && $3 != "PCR" && $3 != "NOT_FOUND"' "$ENA_FULL") > "$FINAL"
awk -F'\t' 'NR>1 {print $1}' "$FINAL" > "$FINAL_IDS"

echo "library_selection distribution:"
awk -F'\t' 'NR>1 {print $3}' "$ENA_FULL" | sort | uniq -c | sort -rn
echo ""
echo "Final non-PCR runs: $(wc -l < $FINAL_IDS)"
