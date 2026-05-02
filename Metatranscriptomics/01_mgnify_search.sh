#!/bin/bash
# Step 1-3: Search MGnify for metatranscriptomic samples with archaea ≥ 1%
# 
# Output:
#   all_analyses.tsv             - All metatranscriptomic analyses (2651)
#   archaea_results.tsv          - With archaea abundance check
#   archaea_kept.tsv             - Records with archaea ≥ 1% (671)
#   archaea_kept_uniq_runs.tsv   - Deduplicated unique runs (656)

set -u

OUTPUT_DIR="./mgnify_search"
mkdir -p "$OUTPUT_DIR"

ALL_ANALYSES="${OUTPUT_DIR}/all_analyses.tsv"
RESULTS="${OUTPUT_DIR}/archaea_results.tsv"
KEPT="${OUTPUT_DIR}/archaea_kept.tsv"
KEPT_UNIQ="${OUTPUT_DIR}/archaea_kept_uniq_runs.tsv"
CHECKPOINT="${OUTPUT_DIR}/checkpoint.txt"

VERSIONS=("1.0" "2.0" "3.0" "4.0" "4.1" "5.0")
THRESHOLD="1.0"

# === Step 1: Collect all metatranscriptomic analyses ===
echo "Step 1: Collecting analyses from MGnify..."
echo -e "analysis_id\trun_id\tstudy_id\tpipeline\tbiome" > "$ALL_ANALYSES"

for VERSION in "${VERSIONS[@]}"; do
    URL="https://www.ebi.ac.uk/metagenomics/api/v1/analyses?experiment_type=metatranscriptomic&pipeline_version=${VERSION}&page_size=100"
    while [ -n "$URL" ] && [ "$URL" != "null" ]; do
        RESP=$(curl -s --max-time 60 --retry 3 "$URL")
        echo "$RESP" | jq empty 2>/dev/null || break
        echo "$RESP" | jq -r --arg ver "$VERSION" '
            .data[] | [
                .id,
                .relationships.run.data.id // "NA",
                .relationships.study.data.id // "NA",
                $ver,
                .attributes."biome-name" // "NA"
            ] | @tsv
        ' >> "$ALL_ANALYSES"
        URL=$(echo "$RESP" | jq -r '.links.next // empty')
        sleep 0.3
    done
done

# === Step 2: Check archaea abundance per analysis ===
echo "Step 2: Checking archaea abundance..."
[ ! -f "$CHECKPOINT" ] && touch "$CHECKPOINT"
[ ! -s "$RESULTS" ] && echo -e "analysis_id\trun_id\tstudy_id\tpipeline\tbiome\tarchaea_pct\tkeep" > "$RESULTS"

get_archaea_pct() {
    local aid=$1
    local pipeline=$2
    
    case "$pipeline" in
        "1.0"|"2.0"|"3.0") local endpoints="taxonomy" ;;
        *) local endpoints="taxonomy/ssu taxonomy/lsu" ;;
    esac
    
    local max_pct="0"
    local has_data=0
    
    for ep in $endpoints; do
        local all_items='[]'
        local url="https://www.ebi.ac.uk/metagenomics/api/v1/analyses/${aid}/${ep}?page_size=500"
        while [ -n "$url" ] && [ "$url" != "null" ]; do
            local resp=$(curl -s --max-time 30 --retry 2 "$url")
            echo "$resp" | jq empty 2>/dev/null || break
            all_items=$(echo "$resp" | jq --argjson prev "$all_items" '$prev + .data')
            url=$(echo "$resp" | jq -r '.links.next // empty')
            sleep 0.1
        done
        
        local count=$(echo "$all_items" | jq 'length')
        if [ "$count" -gt 0 ]; then
            has_data=1
            # Match Archaea across multiple field names (compatibility across versions)
            local pct=$(echo "$all_items" | jq -r '
                [.[].attributes] as $all
                | ($all | map(.count // 0) | add // 0) as $total
                | ($all | map(
                    select(
                        .hierarchy."super kingdom" == "Archaea" or
                        .hierarchy.superkingdom == "Archaea" or
                        .hierarchy.kingdom == "Archaea" or
                        .domain == "Archaea"
                    ) | .count // 0
                  ) | add // 0) as $archaea
                | if ($total > 0) then ($archaea / $total * 100 | tostring) else "0" end
            ')
            awk "BEGIN {exit !($pct > $max_pct)}" 2>/dev/null && max_pct=$pct
        fi
    done
    
    if [ $has_data -eq 0 ]; then echo "NA"; else echo "$max_pct"; fi
}

while IFS=$'\t' read -r AID RUN_ID STUDY PIPELINE BIOME; do
    [ "$AID" = "analysis_id" ] && continue
    grep -qxF "$AID" "$CHECKPOINT" 2>/dev/null && continue
    
    PCT=$(get_archaea_pct "$AID" "$PIPELINE")
    if [ "$PCT" = "NA" ]; then
        KEEP="NO_DATA"; PCT_FMT="NA"
    else
        PCT_FMT=$(printf "%.3f" "$PCT")
        if awk "BEGIN {exit !($PCT >= $THRESHOLD)}" 2>/dev/null; then
            KEEP="YES"
        else
            KEEP="NO"
        fi
    fi
    
    echo -e "${AID}\t${RUN_ID}\t${STUDY}\t${PIPELINE}\t${BIOME}\t${PCT_FMT}\t${KEEP}" >> "$RESULTS"
    echo "$AID" >> "$CHECKPOINT"
done < "$ALL_ANALYSES"

# === Step 3: Filter and deduplicate ===
echo "Step 3: Deduplicating..."
(head -1 "$RESULTS"; awk -F'\t' 'NR>1 && $7=="YES"' "$RESULTS") > "$KEPT"

# Keep one record per run (highest pipeline version)
(head -1 "$KEPT"; \
 awk -F'\t' 'NR>1' "$KEPT" \
    | sort -t$'\t' -k2,2 -k4,4rn \
    | awk -F'\t' '!seen[$2]++') > "$KEPT_UNIQ"

echo "Done. Kept $(($(wc -l < $KEPT_UNIQ) - 1)) unique runs."
