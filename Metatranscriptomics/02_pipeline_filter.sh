#!/bin/bash
# Step 4: Filter to pipeline v3.0+ only
# Reason: v1.0/v2.0 used outdated, archived reference databases
#
# Input:  archaea_kept_uniq_runs.tsv
# Output: runs_v3plus.txt (370 runs)

cd ./mgnify_search

awk -F'\t' 'NR>1 && ($4=="3.0" || $4=="4.0" || $4=="4.1" || $4=="5.0") {print $2}' \
    archaea_kept_uniq_runs.tsv | sort -u > runs_v3plus.txt

echo "v3+ runs: $(wc -l < runs_v3plus.txt)"
