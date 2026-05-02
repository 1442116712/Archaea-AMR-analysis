#!/bin/bash
# Step 6: Annotate biome for each unique study via MGnify API
#
# Output: study_biomes.tsv

cd ./mgnify_search

# Extract unique studies from final_real_runs.tsv
awk -F'\t' 'NR>1 {print $8}' final_real_runs.tsv | sort -u > final_unique_studies.txt

echo -e "study_id\tbiome\tstudy_name" > study_biomes.tsv

while read STUDY; do
    RESP=$(curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/studies/${STUDY}")
    BIOME=$(echo "$RESP" | jq -r '.data.relationships.biomes.data[0].id // "NA"')
    NAME=$(echo "$RESP" | jq -r '.data.attributes."study-name" // "NA"')
    echo -e "${STUDY}\t${BIOME}\t${NAME}" >> study_biomes.tsv
    sleep 0.3
done < final_unique_studies.txt

echo "Biome distribution:"
awk -F'\t' 'NR>1 {print $2}' study_biomes.tsv | sort | uniq -c | sort -rn
