#!/bin/bash
# Step 7: Select top 3 runs per target biome by archaeal abundance
#
# Target biomes: Rumen, Marine, Hydrothermal, Wastewater, Freshwater, Human
# Output: selected_top3_per_category.tsv, selected_run_ids.txt

cd ./mgnify_search

# Merge biome info into archaea_kept_uniq_runs
awk -F'\t' 'NR==FNR {biome[$1]=$2; next} 
            FNR==1 {print $0"\tbiome"; next} 
            {print $0"\t"(biome[$3]?biome[$3]:"NA")}' \
    study_biomes.tsv archaea_kept_uniq_runs.tsv > archaea_with_biome.tsv

# Merge ENA info
awk -F'\t' '
    NR==FNR { ena[$1] = $2"\t"$3"\t"$5"\t"$6; next }
    FNR==1 { print $0"\tlibrary_strategy\tlibrary_selection\tinstrument\tread_count"; next }
    { print $0"\t"(ena[$2] ? ena[$2] : "NA\tNA\tNA\tNA") }
' final_real_runs.tsv archaea_with_biome.tsv > samples_full_info.tsv

python3 << 'EOF'
import csv

def classify_biome(biome):
    if not biome or biome == "NA":
        return None
    if "Hydrothermal vents" in biome: return "Hydrothermal"
    if "Salt crystallizer" in biome:  return "Hypersaline"
    if "Freshwater:Lake" in biome:    return "Freshwater"
    if "Engineered:Wastewater" in biome: return "Wastewater"
    # Restrict Host-associated to Human / Mammalian (incl. Rumen, gut)
    if "Host-associated:Human" in biome or "Host-associated:Mammals" in biome:
        if "Rumen" in biome:
            return "Rumen"
        return "Host-associated"
    if "Marine" in biome:             return "Marine"
    return None

samples = []
with open("samples_full_info.tsv") as f:
    for r in csv.DictReader(f, delimiter='\t'):
        cat = classify_biome(r["biome"])
        if not cat: continue
        try:
            archaea = float(r["archaea_pct"])
            reads = int(r["read_count"]) if r["read_count"] != "NA" else 0
        except: continue
        if reads < 1000: continue
        samples.append({**r, "category": cat, "archaea_pct": archaea, "read_count": reads})

from collections import defaultdict
groups = defaultdict(list)
for s in samples:
    groups[s["category"]].append(s)

with open("selected_top3_per_category.tsv", "w") as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(["category", "run_id", "study_id", "biome", "archaea_pct",
                "read_count", "library_selection", "instrument"])
    for cat in ["Host-associated", "Rumen", "Marine", "Hydrothermal",
                "Wastewater", "Freshwater"]:
        items = sorted(groups.get(cat, []), key=lambda x: -x["archaea_pct"])[:3]
        for s in items:
            w.writerow([cat, s["run_id"], s["study_id"], s["biome"],
                        s["archaea_pct"], s["read_count"],
                        s["library_selection"], s["instrument"]])
EOF

awk -F'\t' 'NR>1 {print $2}' selected_top3_per_category.tsv > selected_run_ids.txt

# Cleanup intermediates
rm -f archaea_with_biome.tsv samples_full_info.tsv

echo "Selected runs: $(wc -l < selected_run_ids.txt)"
