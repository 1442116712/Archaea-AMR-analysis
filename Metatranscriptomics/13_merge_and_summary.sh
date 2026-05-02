#!/bin/bash
# Merge 16 sample TPM tables + cross-sample summary tables
set -euo pipefail

META=/users/40335635/sharedscratch/SRA/downloads/meta
OUT=${META}/contig_arg

# Combine 16 tables
{
  H=0
  for f in ${OUT}/*/*_ARG_TPM_contig.tsv; do
    if [[ $H -eq 0 ]]; then head -1 $f; H=1; fi
    tail -n +2 $f
  done
} > ${OUT}/ALL_16_samples_ARG_TPM_contig.tsv

NROW=$(( $(wc -l < ${OUT}/ALL_16_samples_ARG_TPM_contig.tsv) - 1 ))
echo "Combined: $NROW rows (expect 16×41 = 656)"

# Per-sample summary
echo ""
echo "=== Per-sample summary ==="
awk -F'\t' 'NR>1 {
  total[$1]++
  if($6>0) g[$1]++
  if($9>0) a[$1]++
} END {
  printf "%-15s  %5s  %12s  %15s\n","Sample","ARGs","TPM_Global+","TPM_Archaeal+"
  for(k in total) printf "%-15s  %5d  %12d  %15d\n",
    k, total[k], (k in g ? g[k] : 0), (k in a ? a[k] : 0)
}' ${OUT}/ALL_16_samples_ARG_TPM_contig.tsv | sort

# Cross-sample ARG ranking
echo ""
echo "=== ARG ranking by max TPM_Global across 16 samples ==="
awk -F'\t' 'NR>1 {
  if($6 > max[$2]) {max[$2]=$6; smp[$2]=$1}
  total[$2] += $6
  n_pos[$2] += ($6>0 ? 1 : 0)
} END {
  printf "%-25s  %12s  %15s  %12s  %s\n",
         "ARG","Max_TPM","Sum_TPM_16","N_samples+","Top_sample"
  for(a in max) printf "%-25s  %12.2f  %15.2f  %12d  %s\n",
    a, max[a], total[a], n_pos[a], smp[a]
}' ${OUT}/ALL_16_samples_ARG_TPM_contig.tsv | sort -k2 -gr

# Archaeal-positive (expect empty)
echo ""
echo "=== Archaeal-positive ARG hits (16 samples) ==="
awk -F'\t' 'NR==1 || $9>0' ${OUT}/ALL_16_samples_ARG_TPM_contig.tsv \
  | column -t -s $'\t'

# By habitat
echo ""
echo "=== By habitat ==="
declare -A hab=(
  [ERR1711948]=Marine    [ERR1711870]=Marine    [ERR1711945]=Marine
  [ERR694391]=Hydroth    [ERR2021508]=Hydroth   [ERR694360]=Hydroth
  [DRR066667]=Wastewater [DRR066668]=Wastewater [ERR2088989]=Wastewater
  [ERR3132363]=Freshwater [ERR3132364]=Freshwater [ERR3132366]=Freshwater
  [SRR3313097]=Human_gut
  [ERR747934]=Rumen      [ERR747931]=Rumen      [ERR747932]=Rumen
)
{
  for s in "${!hab[@]}"; do
    f=${OUT}/${s}/${s}_ARG_TPM_contig.tsv
    [[ ! -s $f ]] && continue
    g=$(awk -F'\t' 'NR>1 && $6>0' $f | wc -l)
    a=$(awk -F'\t' 'NR>1 && $9>0' $f | wc -l)
    printf "%s\t%s\t%d\t%d\n" "${hab[$s]}" "$s" "$g" "$a"
  done
} | sort | awk -F'\t' '
  BEGIN { OFS="\t"; print "Habitat","N_samples","ARG_pos","Archaeal_pos" }
  { h=$1; n[h]++; tg[h]+=$3; ta[h]+=$4 }
  END {
    for(h in n) printf "%-12s\t%d\t%d\t%d\n", h, n[h], tg[h], ta[h]
  }
'

echo ""
echo "Output: ${OUT}/ALL_16_samples_ARG_TPM_contig.tsv"
echo "Done at $(date)"
