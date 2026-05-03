#!/bin/bash
set -euo pipefail

META=/users/40335635/sharedscratch/SRA/downloads/meta
OUT=${META}/contig_arg_min300

{
  H=0
  for f in ${OUT}/*/*_ARG_TPM_contig.tsv; do
    if [[ $H -eq 0 ]]; then head -1 $f; H=1; fi
    tail -n +2 $f
  done
} > ${OUT}/ALL_16_samples_ARG_TPM_contig_min300.tsv

NROW=$(( $(wc -l < ${OUT}/ALL_16_samples_ARG_TPM_contig_min300.tsv) - 1 ))
echo "Combined: $NROW rows (expect 656)"

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
}' ${OUT}/ALL_16_samples_ARG_TPM_contig_min300.tsv | sort

echo ""
echo "=== ARG ranking (max TPM_Global across 16) ==="
awk -F'\t' 'NR>1 {
  if($6 > max[$2]) {max[$2]=$6; smp[$2]=$1}
  total[$2] += $6
  n_pos[$2] += ($6>0 ? 1 : 0)
} END {
  printf "%-25s  %12s  %15s  %12s  %s\n",
         "ARG","Max_TPM","Sum_TPM_16","N_samples+","Top_sample"
  for(a in max) printf "%-25s  %12.2f  %15.2f  %12d  %s\n",
    a, max[a], total[a], n_pos[a], smp[a]
}' ${OUT}/ALL_16_samples_ARG_TPM_contig_min300.tsv | sort -k2 -gr

echo ""
echo "=== Archaeal-positive (expect empty) ==="
awk -F'\t' 'NR==1 || $9>0' ${OUT}/ALL_16_samples_ARG_TPM_contig_min300.tsv \
  | column -t -s $'\t'

echo ""
echo "Output: ${OUT}/ALL_16_samples_ARG_TPM_contig_min300.tsv"
