#!/bin/bash
# Extract unclassified ARG-bearing contigs for manual NCBI BLAST validation
set -euo pipefail

module load apps/samtools/1.17/gcc-14.1.0

SAMPLES=(
  ERR1711948 ERR1711870 ERR1711945
  ERR694391  ERR2021508 ERR694360
  DRR066667  DRR066668  ERR2088989
  ERR3132363 ERR3132364 ERR3132366
  SRR3313097
  ERR747934  ERR747931  ERR747932
)

META=/users/40335635/sharedscratch/SRA/downloads/meta
UNCLASS_DIR=${META}/unclassified_arg_contigs
mkdir -p ${UNCLASS_DIR}

> ${UNCLASS_DIR}/all_unclass_arg_contigs.fasta
> ${UNCLASS_DIR}/contig_to_sample.tsv

TOTAL=0
for s in "${SAMPLES[@]}"; do
  bc=${META}/contig_arg/${s}/${s}_arg_bedcov.tsv
  named=${META}/kaiju_out/${s}_contigs/${s}_contigs_kaiju_named.out
  contigs=${META}/assembly/${s}/transcripts.fasta
  
  if [[ ! -s $bc || ! -s $named ]]; then
    echo "  $s: skipped"
    continue
  fi
  
  # Get contigs that are (i) ARG-bearing AND (ii) Kaiju-unclassified
  awk -F'\t' '$1=="U"{print $2}' $named | sort -u > /tmp/_unc_$$
  
  awk -F'\t' 'NR==FNR{u[$1]=1; next} ($1 in u){print $1}' \
    /tmp/_unc_$$ $bc | sort -u > ${UNCLASS_DIR}/${s}_unclass_arg_ids.txt
  
  rm -f /tmp/_unc_$$
  
  N=$(wc -l < ${UNCLASS_DIR}/${s}_unclass_arg_ids.txt)
  if [[ $N -gt 0 ]]; then
    samtools faidx $contigs
    while read cid; do
      samtools faidx $contigs "$cid" \
        | awk -v s="$s" -v c="$cid" 'NR==1{print ">"s"|"c} NR>1{print}' \
        >> ${UNCLASS_DIR}/all_unclass_arg_contigs.fasta
      printf "%s|%s\t%s\t%s\n" "$s" "$cid" "$s" "$cid" \
        >> ${UNCLASS_DIR}/contig_to_sample.tsv
    done < ${UNCLASS_DIR}/${s}_unclass_arg_ids.txt
    echo "  $s: $N unclassified ARG-bearing contigs"
    TOTAL=$((TOTAL + N))
  else
    echo "  $s: 0"
  fi
done

echo ""
echo "Total: $TOTAL unclassified ARG-bearing contigs"
echo "FASTA: ${UNCLASS_DIR}/all_unclass_arg_contigs.fasta"

if [[ $TOTAL -gt 0 ]]; then
  samtools faidx ${UNCLASS_DIR}/all_unclass_arg_contigs.fasta
  echo ""
  echo "=== Length distribution ==="
  awk '{
    n++; sum+=$2
    if($2>=500)c5++; if($2>=1000)c1k++; if($2>=3000)c3k++
    if($2>max)max=$2; if(min==""||$2<min)min=$2
  }
  END {
    printf "  Total: %d\n  Mean: %.0f bp\n  Min/Max: %d/%d bp\n  >=500: %d  >=1000: %d  >=3000: %d\n",
           n, sum/n, min, max, c5+0, c1k+0, c3k+0
  }' ${UNCLASS_DIR}/all_unclass_arg_contigs.fasta.fai
  
  echo ""
  echo "Action: paste FASTA contents into https://blast.ncbi.nlm.nih.gov"
  echo "  Database: nt (Nucleotide collection)"
  echo "  Program:  megablast"
  echo "  Record top hit organism for each contig in Supplementary Table S4"
fi
