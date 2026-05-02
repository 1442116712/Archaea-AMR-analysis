#!/bin/bash
# Pre-flight: verify all references and per-sample inputs exist
set -euo pipefail

META=/users/40335635/sharedscratch/SRA/downloads/meta
ARG_FASTA=/users/40335635/sharedscratch/SRA/downloads/ARGs_41.fasta
ARG_BLASTDB=/users/40335635/sharedscratch/SRA/downloads/blastndb/blastn_db
ARG_LEN=/users/40335635/sharedscratch/SRA/downloads/arg_db.lengths.tsv
KAIJU_DB=/users/40335635/sharedscratch/kaiju_db

module load apps/ncbiblast/2.15.0/gcc-14.1.0
module load apps/samtools/1.17/gcc-14.1.0

echo "=========================================="
echo " ARG REFERENCE CHECK"
echo "=========================================="

if [[ ! -s ${ARG_FASTA}.fai ]]; then samtools faidx ${ARG_FASTA}; fi
cut -f1,2 ${ARG_FASTA}.fai > ${ARG_LEN}

N_FA=$(grep -c '^>' ${ARG_FASTA})
N_DB=$(blastdbcmd -db ${ARG_BLASTDB} -info 2>/dev/null | grep -oP '\d+(?= sequences)')
N_LN=$(wc -l < ${ARG_LEN})

printf "  ARG fasta:     %d sequences\n" $N_FA
printf "  Blast db:      %s sequences\n" "$N_DB"
printf "  Length tbl:    %d entries\n" $N_LN

if [[ $N_FA -ne 41 || $N_LN -ne 41 ]]; then
  echo "  ⚠ FAIL: reference inconsistent"
  exit 1
fi

echo ""
echo "=========================================="
echo " KAIJU DB CHECK"
echo "=========================================="
ls -lh ${KAIJU_DB}/{kaiju_db_refseq_ref.fmi,nodes.dmp,names.dmp} 2>&1

echo ""
echo "=========================================="
echo " PER-SAMPLE INPUT INVENTORY"
echo "=========================================="

SAMPLES=(
  ERR1711948 ERR1711870 ERR1711945
  ERR694391  ERR2021508 ERR694360
  DRR066667  DRR066668  ERR2088989
  ERR3132363 ERR3132364 ERR3132366
  SRR3313097
  ERR747934  ERR747931  ERR747932
)

printf "%-15s %5s %5s %5s\n" "Sample" "mRNA" "asm" "type"
ALL_OK=1
for s in "${SAMPLES[@]}"; do
  R1=${META}/mRNA/${s}_mRNA_fwd.fq.gz
  R2=${META}/mRNA/${s}_mRNA_rev.fq.gz
  SE=${META}/mRNA/${s}_mRNA.fq.gz
  CTG=${META}/assembly/${s}/transcripts.fasta
  
  if [[ -s $R1 && -s $R2 ]]; then mrna="PE"
  elif [[ -s $SE ]]; then mrna="SE"
  else mrna="✗"; ALL_OK=0; fi
  
  if [[ -s $CTG ]]; then asm="✓"
  else asm="✗"; ALL_OK=0; fi
  
  printf "%-15s %5s %5s\n" "$s" "$mrna" "$asm"
done

echo ""
[[ $ALL_OK -eq 1 ]] && echo "  ✓ All inputs ready" || echo "  ⚠ Some inputs missing"
