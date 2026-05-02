META=/users/40335635/sharedscratch/SRA/downloads/meta
SAMPLE_LEN_FILE=${META}/sample_mean_read_length.tsv

SAMPLES=(
  ERR1711948 ERR1711870 ERR1711945
  ERR694391  ERR2021508 ERR694360
  DRR066667  DRR066668  ERR2088989
  ERR3132363 ERR3132364 ERR3132366
  SRR3313097
  ERR747934  ERR747931  ERR747932
)

> ${SAMPLE_LEN_FILE}
echo "Computing mean post-fastp read length per sample (full file scan)..."

for s in "${SAMPLES[@]}"; do
  R1=${META}/mRNA/${s}_mRNA_fwd.fq.gz
  SE=${META}/mRNA/${s}_mRNA.fq.gz
  if   [[ -s $R1 ]]; then F=$R1; type=PE
  elif [[ -s $SE ]]; then F=$SE; type=SE
  else echo "$s: NO FILE"; continue; fi
  
  # 全文件扫描(不是 head),最准
  mean=$(zcat $F | awk 'NR%4==2 {sum+=length($0); n++} END{printf "%.2f", sum/n}')
  printf "%s\t%s\t%s\n" "$s" "$mean" "$type" >> ${SAMPLE_LEN_FILE}
  printf "  %-15s %s bp (%s)\n" "$s" "$mean" "$type"
done

echo ""
echo "Output: ${SAMPLE_LEN_FILE}"
echo ""
cat ${SAMPLE_LEN_FILE}

#!/bin/bash
# Compute TPM_Global and TPM_Archaeal with sample-specific mean read length
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
ARG_LEN=/users/40335635/sharedscratch/SRA/downloads/arg_db.lengths.tsv
SAMPLE_LEN_FILE=${META}/sample_mean_read_length.tsv

# Validate: sample length file exists with all 16 samples
if [[ ! -s ${SAMPLE_LEN_FILE} ]]; then
  echo "ERROR: ${SAMPLE_LEN_FILE} missing. Run step 1 first."
  exit 1
fi

for SAMPLE in "${SAMPLES[@]}"; do
  echo "==========================================="
  echo " ${SAMPLE}"
  echo "==========================================="
  
  # Sample-specific mean read length
  MEAN_READ_LEN=$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2}' ${SAMPLE_LEN_FILE})
  if [[ -z "$MEAN_READ_LEN" ]]; then
    echo "  [SKIP] no mean read length found for $SAMPLE"; continue
  fi
  echo "  mean post-fastp read length: ${MEAN_READ_LEN} bp"
  
  KAIJU_DIR=${META}/kaiju_out/${SAMPLE}_contigs_min300
  ARCH_CTG=${KAIJU_DIR}/${SAMPLE}_archaeal_contigs.txt
  BLAST_DIR=${META}/blastn_arg_min300/${SAMPLE}
  BLAST_FILT=${BLAST_DIR}/${SAMPLE}_blastn_filt.tsv
  ALIGN_DIR=${META}/contig_arg_min300/${SAMPLE}
  IDXSTATS=${ALIGN_DIR}/${SAMPLE}_ctg_idxstats.tsv
  BAM=${ALIGN_DIR}/${SAMPLE}_reads_vs_ctg.bam
  
  ok=1
  for f in $ARCH_CTG $BLAST_FILT $IDXSTATS $BAM; do
    [[ ! -s $f ]] && { echo "  [SKIP] missing $f"; ok=0; }
  done
  [[ $ok -eq 0 ]] && continue
  
  cd ${ALIGN_DIR}
  
  # Step 1: BED + bedcov
  N_HSP=$(wc -l < ${BLAST_FILT})
  if [[ ${N_HSP} -gt 0 ]]; then
    awk -F'\t' 'BEGIN{OFS="\t"}
      { s=($7<$8)?$7:$8; e=($7<$8)?$8:$7; print $1, s-1, e, $2 }
    ' ${BLAST_FILT} | sort -k1,1 -k2,2n \
      > ${SAMPLE}_arg_hits.sorted.bed
    samtools bedcov ${SAMPLE}_arg_hits.sorted.bed ${BAM} \
      > ${SAMPLE}_arg_bedcov.tsv
  else
    : > ${SAMPLE}_arg_hits.sorted.bed
    : > ${SAMPLE}_arg_bedcov.tsv
  fi
  
  # Step 2: archaeal subset
  awk -F'\t' 'NR==FNR{arch[$1]=1; next} ($1 in arch)' \
    ${ARCH_CTG} ${SAMPLE}_arg_bedcov.tsv \
    > ${SAMPLE}_arg_bedcov_archaeal.tsv
  
  # Step 3: ΣRPK denominators (idxstats reads, not affected by read length)
  SUM_RPK_ALL=$(awk -F'\t' 'BEGIN{s=0} $3>0 {s += $3/($2/1000)} END{printf "%.6f", s}' ${IDXSTATS})
  SUM_RPK_ARCH=$(awk -F'\t' 'BEGIN{s=0}
    NR==FNR { arch[$1]=1; next }
    ($1 in arch) && $3>0 { s += $3/($2/1000) }
    END { printf "%.6f", s }
  ' ${ARCH_CTG} ${IDXSTATS})
  
  # Step 4: per-ARG RPK using sample-specific MEAN_READ_LEN
  for SCOPE in global archaeal; do
    if [[ $SCOPE == "global" ]]; then INPUT_BC=${SAMPLE}_arg_bedcov.tsv
    else INPUT_BC=${SAMPLE}_arg_bedcov_archaeal.tsv; fi
    
    awk -v rl=${MEAN_READ_LEN} -F'\t' '
    BEGIN{OFS="\t"}
    NR==FNR { arglen[$1]=$2; next }
    {
      arg=$4; bp_cov=$5
      arg_total[arg] += bp_cov / rl
    }
    END {
      for (a in arg_total) if (a in arglen)
        printf "%s\t%d\t%.4f\t%.6f\n", a, arglen[a], arg_total[a],
               arg_total[a]/(arglen[a]/1000)
    }
    ' ${ARG_LEN} ${INPUT_BC} > ${SAMPLE}_arg_rpk_${SCOPE}.tsv
  done
  
  # Step 5: TPM table
  awk -v sum_all="$SUM_RPK_ALL" \
      -v sum_arch="$SUM_RPK_ARCH" \
      -v sample="$SAMPLE" \
      -v rl="$MEAN_READ_LEN" \
      -F'\t' '
  BEGIN { OFS="\t"
    print "Sample","Gene_ID","Length","Reads_All","RPK_All","TPM_Global",
          "Reads_Archaeal","RPK_Archaeal","TPM_Archaeal","Mean_read_len"
  }
  ARGIND==1 { arch_reads[$1]=$3; arch_rpk[$1]=$4; next }
  ARGIND==2 {
    arg=$1; len=$2; reads=$3; rpk=$4
    tpm_g = (sum_all>0)  ? rpk/sum_all*1e6 : 0
    ar    = (arg in arch_reads) ? arch_reads[arg] : 0
    arpk  = (arg in arch_rpk)   ? arch_rpk[arg]   : 0
    tpm_a = (sum_arch>0) ? arpk/sum_arch*1e6 : 0
    printf "%s\t%s\t%d\t%.2f\t%.4f\t%.4f\t%.2f\t%.4f\t%.4f\t%.2f\n",
           sample, arg, len, reads, rpk, tpm_g, ar, arpk, tpm_a, rl
    seen[arg]=1
  }
  ARGIND==3 && !($1 in seen) {
    printf "%s\t%s\t%d\t0\t0.0000\t0.0000\t0\t0.0000\t0.0000\t%.2f\n", sample, $1, $2, rl
  }
  ' ${SAMPLE}_arg_rpk_archaeal.tsv \
    ${SAMPLE}_arg_rpk_global.tsv \
    ${ARG_LEN} \
    > ${SAMPLE}_ARG_TPM_contig.raw
  
  {
    head -1 ${SAMPLE}_ARG_TPM_contig.raw
    tail -n +2 ${SAMPLE}_ARG_TPM_contig.raw | sort -k6,6 -gr
  } > ${SAMPLE}_ARG_TPM_contig.tsv
  rm ${SAMPLE}_ARG_TPM_contig.raw
  
  ROWS=$(( $(wc -l < ${SAMPLE}_ARG_TPM_contig.tsv) - 1 ))
  GP=$(awk -F'\t' 'NR>1 && $6>0' ${SAMPLE}_ARG_TPM_contig.tsv | wc -l)
  AP=$(awk -F'\t' 'NR>1 && $9>0' ${SAMPLE}_ARG_TPM_contig.tsv | wc -l)
  printf "  >>> %s rows=%d  TPM_Global+=%d  TPM_Archaeal+=%d  read_len=%s\n" \
    "$SAMPLE" "$ROWS" "$GP" "$AP" "$MEAN_READ_LEN"
done

echo ""
echo "All samples processed at $(date)"
