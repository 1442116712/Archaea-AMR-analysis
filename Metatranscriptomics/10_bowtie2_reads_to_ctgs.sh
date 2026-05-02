#!/bin/bash
#SBATCH --job-name=bt2_ctg
#SBATCH --partition=k2-medpri
#SBATCH --array=0-15%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/bt2_ctg_%A_%a.out
#SBATCH --error=logs/bt2_ctg_%A_%a.err
#SBATCH --mail-user=zwu18@qub.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs

SAMPLES=(
  ERR1711948 ERR1711870 ERR1711945
  ERR694391  ERR2021508 ERR694360
  DRR066667  DRR066668  ERR2088989
  ERR3132363 ERR3132364 ERR3132366
  SRR3313097
  ERR747934  ERR747931  ERR747932
)
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}

META=/users/40335635/sharedscratch/SRA/downloads/meta
mRNA_DIR=${META}/mRNA
CONTIGS=${META}/assembly/${SAMPLE}/transcripts.fasta
ALIGN_DIR=${META}/contig_arg/${SAMPLE}
mkdir -p ${ALIGN_DIR}

echo "=========================================="
echo " ${SAMPLE} on $(hostname) at $(date)"
echo "=========================================="

module load apps/bowtie2/2.5.2/gcc-14.1.0
module load apps/samtools/1.17/gcc-14.1.0

# Detect SE/PE
R1=${mRNA_DIR}/${SAMPLE}_mRNA_fwd.fq.gz
R2=${mRNA_DIR}/${SAMPLE}_mRNA_rev.fq.gz
SE=${mRNA_DIR}/${SAMPLE}_mRNA.fq.gz
if [[ -s $R1 && -s $R2 ]]; then MODE=PE
elif [[ -s $SE ]]; then MODE=SE
else echo "ERROR: no mRNA"; exit 1; fi
echo "Mode: $MODE"

CTG_IDX=${ALIGN_DIR}/${SAMPLE}_ctg_idx
BAM=${ALIGN_DIR}/${SAMPLE}_reads_vs_ctg.bam

if [[ -s ${BAM}.bai ]]; then
  echo "[skip] already done"; exit 0
fi

# Build index
if [[ ! -s ${CTG_IDX}.1.bt2 && ! -s ${CTG_IDX}.1.bt2l ]]; then
  bowtie2-build --threads 16 ${CONTIGS} ${CTG_IDX}
fi

cd ${ALIGN_DIR}

if [[ $MODE == "PE" ]]; then
  bowtie2 -x ${CTG_IDX} -1 ${R1} -2 ${R2} \
          --very-sensitive --no-unal \
          -p 16 -X 1000 \
          2> ${SAMPLE}.bt2.log \
    | samtools view -bS -q 10 -F 4 - \
    | samtools sort -@ 4 -o ${BAM} -
else
  bowtie2 -x ${CTG_IDX} -U ${SE} \
          --very-sensitive --no-unal \
          -p 16 \
          2> ${SAMPLE}.bt2.log \
    | samtools view -bS -q 10 -F 4 - \
    | samtools sort -@ 4 -o ${BAM} -
fi
samtools index ${BAM}
samtools idxstats ${BAM} > ${SAMPLE}_ctg_idxstats.tsv

echo ""
cat ${SAMPLE}.bt2.log
echo ""
echo "Done at $(date)"
