#!/bin/bash
#SBATCH --job-name=kaiju_ctg
#SBATCH --partition=k2-medpri
#SBATCH --array=0-15%2
#SBATCH --cpus-per-task=16
#SBATCH --mem=140G
#SBATCH --time=8:00:00
#SBATCH --output=logs/kaiju_ctg_%A_%a.out
#SBATCH --error=logs/kaiju_ctg_%A_%a.err
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
KAIJU_DB=/users/40335635/sharedscratch/kaiju_db
NODES=${KAIJU_DB}/nodes.dmp
NAMES=${KAIJU_DB}/names.dmp
FMI=${KAIJU_DB}/kaiju_db_refseq_ref.fmi

CONTIGS=${META}/assembly/${SAMPLE}/transcripts_min300.fasta
KAIJU_DIR=${META}/kaiju_out/${SAMPLE}_contigs_min300
mkdir -p ${KAIJU_DIR}

echo "=========================================="
echo " ${SAMPLE} on $(hostname) at $(date)"
echo " Input: $(grep -c '^>' ${CONTIGS}) filtered contigs"
echo "=========================================="

source /users/40335635/sharedscratch/miniconda3/etc/profile.d/conda.sh
conda activate kaiju

KOUT=${KAIJU_DIR}/${SAMPLE}_contigs_kaiju.out
KNAMED=${KAIJU_DIR}/${SAMPLE}_contigs_kaiju_named.out
ARCH_CTG=${KAIJU_DIR}/${SAMPLE}_archaeal_contigs.txt

if [[ -s ${ARCH_CTG} ]]; then
  echo "[skip] already done"; exit 0
fi

kaiju -t ${NODES} -f ${FMI} \
      -i ${CONTIGS} \
      -o ${KOUT} \
      -z 16 -a greedy -e 5 -E 1e-5

kaiju-addTaxonNames -t ${NODES} -n ${NAMES} -p \
    -i ${KOUT} -o ${KNAMED}

awk -F'\t' '$1=="C" && $NF ~ /Archaea;/ {print $2}' ${KNAMED} > ${ARCH_CTG}

# 统计
echo ""
echo "Total contigs    : $(wc -l < ${KOUT})"
echo "Archaeal contigs : $(wc -l < ${ARCH_CTG})"
awk -F'\t' '
  $1=="C" && $NF ~ /Archaea;/   {a++}
  $1=="C" && $NF ~ /Bacteria;/  {b++}
  $1=="C" && $NF ~ /Eukaryota;/ {e++}
  $1=="C" && $NF ~ /Viruses;/   {v++}
  $1=="U"                        {u++}
  END {
    printf "  Bacteria:    %d\n  Archaea:     %d\n  Eukaryota:   %d\n  Viruses:     %d\n  Unclassified: %d\n",
           b+0, a+0, e+0, v+0, u+0
  }
' ${KNAMED}
echo ""
echo "Done at $(date)"
