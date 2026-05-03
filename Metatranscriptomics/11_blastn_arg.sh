#!/bin/bash
#SBATCH --job-name=blastn_arg
#SBATCH --partition=k2-medpri
#SBATCH --array=0-15%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/blastn_arg_%A_%a.out
#SBATCH --error=logs/blastn_arg_%A_%a.err
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
ARG_BLASTDB=/users/40335635/sharedscratch/SRA/downloads/blastndb/blastn_db
CONTIGS=${META}/assembly/${SAMPLE}/transcripts_min300.fasta
BLAST_DIR=${META}/blastn_arg_min300/${SAMPLE}
mkdir -p ${BLAST_DIR}

echo "=========================================="
echo " ${SAMPLE} on $(hostname) at $(date)"
echo "=========================================="

module load apps/ncbiblast/2.15.0/gcc-14.1.0

BLAST_RAW=${BLAST_DIR}/${SAMPLE}_blastn_all.tsv
BLAST_FILT=${BLAST_DIR}/${SAMPLE}_blastn_filt.tsv

blastn -query ${CONTIGS} -db ${ARG_BLASTDB} \
  -evalue 1e-6 \
  -perc_identity 70 \
  -num_threads 8 \
  -max_target_seqs 41 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp" \
  -out ${BLAST_RAW}

awk -F'\t' '$16 >= 60' ${BLAST_RAW} > ${BLAST_FILT}

echo "raw HSPs:    $(wc -l < ${BLAST_RAW})"
echo "filt HSPs:   $(wc -l < ${BLAST_FILT})"
echo "unique ARGs: $(cut -f2 ${BLAST_FILT} | sort -u | wc -l)"
echo "Done at $(date)"
