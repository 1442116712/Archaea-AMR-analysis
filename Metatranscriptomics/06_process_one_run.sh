#!/bin/bash
# Step 8: Per-sample pipeline (SLURM array job)
#   download → fastp → SortMeRNA → keep mRNA only
#
# Submit: sbatch 06_process_one_run.sh

#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --job-name=mtx_run
#SBATCH --error=logs/stderr_%A_%a.txt
#SBATCH --output=logs/stdout_%A_%a.txt
#SBATCH --partition=k2-medpri
#SBATCH --time=24:00:00
#SBATCH --array=1-16%3

set -eu

RUN_LIST="./mgnify_search/selected_run_ids.txt"
RUN_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$RUN_LIST")

[ -z "$RUN_ID" ] && { echo "ERROR: no RUN_ID for task $SLURM_ARRAY_TASK_ID"; exit 1; }
echo "[$(date)] Task $SLURM_ARRAY_TASK_ID: $RUN_ID"

# === Paths ===
BASE="/users/40335635/sharedscratch/SRA/downloads/meta"
RAW_DIR="${BASE}/raw_fastq"
FASTP_DIR="${BASE}/fastp_out"
MRNA_DIR="${BASE}/mRNA"
mkdir -p "$RAW_DIR" "$FASTP_DIR" "$MRNA_DIR"

# === Tools ===
FASTP="/users/40335635/sharedscratch/SRA/downloads/fastp"
SMR_DB="/users/40335635/sharedscratch/sortmerna_db/rRNA_databases_v4/smr_v4.3_default_db.fasta"
source /users/40335635/sharedscratch/miniconda3/etc/profile.d/conda.sh

# === Skip if already done ===
ls "${MRNA_DIR}/${RUN_ID}_mRNA"*.fastq.gz 2>/dev/null | grep -q . && \
    { echo "$RUN_ID: mRNA exists, SKIP"; exit 0; }

# === Determine if to skip download+fastp ===
R1_FASTP="${FASTP_DIR}/${RUN_ID}_1.fastp.fastq.gz"
R2_FASTP="${FASTP_DIR}/${RUN_ID}_2.fastp.fastq.gz"
SE_FASTP="${FASTP_DIR}/${RUN_ID}.fastp.fastq.gz"
SKIP_DL=0; MODE=""
if [ -f "$R1_FASTP" ] && [ -f "$R2_FASTP" ]; then SKIP_DL=1; MODE="PE"
elif [ -f "$SE_FASTP" ]; then SKIP_DL=1; MODE="SE"
fi

# === Download + fastp ===
if [ $SKIP_DL -eq 0 ]; then
    URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN_ID}&result=read_run&fields=fastq_ftp&format=tsv" \
        | tail -n 1 | awk -F'\t' '{print $2}')
    [ -z "$URLS" ] && { echo "ERROR: no URLs"; exit 1; }
    
    cd "$RAW_DIR"
    IFS=';' read -ra URL_ARR <<< "$URLS"
    for url in "${URL_ARR[@]}"; do
        aria2c -x 4 -c --auto-file-renaming=false --max-tries=5 "https://${url}"
    done
    
    R1="${RAW_DIR}/${RUN_ID}_1.fastq.gz"
    R2="${RAW_DIR}/${RUN_ID}_2.fastq.gz"
    SE="${RAW_DIR}/${RUN_ID}.fastq.gz"
    if [ -f "$R1" ] && [ -f "$R2" ]; then MODE="PE"
    elif [ -f "$SE" ]; then MODE="SE"
    else echo "ERROR: files not found"; exit 1; fi
    
    if [ "$MODE" = "PE" ]; then
        "$FASTP" -i "$R1" -I "$R2" -o "$R1_FASTP" -O "$R2_FASTP" \
            --json "${FASTP_DIR}/${RUN_ID}.fastp.json" \
            --html "${FASTP_DIR}/${RUN_ID}.fastp.html" \
            --thread $SLURM_CPUS_PER_TASK --detect_adapter_for_pe
    else
        "$FASTP" -i "$SE" -o "$SE_FASTP" \
            --json "${FASTP_DIR}/${RUN_ID}.fastp.json" \
            --html "${FASTP_DIR}/${RUN_ID}.fastp.html" \
            --thread $SLURM_CPUS_PER_TASK
    fi
    rm -f "$R1" "$R2" "$SE"
fi

# === SortMeRNA (per-sample independent workdir) ===
conda activate /users/40335635/sharedscratch/.conda/envs/sortmerna_legacy
WORKDIR="${FASTP_DIR}/smr_workdir_${RUN_ID}"
rm -rf "$WORKDIR"; mkdir -p "$WORKDIR"

if [ "$MODE" = "PE" ]; then
    sortmerna --ref "$SMR_DB" \
        --reads "$R1_FASTP" --reads "$R2_FASTP" \
        --workdir "$WORKDIR" \
        --aligned "${MRNA_DIR}/${RUN_ID}_rRNA" \
        --other "${MRNA_DIR}/${RUN_ID}_mRNA" \
        --fastx --paired_in --out2 \
        --threads $SLURM_CPUS_PER_TASK --num_alignments 1
else
    sortmerna --ref "$SMR_DB" \
        --reads "$SE_FASTP" \
        --workdir "$WORKDIR" \
        --aligned "${MRNA_DIR}/${RUN_ID}_rRNA" \
        --other "${MRNA_DIR}/${RUN_ID}_mRNA" \
        --fastx \
        --threads $SLURM_CPUS_PER_TASK --num_alignments 1
fi

rm -f "$R1_FASTP" "$R2_FASTP" "$SE_FASTP"
rm -rf "$WORKDIR"
rm -f "${MRNA_DIR}/${RUN_ID}_rRNA"*

# === Compress mRNA ===
for f in "${MRNA_DIR}/${RUN_ID}_mRNA"*.fastq; do
    [ -f "$f" ] && pigz -p $SLURM_CPUS_PER_TASK "$f"
done

echo "[$(date)] DONE: $RUN_ID"
