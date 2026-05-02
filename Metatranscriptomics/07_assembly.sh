#!/bin/bash
# Step 7: De novo assembly with rnaSPAdes (per-sample SLURM array)
#
# Input:  mRNA/<RUN_ID>_mRNA*.fq.gz  (from step 06)
# Output: assembly/<RUN_ID>/transcripts.fasta
#
# Submit: sbatch scripts/07_assembly.sh

#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --job-name=mtx_asm
#SBATCH --error=logs/asm_%A_%a.err
#SBATCH --output=logs/asm_%A_%a.out
#SBATCH --partition=k2-medpri
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --array=1-16%3
#SBATCH --mail-user=zwu18@qub.ac.uk
#SBATCH --mail-type=END,FAIL

set -euxo pipefail

# === Get RUN_ID from selected list ===
RUN_LIST="./mgnify_search/selected_run_ids.txt"
RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$RUN_LIST")
[ -z "$RUN" ] && { echo "ERROR: no RUN for task $SLURM_ARRAY_TASK_ID"; exit 1; }

# === Paths ===
BASE=/users/40335635/sharedscratch/SRA/downloads/meta
mRNA_DIR=$BASE/mRNA
ASSEMBLY_DIR=$BASE/assembly
CONDA_BASE=/users/40335635/sharedscratch/miniconda3
CONDA_ENV=rnaspades

mkdir -p "$ASSEMBLY_DIR" logs

echo "[$(date)] === Task $SLURM_ARRAY_TASK_ID: $RUN on $(hostname) ==="
echo "TMPDIR: $TMPDIR"
echo "Free disk on TMPDIR: $(df -h $TMPDIR | tail -1)"

# === Detect PE / SE layout ===
PE_FWD=$mRNA_DIR/${RUN}_mRNA_fwd.fq.gz
PE_REV=$mRNA_DIR/${RUN}_mRNA_rev.fq.gz
SE_R=$mRNA_DIR/${RUN}_mRNA.fq.gz

if [ -f "$PE_FWD" ] && [ -f "$PE_REV" ]; then
    LAYOUT="PE"
elif [ -f "$SE_R" ]; then
    LAYOUT="SE"
else
    echo "ERROR: No mRNA input found for $RUN"
    ls -lah $mRNA_DIR/${RUN}* 2>/dev/null || true
    exit 1
fi
echo "Detected layout: $LAYOUT"

# === Skip if already done ===
FINAL_OUT=$ASSEMBLY_DIR/${RUN}
if [ -f "$FINAL_OUT/transcripts.fasta" ]; then
    echo "Assembly already exists, SKIP"
    exit 0
fi

# === Local workdir on node's TMPDIR (faster IO than NFS) ===
LOCAL_WORK=$TMPDIR/spades_${RUN}
LOCAL_READS=$TMPDIR/reads_${RUN}
rm -rf $LOCAL_WORK $LOCAL_READS
mkdir -p $LOCAL_WORK $LOCAL_READS

# === Copy reads to local TMPDIR ===
echo "[$(date)] Copying reads to TMPDIR..."
if [ "$LAYOUT" = "PE" ]; then
    cp $PE_FWD $PE_REV $LOCAL_READS/
else
    cp $SE_R $LOCAL_READS/
fi
ls -lah $LOCAL_READS/

# === Activate conda env ===
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV
echo "[$(date)] rnaSPAdes version:"
rnaspades.py --version

# === Run rnaSPAdes ===
echo "[$(date)] === Starting rnaSPAdes ($LAYOUT) ==="
if [ "$LAYOUT" = "PE" ]; then
    rnaspades.py \
        -1 $LOCAL_READS/${RUN}_mRNA_fwd.fq.gz \
        -2 $LOCAL_READS/${RUN}_mRNA_rev.fq.gz \
        -o $LOCAL_WORK \
        -t $SLURM_CPUS_PER_TASK \
        -m 100
else
    rnaspades.py \
        -s $LOCAL_READS/${RUN}_mRNA.fq.gz \
        -o $LOCAL_WORK \
        -t $SLURM_CPUS_PER_TASK \
        -m 100
fi
echo "[$(date)] === rnaSPAdes finished ==="

# === Sanity check ===
[ ! -f "$LOCAL_WORK/transcripts.fasta" ] && { echo "ERROR: no transcripts.fasta"; exit 1; }

echo "=== Assembly stats ==="
echo "Total contigs: $(grep -c '^>' $LOCAL_WORK/transcripts.fasta)"
echo "File size: $(ls -lah $LOCAL_WORK/transcripts.fasta | awk '{print $5}')"

# === Copy outputs back to sharedscratch ===
mkdir -p $FINAL_OUT
echo "[$(date)] Copying outputs to $FINAL_OUT..."
cp -v $LOCAL_WORK/transcripts.fasta              $FINAL_OUT/
cp -v $LOCAL_WORK/hard_filtered_transcripts.fasta $FINAL_OUT/ 2>/dev/null || true
cp -v $LOCAL_WORK/soft_filtered_transcripts.fasta $FINAL_OUT/ 2>/dev/null || true
cp -v $LOCAL_WORK/transcripts.paths              $FINAL_OUT/ 2>/dev/null || true
cp -v $LOCAL_WORK/spades.log                     $FINAL_OUT/ 2>/dev/null || true
# K* intermediate files NOT copied (too large)

ls -lah $FINAL_OUT/

# === Cleanup TMPDIR ===
rm -rf $LOCAL_WORK $LOCAL_READS
echo "[$(date)] === Done: $RUN ==="
