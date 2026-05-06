# Archaeal ARG Metatranscriptome Pipeline

> **Detection of archaeal antibiotic resistance gene (ARG) expression across 16 globally distributed metatranscriptomes from 7 habitats.**

A reproducible pipeline for identifying and quantifying ARG transcripts on archaeal contigs in metatranscriptomic data, with strong bacterial positive controls and orthogonal validation.

---

## Overview

This repository contains the full computational pipeline used to:

1. Select archaea-rich metatranscriptome samples from [MGnify](https://www.ebi.ac.uk/metagenomics/) across diverse habitats
2. Process raw reads (QC + rRNA depletion) and assemble *de novo* transcriptomes
3. Classify contigs taxonomically and identify the archaeal subset
4. Detect ARG transcripts via BLASTn against a custom 41-ARG reference panel
5. Quantify ARG expression as TPM in two scopes (archaeal-only and whole-sample)

The pipeline is designed for SLURM-managed HPC clusters but the analysis steps can be adapted to local execution.

---

## Key findings (summary)

Across **16 metatranscriptomes** spanning **7 habitats** (marine, hydrothermal, wastewater, freshwater, human gut, rumen):

- **71,002 high-confidence archaeal contigs** (≥300 bp) recovered across samples
- **72 ARG transcripts** detected on bacterial/unclassified contigs (positive control: 10/16 samples)
- **0 ARG transcripts on archaeal contigs** across all 656 ARG-sample observations
- Triangulating evidence (stringent + relaxed BLASTn + reads-level mapping + manual NCBI BLAST) converges on zero archaeal ARG signal

---

## Repository structure

```
mgnify-archaea-arg-pipeline/
├── README.md
├── scripts/
│   ├── selection/                  # Sample selection from MGnify
│   │   ├── 01_mgnify_search.sh
│   │   ├── 02_pipeline_filter.sh
│   │   ├── 03_ena_verify.sh
│   │   ├── 04_biome_annotation.sh
│   │   └── 05_ecosystem_filter.sh
│   ├── processing/                 # Per-sample data processing
│   │   ├── 06_process_one_run.sh   # download → fastp → SortMeRNA
│   │   └── 07_assembly.sh          # rnaSPAdes
│   └── arg_analysis/               # Contig-level ARG analysis
│       ├── 08_prereq_check.sh
│       ├── 09_kaiju_contigs.sh
│       ├── 10_bowtie2_reads_to_ctgs.sh
│       ├── 11_blastn_arg.sh
│       ├── 12_compute_tpm.sh
│       ├── 13_merge_and_summary.sh
│       └── 14_validate_unclass.sh
├── data/
│   ├── selected_top3_per_category.tsv   # Final 16 samples
│   ├── ARGs_41.fasta                    # Custom ARG reference
│   └── arg_db.lengths.tsv
├── docs/
│   └── pipeline_flowchart.txt
└── .gitignore
```

---

## Pipeline architecture

```
Raw FASTQ (ENA)
    ↓ fastp (v0.23.4) — adapter trim + QC
mRNA + rRNA reads
    ↓ SortMeRNA (v4.3.4) — SILVA + Rfam depletion
mRNA reads
    ↓ rnaSPAdes (v4.0.0)
de novo contigs
    ↓ length filter ≥300 bp
filtered contigs
    ↓
    ├──→ Kaiju (v1.10.1, refseq_ref_2024-08-14, greedy)
    │       └──→ Archaeal contig subset
    │
    ├──→ Bowtie2 (v2.5.2, --very-sensitive) → BAM + idxstats
    │
    └──→ BLASTn (v2.15.0)
            ↓ -evalue 1e-6 -perc_identity 70 -max_target_seqs 41
            ↓ post-filter qcovhsp ≥60%
       ARG-bearing contigs
            ↓ samtools bedcov on ARG HSP coordinates
            ÷ sample-specific mean read length
       per-ARG read counts → RPK → TPM
            ├──→ TPM_Archaeal (denominator: archaeal contigs)
            └──→ TPM_Global (denominator: all contigs)
```

---

## Requirements

### Software (tested versions)

| Tool        | Version           | Purpose                              |
|-------------|-------------------|--------------------------------------|
| fastp       | 0.23.4            | Read QC + adapter trimming           |
| SortMeRNA   | 4.3.4             | rRNA depletion                       |
| rnaSPAdes   | 4.0.0             | De novo transcriptome assembly       |
| Kaiju       | 1.10.1            | Protein-level taxonomic classification |
| Bowtie2     | 2.5.2             | Read alignment to contigs            |
| SAMtools    | 1.17              | BAM processing                       |
| BLAST+      | 2.15.0            | ARG detection (BLASTn)               |
| GNU awk     | ≥4.0              | Data processing                      |

### Reference databases

- **Kaiju RefSeq protein index** (`refseq_ref_2024-08-14`):
  Download from [Kaiju server](https://kaiju.binf.ku.dk/server)
- **Custom 41-ARG nucleotide reference** (`data/ARGs_41.fasta`):
  Compiled by BLASTn search of archaeal genomes against PanRes-derived ARG database at ≥70% identity threshold

### Compute resources (suggested)

- **Per-sample peak RAM**: ~140 GB (Kaiju on full RefSeq index)
- **Per-sample wall time**: 4–10 hours depending on read depth
- **Disk space**: ~500 GB total for 16 samples (intermediate + final outputs)

---

## Quick start

### Step 1 — Sample selection (run once)

```bash
cd scripts/selection
bash 01_mgnify_search.sh
bash 02_pipeline_filter.sh
bash 03_ena_verify.sh
bash 04_biome_annotation.sh
bash 05_ecosystem_filter.sh
# Output: data/selected_top3_per_category.tsv (16 samples)
```

### Step 2 — Per-sample processing (SLURM array, parallel)

```bash
cd scripts/processing
sbatch 06_process_one_run.sh    # FASTQ download + fastp + SortMeRNA
sbatch 07_assembly.sh           # rnaSPAdes assembly
```

### Step 3 — ARG analysis (sequential)

```bash
cd scripts/arg_analysis
bash 08_prereq_check.sh         # Verify all references
sbatch 09_kaiju_contigs.sh      # Kaiju on filtered contigs
sbatch 10_bowtie2_reads_to_ctgs.sh   # Bowtie2 alignment
sbatch 11_blastn_arg.sh         # BLASTn ARG detection
bash 12_compute_tpm.sh          # bedcov + RPK + TPM
bash 13_merge_and_summary.sh    # 16-sample merge + cross-sample stats
bash 14_validate_unclass.sh     # Extract unclassified ARG-bearing contigs
```

---

## Key parameters

### BLASTn (ARG detection)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `-evalue` | `1e-6` | Standard E-value cutoff |
| `-perc_identity` | `70` | Matches threshold used to compile the 41-ARG panel from archaeal genomes |
| `-max_target_seqs` | `41` | Equal to database size; avoids cutoff bias [(Shah et al. 2018)](https://doi.org/10.1093/bioinformatics/bty833) |
| `qcovhsp` (post-filter) | `≥60%` | Excludes partial low-quality alignments |

### Kaiju

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Mode | `greedy` (`-a greedy`) | Higher sensitivity than mem mode |
| Min E-value | `1e-5` (`-E 1e-5`) | Default Kaiju protein matching |
| Max mismatches | `5` (`-e 5`) | Default Kaiju setting |
| Index | `refseq_ref_2024-08-14` | RefSeq representative-genomes |

### Bowtie2

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Preset | `--very-sensitive` | Maximize alignment sensitivity for archaeal mRNA |
| Max insert | `-X 1000` | Standard for paired-end mRNA |
| MAPQ filter | `≥10` | Excludes multi-mapping reads |

### Length filter

- Contigs **≥300 bp** retained (corresponds to ≥100 amino acids when translated, the minimum length for reliable Kaiju protein-level classification)

---

## Output structure

After running the full pipeline, the following per-sample outputs are produced:

```
META=/path/to/working/directory
├── assembly/
│   └── {SAMPLE}/
│       ├── transcripts.fasta              # Original rnaSPAdes output
│       └── transcripts_min300.fasta       # Length-filtered
├── kaiju_out/
│   └── {SAMPLE}_contigs_min300/
│       ├── {SAMPLE}_contigs_kaiju.out
│       ├── {SAMPLE}_contigs_kaiju_named.out
│       └── {SAMPLE}_archaeal_contigs.txt  # Archaeal contig IDs
├── contig_arg_min300/
│   └── {SAMPLE}/
│       ├── {SAMPLE}_reads_vs_ctg.bam
│       ├── {SAMPLE}_ctg_idxstats.tsv
│       ├── {SAMPLE}_arg_bedcov.tsv
│       ├── {SAMPLE}_arg_bedcov_archaeal.tsv
│       └── {SAMPLE}_ARG_TPM_contig.tsv     # Per-sample TPM table (41 rows)
├── blastn_arg_min300/
│   └── {SAMPLE}/
│       ├── {SAMPLE}_blastn_all.tsv
│       └── {SAMPLE}_blastn_filt.tsv        # qcovhsp ≥60%
└── unclassified_arg_contigs_min300/
    ├── all_unclass_arg_contigs.fasta       # For manual NCBI BLAST validation
    └── contig_to_sample.tsv
```

### Final main table

```
contig_arg_min300/ALL_16_samples_ARG_TPM_contig_min300.tsv  (657 rows)
```

Columns: `Sample | Gene_ID | Length | Reads_All | RPK_All | TPM_Global | Reads_Archaeal | RPK_Archaeal | TPM_Archaeal | Mean_read_len`

---

## Sample inventory

| Habitat       | Samples (n) | ENA Run IDs                                        |
|---------------|-------------|----------------------------------------------------|
| Marine        | 3           | ERR1711948, ERR1711870, ERR1711945                 |
| Hydrothermal  | 3           | ERR694391, ERR2021508, ERR694360                   |
| Wastewater    | 3           | DRR066667, DRR066668, ERR2088989                   |
| Freshwater    | 3           | ERR3132363, ERR3132364, ERR3132366                 |
| Human gut     | 1           | SRR3313097                                         |
| Rumen         | 3           | ERR747934, ERR747931, ERR747932                    |

Per-sample metadata (BioProject, platform, depth) in `data/selected_top3_per_category.tsv`.

---

## Tool references

- **fastp**: [Chen et al. 2018](https://doi.org/10.1093/bioinformatics/bty560)
- **SortMeRNA**: [Kopylova et al. 2012](https://doi.org/10.1093/bioinformatics/bts611)
- **rnaSPAdes**: [Bushmanova et al. 2019](https://doi.org/10.1093/gigascience/giz100)
- **Kaiju**: [Menzel et al. 2016](https://doi.org/10.1038/ncomms11257)
- **Bowtie2**: [Langmead & Salzberg 2012](https://doi.org/10.1038/nmeth.1923)
- **SAMtools**: [Danecek et al. 2021](https://doi.org/10.1093/gigascience/giab008)
- **TPM**: [Wagner et al. 2012](https://doi.org/10.1007/s12064-012-0162-3)
- **BLAST max_target_seqs**: [Shah et al. 2018](https://doi.org/10.1093/bioinformatics/bty833)
- **Database bias**: [Pasolli et al. 2019](https://doi.org/10.1016/j.cell.2019.01.001), [Nayfach et al. 2021](https://doi.org/10.1038/s41587-020-0718-6)
