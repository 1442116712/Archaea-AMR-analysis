# Archaea-AMR-analysis

Reproducible analysis code for **"A domain-wide survey of archaeal antimicrobial resistance genes reveals an ecologically structured and compositionally heterogeneous resistome"** (Wu Z. *et al.*, 2026).

This repository contains all R scripts used to detect, quantify, contextualise and visualise antimicrobial resistance genes (ARGs) and mobile genetic elements (MGEs) across 12,477 archaeal genomes from GTDB, together with the codon-usage, metatranscriptomic and PCR/MIC-validation analyses reported in the manuscript.

The scripts are grouped by pipeline stage and produce the input data underlying every figure and Supplementary Table in the paper.

---

## Citation

If you use any of the code or intermediate outputs from this repository, please cite:

> Wu Z., Godoy-Santos F., Sabino Y.N.V., Huws S.A., Oyama L.B. (2026). *A domain-wide survey of archaeal antimicrobial resistance genes reveals an ecologically structured and compositionally heterogeneous resistome*. **[Journal]** (accepted). DOI: [to be added on publication].

---

## Repository structure

```
Archaea-AMR-analysis/
├── Genomics/                              # per-genome ARG + MGE + codon-usage analyses
│   ├── De-duplication_pipeline_ARG.R      # collapse overlapping ARG hits per contig/strand (Supp Fig 6B)
│   ├── De-duplication_pipeline_MGE.R      # collapse overlapping MGE-associated protein hits
│   ├── GC_calculator.R                    # per-genome GC and GC3 content
│   └── codon analysis.R                   # RSCU / MILC / GC3 for prioritised ARGs and comparators
│
├── Metatranscriptomics/                   # metatranscriptomic ARG expression analyses
│   └── (contig taxonomic assignment, ARG-BLASTn, TPM quantification wrappers)
│
# --- Plotting scripts (ordered by figure) ---
├── Archaea_heatmap.R                      # Figure 1D–E ecological/lineage ARG heatmap
├── ARG_sankey_diagram.R                   # Figure 1G & Figure 2A alluvial views of ARG distribution
├── stackbar_mobileOG-db.R                 # Figure 3A MGE-category stacked bars
├── MGE_ARG_co-localization.R              # Figure 3B ARG–MGE co-localisation matrices
├── MGE_genomic_structure_annotaion.R      # Figure 3B/C gene-arrangement plots (gggenes)
├── stackbar_homologouos_ARG.R             # Figure 4A ARG homologue stacked bars (host-associated vs environmental)
├── Metadata_network_homologous_MGEs.R     # Figure 4B ecological-metadata network of homologues
├── MGE-ARG_worldmap.R                     # Figure 4D–G geographic distribution of homologues
├── GC_MILC.R                              # Figure 6H MILC + |ΔGC3| bubble plot
└── README.md
```

---

## Software requirements

### R environment
- **R** ≥ 4.5.0
- **RStudio** (optional but recommended)

### R packages
Install with:
```r
install.packages(c(
  "tidyverse", "ggplot2", "ggrepel", "ggalluvial", "ggraph", "aplot",
  "patchwork", "stringr", "dplyr", "tidyr", "readxl", "openxlsx",
  "seqinr", "maps", "mapdata", "rentrez", "plot3D", "gggenes"
))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "coRdon", "ComplexHeatmap"))
```

### Command-line tools (called externally by the pipeline)
Used to generate the intermediate inputs on which the R scripts operate:

| Tool | Version | Purpose |
|------|---------|---------|
| BLASTn / BLAST+ | v2.15.0 | ARG and MGE homology searches |
| Prokka | v1.14.6 | Genome annotation |
| Roary | v3.13.0 | Pan-genome / core-gene reconstruction |
| CheckM2 | — (as reported in GTDB release 220) | Genome completeness/contamination QC |
| PHASTEST | against PHAST-BSD bacterial database (release 22 Dec 2020) | Prophage annotation |
| ISfinder BLAST | against ISfinder database (release Oct 2020) | IS element annotation |
| ICEberg | v3.0 | Integrative and conjugative element annotation |
| mobileOG-db (reference pipeline) | v2.0 (internal Prodigal v2.6.3 + DIAMOND blastp) | MGE-associated protein detection |
| PanRes | v1.0.2 | ARG reference database |
| fastp | v0.23.4 | Metatranscriptomic read QC |
| SortMeRNA | v4.2.6 | rRNA depletion |
| rnaSPAdes | v4.0.0 | Metatranscriptome assembly |
| Kaiju | v1.10.1 | Contig protein-level taxonomy |
| Bowtie2 | v2.5.2 | Read mapping for TPM quantification |
| samtools | v1.19 | BAM handling |
| MAFFT | v7.525 | Multiple sequence alignment |
| TrimAl | v1.5.0 | Alignment trimming |
| RAxML | v8.2.12 | Maximum-likelihood phylogeny |

---

## Input data (all publicly available)

| Source | Access |
|---|---|
| **GTDB release 220** (12,477 archaeal genomes) | https://data.gtdb.aau.ecogenomic.org/releases/ |
| **NCBI prokaryotic nucleotide DB v6** (bacterial homologues) | https://ftp.ncbi.nlm.nih.gov/blast/db |
| **PanRes v1.0.2** (ARG reference DB) | https://doi.org/10.5281/zenodo.8055115 |
| **MGnify** (metatranscriptomic runs, n = 16) | https://www.ebi.ac.uk/metagenomics/ (accessions in Supplementary Table 20) |
| **NCBI SRA** (pure-culture RNA-seq) | Accessions in Supplementary Table 20 |
| **mobileOG-db, ISfinder, PHAST-BSD, ICEberg 3.0** | See tool references in Methods |

Intermediate data tables (Supplementary Tables 1–22) accompany the manuscript and can be used to jump into the pipeline at any stage.

---

## Reproduction workflow

The pipeline follows the manuscript flow. **Each block lists the code, its inputs, and the manuscript output(s) it produces.**

### 1. Genome retrieval and quality control
- **Input**: GTDB release 220 metadata; download archaeal genome FASTAs.
- **Output**: Supplementary Table 1 (12,477 rows × metadata).
- Filtering criteria (GTDB default): CheckM2 completeness ≥ 50 %, contamination ≤ 10 %.

### 2. ARG detection (BLASTn against PanRes)
- **CLI**: BLASTn each genome against PanRes v1.0.2, retain hits ≥ 70 % identity and ≥ 70 % coverage.
- **Script**: `Genomics/De-duplication_pipeline_ARG.R`
  - Input: raw BLASTn hits (per genome / per contig).
  - Method: hits assigned to the same ARG class on the same contig and strand are considered overlapping when the distance between their midpoints < ½ (sum of their lengths) + 3 nt (biological overlap tolerance for start/stop codon overlaps).
  - Output: de-duplicated ARG hit table used as input for Supplementary Tables 2–4.
- **Figures**: contributes to Figure 1E–G; the alluvial views (Figure 1G, Figure 2A) are built with `ARG_sankey_diagram.R`.
- **Heatmap**: `Archaea_heatmap.R` → Figure 1D–E (class × phylum ARG heatmap).

### 3. MGE-associated protein detection
- **Pipeline**: all genomes (nucleotide FASTA) were processed through the **mobileOG-db v2.0 reference pipeline** (Brown *et al.*, 2022; https://github.com/clb21565/mobileOG-db). The pipeline internally predicts ORFs with **Prodigal v2.6.3** and queries the translated proteins against the mobileOG-db amino-acid reference database with **DIAMOND blastp**. Retention filters: ≥ 70 % amino-acid identity, ≥ 90 % coverage of the mobileOG-db reference (subject) sequence, e-value ≤ 1 × 10⁻²⁰.
- **Downstream validation** of specific MGE candidates (Supplementary Tables 8–10):
  - Plasmid RefSeq (manual inspection)
  - **PHASTEST** against the PHAST-BSD bacterial database (release 22 Dec 2020) — prophage
  - **ISfinder** BLAST (release Oct 2020) — insertion sequences
  - **ICEberg v3.0** — integrative and conjugative elements
- **De-duplication**: `Genomics/De-duplication_pipeline_MGE.R` (same overlap logic as ARG dedup; used to collapse redundant mobileOG hits per contig/strand).
- **Figures**:
  - `stackbar_mobileOG-db.R` → Figure 3A MGE-category composition (per class, per genome category).
  - `MGE_ARG_co-localization.R` → Figure 3B ARG–MGE co-localisation on the same contig (uses `gggenes`).
  - `MGE_genomic_structure_annotaion.R` → Figure 3B/C per-contig gene-arrangement plots.

### 4. Homology search of prioritised ARGs and ARG-MGE regions
- **CLI**: BLASTn of 16 prioritised ARGs and 3 linked ARG-MGE regions against the NCBI prokaryotic nucleotide database v6. Retention filters: ≥ 90 % nucleotide identity, ≥ 70 % subject coverage, e-value ≤ 1 × 10⁻⁶.
- **Post-processing scripts**:
  - `Metadata_network_homologous_MGEs.R` → Figure 4B ecological-metadata network.
  - `stackbar_homologouos_ARG.R` → Figure 4A source-environment composition of bacterial homologues.
  - `MGE-ARG_worldmap.R` → Figure 4D–G global geographic map of homologue distributions.

### 5. Phylogenetic reconstruction
- **CLI**: MAFFT alignment → TrimAl → RAxML maximum-likelihood trees (Figure 5).
- No R helper needed; trees rendered externally in iTOL.

### 6. Codon-usage compositional analysis
- **Prerequisite**: Prokka + Roary per-species pan-genomes for 5 archaeal hosts + 21 bacterial comparators (metadata in Supplementary Tables 14–15).
- **Scripts**:
  - `Genomics/GC_calculator.R` → per-sequence and per-genome GC / GC3.
  - `Genomics/codon analysis.R` → RSCU, weighted RSCU, MILC, per-locus PCA, positive-control fragment calibration (Ma > 0.8, |ΔGC3| > 6 percentage points).
- **Global rank-based robustness test** (Supplementary Table 19): weighted RSCU Euclidean distance between each of 12 prioritised targets and 26 core-genome references; ranks + gap statistic + z-score.
- **Figure 6**:
  - Panels A–F: per-locus PCA (from `codon analysis.R`).
  - Panel G: amino-acid-level Ma heatmap (Supplementary Table 17).
  - Panel H: MILC + |ΔGC3| bubble plot from `GC_MILC.R`.
  - Supplementary Figure 4.x: nucleotide-level weighted RSCU stacked bar plots.

### 7. Metatranscriptomic ARG expression
- **Directory**: `Metatranscriptomics/`
- **CLI pipeline**: fastp (QC) → SortMeRNA (rRNA depletion) → rnaSPAdes (assembly) → Kaiju (contig taxonomy) → BLASTn (ARG detection on contigs) → Bowtie2 + samtools (per-ARG read mapping) → TPM quantification.
- **Output**: Figure 7A (rRNA vs. mRNA-contig archaeal fraction), Figure 7B (ARG TPM), Supplementary Table 20.

### 8. Pure-culture transcriptomic validation
- Bowtie2 mapping of SRA RNA-seq to a strain-specific ARG reference set (5 archaeal genomes).
- TPM quantification; results reported in Figure 7B (right panel).

### 9. PCR and Sanger validation of selected ARGs
- Primers listed in Supplementary Table 21.
- Sanger alignments (CLUSTAL Omega) in Supplementary Notes.

### 10. Antimicrobial susceptibility testing
- Broth microdilution MIC assay of *Halobacterium salinarum* NRC-1 against 6 antibiotics; adapted EUCAST methodology.
- Results in Supplementary Table 22.

---

## Data + code archive

- **GitHub**: https://github.com/1442116712/Archaea-AMR-analysis

---

## Contact

- **First author / code maintainer**: Ziming Wu — zwu18@qub.ac.uk
- **Corresponding senior author**: Dr Linda B. Oyama — l.oyama@qub.ac.uk

School of Biological Sciences, Institute for Global Food Security
Queen's University Belfast, Belfast BT9 5DL, United Kingdom

---
