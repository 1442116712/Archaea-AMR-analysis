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
