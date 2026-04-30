# Install packages (if not already installed)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "coRdon"))
install.packages("aplot")

# Load libraries
library(Biostrings)
library(coRdon)
library(ggplot2)
library(stringr)
library(ggrepel)

# MILC Analysis ####
# Read host genome
host_genome <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_CUP_analysis/Klebsiella_pneumoniae_core_nogap.fasta")
# Create codon table for host genome
cThost <- codonTable(host_genome)
# Count codons and lengths for host genome
counts_host <- codonCounts(cThost)

# Read target ARG (Antibiotic Resistance Gene) sequences
ARG_seq <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/ARG_MGE_CDS16x3.ffn")
names(ARG_seq)

# Create codon table and check matrix status
cTARG <- codonTable(ARG_seq)  
if (!is.matrix(cTARG@counts)) {
  cTARG@counts <- matrix(cTARG@counts, nrow = 1)
}
class(cTARG)            # "codonTable"
is.matrix(cTARG@counts) # TRUE
dim(cTARG@counts)       # 1 64

# Count codons of ARG
counts_ARG <- codonCounts(cTARG)

# MILC (Measure Independent of Length and Composition)
MILC(cTARG, subsets = list(cThost), self = FALSE, ribosomal = FALSE,
     id_or_name2 = "1", alt.init = TRUE, stop.rm = TRUE,
     filtering = "hard", len.threshold = 80)



# RSCU Analysis ####
# Load required packages (assuming Biostrings and coRdon are used)
library(Biostrings)

# 1. Basic data preparation: Codon-Amino Acid mapping table (placed outside the loop to avoid redundant execution)
std_codons <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
                "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
                "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
                "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
                "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
                "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")

std_aas <- c("K", "N", "K", "N", "T", "T", "T", "T",
             "R", "S", "R", "S", "I", "I", "M", "I",
             "Q", "H", "Q", "H", "P", "P", "P", "P",
             "R", "R", "R", "R", "L", "L", "L", "L",
             "E", "D", "E", "D", "A", "A", "A", "A",
             "G", "G", "G", "G", "V", "V", "V", "V",
             "*", "Y", "*", "Y", "S", "S", "S", "S",
             "*", "C", "W", "C", "L", "F", "L", "F")

std_aas_string <- paste(std_aas, collapse = "")
codon_to_aa <- setNames(strsplit(std_aas_string, "")[[1]], std_codons)

aa_codons <- list()
for (codon in std_codons) {
  aa <- codon_to_aa[codon]
  if (!is.na(aa) && aa != "*") {  # Exclude NA and stop codons
    if (!aa %in% names(aa_codons)) {
      aa_codons[[aa]] <- c()
    }
    aa_codons[[aa]] <- c(aa_codons[[aa]], codon)
  }
}

stop_codons <- c("TAA", "TAG", "TGA")

# 2. Define function to calculate weighted average RSCU
calculate_weighted_avg_rscu <- function(rscu_values, total_codon_number) {
  weighted_avg <- numeric(ncol(rscu_values))
  names(weighted_avg) <- colnames(rscu_values)
  # Weights: Total codon count per gene
  gene_weights <- rowSums(total_codon_number)
  
  for (codon in colnames(rscu_values)) {
    rscu_codon <- rscu_values[, codon]
    # Avoid division by zero if total weight is 0
    if(sum(gene_weights) > 0) {
      weighted_avg[codon] <- sum(rscu_codon * gene_weights) / sum(gene_weights)
    } else {
      weighted_avg[codon] <- 0 
    }
  }
  return(weighted_avg)
}

# 3. Set directory and retrieve all FASTA file paths
input_dir <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_CUP_analysis/"
# Automatically match all files ending in .fasta in the directory
fasta_files <- list.files(path = input_dir, pattern = "\\.fasta$", full.names = TRUE)

# Create an empty list to store the weighted average RSCU results for each file
all_results_list <- list()

# 4. Loop through each FASTA file
for (file_path in fasta_files) {
  # Extract file name (used as row identifier in the final table)
  file_name <- basename(file_path)
  cat("Processing:", file_name, "...\n")
  
  # Read genome
  host_genome <- readDNAStringSet(file_path)
  
  # Generate codon table and calculate counts
  cThost <- codonTable(host_genome)
  counts_host <- codonCounts(cThost)
  
  # Remove stop codons
  counts_final <- counts_host[, !colnames(counts_host) %in% stop_codons, drop = FALSE]
  
  # Add pseudocounts
  pseudocount <- 0.001
  counts_pseudo <- counts_final + pseudocount
  
  # Calculate RSCU matrix for this file
  rscu_mat <- counts_pseudo
  for (aa in names(aa_codons)) {
    codons <- aa_codons[[aa]]
    n_codons <- length(codons)
    # Calculation based on RSCU formula
    rscu_mat[, codons] <- counts_final[, codons] / (rowSums(counts_pseudo[, codons, drop=FALSE]) / n_codons)
  }
  
  # Calculate weighted average RSCU for the file
  weighted_rscu <- calculate_weighted_avg_rscu(rscu_mat, counts_final)
  
  # Save the results of 61 codons to the list, indexed by file name
  all_results_list[[file_name]] <- weighted_rscu
}

# 5. Integrate results and output
# Convert list to data frame (Rows = Fasta files, Columns = 61 codons)
final_rscu_df <- do.call(rbind, all_results_list)

# Round results to three decimal places for readability
final_rscu_df <- round(final_rscu_df, 3)

# Print a summary of the first few files in the console
print("All files processed! Overview of the first few files:")
print(head(final_rscu_df))

# Export to CSV file in the same directory
output_csv <- paste0(input_dir, "All_Genomes_Weighted_RSCU.csv")
write.csv(final_rscu_df, file = output_csv, row.names = TRUE)
cat("\nWeighted average RSCU results exported successfully to:\n", output_csv, "\n")


# Calculate RSCU for each sequence in a specific file - For ARG ####
# 1. Set file path and read sequences
ribo_file <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/ARG_MGE_CDS16x3.ffn"
ribo_seq <- readDNAStringSet(ribo_file)

# 2. Generate codon table and check matrix status
cTribo <- codonTable(ribo_seq)
if (!is.matrix(cTribo@counts)) {
  cTribo@counts <- matrix(cTribo@counts, nrow = 1)
}

# 3. Extract codon counts
counts_ribo <- codonCounts(cTribo)

# 4. Remove stop codons
counts_ribo_final <- counts_ribo[, !colnames(counts_ribo) %in% stop_codons, drop = FALSE]

# 5. Add pseudocounts to avoid division by zero
pseudocount <- 0.001
counts_ribo_pseudo <- counts_ribo_final + pseudocount

# 6. Calculate RSCU matrix for each ribosome protein sequence
rscu_ribo <- counts_ribo_pseudo
for (aa in names(aa_codons)) {
  codons <- aa_codons[[aa]]
  n_codons <- length(codons)
  rscu_ribo[, codons] <- counts_ribo_final[, codons] / (rowSums(counts_ribo_pseudo[, codons, drop=FALSE]) / n_codons)
}

# 7. Organize row names and round to 3 decimal places
rownames(rscu_ribo) <- names(ribo_seq)
rscu_ribo_df <- as.data.frame(round(rscu_ribo, 3))

# 8. Preview data in the console
cat("\nIndividual sequence RSCU calculation for ribosome proteins complete! Preview:\n")
print(head(rscu_ribo_df[, 1:5])) 

# 9. Export results to CSV
output_ribo_csv <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/ARG_RSCU.csv"
write.csv(rscu_ribo_df, file = output_ribo_csv, row.names = TRUE)

cat("\nIndividual sequence RSCU results saved to:\n", output_ribo_csv, "\n")


# Codon Usage Bias Visualization Plot - CORE GENOME #### 
library(ggplot2)
library(aplot)

# 1. Basic settings and read external annotation table
out_dir <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/RSCU/"
cols <- c("#BDB76B", "#CAB2D6", "#FF69B4", "#B2DF8A", "#FF9D1E", "#6495ED")

aa_table_path <- paste0(out_dir, "AA_table.txt")
df_base <- read.table(aa_table_path, header = F, stringsAsFactors = F)

# 2. Begin Loop
for (i in 1:nrow(final_rscu_df)) {
  
  # Extract original file name
  raw_name <- rownames(final_rscu_df)[i]
  
  # --- Enhanced Name Cleaning ---
  clean_name <- raw_name
  # Remove common file extensions
  clean_name <- gsub("_CDS_merged\\.fasta|\\.fasta", "", clean_name) 
  # Remove 'core' and 'nogap' (case-insensitive)
  clean_name <- gsub("_core|core|_nogap|nogap", "", clean_name, ignore.case = TRUE)
  # Remove redundant underscores or leading/trailing underscores
  clean_name <- gsub("_{2,}", "_", clean_name)
  clean_name <- sub("^_|_$", "", clean_name) 
  
  # Replace underscores with spaces for chart titles
  title_name <- gsub("_", " ", clean_name)
  
  cat("Plotting:", title_name, "...\n")
  
  # Extract RSCU values for the 61 codons of the current row
  current_rscu <- as.numeric(final_rscu_df[i, ])
  names(current_rscu) <- colnames(final_rscu_df)
  
  # 3. Integrate plotting data
  rscu_df <- data.frame(
    V1 = names(current_rscu), 
    V4 = round(current_rscu, 3)
  )
  
  df <- merge(df_base, rscu_df, by = "V1", all.x = TRUE)
  
  # 4. Draw main bar chart (p1)
  p1 <- ggplot(df, aes(x = V2, y = V4, fill = as.character(V3))) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = c(cols)) +
    geom_text(aes(label = ifelse(V4 > 0.1, paste0(V4, ""), "")),
              position = position_stack(vjust = 0.5), size = 2.5, color = "black", fontface = "plain") +
    labs(x = NULL, y = "RSCU", 
         title = paste0("RSCU of ", title_name)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0,6.2)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),  
          axis.title.y = element_text(size = 15, face = "bold"),  
          plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  # 5. Draw bottom label plot (p2)
  p2 <- ggplot(df, aes(x = V2, y = V3))+
    geom_label(aes(label=V1, fill=as.character(V3)),size=2)+
    scale_fill_manual(values = c(cols)) +
    labs(x="",y="") + ylim(3.4,6.3)+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  
  # 6. Combine and save using aplot
  p <- insert_bottom(p1, p2, height = 0.2) 
  
  save_path <- file.path(out_dir, paste0(clean_name, ".png"))
  
  ggsave(
    filename = save_path,     
    plot = p,      
    width = 8,         
    height = 8,         
    dpi = 300,         
    units = "in",
    bg = "white" 
  )
}


# Codon Usage Bias Visualization Plot - ARG/MEG-ARG #### 
library(ggplot2)
library(aplot)

# 1. Basic settings
out_dir <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/RSCU/"
cols <- c("#BDB76B", "#CAB2D6", "#FF69B4", "#B2DF8A", "#FF9D1E", "#6495ED")

aa_table_path <- paste0(out_dir, "AA_table.txt")
df_base <- read.table(aa_table_path, header = F, stringsAsFactors = F)

cat("\nStarting to draw detailed RSCU plots for each individual sequence...\n")

# 2. Loop through each sequence in rscu_ribo_df
for (i in 1:nrow(rscu_ribo_df)) {
  
  # Extract original sequence name
  raw_name <- rownames(rscu_ribo_df)[i]
  
  # === File name and title cleaning ===
  clean_file_name <- gsub("[<>:/\\\\|?*]", "_", raw_name)
  title_name <- gsub("_", " ", raw_name)
  
  cat(sprintf("Plotting sequence %d/%d: %s ...\n", i, nrow(rscu_ribo_df), raw_name))
  
  # Extract RSCU values
  current_rscu <- as.numeric(rscu_ribo_df[i, ])
  names(current_rscu) <- colnames(rscu_ribo_df)
  
  # 3. Integrate data
  temp_rscu_df <- data.frame(
    V1 = names(current_rscu), 
    V4 = round(current_rscu, 3) 
  )
  
  df <- merge(df_base, temp_rscu_df, by = "V1", all.x = TRUE)
  
  # 4. Main bar chart (p1)
  p1 <- ggplot(df, aes(x = V2, y = V4, fill = as.character(V3))) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = c(cols)) +
    geom_text(aes(label = ifelse(V4 > 0.1, paste0(V4, ""), "")),
              position = position_stack(vjust = 0.5), size = 2.5, color = "black", fontface = "plain") +
    labs(x = NULL, y = "RSCU", 
         title = paste0("RSCU of ", title_name)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0,6.2)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),  
          axis.title.y = element_text(size = 15, face = "bold"),  
          plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  # 5. Bottom label plot (p2)
  p2 <- ggplot(df, aes(x = V2, y = V3))+
    geom_label(aes(label=V1, fill=as.character(V3)),size=2)+
    scale_fill_manual(values = c(cols)) +
    labs(x="",y="") + ylim(3.4,6.3)+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  
  # 6. Combine and save
  p <- insert_bottom(p1, p2, height = 0.2) 
  
  save_path <- file.path(out_dir, paste0(clean_file_name, ".png"))
  
  ggsave(
    filename = save_path,     
    plot = p,      
    width = 8,         
    height = 8,         
    dpi = 300,         
    units = "in",
    bg = "white" 
  )
}


# Core genome validation ####
library(Biostrings)
library(coRdon)

# ==========================================
# 1. Parameter Settings
# ==========================================
input_dir <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_CUP_analysis/"
target_file <- paste0(input_dir, "Escherichia_coli_core_nogap.fasta")

# Set the length of the head and tail segments you wish to extract
clip_length <- 1800 
# ==========================================


# 2. Preparation Phase
stop_codons <- c("TAA", "TAG", "TGA")

# 3. Read specific file and merge into a long sequence
cat("Reading large file:", target_file, "...\n")
raw_seqs <- readDNAStringSet(target_file)
full_genome_seq <- unlist(raw_seqs) 
total_len <- length(full_genome_seq)

# Safety check
if (total_len < clip_length * 2) {
  warning("Warning: Total sequence length is insufficient to extract non-overlapping head and tail segments!")
}

# 4. Dynamically extract head and tail fragments
cat("Total length:", total_len, "bp. Extracting head and tail fragments of", clip_length, "bp each...\n")

# First clip_length bases
head_seq <- subseq(full_genome_seq, start = 1, end = clip_length)
# Last clip_length bases
tail_seq <- subseq(full_genome_seq, start = total_len - (clip_length - 1), end = total_len)

# 5. Core RSCU calculation function
get_weighted_rscu <- function(dna_set) {
  cT <- codonTable(dna_set)
  counts <- codonCounts(cT)
  
  counts_final <- counts[, !colnames(counts) %in% stop_codons, drop = FALSE]
  
  pseudocount <- 0.001
  counts_pseudo <- counts_final + pseudocount
  
  rscu_mat <- counts_pseudo
  for (aa in names(aa_codons)) {
    codons <- aa_codons[[aa]]
    codons <- codons[!codons %in% stop_codons]
    if(length(codons) == 0) next
    
    n_codons <- length(codons)
    rscu_mat[, codons] <- counts_final[, codons] / (rowSums(counts_pseudo[, codons, drop=FALSE]) / n_codons)
  }
  
  res <- calculate_weighted_avg_rscu(rscu_mat, counts_final)
  return(res)
}

# 6. Execute calculation and dynamic naming
results_segments <- list()

head_name <- paste0("Head_", clip_length)
tail_name <- paste0("Tail_", clip_length)

results_segments[[head_name]] <- get_weighted_rscu(DNAStringSet(head_seq))
results_segments[[tail_name]] <- get_weighted_rscu(DNAStringSet(tail_seq))

# 7. Integrate and output results
final_segment_df <- do.call(rbind, results_segments)
final_segment_df <- round(final_segment_df, 3)

cat("\nE. coli core genome head/tail fragment (", clip_length, "bp) RSCU calculation complete:\n")
print(final_segment_df)

output_segment_csv <- paste0(input_dir, "Ecoli_Head_Tail_", clip_length, "bp_RSCU.csv")
write.csv(final_segment_df, file = output_segment_csv)
cat("\nSegment results saved to:", output_segment_csv, "\n")


# catA species 3D RSCU PCA ####
library(Biostrings)
library(coRdon)
library(scatterplot3d)
library(plot3D)

# ============================================================
# 0. Basic data preparation
# ============================================================
# (Variables defined previously in the script)

# General RSCU calculation function (returns RSCU matrix for each sequence)
calc_rscu_matrix <- function(seqs, pseudocount = 0.001) {
  cT <- codonTable(seqs)
  if (!is.matrix(cT@counts)) cT@counts <- matrix(cT@counts, nrow = 1)
  
  counts_all   <- codonCounts(cT)
  counts_final <- counts_all[, !colnames(counts_all) %in% stop_codons, drop = FALSE]
  counts_pseudo <- counts_final + pseudocount
  
  rscu_mat <- counts_pseudo
  for (aa in names(aa_codons)) {
    codons   <- aa_codons[[aa]]
    n_codons <- length(codons)
    rscu_mat[, codons] <- counts_final[, codons] / 
      (rowSums(counts_pseudo[, codons, drop = FALSE]) / n_codons)
  }
  list(rscu = rscu_mat, counts_final = counts_final)
}

# ============================================================
# 1. Set directory and automatically identify ARG / Host files
# ============================================================
input_dir <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_CUP_analysis/catA/"
species_name <- basename(normalizePath(input_dir))   # "catA"

all_fastas <- list.files(input_dir, pattern = "\\.(fasta|ffn|fna|fa)$", 
                         full.names = TRUE, ignore.case = TRUE)

# Identify ARG files: Filename contains "ARG" (case-insensitive)
arg_idx_file  <- grep("ARG", basename(all_fastas), ignore.case = TRUE)
arg_files     <- all_fastas[arg_idx_file]
host_files    <- all_fastas[-arg_idx_file]

if (length(arg_files) == 0) stop("No ARG.fasta file found in directory!")
if (length(host_files) == 0) stop("No host fasta files found in directory!")

cat("Found ARG file:", basename(arg_files), "\n")
cat("Found", length(host_files), "host files\n\n")

# ============================================================
# 2. Calculate weighted average RSCU for each host file
# ============================================================
host_results_list <- list()
for (file_path in host_files) {
  file_name <- basename(file_path)
  cat("Processing host:", file_name, "\n")
  
  host_seq <- readDNAStringSet(file_path)
  res      <- calc_rscu_matrix(host_seq)
  weighted <- calculate_weighted_avg_rscu(res$rscu, res$counts_final)
  
  host_results_list[[file_name]] <- weighted
}
final_rscu_df <- do.call(rbind, host_results_list)

# ============================================================
# 3. Calculate RSCU for each ARG sequence
# ============================================================
cat("\nProcessing ARG:", basename(arg_files), "\n")
ARG_seq  <- readDNAStringSet(arg_files[1])
res_ARG  <- calc_rscu_matrix(ARG_seq)
rscu_ARG <- res_ARG$rscu
rownames(rscu_ARG) <- names(ARG_seq)


# ============================================================
# 4. Merge matrices and perform PCA
# ============================================================
rscu_ARG <- rscu_ARG[, colnames(final_rscu_df), drop = FALSE]
combined_rscu_mat <- rbind(final_rscu_df, rscu_ARG)

host_groups <- gsub("_CDS_merged\\.fasta|\\.fasta|\\.ffn|\\.fna|\\.fa", "",
                    rownames(final_rscu_df), ignore.case = TRUE)
arg_groups  <- rep("ARG/MGE-ARG", nrow(rscu_ARG))
group       <- c(host_groups, arg_groups)
arg_labels  <- names(ARG_seq)

pca_res <- prcomp(combined_rscu_mat, scale. = TRUE)
pca_df  <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3],
  group = group,
  stringsAsFactors = FALSE
)

pc1_var <- round(summary(pca_res)$importance[2, "PC1"] * 100, 1)
pc2_var <- round(summary(pca_res)$importance[2, "PC2"] * 100, 1)
pc3_var <- round(summary(pca_res)$importance[2, "PC3"] * 100, 1)

# ============================================================
# 5. Color and Plotting
# ============================================================
unique_groups <- c("ARG/MGE-ARG", sort(unique(pca_df$group[pca_df$group != "ARG/MGE-ARG"])))

colours <- c(
  "#D62728", "#1F77B4", "#2CA02C", "#FF7F0E", 
  "#9467BD", "#E377C2", "#17BECF", "#BCBD22", 
  "#8C564B", "#FF9896", "#98DF8A", "#C49C94", 
  "#C5B0D5", "#F7B6D2", "#DBDB8D", "#393B79",
  "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#5254A3", "#9C9EDE", "#CEDB9C", "#E7BA52",
  "#AD494A", "#A55194", "#606060"
)
group_colour_map <- setNames(colours[seq_along(unique_groups)], unique_groups)

pca_df$colour <- group_colour_map[pca_df$group]
pca_df$pch    <- ifelse(pca_df$group == "ARG/MGE-ARG", 17, 16)

# Output Path
output_dir  <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/RSCU/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_path <- paste0(output_dir, species_name, "_3D_PCA.png")

png(filename = output_path, width = 12, height = 7.5, 
    units = "in", res = 600, bg = "white")

par(mar = c(2, 2, 3, 8), xpd = TRUE)

scatter3D(
  x = pca_df$PC1, y = pca_df$PC2, z = pca_df$PC3,
  colvar    = NULL,            # Turn off default continuous color mapping
  col       = pca_df$colour,
  pch       = 21,              # Filled circle with border
  bg        = pca_df$colour,   # Sphere fill color
  cex       = 1.8,
  
  # —— Parameters for 3D depth —— 
  bty       = "b2",            # Shaded box for depth perception
  theta     = 40,              # Azimuthal rotation
  phi       = 20,              # Elevation angle
  d         = 2,               # Perspective intensity
  expand    = 0.8,             # Z-axis stretch
  
  ticktype  = "detailed",      # Detailed scale
  nticks    = 5,
  
  # —— Projection lines to enhance depth ——
  type      = "h",             # Vertical lines to the XY plane
  lwd       = 0.6,
  
  xlab      = paste0("PC1 (", pc1_var, "%)"),
  ylab      = paste0("PC2 (", pc2_var, "%)"),
  zlab      = paste0("PC3 (", pc3_var, "%)"),
  main      = paste0("RSCU 3D PCA - ", species_name),
  
  colkey    = FALSE            # No default color key
)

# ARG text labels
arg_idx <- which(pca_df$group == "ARG/MGE-ARG")

# Calculate offset for labels
x_offset <- diff(range(pca_df$PC1)) * 0.015

text3D(
  x = pca_df$PC1[arg_idx] + x_offset, 
  y = pca_df$PC2[arg_idx], 
  z = pca_df$PC3[arg_idx],
  labels = arg_labels,
  add    = TRUE,
  cex    = 0.75, 
  font   = 2,
  col    = "black",
  adj    = 0                   # Left alignment
)

# Legend (placed outside the plot area)
legend("right",
       inset     = c(-0.15, 0),
       legend    = unique_groups,
       col       = "black",
       pt.bg     = group_colour_map[unique_groups],
       pch       = 21,
       bty       = "n",
       cex       = 0.85,
       pt.cex    = 1.5,
       title     = "Group",
       title.adj = 0)

dev.off()

cat("\n3D PCA plot exported to:\n", output_path, "\n")

# ============================================================
# 6. Additional generation of 2D PCA (PC1 vs PC2)
# ============================================================
library(ggplot2)
library(ggrepel)

output_path_2d <- paste0(output_dir, species_name, "_2D_PCA.png")

# Adjust legend order: ARG/MGE-ARG first, others alphabetical
pca_df$group <- factor(pca_df$group, levels = unique_groups)

# Label vector: Only ARG points have labels
pca_df$label <- ifelse(pca_df$group == "ARG/MGE-ARG", 
                       c(rep("", nrow(final_rscu_df)), arg_labels), 
                       "")

PCA_2D <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(
    aes(shape = group), 
    size = 3.5, 
    alpha = 0.85,
    stroke = 0.5
  ) +
  geom_text_repel(
    data        = subset(pca_df, group == "ARG/MGE-ARG"),
    aes(label = label),
    size        = 3.2,
    fontface    = "bold",
    color       = "black",
    max.overlaps = Inf,
    box.padding  = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size  = 0.3
  ) +
  scale_color_manual(values = group_colour_map) +
  scale_shape_manual(values = setNames(
    ifelse(unique_groups == "ARG/MGE-ARG", 17, 16),
    unique_groups
  )) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  ggtitle(paste0("RSCU 2D PCA - ", species_name)) +
  labs(color = "Group", shape = "Group") +
  theme_bw() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.key.size  = unit(0.6, "cm"),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(PCA_2D)

ggsave(
  filename = output_path_2d,
  plot     = PCA_2D,
  width    = 9.5,
  height   = 5.5,
  dpi      = 600,
  units    = "in",
  bg       = "white"
)

cat("\n2D PCA plot exported to:\n", output_path_2d, "\n")

# ============================================================
# 7. Calculate 3D variance-weighted distance between ARGs and Hosts
# ============================================================

# Split PCA coordinates
host_idx <- which(pca_df$group != "ARG/MGE-ARG")
arg_idx  <- which(pca_df$group == "ARG/MGE-ARG")

host_coords <- pca_df[host_idx, c("PC1", "PC2", "PC3")]
arg_coords  <- pca_df[arg_idx,  c("PC1", "PC2", "PC3")]
rownames(host_coords) <- pca_df$group[host_idx]
rownames(arg_coords)  <- arg_labels

# Use variance proportions as weights for the three axes
var_props <- summary(pca_res)$importance[2, ]
w1 <- var_props["PC1"]; w2 <- var_props["PC2"]; w3 <- var_props["PC3"]

# Distance matrix: Rows=ARG, Columns=Host
n_arg  <- nrow(arg_coords)
n_host <- nrow(host_coords)
dist_mat <- matrix(NA, n_arg, n_host, 
                    dimnames = list(arg_labels, rownames(host_coords)))

for (i in seq_len(n_arg)) {
  for (j in seq_len(n_host)) {
    dx <- arg_coords$PC1[i] - host_coords$PC1[j]
    dy <- arg_coords$PC2[i] - host_coords$PC2[j]
    dz <- arg_coords$PC3[i] - host_coords$PC3[j]
    dist_mat[i, j] <- sqrt(w1*dx^2 + w2*dy^2 + w3*dz^2)
  }
}

# ---- Nearest host for each ARG ----
nearest_per_ARG <- data.frame(
  ARG          = arg_labels,
  Nearest_Host = colnames(dist_mat)[apply(dist_mat, 1, which.min)],
  Distance     = round(apply(dist_mat, 1, min), 3),
  stringsAsFactors = FALSE
)

cat("\n========== Nearest Host for each ARG (3D variance-weighted distance) ==========\n")
print(nearest_per_ARG, row.names = FALSE)

# ---- Top-3 nearest hosts for each ARG ----
top3_combined <- do.call(rbind, lapply(seq_len(n_arg), function(i) {
  ranks <- order(dist_mat[i, ])[1:min(3, n_host)]
  data.frame(
    ARG      = arg_labels[i],
    Rank     = seq_along(ranks),
    Host     = colnames(dist_mat)[ranks],
    Distance = round(dist_mat[i, ranks], 3),
    stringsAsFactors = FALSE
  )
}))

cat("\n========== Top-3 Nearest Hosts for each ARG ==========\n")
print(top3_combined, row.names = FALSE)


# ============================================================
# 8. Draw "ARG -> Nearest Host" connection lines on 2D PCA
# ============================================================
nearest_host_xy <- host_coords[nearest_per_ARG$Nearest_Host, c("PC1", "PC2")]

segment_df <- data.frame(
  x_arg  = arg_coords$PC1,
  y_arg  = arg_coords$PC2,
  x_host = nearest_host_xy$PC1,
  y_host = nearest_host_xy$PC2,
  ARG    = arg_labels
)

PCA_2D_with_lines <- PCA_2D +
  geom_segment(
    data = segment_df,
    aes(x = x_arg, y = y_arg, xend = x_host, yend = y_host),
    inherit.aes = FALSE,
    linetype    = "dashed",
    color       = "grey40",
    alpha       = 0.7,
    arrow       = arrow(length = unit(0.15, "cm"), ends = "first")
  )

ggsave(
  filename = paste0(output_dir, species_name, "_2D_PCA_with_nearest.png"),
  plot     = PCA_2D_with_lines,
  width    = 9.5, height = 5.5, dpi = 600, units = "in", bg = "white"
)

cat("\n2D PCA with nearest host connections exported.\n")
