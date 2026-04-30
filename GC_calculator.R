library(Biostrings)
# Using Biostrings alone is sufficient and more efficient than using seqinr for this workflow

# Read sequences (DNAStringSet objects store entries separately, preventing memory overflow issues)
fasta_path <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_CUP_analysis/Klebsiella_pneumoniae_core_nogap.fasta"
sequences <- readDNAStringSet(fasta_path)

# ---------------------------------------------------------
# Calculate the overall weighted average GC3 content
# ---------------------------------------------------------
# Iterate through each sequence to count the third-position base statistics
gc3_info <- lapply(sequences, function(seq) {
  l <- length(seq)
  
  # Exclude sequence anomalies that are shorter than a full codon
  if (l < 3) {
    return(c(gc3_count = 0, pos3_length = 0))
  }
  
  # Get indices for the third positions (3, 6, 9...)
  third_pos <- seq(3, l, by = 3)
  
  # Extract third-position bases and count the frequency of G and C
  sub_seq <- seq[third_pos]
  gc_num <- sum(letterFrequency(sub_seq, "GC"))
  
  # Return the GC3 count and the total number of third positions in this sequence
  return(c(gc3_count = gc_num, pos3_length = length(third_pos)))
})

# Merge the statistics of all sequences into a matrix
gc3_matrix <- do.call(rbind, gc3_info)

# Sum the GC3 counts and total third-position base counts for all sequences
total_gc3_count <- sum(gc3_matrix[, "gc3_count"])
total_pos3_length <- sum(gc3_matrix[, "pos3_length"])

# Calculate weighted average GC3 (Total GC3 count / Total third-position base count)
genome_gc3 <- (total_gc3_count / total_pos3_length) * 100

# Print final result
cat("total GC3 content:", round(genome_gc3, 2), "%\n")
