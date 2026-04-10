library(readxl)
library(openxlsx)
library(dplyr)
df <- read_excel('C:/Users/CFL/Desktop/1.xlsx', sheet = "Sheet2")

df$gene_length <- as.numeric(df$gene_length)
# filter subject coverage < 60%
df$coverage <- (df$length) / df$gene_length 
df <- df[df$coverage >= 0.6, ]
df$coverage <- ifelse(df$coverage > 1, 1, df$coverage)

# determine strand and centre
df$strand <- ifelse(df$send - df$sstart > 0, "F", "R")
df$center_ARG <- df$qstart + (df$qend - df$qstart) / 2

remove_overlapping_genes <- function(df) {
  # Pre-calculate original indices to avoid using the expensive which() function inside loops
  df$orig_idx <- seq_len(nrow(df))
  # Initialize a logical vector to track which genes to keep
  keep <- rep(TRUE, nrow(df))
  
  # Identify all groups (Contig + ARG_class + strand) that have potential overlaps
  groups <- df %>% 
    group_by(Contig, ARG_class, strand) %>% 
    summarise(n = n(), .groups = 'drop') %>% 
    filter(n >= 2) # Only process groups with at least 2 genes
  
  # Iterate through each candidate group
  for (g in seq_len(nrow(groups))) {
    # Extract original row indices for the current group
    sub_indices <- which(df$Contig == groups$Contig[g] & 
                         df$ARG_class == groups$ARG_class[g] & 
                         df$strand == groups$strand[g])
    
    # Create a temporary subset for faster calculation
    sub_df <- df[sub_indices, ]
    n_sub <- nrow(sub_df)
    
    # Pairwise comparison within the group
    for (i in 1:(n_sub - 1)) {
      for (j in (i + 1):n_sub) {
        
        # Map back to the global indices in the original data frame
        idx_i <- sub_indices[i]
        idx_j <- sub_indices[j]
        
        # CRITICAL: If either gene has already been marked for removal, skip this comparison
        if (!keep[idx_i] | !keep[idx_j]) next
        
        # Calculate the physical distance between gene centers
        dist <- abs(sub_df$center_ARG[i] - sub_df$center_ARG[j])
        # Determine the threshold for overlap (allowing a 4bp buffer)
        allowed_dist <- (sub_df$length[i] + sub_df$length[j]) / 2 - 4
        
        # If the actual distance is less than the threshold, an overlap is confirmed
        if (!is.na(dist) && dist < allowed_dist) {
          # Competition: Retain the gene with the higher bitscore
          if (sub_df$bitscore[i] < sub_df$bitscore[j]) {
            keep[idx_i] <- FALSE # Mark gene i for removal
          } else {
            keep[idx_j] <- FALSE # Mark gene j for removal
          }
        }
      }
    }
  }
  
  # Return the filtered data frame and remove the helper index column
  return(df[keep, -which(names(df) == "orig_idx")]) 
}

# Output the filter result
filtered_df <- remove_overlapping_genes(df)

# Print results
write.xlsx(filtered_df, 'C:/Users/CFL/Desktop/UK/PhD/AMR/output_file.xlsx')
