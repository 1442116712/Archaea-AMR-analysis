library(readxl)
library(openxlsx)
library(dplyr)

# Overlap gene de-duplication ####
df <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "All_ARG-carrying_MGE")

# Calculate coverage and filter (Threshold: 90%)
df$coverage <- df$Query_Sequence_Length / df$Subject_Sequence_Length 
df <- df[df$coverage >= 0.9, ]
df$coverage <- ifelse(df$coverage > 1, 1, df$coverage)

# Calculate the midpoint location for MGEs
# Logic: Start position + half of the length
df$midpoint_MGE <- df$ORF_Start + (df$ORF_End - df$ORF_Start) / 2

# De-duplication Function
remove_overlapping_MGEs <- function(df) {
  # Pre-calculate original indices to map back from sub-groups to the main data frame
  df$orig_idx <- seq_len(nrow(df))
  # Initialize a logical vector: TRUE means keep, FALSE means remove
  keep <- rep(TRUE, nrow(df))
  
  # Optimization: Identify groups (Contig + Strand) that have at least 2 entries
  # Note: MGEs are grouped by Contig and Strand to resolve spatial overlaps
  groups <- df %>% 
    group_by(Specific_Contig, Sense_or_Antisense_Strand) %>% 
    summarise(n = n(), .groups = 'drop') %>% 
    filter(n >= 2)
  
  # Iterate through groups with potential overlaps
  for (g in seq_len(nrow(groups))) {
    # Extract global indices belonging to the current group
    idx <- which(df$Specific_Contig == groups$Specific_Contig[g] & 
                   df$Sense_or_Antisense_Strand == groups$Sense_or_Antisense_Strand[g])
    
    # Create a temporary subset for the current group
    df_sub <- df[idx, ]
    n_sub <- nrow(df_sub)
    
    # Pairwise comparison within the group
    for (i in 1:(n_sub - 1)) {
      for (j in (i + 1):n_sub) {
        
        # Get the global row indices for the two MGEs being compared
        idx_i <- idx[i]
        idx_j <- idx[j]
        
        # CRITICAL: If either MGE has already been marked for removal, skip this pair
        if (!keep[idx_i] | !keep[idx_j]) next
        
        # Calculate the physical distance between midpoints
        dist <- abs(df_sub$midpoint_MGE[i] - df_sub$midpoint_MGE[j])
        
        # Calculate the allowed distance (Threshold for overlap: 2 bp buffer as per your setting)
        # Formula: (Length_i + Length_j) / 2 - Buffer
        allowed <- (df_sub$Query_Sequence_Length[i] + df_sub$Query_Sequence_Length[j]) / 2 - 2
        
        # If an overlap is detected (distance < allowed)
        if (!is.na(dist) && dist < allowed) {
          # Competition: Retain the hit with the higher Bitscore
          if (df_sub$Bitscore[i] < df_sub$Bitscore[j]) {
            keep[idx_i] <- FALSE # Mark i for removal
          } else {
            keep[idx_j] <- FALSE # Mark j for removal
          }
        }
      }
    }
  }
  
  # Return the filtered data frame, removing the helper index column
  return(df[keep, -which(names(df) == "orig_idx")])
}

# Execute Filtering
df_final <- remove_overlapping_MGEs(df)

# Final cleanup: Remove any entirely empty or NA rows
df_final <- df_final %>% filter_all(any_vars(. !="" & !is.na(.)))

write.xlsx(df_final, 'C:/Users/CFL/Desktop/UK/PhD/output_file_MGE.xlsx')
