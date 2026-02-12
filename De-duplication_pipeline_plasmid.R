library(readxl)
library(openxlsx)
library(dplyr)
# Overlap gene de-duplication
df <- read_excel('C:/Users/KEN/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "Phage")

# determine strand and centre
df$strand <- ifelse(df$qend - df$qstart > 0, "forward", "reverse")
df$center <- df$qstart + (df$qend - df$qstart) / 2

# Backup the data
df_backup <- df

# De-duplication step Two
remove_overlapping_genes <- function(df) {
  
  # Initialise the logic vector to indicate retention of all genes
  keep <- rep(TRUE, nrow(df))
  
  # Group comparable vectors based on the same contig + class + strand
  for (contig in unique(df$Contig)) {
      for (strand in unique(df$strand)) {  
        df_contig_strand <- df[df$Contig == contig & df$strand == strand, ]
        
        # Skip groups with only one vector
        if (nrow(df_contig_strand) < 2) {
          next
        }
        
        # The vectors of the group compare each other in turn 
        for (i in 1:(nrow(df_contig_strand) - 1)) {
          for (j in (i + 1):nrow(df_contig_strand)) {
            # Get the center and length of all genes
            center_i <- df_contig_strand$center[i]
            center_j <- df_contig_strand$center[j]
            length_i <- df_contig_strand$length[i]
            length_j <- df_contig_strand$length[j]
            
            # Check for NA values and skip 
            if (is.na(center_i) || is.na(center_j) || is.na(length_i) || is.na(length_j) || 
                is.na(df_contig_strand$bitscore[i]) || is.na(df_contig_strand$bitscore[j])) {
              next
            }
            
            # Calculate the distance from the centers
            distance <- abs(center_i - center_j)
            # Calculate the allowed distance from length (less than 4 units are allowed overlapping, can change to 12 if using nucleotide data)
            allowed_distance <- (length_i + length_j) / 2 - 12
            
            # If the overlap area is greater than or equal to 4 units, it is considered to be an overlap gene
            if (!is.na(distance) && distance < allowed_distance) {
              # label and delete overlap genes with smaller bitscore, the labelled objects will not be compared next turn
              if (df_contig_strand$bitscore[i] < df_contig_strand$bitscore[j]) {
                keep[which(df$Contig == contig & df$strand == strand)[i]] <- FALSE
              } else {
                keep[which(df$Contig == contig & df$strand == strand)[j]] <- FALSE
              }
            }
          }
        }
      }
  }
  
  # Return of retained genes
  return(df[keep, ])
}

# Output the filter result
filtered_df <- remove_overlapping_genes(df)

# Print results
write.xlsx(filtered_df, 'output_file.xlsx')

