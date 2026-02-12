library(readxl)
library(openxlsx)
library(dplyr)
# Overlap gene de-duplication
df <- read_excel('C:/Users/KEN/Desktop/UK/PhD/AMR/Phenotypic_test/Dorea/Dorea.xlsx', sheet = "Sheet2")

df$gene_length <- as.numeric(df$gene_length)
# filter ARG coverage < 60% (no need for '* 3' if nucleotide sequence)
df$coverage <- (df$length) / df$gene_length 
df <- df[df$coverage >= 0.6, ]
df$coverage <- ifelse(df$coverage > 1, 1, df$coverage)

# determine strand and centre
df$strand <- ifelse(df$qend - df$qstart > 0, "F", "R")
df$q_direction <- NULL
df$s_direction <- NULL
df$center_ARG <- df$qstart + (df$qend - df$qstart) / 2

# Backup the data
df_backup <- df
# remove unsure class
df <- df[df$class != "other", ]

# De-duplication step Two
remove_overlapping_genes <- function(df) {
  
  # Initialise the logic vector to indicate retention of all genes
  keep <- rep(TRUE, nrow(df))
  
  # Group comparable vectors based on the same contig + class + strand
  for (contig in unique(df$Contig)) {
    for (class in unique(df$class)) {
      for (strand in unique(df$strand)) {  
        df_contig_class_strand <- df[df$Contig == contig & df$class == class & df$strand == strand, ]
        
        # Skip groups with only one vector
        if (nrow(df_contig_class_strand) < 2) {
          next
        }
        
        # The vectors of the group compare each other in turn 
        for (i in 1:(nrow(df_contig_class_strand) - 1)) {
          for (j in (i + 1):nrow(df_contig_class_strand)) {
            # Get the center and length of all genes
            center_ARG_i <- df_contig_class_strand$center_ARG[i]
            center_ARG_j <- df_contig_class_strand$center_ARG[j]
            length_i <- df_contig_class_strand$length[i]
            length_j <- df_contig_class_strand$length[j]
            
            # Check for NA values and skip 
            if (is.na(center_ARG_i) || is.na(center_ARG_j) || is.na(length_i) || is.na(length_j) || 
                is.na(df_contig_class_strand$bitscore[i]) || is.na(df_contig_class_strand$bitscore[j])) {
              next
            }
            
            # Calculate the distance from the centers
            distance <- abs(center_ARG_i - center_ARG_j)
            # Calculate the allowed distance from length (less than 4 units are allowed overlapping, can change to 12 if using nucleotide data)
            allowed_distance <- (length_i + length_j) / 2 - 4
            
            # If the overlap area is greater than or equal to 4 units, it is considered to be an overlap gene
            if (!is.na(distance) && distance < allowed_distance) {
              # label and delete overlap genes with smaller bitscore, the labelled objects will not be compared next turn
              if (df_contig_class_strand$bitscore[i] < df_contig_class_strand$bitscore[j]) {
                keep[which(df$Contig == contig & df$class == class & df$strand == strand)[i]] <- FALSE
              } else {
                keep[which(df$Contig == contig & df$class == class & df$strand == strand)[j]] <- FALSE
              }
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

# Extract rows with class "other" from df_backup and add them back to the filtered data
df_other <- df_backup[df_backup$class == "other", ]
df_combined <- rbind(filtered_df, df_other)


df_final <- df_combined %>% filter_all(any_vars(. !="" & !is.na(.)))

# Print results
write.xlsx(df_final, 'output_file.xlsx')



# De-duplication step One
# df <- df[!duplicated(df[, c('Contig', 'cluster_representative')]), ] # can add more variable 
# df <- df[!duplicated(df[, c('Contig', 'class','qstart')]), ]
# df <- df[!duplicated(df[, c('Contig', 'class','qend')]), ]