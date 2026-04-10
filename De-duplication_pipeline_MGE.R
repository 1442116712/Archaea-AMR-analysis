library(readxl)
library(openxlsx)
library(dplyr)

# Overlap MGE de-duplication ####
df <- read_excel('C:/Users/CFL/Desktop/UK/PhD/1.xlsx', sheet = "Sheet1")

# Calculate coverage and filter (Threshold: 90%)
# Updated variable names to lowercase: query_sequence_length, subject_sequence_length
df$coverage <- df$query_sequence_length / df$subject_sequence_length 
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
  # Updated variable names: specific_contig, strand
  groups <- df %>% 
    group_by(specific_contig, strand) %>% 
    summarise(n = n(), .groups = 'drop') %>% 
    filter(n >= 2)
  
  # Iterate through groups with potential overlaps
  for (g in seq_len(nrow(groups))) {
    # Extract global indices belonging to the current group
    idx <- which(df$specific_contig == groups$specific_contig[g] & 
                   df$strand == groups$strand[g])
    
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
        
        # Calculate the allowed distance (Threshold for overlap: 2 bp buffer)
        # Updated variable names: query_sequence_length
        allowed <- (df_sub$query_sequence_length[i] + df_sub$query_sequence_length[j]) / 2 - 2
        
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





# ARG - MGE link ####
df1 <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "ARG")
df_final <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "All_unique_ARG-carrying_MGE")
link_MGE_ARG <- function(mge_df, arg_df) {
  
  # 1. 提取 ARG 表中需要的信息 (注意：使用我们改好的 midpoint_ARG)
  arg_subset <- arg_df %>%
    select(Contig, class, fa_name, midpoint_ARG, Kingdom, Phylum, Class, Order, Family, Genera, Species, Genome) %>%
    # 重命名列以防混淆
    rename(ARG_class = class, ARG_fa_name = fa_name)
  
  # 2. 根据 Contig 建立所有的潜在组合 (many-to-many 关系)
  linked_df <- mge_df %>%
    inner_join(arg_subset, by = c("specific_contig" = "Contig"), relationship = "many-to-many")
  
  # 3. 计算距离并执行 10kb (10000 bp) 过滤
  filtered_links <- linked_df %>%
    # 注意：使用 midpoint_MGE
    mutate(gap_distance = abs(midpoint_MGE - midpoint_ARG)) %>%
    filter(gap_distance <= 10000)
  
  return(filtered_links)
}

# 运行函数
df_final_linked <- link_MGE_ARG(df_final, df1)

# Print results
write.xlsx(df_final_linked, 'C:/Users/CFL/Desktop/UK/PhD/output_file_MGE.xlsx')
