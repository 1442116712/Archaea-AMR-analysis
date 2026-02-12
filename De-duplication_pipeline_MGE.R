library(readxl)
library(openxlsx)
library(dplyr)

#Necessary information in excel
#MGE - Specific_Contig,Gene_Name,ORF_End,ORF_Start,End/Start_of_Alignment_in_Subject/Query,Query_Sequence_Length,Subject_Sequence_Length,Bitscore
#ARG - Contig,class,fa_name,qend,qstart,send,start

# Overlap gene de-duplication ####
df <- read_excel('C:/Users/KEN/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "MGE")

# filter ARG coverage < 90%
df$coverage <- df$Query_Sequence_Length / df$Subject_Sequence_Length 
df <- df[df$coverage >= 0.9, ]
df$coverage <- ifelse(df$coverage > 1, 1, df$coverage)

# determine strand and midpoint
df$strand <- ifelse(df$End_of_Alignment_in_Query - df$Start_of_Alignment_in_Query > 0, "F", "R")
df$q_direction <- NULL
df$s_direction <- NULL
df$center_MGE <- df$ORF_Start + (df$ORF_End - df$ORF_Start) / 2

# De-duplication
remove_overlapping_MGEs <- function(df) {
  
  # Initialise the logic vector to indicate retention of all genes
  keep <- rep(TRUE, nrow(df))
  
  # Group comparable vectors based on the same contig + class + strand
  for (contig in unique(df$Specific_Contig)) {
    for (gene in unique(df$Gene_Name)) {
      for (strand in unique(df$strand)) {  
        df_contig_gene_strand <- df[df$Specific_Contig == contig & df$Gene_Name == gene & df$strand == strand, ]
        
        # Skip groups with only one vector
        if (nrow(df_contig_gene_strand) < 2) {
          next
        }
        
        # The vectors of the group compare each other in turn 
        for (i in 1:(nrow(df_contig_gene_strand) - 1)) {
          for (j in (i + 1):nrow(df_contig_gene_strand)) {
            # Get the center and length of all genes
            center_i <- df_contig_gene_strand$center_MGE[i]
            center_j <- df_contig_gene_strand$center_MGE[j]
            length_i <- df_contig_gene_strand$Query_Sequence_Length[i]
            length_j <- df_contig_gene_strand$Query_Sequence_Length[j]
            
            # Check for NA values and skip 
            if (is.na(center_i) || is.na(center_j) || is.na(length_i) || is.na(length_j) || 
                is.na(df_contig_gene_strand$Bitscore[i]) || is.na(df_contig_gene_strand$Bitscore[j])) {
              next
            }
            
            # Calculate the distance from the centers
            actual_distance <- abs(center_i - center_j)
            # Calculate the allowed distance from length (less than 4 units are allowed overlapping, can change to 12 if using nucleotide data)
            allowed_distance <- (length_i + length_j) / 2 - 12
            
            # If the overlap area is greater than or equal to 12 units, it is considered to be an overlap gene
            if (!is.na(actual_distance) && actual_distance < allowed_distance) {
              # label and delete overlap genes with smaller bitscore, the labelled objects will not be compared next turn
              if (df_contig_gene_strand$Bitscore[i] < df_contig_gene_strand$Bitscore[j]) {
                keep[which(df$Specific_Contig == contig & df$Gene_Name == gene & df$strand == strand)[i]] <- FALSE
              } else {
                keep[which(df$Specific_Contig == contig & df$Gene_Name == gene & df$strand == strand)[j]] <- FALSE
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
filtered_df <- remove_overlapping_MGEs(df)
df_final <- filtered_df %>% filter_all(any_vars(. !="" & !is.na(.)))
# Backup the data
df_backup <- df_final

# ARG Grouping ####
df1 <- read_excel('C:/Users/KEN/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "ARG")

# 定义分组函数
group_ARGs <- function(df1) {
  df1$Group <- NA  # 创建 Group 列
  
  # 按 Contig 分组
  for (contig in unique(df1$Contig)) {
    subset_df <- df1[df1$Contig == contig, ]  # 选出当前 Contig 的数据
    
    # 按 class 分组
    for (class in unique(subset_df$class)) {
      class_rows <- subset_df[subset_df$class == class, ]  # 选出当前 class 的数据
      
      # 按行分配组号
      if (nrow(class_rows) > 0) {
        df1$Group[df1$Contig == contig & df1$class == class] <- 1:nrow(class_rows)
      }
    }
  }
  
  return(df1)
}

df1_grouped <- group_ARGs(df1)

# ARG与MGE匹配,并筛选中点间距小于10kb的结果，填入对应的ARG信息 ####
MGE_compilation <- function(df_final, df1_grouped) {
  # 获取所有唯一的class值
  unique_classes <- unique(df1_grouped$class)
  
  # 在df_final中为每个class创建一列，列名加上 _gap
  df_final[, paste0(unique_classes, "_gap")] <- NA
  
  # 为每个class创建 fa_name 列
  for (class_name in unique_classes) {
    fa_name_col <- paste0(class_name, "_fa_name")
    df_final[[fa_name_col]] <- NA
  }
  
  # 逐步匹配并填充数据
  for (group in unique(df1_grouped$Group)) {
    df_grouped_1 <- df1_grouped %>% filter(Group == group)
    
    for (class_name in unique(df_grouped_1$class)) {
      df_class <- df_grouped_1 %>% filter(class == class_name)
      
      # 使用 center_ARG 列进行匹配，顺便匹配 fa_name
      matched_data <- df_class %>% select(Contig, center_ARG, fa_name)
      
      # 执行左连接，将 center_ARG 和 fa_name 填充到 df_final 中
      df_final <- df_final %>%
        left_join(matched_data, by = c("Specific_Contig" = "Contig"))
      
      if ("center_ARG" %in% colnames(df_final)) {  # 确保 center_ARG 列存在
        # 填充 center_ARG 和 fa_name
        df_final <- df_final %>%
          mutate(
            !!paste0(class_name, "_gap") := ifelse(!is.na(center_ARG),
                                                   ifelse(is.na(!!sym(paste0(class_name, "_gap"))), center_ARG, paste(!!sym(paste0(class_name, "_gap")), center_ARG, sep = ";")),
                                                   !!sym(paste0(class_name, "_gap"))),
            !!paste0(class_name, "_fa_name") := ifelse(!is.na(fa_name),
                                                       ifelse(is.na(!!sym(paste0(class_name, "_fa_name"))), fa_name, paste(!!sym(paste0(class_name, "_fa_name")), fa_name, sep = ";")),
                                                       !!sym(paste0(class_name, "_fa_name")))
          ) %>%
          select(-center_ARG, -fa_name)  # 删除临时匹配的 center_ARG 和 fa_name 列
      }
    }
  }
  
  # 计算 center_ARG 与 center_MGE 之间的差的绝对值，并更新 center_ARG 列
  df_final <- df_final %>%
    mutate(across(all_of(paste0(unique_classes, "_gap")), 
                  ~ sapply(seq_along(.), function(i) {
                    val <- .[i]
                    if (!is.na(val)) {  # 只有非NA值才进行处理
                      # 确保 val 是字符类型
                      val <- as.character(val)
                      
                      # 检查是否有分号分隔的多个值
                      values <- unlist(strsplit(val, ";"))
                      
                      # 获取该行对应的 center_MGE 值
                      center_MGE_value <- center_MGE[i]
                      
                      # 对每个值计算与该行对应的 center_MGE 的差的绝对值
                      abs_diffs <- sapply(values, function(v) abs(as.numeric(v) - center_MGE_value))
                      
                      # 只保留差值小于等于 10000 的值
                      filtered_values <- values[abs_diffs <= 10000]
                      filtered_diffs <- abs_diffs[abs_diffs <= 10000]
                      
                      # 如果有多个差值，合并为一个字符串
                      if (length(filtered_diffs) == 1) {
                        return(filtered_diffs[1])  # 如果只有一个差值，直接返回
                      } else if (length(filtered_diffs) > 1) {
                        return(paste(filtered_diffs, collapse = ";"))  # 如果有多个差值，合并为一个字符串
                      } else {
                        return(NA)  # 如果没有符合条件的值，返回 NA
                      }
                    } else {
                      return(NA)  # 如果是NA，直接返回NA
                    }
                  })))
  
  # 过滤掉所有差值大于 10000 的行对应的fa_name
  df_final <- df_final %>%
    filter(if_any(all_of(paste0(unique_classes, "_gap")), ~ !is.na(.))) %>%
    mutate(across(ends_with("_fa_name"), ~ ifelse(is.na(.), "NA", .))) %>%
    mutate(across(ends_with("_gap"), ~ ifelse(is.na(.), "NA", .))) %>% # 空白值全部换成NA字符
    select_if(~ !all(. == "NA")) # 删除没有结果的列
  
  df_final <- df_final %>%
    left_join(
      df1 %>%
        select(Contig, Kingdom, Phylum, Class, Order, Family, Genera, Species, Genome) %>%
        distinct(Contig, .keep_all = TRUE),  # 去重，保留每个 Contig 的第一行
      by = c("Specific_Contig" = "Contig")
    )
  
  # 返回最终的df_final
  return(df_final)
}

# 调用函数并将结果输出到文件
df_final <- MGE_compilation(df_final, df1_grouped)

# Print results
write.xlsx(df_final, 'C:/Users/KEN/Desktop/UK/PhD/output_file_MGE.xlsx')
