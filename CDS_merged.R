# 加载必要的包
library(Biostrings)

merge_cds <- function(input_file, output_file) {
  # 1. 读取FASTA
  fasta <- readDNAStringSet(input_file)
  
  # 2. 提取GCA号
  headers <- names(fasta)
  sequences <- as.character(fasta)
  
  # 存储GCA号
  gca_ids <- character(length(headers))
  
  for(i in 1:length(headers)) {
    header <- headers[i]
    # 去掉>
    header <- gsub("^>", "", header)
    # 取第一个空格前的内容
    header <- strsplit(header, " ")[[1]][1]
    # 取GCA_数字.数字
    if(grepl("GCA_", header)) {
      # 找到GCA_的位置并提取到第二个下划线
      parts <- strsplit(header, "_")[[1]]
      if(length(parts) >= 2) {
        gca_ids[i] <- paste(parts[1], parts[2], sep = "_")
      } else {
        gca_ids[i] <- header
      }
    } else {
      gca_ids[i] <- header
    }
  }
  
  # 3. 按GCA分组并直接拼接
  unique_gcas <- unique(gca_ids)
  
  # 存储结果
  merged_seqs <- character(length(unique_gcas))
  merged_headers <- character(length(unique_gcas))
  fragment_counts <- integer(length(unique_gcas))
  
  for(i in 1:length(unique_gcas)) {
    gca <- unique_gcas[i]
    
    # 找到所有属于这个GCA的序列
    idx <- which(gca_ids == gca)
    
    # 直接拼接（完全按照原始顺序）
    merged_seq <- paste(sequences[idx], collapse = "")
    merged_seqs[i] <- merged_seq
    fragment_counts[i] <- length(idx)
    merged_headers[i] <- paste0(">", gca, "_coregene_", length(idx), "CDS")
  }
  # 创建DNAStringSet
  result_seqs <- DNAStringSet(merged_seqs)
  names(result_seqs) <- merged_headers
  
  # 写入文件
  writeXStringSet(result_seqs, output_file, width = 60)
  
  # 5. 简单统计
  cat(paste(rep("#", 60), collapse = ""), "\n")
  cat("                   完成！\n")
  cat(paste(rep("#", 60), collapse = ""), "\n")
  cat("\n统计信息:\n")
  cat("  输入片段数:", length(fasta), "\n")
  cat("  输出基因数:", length(result_seqs), "\n")
  cat("  平均每个基因的片段数:", round(mean(fragment_counts), 2), "\n")
  cat("  最大片段数:", max(fragment_counts), "\n")
  cat("  总碱基数:", sum(width(result_seqs)), "bp\n")
  cat("  平均基因长度:", round(mean(width(result_seqs))), "bp\n")
  
  # 返回结果
  return(list(
    gca_ids = unique_gcas,
    sequences = result_seqs,
    fragment_counts = fragment_counts,
    merged_seqs = merged_seqs
  ))
}

batch_merge_fna <- function(input_folder, output_folder) {
  
  if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  fna_files <- list.files(input_folder, pattern = "\\.fna$", full.names = TRUE, ignore.case = TRUE)
  
  for(f in 1:length(fna_files)) {
    input_file <- fna_files[f]
    file_name <- basename(input_file)
    file_name_noext <- tools::file_path_sans_ext(file_name)
    output_file <- file.path(output_folder, paste0(file_name_noext, "_merged.fasta"))
    
    tryCatch({
      merge_cds(input_file, output_file)
    }, error = function(e) {
    })
  }
}

input_folder <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/"
output_folder <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/Merged/"
batch_merge_fna(input_folder, output_folder)

