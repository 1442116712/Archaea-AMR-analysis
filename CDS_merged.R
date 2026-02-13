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
    # 提取到第二个下划线
    parts <- strsplit(header, "_")[[1]]
    if(length(parts) >= 2) {
      gca_ids[i] <- paste(parts[1], parts[2], sep = "_")
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



# For host genome, group each 20 CDS into one to reduce the number of CDS
merge_cds <- function(input_file, output_file, group_size = 20) {
  # 1. 读取FASTA
  fasta <- readDNAStringSet(input_file)
  
  # 2. 提取序列信息
  headers <- names(fasta)
  sequences <- as.character(fasta)
  
  # 3. 按每group_size个序列直接合并
  n_seq <- length(sequences)
  n_groups <- ceiling(n_seq / group_size)
  
  # 存储结果
  merged_seqs <- character(n_groups)
  merged_headers <- character(n_groups)
  fragment_counts <- integer(n_groups)
  
  for(i in 1:n_groups) {
    # 确定当前组的序列索引
    start_idx <- (i-1) * group_size + 1
    end_idx <- min(i * group_size, n_seq)
    
    # 拼接当前组的序列
    merged_seq <- paste(sequences[start_idx:end_idx], collapse = "")
    merged_seqs[i] <- merged_seq
    
    # 记录片段数
    fragment_counts[i] <- end_idx - start_idx + 1
    
    # 生成header
    merged_headers[i] <- paste0(">merged_group_", i, "_", fragment_counts[i], "CDS")
  }
  
  # 4. 创建DNAStringSet
  result_seqs <- DNAStringSet(merged_seqs)
  names(result_seqs) <- merged_headers
  
  # 5. 写入文件
  writeXStringSet(result_seqs, output_file, width = 60)
  
  # 6. 统计信息
  cat(paste(rep("#", 60), collapse = ""), "\n")
  cat("                   完成！\n")
  cat(paste(rep("#", 60), collapse = ""), "\n")
  cat("\n统计信息:\n")
  cat("  输入片段数:", length(fasta), "\n")
  cat("  输出基因数:", length(result_seqs), "\n")
  cat("  分组大小:", group_size, "个序列/组\n")
  cat("  平均每个输出的片段数:", round(mean(fragment_counts), 2), "\n")
  cat("  最大片段数:", max(fragment_counts), "\n")
  cat("  最小片段数:", min(fragment_counts), "\n")
  cat("  总碱基数:", sum(width(result_seqs)), "bp\n")
  cat("  平均基因长度:", round(mean(width(result_seqs))), "bp\n")
  
  # 返回结果
  return(list(
    sequences = result_seqs,
    fragment_counts = fragment_counts,
    merged_seqs = merged_seqs,
    group_size = group_size
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
