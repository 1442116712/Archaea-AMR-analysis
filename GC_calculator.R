# 如果尚未安装BiocManager，先安装它
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

install.packages("seqinr")

library(Biostrings)
library(seqinr)

# 方法1：使用Biostrings读取FASTA文件
sequences <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Pan_genome/Methanosarcina_mazei_core.aln")

# 计算每个序列的GC含量
gc_contents <- letterFrequency(sequences, "GC", as.prob = TRUE) * 100
sequence_lengths <- width(sequences)

# 创建详细的结果表格
results <- data.frame(
  Sequence_ID = names(sequences),
  GC_Content = as.numeric(gc_contents),
  Length = sequence_lengths
)

# 计算加权平均GC含量（按序列长度加权）
weighted_mean_gc <- sum(gc_contents * sequence_lengths) / sum(sequence_lengths)

# 计算简单平均GC含量
simple_mean_gc <- mean(gc_contents)

# 计算总GC含量（将所有序列合并后计算）
all_sequences_combined <- paste(as.character(sequences), collapse = "")
total_gc <- letterFrequency(DNAString(all_sequences_combined), "GC", as.prob = TRUE) * 100

cat("\n=== GC含量统计结果 ===\n")
cat("序列数量:", length(sequences), "\n")
cat("总碱基数:", sum(sequence_lengths), "\n")
cat("加权平均GC含量:", round(weighted_mean_gc, 2), "%\n")
cat("简单平均GC含量:", round(simple_mean_gc, 2), "%\n")
cat("总体GC含量:", round(total_gc, 2), "%\n")

