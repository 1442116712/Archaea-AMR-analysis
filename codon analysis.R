# 安装包（如果尚未安装）
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "coRdon"))
install.packages("aplot")
# 加载包
library(Biostrings)
library(coRdon)
library(ggplot2)
library(stringr)
help(MILC)

# read host genome ####
host_genome <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/Forterrea_sp023230325_CDS.fna")
# create codon table for host genome
cThost <- codonTable(host_genome)
# count codons and lengths for host genome
counts_host <- codonCounts(cThost)
gene_lengths <- getlen(cThost)
# MILC (Measure Independent of Length and Composition) calculation
milc_host <- MILC(cThost, subsets = list(), self = TRUE, ribosomal = FALSE,
                 id_or_name2 = "1", alt.init = TRUE, stop.rm = TRUE,
                 filtering = "none", len.threshold = 80)
milc_host
# weighted mean MILC of host genome 
weighted_milc <- sum(milc_host * gene_lengths) / sum(gene_lengths)
weighted_milc

# read target ARG ####
ARG_seq <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/ARG_MGE_Forterrea_sp023230325.fna")
names(ARG_seq)
# Ensure the length is a multiple of three (essential!)
ARG_seq <- narrow(ARG_seq, start = 1, end = width(ARG_seq) - (width(ARG_seq) %% 3))
# create codon table and check matrix status
cTARG <- codonTable(ARG_seq)  
if (!is.matrix(cTARG@counts)) {
  cTARG@counts <- matrix(cTARG@counts, nrow = 1)
}
class(cTARG)          # "codonTable"
is.matrix(cTARG@counts) # TRUE
dim(cTARG@counts)      # 1 64
# count codons of ARG
counts_ARG <- codonCounts(cTARG)
# MILC (Measure Independent of Length and Composition)
MILC(cTARG, subsets = list(cThost), self = TRUE, ribosomal = FALSE,
             id_or_name2 = "1", alt.init = TRUE, stop.rm = TRUE,
             filtering = "hard", len.threshold = 80)
# CAI (Codon Adaptation Index)
CAI(cTARG, subsets = list(cThost), ribosomal = FALSE, id_or_name2 = "1",
           alt.init = TRUE, stop.rm = TRUE, filtering = "hard",
           len.threshold = 80)
# Fop (frequency of optimal codons)
Fop(cTARG, subsets = list(cThost), ribosomal = FALSE, id_or_name2 = "1",
           alt.init = TRUE, stop.rm = TRUE, filtering = "hard",
           len.threshold = 80)





# References input for RSCU analysis####
host_genome1 <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/Enterococcus_faecalis_CDS.fna")
# create codon table for host genome
cThost1 <- codonTable(host_genome1)
# count codons and lengths for host genome
counts_host1 <- codonCounts(cThost1)

host_genome2 <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/Mammaliicoccus_lentus_CDS.fna")
# create codon table for host genome
cThost2 <- codonTable(host_genome2)
# count codons and lengths for host genome
counts_host2 <- codonCounts(cThost2)

host_genome3 <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/Staphylococcus_aureus_CDS.fna")
# create codon table for host genome
cThost3 <- codonTable(host_genome3)
# count codons and lengths for host genome
counts_host3 <- codonCounts(cThost3)

host_genome4 <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/Staphylococcus_equorum_CDS.fna")
# create codon table for host genome
cThost4 <- codonTable(host_genome4)
# count codons and lengths for host genome
counts_host4 <- codonCounts(cThost4)

host_genome5 <- readDNAStringSet("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Pan_genome/CDS_codon_analysis/GCA_023230325.1_CDS.fna")
# create codon table for host genome
cThost5 <- codonTable(host_genome5)
# count codons and lengths for host genome
counts_host5 <- codonCounts(cThost5)


# combine reference genomes and target gene
codon_counts <- rbind(counts_host, counts_host1, counts_host2, counts_host3, counts_host4, counts_host5, counts_ARG)
group <- c(rep("Forterrea p023230325 core", nrow(counts_host)), 
           rep("Enterococcus faecalis core", nrow(counts_host1)), 
           rep("Mammaliicoccus lentus core", nrow(counts_host2)), 
           rep("Staphylococcus aureus core", nrow(counts_host3)), 
           rep("Staphylococcus equorum core", nrow(counts_host4)),  
           rep("Forterrea p023230325 host", nrow(counts_host5)),  
           rep("ARG/MGE", nrow(counts_ARG)))
# create amino acid table
std_codons <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
                "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
                "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
                "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
                "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
                "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")

std_aas <- c("K", "N", "K", "N", "T", "T", "T", "T",
                     "R", "S", "R", "S", "I", "I", "M", "I",
                     "Q", "H", "Q", "H", "P", "P", "P", "P",
                     "R", "R", "R", "R", "L", "L", "L", "L",
                     "E", "D", "E", "D", "A", "A", "A", "A",
                     "G", "G", "G", "G", "V", "V", "V", "V",
                     "*", "Y", "*", "Y", "S", "S", "S", "S",
                     "*", "C", "W", "C", "L", "F", "L", "F")

std_aas_string <- paste(std_aas, collapse = "")
codon_to_aa <- setNames(strsplit(std_aas_string, "")[[1]], std_codons)
aa_codons <- list()
for (codon in colnames(codon_counts)) {
  aa <- codon_to_aa[codon]
  if (!is.na(aa) && aa != "*") {  # 排除NA和终止密码子
    if (!aa %in% names(aa_codons)) {
      aa_codons[[aa]] <- c()
    }
    aa_codons[[aa]] <- c(aa_codons[[aa]], codon)
  }
}

# remove stop codon
stop_codons <- c("TAA", "TAG", "TGA")
codon_counts_final <- codon_counts[, !colnames(codon_counts) %in% stop_codons, drop = FALSE]

# pseudo-count to avoid 0 value in certain amino acid
pseudocount <- 0.001
codon_counts_pseudo <- codon_counts_final + pseudocount


# RSCU (Relative Synonymous Codon Usage) calculation
rscu_mat <- codon_counts_pseudo
for (aa in names(aa_codons)) {
  codons <- aa_codons[[aa]]
  n_codons <- length(codons)
  rscu_mat[, codons] <- codon_counts[, codons] / (rowSums(codon_counts_pseudo[, codons, drop=FALSE]) / n_codons)
}

# host gene name setting
if (is.null(rownames(rscu_mat))) {
  rownames(rscu_mat) <- paste0("gene", seq_len(nrow(rscu_mat)))
}

# PCA calculation and create dataset
pca_res <- prcomp(rscu_mat, scale.=TRUE)
pca_df <- data.frame(PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     gene = rownames(rscu_mat),
                     group = group)
# Setting of the ARG names
pca_df$gene[pca_df$group == "ARG/MGE"] <- c("str", "plasmid_str")
# PCA visualisation
library(ggrepel)
colours <- c("#D62728", "#1F77B4", "#2CA02C", "#BA55D3", "#FF7F0E", "#FFD700","#FF69B4") 

PCA <- ggplot(pca_df, aes(x=PC1, y=PC2, label=gene, color=group)) +
  geom_point(size=3, alpha=0.6) +
  geom_text_repel(data=subset(pca_df, group=="ARG/MGE"), 
                  aes(label=gene),
                  size=3) +
  scale_color_manual(values = c(colours, "60")) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  ggtitle("RSCU PCA of Forterrea sp023230325") + 
  labs(color="") + 
  coord_cartesian(xlim=c(-20, 20), ylim=c(-20, 20)) +
  guides(color = guide_legend(override.aes = list(label=""))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12)) #WIDTH 1000 HEIGHT 750

ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/RSCU/Forterrea_p023230325_PCA.png",     
  plot = PCA,      
  width = 8,       
  height = 8,       
  dpi = 300,       
  units = "in",
  bg = "white" 
)

write.csv(rscu_mat, "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/rscu_mat.csv")



# 宿主的加权平均RSCU #### 
start <- 283635
end <- 283635
rscu_host <- rscu_mat[start:end, ] #提取RSCU
host_counts <- codon_counts_final[start:end, ] #提取总密码子数（无终止子）
# 设计计算加权平均的函数
calculate_weighted_avg_rscu <- function(rscu_values, total_codon_number, aa_codons) {
  # 初始化结果矩阵
  weighted_avg <- numeric(ncol(rscu_values))
  names(weighted_avg) <- colnames(rscu_values)
  # 计算目标基因群中每个基因的总密码子数（无终止子）作为权重
  gene_weights <- rowSums(total_codon_number)
  # 对每个密码子计算加权平均
  for (codon in colnames(rscu_values)) {
    # 获取该密码子的RSCU值
    rscu_codon <- rscu_values[, codon]
    # 计算加权平均：权重是基因的总密码子数（无终止子）
    weighted_avg[codon] <- sum(rscu_codon * gene_weights) / sum(gene_weights)
  }
  return(list(
    weighted_avg_all = weighted_avg,
    gene_weights = gene_weights
  ))
}

# 计算counts_host3的加权平均RSCU并显示
results_host <- calculate_weighted_avg_rscu(rscu_host, host_counts, aa_codons)
AA61_weighted_avg_rscu <- head(sort(results_host$weighted_avg_all, decreasing = TRUE), 61)
print(round(AA61_weighted_avg_rscu, 3))



# 密码子使用偏好详解图 #### 
library(aplot)
cols <- c("#BDB76B", "#CAB2D6", "#FF69B4", "#B2DF8A", "#FF9D1E", "#6495ED")
df <- read.table("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/RSCU/AA_table.txt",header = F, stringsAsFactors = F)
rscu_df <- data.frame(
  V1 = names(AA61_weighted_avg_rscu),
  V4 = round(AA61_weighted_avg_rscu, 3)
)
df <- merge(df, rscu_df, by = "V1", all.x = TRUE)

p1 <- ggplot(df, aes(x = V2, y = V4, fill = as.character(V3))) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = c(cols)) +
  geom_text(aes(label = ifelse(V4 > 0.1, paste0(V4, ""), "")),
            position = position_stack(vjust = 0.5), size = 2.5, color = "black", fontface = "plain") +
  labs(x = NULL, y = "RSCU", 
       title = "RSCU of plasmid_str") +
  scale_y_continuous(expand = c(0,0), limits = c(0,6.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),  # x轴标题大小和样式
        axis.title.y = element_text(size = 15, face = "bold"),  # y轴标题大小和样式
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "none")
p1
p2 <- ggplot(df, aes(x = V2, y = V3))+
  geom_label(aes(label=V1, fill=as.character(V3)),size=2)+
  scale_fill_manual(values = c(cols)) +
  labs(x="",y="") + ylim(3.4,6.3)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
p2
p <- insert_bottom(p1,p2,height = 0.2) # width 700 height 650
  
p
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/RSCU/plasmid_str.png",     
  plot = p,      
  width = 8,       
  height = 8,       
  dpi = 300,       
  units = "in",
  bg = "white" 
)
# NOTE!!!
#MILC → Antibiotic Resistance Genes and Host Deviation
#Very Low	0 ~ 0.1	Gene codon usage closely resembles host
#Moderate Deviation	0.1 ~ 0.3	Some divergence from host, but not extreme
#Significant Deviation	>0.3	Marked difference from host codon usage → High likelihood of HGT

#CAI → Adaptation level of the drug resistance gene relative to the host
#Closer to 1 → Greater alignment with reference set preferences (highly expressed genes or host-preferred codons)
#Closer to 0 → Deviation from reference set (potentially from exogenous sources or low expression)

#Fop → Frequency of host-optimal codon usage by drug resistance genes
#0–0.3	Low → Codon usage deviates from host-optimal preferences
#0.3–0.6	Moderate → Partial alignment with host preferences
#0.6–1	High → Codon usage aligns with host-optimal codons, likely adapted to host

#ENC/MCB/SCUO → Resistance gene preference intensity, compared to host average
#Resistance gene with significant deviation from host average SCUO → Likely of exogenous origin

#RSCU, Relative Synonymous Codon Usage 相对同义密码子使用偏好
#  | AAA | AAG | Total | n_codons |   RSCU(AAA)  |  RSCU(AAG)   |
#  | --- | --- |  ---  | -------- | ------------ | ------------ |
#  | 8   |  2  |   10  |    2     | 8/(10/2)=1.6 | 2/(10/2)=0.4 |
# > 1.6	明显偏好使用（over-represented）
# 1.2–1.6	偏好使用（preferentially used）
# 0.8–1.2	无明显偏好（neutral）
# 0.5–0.8	偏低使用（under-represented）
# < 0.5	明显不使用（under-represented）
