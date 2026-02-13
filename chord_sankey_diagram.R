install.packages("readxl")
install.packages("tidyverse")
install.packages("devtools")
install.packages("ggrepel")
devtools::install_github("r-lib/conflicted")
install.packages("ggalluvial")

library(readxl)
library(tidyverse)
library(ggalluvial)


data <- read_excel("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx", sheet = "co-relation_figure_MAG")

# sankey diagram ####
sankey_data <- data %>%
  group_by(Country_Archaea, Sample_type_Archaea, Family_Archaea, ARG, ARG_host_PanRes, Sample_type_PanRes, Country_PanRes) %>%
  summarise(Count = n()) %>%
  ungroup()

group_colors <- c(
  "aminoglycosides_ant6" = "#89CFF0",
  "aminoglycosides_ant9" = "#4169E1",
  "polymyxin" = "#6A5ACD",
  "beta_lactams_amp" = "#E74C3C",
  "beta_lactams_tem" = "#fc9a99",
  "tetracyclines_tet44" = "#27AE60",
  "tetracyclines_tetc" = "#2ECC71",
  "chloramphenicol_cata" = "#16A085",
  "mls_mefa" = "#C5B0D5",
  "mls_mefEn2" = "#8E44AD",
  "lincosamides_lnuAn2" = "#756BB1",
  "lincosamides_lnuf" = "#E377C2",
  "rifampicin_rph" = "#F39C12",
  "ef_tu_inhibition_tufab" = "#F1C40F",
  
  # 重新设计的铜抗性颜色 - 增加区分度
  "copper_resistance_copb" = "#8B4513",  # 深棕色
  "copper_resistance_copr" = "#B8860B",  # 深金色
  "copper_resistance_dnak" = "#708090",  # 石板灰
  "copper_resistance_tcrb" = "#A0522D",  # 赭色
  "copper_resistance_tcry" = "#DAA520",  # 金色
  "copper_resistance_tcrz" = "#BC8F8F",  # 玫瑰棕
  
  # 调整其他灰色系颜色
  "multi_biocide_hsmr" = "#2C3E50",      # 深灰蓝色
  "multi_metal_copam" = "#5D6D7E",       # 保持原色
  "multi_metal_wtpa" = "#7D7D7D",        # 中灰色
  "multi_metal_wtpb" = "#A0A0A0",        # 浅灰色
  "multi_metal_wtpc" = "#C8C8C8",         # 更浅的灰色
  
  "Australia" = "#1F77B4",
  "Brazil" = "#2CA02C",
  "Tunisia" = "#9467BD",
  "Italy" = "#FF7F0E",
  "Japan" = "#D62728",
  "Norway" = "#8C564B",
  "Spain" = "#E377C2",
  "Germany" = "#BCBD22",
  "Mexico" = "#17BECF",
  "Finland" = "#FF9896",
  "Iceland" = "#98DF8A",
  "Bangladesh" = "#C49C94",
  "Pacific Ocean" = "#C5B0D5",
  "UK" = "#F7B6D2",
  "India" = "#DBDB8D",
  "China" = "#9EDAE5",
  "USA" = "#393B79",
  "Kazakhstan" = "#8C6D31",
  "Netherlands" = "#843C39",
  "Canada" = "#3182BD",
  "Philippines" = "#31A354",
  "Fiji" = "#756BB1",
  "Bolivia" = "#E6550D",
  "Russia" = "#6BAED6",
  "Ghana" = "#74C476",
  "Indonesia" = "#FD8D3C",
  "Europe" = "#9E9AC8",
  "Switzerland" = "#FDBC85",
  "Venezuela" = "#E7969C",
  "Denmark" = "#7B4173",
  "Poland" = "#A1D99B",
  "Thailand" = "#CEDB9C",
  
  "oil reservoir" = "#E41A1C",
  "pickled food" = "#377EB8",
  "salt lake" = "#4DAF4A",
  "human gut" = "#984EA3",
  "river sediment" = "#FF7F00",
  "hotspring" = "#FFFF33",
  "saltern" = "#A65628",
  "Others" = "#F781BF",
  "animal gut" = "#B15928",   
  "marine sediment" = "#66C2A5",
  "crustal fluids" = "#FC8D62",
  "soil" = "#8DA0CB",
  "wastewater" = "#E78AC3",
  "fresh water" = "#FFD92F",
  "clinical sample" = "#E5C494",
  "food" = "#B3DE69", 
  
  "Archaeoglobaceae" = "#8DD3C7",
  "Halobacteriaceae" = "#FFFFB3",
  "Methanobacteriaceae" = "#BEBADA",
  "Methanosarcinaceae" = "#FB8072",
  "Sulfolobaceae" = "#80B1D3",
  "Thermococcaceae" = "#FDB462",
  "Natronoarchaeaceae" = "#B3DE69",
  "Haloarculaceae" = "#FCCDE5",
  "Haloferacaceae" = "#D9D9D9",
  "Methanomassiliicoccaceae" = "#BC80BD",
  "Acidilobaceae" = "#CCEBC5",
  "RBG-16-68-12" = "#FFED6F",
  "Methanospirillaceae" = "#E5C494",
  "Methanotrichaceae" = "#7FC97F",
  "Methanoculleaceae" = "#BEAED4",
  "CSSED10-239" = "#FDC086",
  "Methanocorpusculaceae" = "#FFFF99",
  
  "Archaeoglobus fulgidus" = "#1F77B4",
  "Halobacterium salinarum" = "#FF7F0E",
  "Klebsiella pneumoniae" = "#2CA02C",
  "Escherichia coli" = "#D62728",
  "Sulfolobus solfataricus" = "#9467BD",
  "Aeromonas salmonicida" = "#8C564B",
  "Gallid alphaherpesvirus" = "#E377C2",
  "Providencia rettgeri" = "#7F7F7F",
  "Mycobacterium tuberculosis" = "#BCBD22",
  "Clostridioides difficile" = "#17BECF",
  "Pyrococcus furiosus" = "#FF9896",
  "Escherichia phage" = "#98DF8A",
  "Clostridium perfringens" = "#C49C94",
  "Pseudomonas aeruginosa" = "#C5B0D5",
  "Planobispora rosea" = "#F7B6D2",
  "Enterococcus faecium" = "#DBDB8D",
  "Proteus mirabilis" = "#9EDAE5",
  "Salmonella enterica" = "#393B79",
  "Clostridium kluyveri" = "#8C6D31",
  "Lactobacillus johnsonii" = "#843C39",
  "Lactococcus lactis" = "#3182BD",
  "Brevibacillus brevis" = "#31A354",
  "Streptomyces cinnamoneus" = "#756BB1",
  "Bacteroides fragilis" = "#E6550D",
  "Campylobacter fetus" = "#6BAED6",
  "Enterobacter asburiae" = "#74C476",
  "Serratia marcescens" = "#FD8D3C",
  "Listeria monocytogenes" = "#9E9AC8",
  
  "unknown" = "#7F7F7F"
)


p <- ggplot(sankey_data,
       aes(axis1 = Country_Archaea, axis2 = Sample_type_Archaea, axis3 = Family_Archaea, axis4 = ARG, 
           axis5 = ARG_host_PanRes, axis6 = Sample_type_PanRes, axis7 = Country_PanRes, y = Count)) +
  # 流动层：按Mobile_ARG填充颜色
  geom_alluvium(aes(fill = ARG), width = 1/24, alpha = 0.65) +
  # 节点层：按group_colors填充颜色
  geom_stratum(
    width = 1/24,
    fill = NA,      # 不填充颜色
    color = "black",     # 保留边框（可选）
    size = 0.6
  ) +
  scale_fill_manual(
    values = group_colors,
    guide = "none"
  ) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum),
                hjust = after_stat(case_when(
                  x < 4 ~ 1,      # 左边3列：右对齐
                  x == 4 ~ 0.5,   # 中间列：居中
                  TRUE ~ 0        # 右边3列：左对齐
                )),
                nudge_x = after_stat(case_when(
                  x < 4 ~ -0.05,  # 左边：向左偏移
                  x == 4 ~ 0,     # 中间：不偏移
                  TRUE ~ 0.05     # 右边：向右偏移
                ))),
            size = 4,
            color = "black",  
            fontface = "bold",
            angle = 0,
            vjust = 0.5,
            check_overlap = TRUE) +
  scale_x_discrete(limits = c("Country-Archaea", "Sample type-Archaea", "Family-Archaea","ARG", 
                              "ARG host-PanRes","Sample type-PanRes", "Country-PanRes"), expansion(add = c(2, 2))) + # 安排列的顺序并左右各扩展1个单位空间
  coord_cartesian(clip = "off") + # 允许绘图溢出面板
  theme_minimal() +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),  
    axis.text.x = element_text(size = 15, face = "bold", color = "black", vjust = 5), 
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
  )
    
print(p)

ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Co-related_MAG.png",       # 文件名
  plot = p + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),# 图形对象
  width = 30,        # 宽度（英寸）
  height = 12,       # 高度（英寸）
  dpi = 300,        # 分辨率（默认 300，可调高至 600+）
  units = "in"      # 单位（英寸）
)

