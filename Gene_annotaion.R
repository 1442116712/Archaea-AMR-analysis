library(readxl)
library(openxlsx)
library(ggplot2)
library(gggenes)
library(ggrepel)
library(dplyr)
# Image size : width-800, Height-400
# Methanocorpusculum petauri ####
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = 'JAQATW010000088.1')

restriction_sites <- data.frame(
  enzyme = c("EcoRI", "HindIII", "HincII", "ClaI", "BglI", "EaeI", "BmgBI"),
  position = c(2166, 3055, 499, 2300, 1121, 327, 4112),        
  color = c("red", "blue", "darkgreen", "purple", "brown", "orange", "black") 
)

p <- ggplot(data, aes(xmin = start, xmax = end, y = "JAQATW010000088.1", fill = type, forward = strand == "Forward")) +
  geom_gene_arrow() +
  scale_fill_manual(values = c("Structural gene" = "#9932CC", "Regulatory gene" = "#00CED1", "Mobile element" = "#EE82EE", 
                               "Resistance gene" = "red", "Recombination gene" = "#3CB371", "Metabolic gene" = "#1E90FF",
                               "Transport gene" = "#FFD700", "Target site" = "#C0C0C0")) +
  theme_genes() + 
  geom_text(aes(x = (start + end) / 2, y = "JAQATW010000088.1", label = gene), angle = 0, vjust = 2.5, hjust = 0.5, size = 3, fontface = "bold") +
  labs(fill = "Function Categories", x = "Gene Position") 

p

p_simple <- p +
  
  geom_text_repel(data = restriction_sites,
                  aes(x = position, y = 1.02,  # 从基因顶部开始
                      label = paste0(enzyme, "(", position, ")"),
                      color = enzyme),
                  inherit.aes = FALSE,
                  direction = "y",
                  nudge_y = 0.55,       # 总偏移量，从1.02到1.08
                  segment.size = 0.6,   # 连接线粗细
                  segment.alpha = 0.5,  # 连接线透明度
                  segment.color = "gray50",  # 连接线颜色
                  segment.linetype = "solid", # 实线
                  size = 3.0,
                  fontface = "bold",
                  max.overlaps = Inf,
                  show.legend = FALSE,
                  point.padding = 0) +  # 减少点周围的填充
  
  scale_color_manual(values = setNames(restriction_sites$color, 
                                       restriction_sites$enzyme)) +
  
  coord_cartesian(ylim = c(5, 1)) +
  
  labs(title = "",
       x = "Gene Position") +
  
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_simple)
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/JAQATW010000088.1.png",
  plot = p_simple + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),
  width = 8,
  height = 5,
  dpi = 300,
  units = "in"
)


# NBU2 ####
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = 'NBU2')
restriction_sites <- data.frame(
  enzyme = c("EcoRI(135)", "HindIII(10418)", "HincII(1803)", "ClaI(1)", "BglI(1180)", "EaeI(1973)", "BmgBI(9362)","SmaI(9029)"),
  position = c(2166, 3055, 499, 2300, 1121, 327, 4112, 4445),        
  color = c("red", "blue", "darkgreen", "purple", "brown", "orange", "black", "#FF1493") 
)

p <- ggplot(data, aes(xmin = start, xmax = end, y = "NBU2", fill = type)) +
  geom_gene_arrow(
    data = subset(data, type != "Target_site"), aes(forward = strand == "Forward")
    ) +
  geom_gene_arrow(
    data = subset(data, type == "Target_site"),
    arrowhead_width = unit(0, "mm"),  
    arrow_body_height = unit(3, "mm")  
  ) +
  scale_fill_manual(values = c("Structural gene" = "#9932CC", "Regulatory gene" = "#00CED1", "Mobile element" = "#EE82EE", 
                               "Resistance gene" = "red", "Recombination gene" = "#3CB371", "Metabolic gene" = "#1E90FF",
                               "oriT" = "#FFD700", "Target site" = "#C0C0C0")) +
  theme_genes() + 
  geom_text(aes(x = (start + end) / 2, y = "NBU2", label = gene), angle = 0, vjust = 2.5, hjust = 0.5, size = 3, fontface = "bold") +
  labs(fill = "Function Categories", x = "Gene Position") 

p

p_simple <- p +
  
  geom_text_repel(data = restriction_sites,
                  aes(x = position, y = 1.02,  # 从基因顶部开始
                      label = paste0(enzyme),
                      color = enzyme),
                  inherit.aes = FALSE,
                  direction = "y",
                  nudge_y = 0.55,       # 总偏移量，从1.02到1.08
                  segment.size = 0.6,   # 连接线粗细
                  segment.alpha = 0.5,  # 连接线透明度
                  segment.color = "gray50",  # 连接线颜色
                  segment.linetype = "solid", # 实线
                  size = 3.0,
                  fontface = "bold",
                  max.overlaps = Inf,
                  show.legend = FALSE,
                  point.padding = 0) +  # 减少点周围的填充
  
  scale_color_manual(values = setNames(restriction_sites$color, 
                                       restriction_sites$enzyme)) +
  
  coord_cartesian(ylim = c(5, 1)) +
  
  labs(title = "",
       x = "Gene Position") +
  
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_simple)
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/NBU2.png",
  plot = p_simple + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),
  width = 8,
  height = 5,
  dpi = 300,
  units = "in"
)




# Methanosarcina mazei 1.H.A.2.8 ####
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = 'NZ_JJQP01000060.1')

p <- ggplot(data, aes(xmin = start, xmax = end, y = "NZ_JJQP01000060.1", fill = type, forward = strand == "Forward")) +
  geom_gene_arrow() +
  scale_fill_manual(values = c("Structural gene" = "#9932CC", "Regulatory gene" = "#00CED1", "Mobile element" = "#EE82EE", 
                               "Resistance gene" = "red", "Recombination gene" = "#3CB371", "Metabolic gene" = "#1E90FF",
                               "Transport gene" = "#FFD700", "Target site" = "#C0C0C0")) +
  theme_genes() + 
  geom_text(aes(x = (start + end) / 2, y = "NZ_JJQP01000060.1", label = gene), angle = 0, vjust = 2.5, hjust = 0.5, size = 3, fontface = "bold") +
  labs(fill = "Function Categories", x = "Gene Position") +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        plot.margin = margin(250, 10, -5, 10, unit = "pt"))

p
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/NZ_JJQP01000060.1.png",
  plot = p + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),
  width = 5,
  height = 5,
  dpi = 300,
  units = "in"
)




# plasmid 2964TF ####
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = '2964TF')

p <- ggplot(data, aes(xmin = start, xmax = end, y = 1, fill = type, forward = strand == "Forward")) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.98, ymax = 1.02), fill = "#C0C0C0", color = NA) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.0001, ymax = 0.9), fill = "white", color = NA) +
  geom_gene_arrow(data = subset(data, type != "oriT"), aes(forward = strand == "Forward")
  ) +
  geom_gene_arrow(
    data = subset(data, type == "oriT"),
    arrowhead_width = unit(0, "mm"),  
    arrow_body_height = unit(3, "mm")  
  ) +
  scale_fill_manual(values = c("IS" = "#9932CC", 
                               "T4SS" = "#00CED1", 
                               "Transposon" = "#EE82EE", 
                               "Resistance gene" = "red", 
                               "Efflux pump" = "#3CB371", 
                               "Relaxase" = "#1E90FF",
                               "oriT" = "#FFD700")) +
  theme_void() +  # geom_gene_label(aes(label = gene), align = "centre") +
  coord_polar(theta = "x") +
  geom_text_repel(aes(x = (start + end) / 2, y = 1, label = gene), 
                  size = 3, angle = 0, fontface = "bold", color = "black",
                  max.overlaps = Inf) +
  labs(fill = "Genetic Element") + 
  
  geom_segment(data = data.frame(pos = seq(0, max(data$end), by = 10000)),  # 每5000bp一个刻度
             aes(x = pos, xend = pos, y = 1.15, yend = 1.2),  # 刻度线位置
             color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  # 添加刻度标签
  geom_text(data = data.frame(pos = seq(0, max(data$end), by = 10000)),
            aes(x = pos, y = 1.25, label = ifelse(pos == 0, "0", paste0(pos/1000, "kb"))),
            size = 3, angle = 0, color = "black", inherit.aes = FALSE)

p
restriction_sites <- data.frame(
  enzyme = c("BciVI(80030)", "Ecil(80129)"),
  position = c(80030, 80129),        
  color = c("black","black") 
)

p_simple <- p +
  
  geom_text_repel(data = restriction_sites,
                  aes(x = position, y = 1,  # 从基因顶部开始
                      label = paste0(enzyme),
                      color = "black"),
                  inherit.aes = FALSE,
                  nudge_y = -0.2,       # 总偏移量，从1.02到1.08
                  segment.size = 0.6,   # 连接线粗细
                  segment.alpha = 0.5,  # 连接线透明度
                  segment.color = "gray50",  # 连接线颜色
                  segment.linetype = "solid", # 实线
                  size = 3.0,
                  fontface = "bold",
                  max.overlaps = Inf,
                  show.legend = FALSE) +  # 减少点周围的填充
  
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_simple)
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/plasmid_2964TF.png",
  plot = p_simple + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),
  width = 8,
  height = 8,
  dpi = 300,
  units = "in"
)







# plasmid p142_A-OXA181 ####
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = 'p142_A-OXA181')

p <- ggplot(data, aes(xmin = start, xmax = end, y = 1, fill = type, forward = strand == "Forward")) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.98, ymax = 1.02), fill = "#C0C0C0", color = NA) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.0001, ymax = 0.9), fill = "white", color = NA) +
  geom_gene_arrow(data = subset(data, type != "oriT"), aes(forward = strand == "Forward")
  ) +
  geom_gene_arrow(
    data = subset(data, type == "oriT"),
    arrowhead_width = unit(0, "mm"),  
    arrow_body_height = unit(3, "mm")  
  ) +
  scale_fill_manual(values = c("IS" = "#9932CC", 
                               "T4SS" = "#00CED1", 
                               "Transposon" = "#EE82EE", 
                               "Resistance gene" = "red", 
                               "Efflux pump" = "#3CB371", 
                               "Relaxase" = "#1E90FF",
                               "oriT" = "#FFD700")) +
  theme_void() +
  coord_polar(theta = "x") +
  geom_text_repel(aes(x = (start + end) / 2, y = 1, label = gene), 
                  size = 3, angle = 0, fontface = "bold", color = "black",
                  max.overlaps = Inf) +
  labs(fill = "Genetic Element") + 
  
  geom_segment(data = data.frame(pos = seq(0, max(data$end), by = 10000)),
               aes(x = pos, xend = pos, y = 1.15, yend = 1.2),
               color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  # 添加刻度标签
  geom_text(data = data.frame(pos = seq(0, max(data$end), by = 10000)),
            aes(x = pos, y = 1.25, label = ifelse(pos == 0, "0", paste0(pos/1000, "kb"))),
            size = 3, angle = 0, color = "black", inherit.aes = FALSE) +
  theme(legend.position = "none") 


print(p)
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/plasmid_p142_A-OXA181.png",
  plot = p + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),
  width = 8,
  height = 8,
  dpi = 300,
  units = "in"
)





# plasmid PU25001 ####
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = 'PU25001')

p <- ggplot(data, aes(xmin = start, xmax = end, y = 1, fill = type, forward = strand == "Forward")) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.98, ymax = 1.02), fill = "#C0C0C0", color = NA) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.0001, ymax = 0.9), fill = "white", color = NA) +
  geom_gene_arrow(data = subset(data, type != "oriT"), aes(forward = strand == "Forward")
  ) +
  geom_gene_arrow(
    data = subset(data, type == "oriT"),
    arrowhead_width = unit(0, "mm"),  
    arrow_body_height = unit(3, "mm")  
  ) +
  scale_fill_manual(values = c("IS" = "#9932CC", 
                               "T4SS" = "#00CED1", 
                               "Transposon" = "#EE82EE", 
                               "Resistance gene" = "red", 
                               "Efflux pump" = "#3CB371", 
                               "Relaxase" = "#1E90FF",
                               "oriT" = "#FFD700")) +
  theme_void() +  # geom_gene_label(aes(label = gene), align = "centre") +
  coord_polar(theta = "x") +
  geom_text_repel(aes(x = (start + end) / 2, y = 1, label = gene), 
                  size = 3, angle = 0, fontface = "bold", color = "black",
                  max.overlaps = Inf) +
  labs(fill = "Genetic Element") + 
  
  geom_segment(data = data.frame(pos = seq(0, max(data$end), by = 10000)),  # 每5000bp一个刻度
               aes(x = pos, xend = pos, y = 1.15, yend = 1.2),  # 刻度线位置
               color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  # 添加刻度标签
  geom_text(data = data.frame(pos = seq(0, max(data$end), by = 10000)),
            aes(x = pos, y = 1.25, label = ifelse(pos == 0, "0", paste0(pos/1000, "kb"))),
            size = 3, angle = 0, color = "black", inherit.aes = FALSE)

p
restriction_sites <- data.frame(
  enzyme = c("BciVI(6506)", "Ecil(6407)"),
  position = c(6506, 6407),        
  color = c("black","black") 
)

p_simple <- p +
  
  geom_text_repel(data = restriction_sites,
                  aes(x = position, y = 1,  # 从基因顶部开始
                      label = paste0(enzyme),
                      color = "black"),
                  inherit.aes = FALSE,
                  nudge_y = -0.2,       # 总偏移量，从1.02到1.08
                  segment.size = 0.6,   # 连接线粗细
                  segment.alpha = 0.5,  # 连接线透明度
                  segment.color = "gray50",  # 连接线颜色
                  segment.linetype = "solid", # 实线
                  size = 3.0,
                  fontface = "bold",
                  max.overlaps = Inf,
                  show.legend = FALSE) +  # 减少点周围的填充
  
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_simple)
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/plasmid_PU25001.png",
  plot = p_simple + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),
  width = 8,
  height = 8,
  dpi = 300,
  units = "in"
)




# plasmid pNDM15-1091 ####
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = 'pNDM15-1091')

p <- ggplot(data, aes(xmin = start, xmax = end, y = 1, fill = type, forward = strand == "Forward")) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.98, ymax = 1.02), fill = "#C0C0C0", color = NA) +
  geom_rect(aes(xmin = 0, xmax = max(end), ymin = 0.0001, ymax = 0.9), fill = "white", color = NA) +
  geom_gene_arrow(data = subset(data, type != "oriT"), aes(forward = strand == "Forward")
  ) +
  geom_gene_arrow(
    data = subset(data, type == "oriT"),
    arrowhead_width = unit(0, "mm"),  
    arrow_body_height = unit(3, "mm")  
  ) +
  scale_fill_manual(values = c("IS" = "#9932CC", 
                               "T4SS" = "#00CED1", 
                               "Transposon" = "#EE82EE", 
                               "Resistance gene" = "red", 
                               "Efflux pump" = "#3CB371", 
                               "Relaxase" = "#1E90FF",
                               "oriT" = "#FFD700")) +
  theme_void() +  # geom_gene_label(aes(label = gene), align = "centre") +
  coord_polar(theta = "x") +
  geom_text_repel(aes(x = (start + end) / 2, y = 1, label = gene), 
                  size = 3, angle = 0, fontface = "bold", color = "black",
                  max.overlaps = Inf) +
  labs(fill = "Genetic Element") + 
  
  geom_segment(data = data.frame(pos = seq(0, max(data$end), by = 10000)),  # 每5000bp一个刻度
               aes(x = pos, xend = pos, y = 1.15, yend = 1.2),  # 刻度线位置
               color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  # 添加刻度标签
  geom_text(data = data.frame(pos = seq(0, max(data$end), by = 10000)),
            aes(x = pos, y = 1.25, label = ifelse(pos == 0, "0", paste0(pos/1000, "kb"))),
            size = 3, angle = 0, color = "black", inherit.aes = FALSE)

p
restriction_sites <- data.frame(
  enzyme = c("BciVI(56169)", "Ecil(56268)"),
  position = c(56169, 56268),        
  color = c("black","black") 
)

p_simple <- p +
  
  geom_text_repel(data = restriction_sites,
                  aes(x = position, y = 1,  # 从基因顶部开始
                      label = paste0(enzyme),
                      color = "black"),
                  inherit.aes = FALSE,
                  nudge_y = -0.2,       # 总偏移量，从1.02到1.08
                  segment.size = 0.6,   # 连接线粗细
                  segment.alpha = 0.5,  # 连接线透明度
                  segment.color = "gray50",  # 连接线颜色
                  segment.linetype = "solid", # 实线
                  size = 3.0,
                  fontface = "bold",
                  max.overlaps = Inf,
                  show.legend = FALSE) +  # 减少点周围的填充
  
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_simple)
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/plasmid_pNDM15-1091.png",
  plot = p_simple + 
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # 白色背景，无边框
      plot.background = element_rect(fill = "white", color = NA),   # 图形背景，无边框
      panel.border = element_blank(),                               # 移除面板边框
      panel.spacing = unit(0, "cm")                                 # 移除面板间距
    ),
  width = 8,
  height = 8,
  dpi = 300,
  units = "in"
)
