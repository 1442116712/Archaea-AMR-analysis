library(ggplot2)
library(tidyr)
library(readxl)
library(openxlsx)

data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Archaea_MGE_nt_prok.xlsx', sheet = "Phyla")

# 转换为长格式
long_data <- data %>%
  pivot_longer(cols = -ID, names_to = "Category", values_to = "Count")

library(dplyr)
long_data_percent <- long_data %>%
  group_by(ID) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()


extended_colors <- c( "#229954","#3C5488",  "#B09C85", "#EF8A47")
#绘制百分比条形图
p <- ggplot(long_data_percent, aes(x = factor(1), y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = extended_colors) +
  geom_text(aes(label = ifelse(Percentage > 5, paste0(round(Percentage), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black", fontface = "plain") +
  facet_wrap(~ ID, ncol = 9, strip.position = "bottom") +  # 按ID分面
  labs(x = "Accession_MGE_ARG", y = "Percentage (%)", 
       title = "Phyla distribution by Accession_MGE_ARG",
       fill = "Phyla") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 隐藏x轴文本，因为每个分面只有一个条形
        axis.ticks.x = element_blank(),  # 隐藏x轴刻度
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),  # x轴标题大小和样式
        axis.title.y = element_text(size = 15, face = "bold"),  # y轴标题大小和样式
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 15, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size =15, face = "bold"),
        legend.key.size = unit(0.8, "cm"),                    # 图例键（颜色方块）大小
        legend.spacing.y = unit(0.3, "cm"),                   # 图例项之间的垂直间距
        strip.text.x = element_text(
          size = 15, 
          face = "bold",
          angle = 90,           # 倾斜45度
          hjust = 1,          # 水平居中
          vjust = 0.5           # 垂直上对齐
        ),
        strip.placement = "outside",  # 将标签放在绘图区域外部
        panel.spacing = unit(0.1, "lines")) +  # 调整分面间距
  guides(fill = guide_legend(ncol = 1))

#绘制堆叠条形图
p <- ggplot(long_data_percent, aes(x = factor(1), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = extended_colors) +
  geom_text(aes(label = ifelse(Count > 10, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black", fontface = "plain") +
  facet_wrap(~ ID, ncol = 9, strip.position = "bottom") +  # 按ID分面
  labs(x = "Species_MGE_ARG", y = "Homogenous sequence counts", 
       title = "Phyla distribution by Species_MGE_ARG",
       fill = "Phyla") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 隐藏x轴文本，因为每个分面只有一个条形
        axis.ticks.x = element_blank(),  # 隐藏x轴刻度
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),  # x轴标题大小和样式
        axis.title.y = element_text(size = 15, face = "bold"),  # y轴标题大小和样式
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_text(size = 15, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size =15, face = "bold"),
        legend.key.size = unit(0.8, "cm"),                    # 图例键（颜色方块）大小
        legend.spacing.y = unit(0.3, "cm"),                   # 图例项之间的垂直间距
        strip.text.x = element_text(
          size = 15, 
          face = "bold",
          angle = 90,           # 倾斜45度
          hjust = 1,          # 水平居中
          vjust = 0.5           # 垂直上对齐
        ),
        strip.placement = "outside",  # 将标签放在绘图区域外部
        panel.spacing = unit(0.1, "lines")) +  # 调整分面间距
  guides(fill = guide_legend(ncol = 1))

p

ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/plot.png",       # 文件名
  plot = p,         # 图形对象
  width = 7,        # 宽度（英寸）
  height = 12,       # 高度（英寸）
  dpi = 300,        # 分辨率（默认 300，可调高至 600+）
  units = "in",
  bg = "white" # 单位（英寸）
)





