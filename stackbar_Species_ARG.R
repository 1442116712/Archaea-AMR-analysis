library(ggplot2)
library(tidyr)
library(readxl)
library(openxlsx)

data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Archaea_ARG_nt_prok.xlsx', sheet = "Species")

# 转换为长格式
long_data <- data %>%
  pivot_longer(cols = -ID, names_to = "Category", values_to = "Count")

library(dplyr)
long_data_percent <- long_data %>%
  group_by(ID) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 在绘图前设置ID的因子水平，保持原始顺序
long_data_percent$ID <- factor(long_data_percent$ID, 
                               levels = unique(data$ID))
# 按视觉光谱顺序重新排序（ROYGBIV + 中性色）
colors_30_enhanced <- c(
  # 红色系 (6色) - 从深到浅
  "#330000",      # 极深红 - 新增
  "#8B0000",      # 深红
  "#E74C3C",      # 亮红
  "#FF6B6B",      # 珊瑚红
  "#FF1493",      # 亮粉色
  "#FFB6C1",      # 浅粉
  
  # 橙色系 (5色) - 橙红到橙黄
  "#D84315",      # 深橙红 - 新增
  "#FF4500",      # 橙红
  "#FF8C00",      # 深橙
  "#DAA520",      # 标准橙
  "#FFD700",      # 金黄
  
  # 黄色系 (4色) - 金黄到黄绿
  "#F9A825",      # 琥珀色 - 新增
  "#F0E68C",      # 卡其黄
  "#BDB76B",      # 暗卡其黄
  "#ADFF2F",      # 绿黄
  
  # 绿色系 (5色) - 黄绿到蓝绿
  "#9ACD32",      # 黄绿
  "#32CD32",      # 酸橙绿
  "#2E8B57",      # 海绿
  "#008080",      # 青色
  "#20B2AA",      # 浅海绿
  
  # 蓝色系 (6色) - 浅蓝到深蓝
  "#E3F2FD",      # 极浅蓝 - 新增
  "#87CEEB",      # 天蓝
  "#1E90FF",      # 道奇蓝
  "#4169E1",      # 宝蓝
  "#0000CD",      # 中蓝
  "#00008B",      # 深蓝
  
  # 紫色系 (4色) - 蓝紫到红紫
  "#8A2BE2",      # 蓝紫
  "#9370DB",      # 中紫
  "#DA70D6",      # 兰紫
  "#C71585",      # 中紫红
  
  # 新增：特殊色系 (保持总数30)
  "#4B0082",      # 靛青 - 深紫
  "#8B008B",      # 深洋红
  "#00CED1",      # 深青色
  "#A0522D"       # 赭石色 - 大地色
)
#绘制百分比条形图
p <- ggplot(long_data_percent, aes(x = factor(1), y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_30_enhanced) +
  geom_text(aes(label = ifelse(Count > 5, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black", fontface = "plain") +
  facet_wrap(~ ID, ncol = 9, strip.position = "bottom") +  # 按ID分面
  labs(x = "ARG", y = "Percentage (%)", 
       title = "Species distribution by ARG",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 隐藏x轴文本，因为每个分面只有一个条形
        axis.ticks.x = element_blank(),  # 隐藏x轴刻度
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),  # x轴标题大小和样式
        axis.title.y = element_text(size = 15, face = "bold"),  # y轴标题大小和样式
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 15, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size =15, face = "bold"),
        legend.key.size = unit(0.8,"cm"),                    # 图例键（颜色方块）大小
        legend.spacing.y = unit(0.3, "cm"),                   # 图例项之间的垂直间距
        strip.text.x = element_text(
          size = 15, 
          face = "bold",
          angle = 90,           # 倾斜60度
          hjust = 1,          # 水平居中
          vjust = 0.5           # 垂直上对齐
        ),
        strip.placement = "outside",  # 将标签放在绘图区域外部
        panel.spacing = unit(0.1, "lines")) +  # 调整分面间距
  guides(fill = guide_legend(ncol = 1))
p
#绘制堆叠条形图
p <- ggplot(long_data_percent, aes(x = factor(1), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = c(colors_27_sequential)) +
  geom_text(aes(label = ifelse(Count > 10, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black", fontface = "plain") +
  facet_wrap(~ ID, ncol = 9, strip.position = "bottom") +  # 按ID分面
  labs(x = "ARG", y = "Homogenous ARG counts", 
       title = "Biosample types by ARG",
       fill = "Biosample type") +
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

p

ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/plot.png",       # 文件名
  plot = p,         # 图形对象
  width = 14.5,        # 宽度（英寸）
  height = 13,       # 高度（英寸）
  dpi = 300,        # 分辨率（默认 300，可调高至 600+）
  units = "in",
  bg = "white" # 单位（英寸）
)





