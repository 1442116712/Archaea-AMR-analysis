# 读取数据（制表符或空格分隔）
df <- read.table("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/1.txt",
                 header = TRUE, sep = "\t", 
                 stringsAsFactors = FALSE,  # 确保字符串不转换为因子
                 check.names = FALSE)       # 保持列名原样，不加点或括号

# 转换为长表
library(tidyr)
df_long <- df %>%
  pivot_longer(cols = -ARG_ID, names_to = "Metric", values_to = "Value")

# 绘图
library(ggplot2)
ggplot(df_long, aes(x = ARG_ID, y = Value, color = Metric, group = Metric)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.margin = unit(c(1, 1, 1, 2), "lines")  # 上右下左的边距
  ) +
  labs(x = "ID", y = "Value", color = "Metric", title = "") +
  coord_cartesian(ylim = c(0, 1.5))

