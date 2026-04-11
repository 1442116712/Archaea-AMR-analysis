# ==========================================
# 1. Load Required Libraries
# ==========================================
library(ggplot2)
library(tidyr)
library(readxl)
library(openxlsx)
library(dplyr)

# ==========================================
# 2. Data Preparation
# ==========================================
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Archaea_ARG_nt_prok.xlsx', sheet = "Biosample")

# Convert to long format and calculate percentages
long_data_percent <- data %>%
  pivot_longer(cols = -ID, names_to = "Category", values_to = "Count") %>%
  group_by(ID) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# Set factor levels to maintain original data order
long_data_percent$ID <- factor(long_data_percent$ID, levels = unique(data$ID))

# ==========================================
# 3. Data Visualization
# ==========================================
p <- ggplot(long_data_percent, aes(x = factor(1), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = c("#F1A7A1", "#F4A582", "#F5D06D", "#A6D854", "#94C5CC", 
                               "#00CED1", "#EED0C1", "#457B9D", "#8E44AD", "#D4A5C5", 
                               "#E6A8D7", "#C3B091", "#808000")) +
  geom_text(aes(label = ifelse(Count > 10, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black", fontface = "plain") +
  facet_wrap(~ ID, ncol = 10, strip.position = "bottom") + 
  labs(x = "ARG", y = "Homogenous ARG counts", 
       title = "Biosample types by ARG",
       fill = "Biosample type") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 15, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size = 15, face = "bold"),
        legend.key.size = unit(0.8, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        strip.text.x = element_text(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
        strip.placement = "outside",
        panel.spacing = unit(0.1, "lines")) +
  guides(fill = guide_legend(ncol = 1))

# Display plot
print(p)

# ==========================================
# 4. Export Figure
# ==========================================
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/plot.png",
  plot = p,
  width = 11,
  height = 10,
  dpi = 300,
  units = "in",
  bg = "white"
)





