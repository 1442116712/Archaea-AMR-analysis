# ==========================================
# 1. Load Required Libraries
# ==========================================
library(ggplot2)
library(tidyr)
library(readxl)
library(openxlsx)
library(dplyr)

# ==========================================
# 2. Data Loading and Preprocessing
# ==========================================
data <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Archaea_ARG_nt_prok.xlsx', sheet = "Species")

# Convert to long format and calculate percentages per ID
long_data_percent <- data %>%
  pivot_longer(cols = -ID, names_to = "Category", values_to = "Count") %>%
  group_by(ID) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# Set factor levels to maintain original data order
long_data_percent$ID <- factor(long_data_percent$ID, levels = unique(data$ID))

# ==========================================
# 3. Define Enhanced Color Palette (34 colors)
# ==========================================
colors_34_enhanced <- c(
  # Red tones
  "#330000", "#8B0000", "#E74C3C", "#FF6B6B", "#FF1493", "#FFB6C1",
  # Orange tones
  "#D84315", "#FF4500", "#FF8C00", "#DAA520", "#FFD700",
  # Yellow tones
  "#F9A825", "#F0E68C", "#BDB76B", "#ADFF2F",
  # Green tones
  "#9ACD32", "#32CD32", "#2E8B57", "#008080", "#20B2AA",
  # Blue tones
  "#E3F2FD", "#87CEEB", "#1E90FF", "#4169E1", "#0000CD", "#00008B",
  # Purple tones
  "#8A2BE2", "#9370DB", "#DA70D6", "#C71585",
  # Special tones
  "#4B0082", "#8B008B", "#00CED1", "#A0522D"
)

# ==========================================
# 4. Create Visualization
# ==========================================
p <- ggplot(long_data_percent, aes(x = factor(1), y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_34_enhanced) +
  # Add labels for counts greater than 5
  geom_text(aes(label = ifelse(Count > 5, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black", fontface = "plain") +
  facet_wrap(~ ID, ncol = 9, strip.position = "bottom") + 
  labs(x = "ARG", y = "Percentage (%)", 
       title = "Species distribution by ARG",
       fill = "Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
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
    panel.spacing = unit(0.1, "lines")
  ) +
  guides(fill = guide_legend(ncol = 1))

# Preview plot
print(p)

# ==========================================
# 5. Export Figure
# ==========================================
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/plot.png",
  plot = p,
  width = 14.5,
  height = 13,
  dpi = 300,
  units = "in",
  bg = "white"
)





