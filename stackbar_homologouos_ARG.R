# ==========================================
# 1. Load Required Libraries
# ==========================================
library(ggplot2)
library(tidyr)
library(readxl)
library(openxlsx)
library(dplyr)

# ==========================================
# 2. Define Common File Paths and Color Palettes
# ==========================================
input_path <- 'C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Archaea_ARG_nt_prok.xlsx'
output_dir <- 'C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/'

# Biosample Palette
colors_biosample <- c("#F1A7A1", "#F4A582", "#F5D06D", "#A6D854", "#94C5CC", 
                      "#00CED1", "#EED0C1", "#457B9D", "#8491B4", "#D4A5C5", 
                      "#E6A8D7", "#C3B091", "#808000")

# Phyla Palette
colors_phyla <- c("#E64B35", "#FFD700", "#229954", "#4DBBD5", "#3C5488", 
                  "#F39B7F", "#8491B4", "#91D1C2", "#B09C85", "#631879", "#EF8A47")

# Species Palette (Enhanced 34 Colors)
colors_species <- c(
  "#330000", "#8B0000", "#E74C3C", "#FF6B6B", "#FF1493", "#FFB6C1", # Reds
  "#D84315", "#FF4500", "#FF8C00", "#DAA520", "#FFD700",           # Oranges
  "#F9A825", "#F0E68C", "#BDB76B", "#ADFF2F",                     # Yellows
  "#9ACD32", "#32CD32", "#2E8B57", "#008080", "#20B2AA",           # Greens
  "#E3F2FD", "#87CEEB", "#1E90FF", "#4169E1", "#0000CD", "#00008B", # Blues
  "#8A2BE2", "#9370DB", "#DA70D6", "#C71585",                     # Purples
  "#4B0082", "#8B008B", "#00CED1", "#A0522D"                      # Specials
)

# ==========================================
# 3. Generate Biosample Plot (Counts)
# ==========================================
data_bio <- read_excel(input_path, sheet = "Biosample")
long_bio <- data_bio %>%
  pivot_longer(cols = -ID, names_to = "Category", values_to = "Count") %>%
  mutate(ID = factor(ID, levels = unique(data_bio$ID)))

p_biosample <- ggplot(long_bio, aes(x = factor(1), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_biosample) +
  geom_text(aes(label = ifelse(Count > 10, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  facet_wrap(~ ID, ncol = 10, strip.position = "bottom") + 
  labs(x = "ARG", y = "Homogenous ARG counts", title = "Biosample types by ARG", fill = "Biosample type") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold", angle = 90, vjust = 0.5, hjust = 1),
        strip.placement = "outside")

# ==========================================
# 4. Generate Phyla Plot (Counts)
# ==========================================
data_phyla <- read_excel(input_path, sheet = "Phyla")
long_phyla <- data_phyla %>%
  pivot_longer(cols = -ID, names_to = "Category", values_to = "Count") %>%
  mutate(ID = factor(ID, levels = unique(data_phyla$ID)))

p_phyla <- ggplot(long_phyla, aes(x = factor(1), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_phyla) +
  geom_text(aes(label = ifelse(Count > 10, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  facet_wrap(~ ID, ncol = 9, strip.position = "bottom") + 
  labs(x = "ARG", y = "Homogenous ARG counts", title = "Phyla distribution by ARG", fill = "Phyla") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold", angle = 90, vjust = 0.5, hjust = 1),
        strip.placement = "outside")

# ==========================================
# 5. Generate Species Plot (Percentage)
# ==========================================
data_species <- read_excel(input_path, sheet = "Species")
long_species <- data_species %>%
  pivot_longer(cols = -ID, names_to = "Category", values_to = "Count") %>%
  group_by(ID) %>%
  mutate(Percentage = Count / sum(Count) * 100,
         ID = factor(ID, levels = unique(data_species$ID))) %>%
  ungroup()

p_species <- ggplot(long_species, aes(x = factor(1), y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_species) +
  geom_text(aes(label = ifelse(Count > 5, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  facet_wrap(~ ID, ncol = 9, strip.position = "bottom") + 
  labs(x = "ARG", y = "Percentage (%)", title = "Species distribution by ARG", fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold", angle = 90, vjust = 0.5, hjust = 1),
        strip.placement = "outside")

# ==========================================
# 6. Save All Plots
# ==========================================
ggsave(paste0(output_dir, "plot_Biosample.png"), plot = p_biosample, width = 11, height = 10, dpi = 300, bg = "white")
ggsave(paste0(output_dir, "plot_Phyla.png"), plot = p_phyla, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(paste0(output_dir, "plot_Species.png"), plot = p_species, width = 14.5, height = 13, dpi = 300, bg = "white")
