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
colors_biosample <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", 
                      "#8491B4", "#91D1C2", "#7E6148", "#B09C85",
                      "#F39C12", "#5F559B", "#B24745", "#0073C2", "#EFC000")

# Phyla Palette
colors_phyla <- c("#E64B35", "#FFD700", "#229954", "#4DBBD5", "#3C5488", 
                  "#F39B7F", "#8491B4", "#91D1C2", "#B09C85", "#631879", "#EF8A47")

# Species Palette (Enhanced 34 Colors)
colors_species <- c(
  # Red Tones (5 colors)
  "#B2182B",      # Deep Brick Red
  "#D6604D",      # Warm Red
  "#E41A1C",      # Bright Red
  "#F4A582",      # Light Coral
  "#FDDBC7",      # Very Light Pink
  
  # Orange Tones (4 colors)
  "#E6550D",      # Bright Orange
  "#FD8D3C",      # Orange-Yellow
  "#FEB24C",      # Light Orange
  "#FEE6CE",      # Cream Orange
  
  # Yellow Tones (4 colors)
  "#FFD700",      # Gold
  "#FED976",      # Pale Yellow
  "#F0E68C",      # Khaki Yellow
  "#BDB76B",      # Dark Khaki
  
  # Green Tones (5 colors)
  "#006837",      # Deep Green
  "#31A354",      # Medium Green
  "#78C679",      # Bright Green
  "#A6D854",      # Yellow-Green
  "#C2E699",      # Light Green
  
  # Blue Tones (3 colors)
  "#08519C",      # Deep Blue
  "#3182BD",      # Medium Blue
  "#9ECAE1",      # Light Blue
  
  # Purple Tones (4 colors)
  "#54278F",      # Deep Purple
  "#756BB1",      # Medium Purple
  "#9E9AC8",      # Light Purple
  "#DADAEB",      # Very Light Purple
  
  # Cyan/Teal Tones (3 colors)
  "#008080",      # Deep Teal
  "#20B2AA",      # Light Sea Green
  "#7FFFD4"       # Aquamarine
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
  facet_wrap(~ ID, ncol = 16, strip.position = "bottom") + 
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
  facet_wrap(~ ID, ncol = 16, strip.position = "bottom") + 
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
  facet_wrap(~ ID, ncol = 16, strip.position = "bottom") + 
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
