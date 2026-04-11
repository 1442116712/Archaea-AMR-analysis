# ==========================================
# 1. Load required R packages
# ==========================================
library(ggplot2)
library(tidyr)
library(readxl)
library(openxlsx)
library(dplyr)
library(ggtext)
library(patchwork)
library(grid)

# ==========================================
# 2. Load data
# ==========================================
file_path <- 'C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx'
data1 <- read_excel(file_path, sheet = "All_MGE_summary")
data2 <- read_excel(file_path, sheet = "ARG_MGE_summary")

# ==========================================
# 3. Data preprocessing (Convert to long format & handle italics)
# ==========================================
# Convert to long format
long_data_percent1 <- data1 %>% pivot_longer(cols = -MGE_type, names_to = "Category", values_to = "Count")
long_data_percent2 <- data2 %>% pivot_longer(cols = -MGE_type, names_to = "Category", values_to = "Count")

# List of taxonomic categories that require italicization
italic_list <- c("Archaeoglobi", "Halobacteria", "Iainarchaeia", "Methanobacteria", 
                 "Methanomicrobia", "Methanosarcinia", "Thermococci", 
                 "Thermoplasmata", "Thermoprotei_A")

# Add <i> tags to matching categories for element_markdown processing
add_italics <- function(x) {
  ifelse(x %in% italic_list, paste0("<i>", x, "</i>"), as.character(x))
}

long_data_percent1 <- long_data_percent1 %>% mutate(MGE_type_plot = add_italics(MGE_type))
long_data_percent2 <- long_data_percent2 %>% mutate(MGE_type_plot = add_italics(MGE_type))

# ==========================================
# 4. Reorder factor levels (Lock plotting order)
# ==========================================
original_levels1 <- unique(data1$MGE_type)
original_levels2 <- unique(data2$MGE_type)

new_levels1 <- add_italics(original_levels1)
new_levels2 <- add_italics(original_levels2)

long_data_percent1$MGE_type_plot <- factor(long_data_percent1$MGE_type_plot, levels = new_levels1)
long_data_percent2$MGE_type_plot <- factor(long_data_percent2$MGE_type_plot, levels = new_levels2)

# ==========================================
# 5. Calculate custom subplot spacing (Widen gap after SCGs)
# ==========================================
total_columns <- length(levels(long_data_percent1$MGE_type_plot))
my_gaps <- rep(0.1, total_columns - 1) # Default gap is 0.1

# Locate the SCGs column and set the gap to its right to 2.0
scgs_col <- which(levels(long_data_percent1$MGE_type_plot) == "SCGs")
if(length(scgs_col) > 0) {
  my_gaps[scgs_col] <- 2.0 
}

# ==========================================
# 6. Define shared theme settings
# ==========================================
my_colors <- c("#E64B35", "#4DBBD5", "#8491B4", "#91D1C2", "#B09C85")

shared_theme <- theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    strip.text.x = element_markdown(size = 15, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
    strip.placement = "outside",
    panel.spacing.y = unit(0.1, "lines"),
    panel.spacing.x = unit(my_gaps, "lines"), # Apply custom gaps
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# ==========================================
# 7. Create plots (y = NULL to remove individual axis titles)
# ==========================================
# Plot P1
p1 <- ggplot(long_data_percent1, aes(x = factor(1), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = my_colors) +
  geom_text(aes(label = ifelse(Count > 20000, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black") +
  facet_wrap(~ MGE_type_plot, ncol = 16, strip.position = "bottom") +  
  # Set y = NULL to allow for a single combined Y-axis label later
  labs(x = "", y = NULL, fill = "MGE types", 
       title = "12,477 Archaeal Genomes") +
  shared_theme +
  guides(fill = guide_legend(ncol = 1))

# Plot P2
p2 <- ggplot(long_data_percent2, aes(x = factor(1), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = my_colors) +
  geom_text(aes(label = ifelse(Count > 20000, paste0(round(Count), ""), "")),
            position = position_stack(vjust = 0.5), size = 3.5, color = "black") +
  facet_wrap(~ MGE_type_plot, ncol = 16, strip.position = "bottom") +  
  labs(x = "", y = NULL, fill = "MGE types", 
       title = "125 ARG-Carrying Genomes") +
  shared_theme +
  guides(fill = guide_legend(ncol = 1))

# ==========================================
# 8. Combine plots, add centered Y-axis title, and save
# ==========================================
# Remove p1 bottom facet labels to make the layout more compact
p1_clean <- p1 + theme(strip.text.x = element_blank())

# Stack the two plots vertically
plots_combined <- p1_clean / p2 + plot_layout(guides = "collect")

# Create a global Y-axis text object (rot = 90 for vertical orientation)
y_title <- textGrob(
  "MGE-associated protein counts", 
  rot = 90, 
  gp = gpar(fontsize = 15, fontface = "bold")
)

# Use "|" to place the Y-axis title on the left and the combined plots on the right
# widths = c(1, 35) allocates the proportional space for text vs graphics
final_plot <- wrap_elements(panel = y_title) | plots_combined
final_plot <- final_plot + plot_layout(widths = c(1, 35))

# Preview the final result
print(final_plot)

# Save the high-resolution image
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/dataset_MGE.png", 
  plot = final_plot,          
  width = 10,        
  height = 16,       
  dpi = 300,          
  units = "in",
  bg = "white"
)
