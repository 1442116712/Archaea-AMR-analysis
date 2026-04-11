# ==========================================
# 1. Load Required Libraries
# ==========================================
library(readxl)
library(tidyverse)
library(ggalluvial)

# ==========================================
# 2. Global Settings and Color Palette
# ==========================================
file_path <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx"
output_dir <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/"

group_colors <- c(
  "aminoglycosides_ANT6-Ia"   = "#1F77B4", "aminoglycosides_str"       = "#38A6D8",  
  "aminoglycosides_ANT9-Ia"   = "#4FA3D1", "polymyxin"                 = "#4B4E9E", 
  "beta_lactams_amp"          = "#D62728", "beta_lactams_TEM"          = "#FF6B6B", 
  "tetracyclines_tet44"       = "#2CA02C", "tetracyclines_tetC"        = "#66C266", 
  "phenicols_catA"            = "#138D75", "MLS_mefA"                  = "#9467BD",  
  "MLS_mefEN2"                = "#6A3D9A", "lincosamides_lnuAN2"       = "#B07CC6", 
  "lincosamides_lnuF"         = "#E377C2", "rifampicin_rph"            = "#FF8C00",  
  "elfamycins_tufab"          = "#D4AC0D", "multi_biocide_hsmR"        = "#1F3A93",
  "copper_resistance_copR"    = "#B03A2E", "copper_resistance_copB"    = "#7B3F00",
  "copper_resistance_dnaK"    = "#A0522D", "copper_resistance_tcrB"    = "#8B5A2B",
  "copper_resistance_tcrY"    = "#CD853F", "copper_resistance_tcrZ"    = "#C39BD3",
  "multi_metal_copAM"         = "#34495E", "multi_metal_wtpA"         = "#566573",
  "multi_metal_wtpB"          = "#85929E", "multi_metal_wtpC"         = "#BFC9CA"        
)

# ==========================================
# 3. Define Plotting Function
# ==========================================
generate_sankey <- function(data_input, save_name) {
  
  # Data aggregation
  sankey_data <- data_input %>%
    group_by(Country_Archaea, Sample_type_Archaea, Class_Archaea, ARG, 
             ARG_host_PanRes, Sample_type_PanRes, Country_PanRes) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Create plot
  p <- ggplot(sankey_data,
              aes(axis1 = Country_Archaea, axis2 = Sample_type_Archaea, axis3 = Class_Archaea, axis4 = ARG, 
                  axis5 = ARG_host_PanRes, axis6 = Sample_type_PanRes, axis7 = Country_PanRes, y = Count)) +
    geom_alluvium(aes(fill = ARG), width = 1/24, alpha = 0.65) +
    geom_stratum(width = 1/24, fill = NA, color = "black", size = 1) +
    scale_fill_manual(values = group_colors, guide = "none") +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum),
                  hjust = after_stat(case_when(x < 4 ~ 1, x == 4 ~ 0.5, TRUE ~ 0)),
                  nudge_x = after_stat(case_when(x < 4 ~ -0.05, x == 4 ~ 0, TRUE ~ 0.05))),
              size = 5, color = "black", fontface = "bold", check_overlap = TRUE) +
    scale_x_discrete(limits = c("Country-Archaea", "Sample type-Archaea", "Class-Archaea", "ARG", 
                                "ARG host-PanRes", "Sample type-PanRes", "Country-PanRes"), 
                     expansion(add = c(2, 2))) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),  
      axis.text.x = element_text(size = 18, face = "bold", color = "black", vjust = 5), 
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Save high-resolution figure
  ggsave(filename = paste0(output_dir, save_name), plot = p, 
         width = 20, height = 7, dpi = 300, units = "in")
}

# ==========================================
# 4. Process Sheets and Output Figures
# ==========================================

# Process Antibiotic Sheet
data_antibiotic <- read_excel(file_path, sheet = "Co-related_for_antibiotic")
generate_sankey(data_antibiotic, "Sankey_Antibiotic_ARG.png")

# Process Biocide/Metal Sheet
data_biocide <- read_excel(file_path, sheet = "Co-related_for_biocide_metal")
generate_sankey(data_biocide, "Sankey_Biocide_Metal_ARG.png")

message("Processing Complete. Plots saved to: ", output_dir)

