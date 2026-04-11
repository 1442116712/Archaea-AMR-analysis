# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(readxl)
library(openxlsx)
library(ggplot2)
library(gggenes)
library(ggrepel)
library(dplyr)

# Global output path
output_dir <- 'C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/'
input_file <- 'C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx'

# ==============================================================================
# 2. LINEAR GENE MAPS
# ==============================================================================

# ---- Methanocorpusculum petauri (JAQATW010000088.1) ----
data_88 <- read_excel(input_file, sheet = 'JAQATW010000088.1')

restr_88 <- data.frame(
  enzyme = c("EcoRI", "HindIII", "HincII", "ClaI", "BglI", "EaeI", "BmgBI"),
  position = c(2166, 3055, 499, 2300, 1121, 327, 4112),        
  color = c("red", "blue", "darkgreen", "purple", "brown", "orange", "black") 
)

p88 <- ggplot(data_88, aes(xmin = start, xmax = end, y = "JAQATW010000088.1", 
                           fill = type, forward = (strand == "Forward"))) +
  geom_gene_arrow() +
  scale_fill_manual(values = c("Structural gene" = "#9932CC", "Regulatory gene" = "#00CED1", 
                               "Mobile element" = "#EE82EE", "Resistance gene" = "red", 
                               "Recombination gene" = "#3CB371", "Metabolic gene" = "#1E90FF",
                               "Transport gene" = "#FFD700", "Target site" = "#C0C0C0")) +
  geom_text(aes(x = (start + end) / 2, label = gene), vjust = 2.5, size = 3, fontface = "bold") +
  geom_text_repel(data = restr_88, aes(x = position, y = 1.02, label = paste0(enzyme, "(", position, ")"), color = enzyme),
                  inherit.aes = FALSE, direction = "y", nudge_y = 0.55, segment.size = 0.6, fontface = "bold", size = 3) +
  scale_color_manual(values = setNames(restr_88$color, restr_88$enzyme)) +
  coord_cartesian(ylim = c(5, 1)) +
  theme_genes() +
  theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA))

ggsave(paste0(output_dir, "JAQATW010000088.1.png"), plot = p88, width = 8, height = 5, dpi = 300)

# ---- NBU2 ----
data_nbu2 <- read_excel(input_file, sheet = 'NBU2')

restr_nbu2 <- data.frame(
  enzyme = c("EcoRI(135)", "HindIII(10418)", "HincII(1803)", "ClaI(1)", "BglI(1180)", "EaeI(1973)", "BmgBI(9362)","SmaI(9029)"),
  position = c(2166, 3055, 499, 2300, 1121, 327, 4112, 4445),        
  color = c("red", "blue", "darkgreen", "purple", "brown", "orange", "black", "#FF1493") 
)

pnbu2 <- ggplot(data_nbu2, aes(xmin = start, xmax = end, y = "NBU2", fill = type)) +
  geom_gene_arrow(data = subset(data_nbu2, type != "Target_site"), aes(forward = (strand == "Forward"))) +
  geom_gene_arrow(data = subset(data_nbu2, type == "Target_site"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(3, "mm")) +
  scale_fill_manual(values = c("Structural gene" = "#9932CC", "Regulatory gene" = "#00CED1", 
                               "Mobile element" = "#EE82EE", "Resistance gene" = "red", 
                               "Recombination gene" = "#3CB371", "Metabolic gene" = "#1E90FF",
                               "oriT" = "#FFD700", "Target site" = "#C0C0C0")) +
  geom_text(aes(x = (start + end) / 2, label = gene), vjust = 2.5, size = 3, fontface = "bold") +
  geom_text_repel(data = restr_nbu2, aes(x = position, y = 1.02, label = enzyme, color = enzyme),
                  inherit.aes = FALSE, direction = "y", nudge_y = 0.55, segment.size = 0.6, fontface = "bold", size = 3) +
  scale_color_manual(values = setNames(restr_nbu2$color, restr_nbu2$enzyme)) +
  coord_cartesian(ylim = c(5, 1)) +
  theme_genes() +
  theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA))

ggsave(paste0(output_dir, "NBU2.png"), plot = pnbu2, width = 8, height = 5, dpi = 300)

# ---- Methanosarcina mazei (NZ_JJQP01000060.1) ----
data_mazei <- read_excel(input_file, sheet = 'NZ_JJQP01000060.1')

p_mazei <- ggplot(data_mazei, aes(xmin = start, xmax = end, y = "NZ_JJQP01000060.1", 
                                  fill = type, forward = (strand == "Forward"))) +
  geom_gene_arrow() +
  scale_fill_manual(values = c("Structural gene" = "#9932CC", "Regulatory gene" = "#00CED1", 
                               "Mobile element" = "#EE82EE", "Resistance gene" = "red", 
                               "Recombination gene" = "#3CB371", "Metabolic gene" = "#1E90FF",
                               "Transport gene" = "#FFD700", "Target site" = "#C0C0C0")) +
  geom_text(aes(x = (start + end) / 2, label = gene), vjust = 2.5, size = 3, fontface = "bold") +
  theme_genes() +
  theme(legend.position = "none", plot.margin = margin(250, 10, -5, 10, unit = "pt"))

ggsave(paste0(output_dir, "NZ_JJQP01000060.1.png"), plot = p_mazei, width = 5, height = 5, dpi = 300)


# ==============================================================================
# 3. CIRCULAR PLASMID MAPS
# ==============================================================================

# Function to generate circular plasmid plots
plot_circular_plasmid <- function(sheet_name, restr_data = NULL) {
  data <- read_excel(input_file, sheet = sheet_name)
  max_pos <- max(data$end)
  
  p <- ggplot(data, aes(xmin = start, xmax = end, y = 1, fill = type)) +
    # Backbone
    geom_rect(aes(xmin = 0, xmax = max_pos, ymin = 0.98, ymax = 1.02), fill = "#C0C0C0", color = NA) +
    # Inner hole
    geom_rect(aes(xmin = 0, xmax = max_pos, ymin = 0.0001, ymax = 0.9), fill = "white", color = NA) +
    # Genes
    geom_gene_arrow(data = subset(data, type != "oriT"), aes(forward = (strand == "Forward"))) +
    geom_gene_arrow(data = subset(data, type == "oriT"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(3, "mm")) +
    scale_fill_manual(values = c("IS" = "#9932CC", "T4SS" = "#00CED1", "Transposon" = "#EE82EE", 
                                 "Resistance gene" = "red", "Efflux pump" = "#3CB371", 
                                 "Relaxase" = "#1E90FF", "oriT" = "#FFD700")) +
    coord_polar(theta = "x") +
    # Gene Labels
    geom_text_repel(aes(x = (start + end) / 2, label = gene), size = 3, fontface = "bold", max.overlaps = Inf) +
    # Scale Marks (every 10kb)
    geom_segment(data = data.frame(pos = seq(0, max_pos, by = 10000)), 
                 aes(x = pos, xend = pos, y = 1.15, yend = 1.2), color = "black", linewidth = 0.3, inherit.aes = FALSE) +
    geom_text(data = data.frame(pos = seq(0, max_pos, by = 10000)),
              aes(x = pos, y = 1.25, label = ifelse(pos == 0, "0", paste0(pos/1000, "kb"))), size = 3, inherit.aes = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA))

  # Add Restriction Sites if provided
  if (!is.null(restr_data)) {
    p <- p + geom_text_repel(data = restr_data, aes(x = position, y = 1, label = enzyme),
                             inherit.aes = FALSE, nudge_y = -0.2, segment.size = 0.6, 
                             segment.alpha = 0.5, fontface = "bold", size = 3)
  }
  return(p)
}

# ---- Plotting and Saving Plasmids ----

# 2964TF
restr_2964 <- data.frame(enzyme = c("BciVI(80030)", "Ecil(80129)"), position = c(80030, 80129))
p_2964 <- plot_circular_plasmid('2964TF', restr_2964)
ggsave(paste0(output_dir, "plasmid_2964TF.png"), plot = p_2964, width = 8, height = 8)

# p142_A-OXA181
p_p142 <- plot_circular_plasmid('p142_A-OXA181')
ggsave(paste0(output_dir, "plasmid_p142_A-OXA181.png"), plot = p_p142, width = 8, height = 8)

# PU25001
restr_PU25001 <- data.frame(enzyme = c("BciVI(6506)", "Ecil(6407)"), position = c(6506, 6407))
p_PU25001 <- plot_circular_plasmid('PU25001', restr_PU25001)
ggsave(paste0(output_dir, "plasmid_PU25001.png"), plot = p_PU25001, width = 8, height = 8)

# pNDM15-1091
restr_pNDM15 <- data.frame(enzyme = c("BciVI(56169)", "Ecil(56268)"), position = c(56169, 56268))
p_pNDM15 <- plot_circular_plasmid('pNDM15-1091', restr_pNDM15)
ggsave(paste0(output_dir, "plasmid_pNDM15-1091.png"), plot = p_pNDM15, width = 8, height = 8)
