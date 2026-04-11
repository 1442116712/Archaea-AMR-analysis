# ==========================================
# 1. Load Required Libraries
# ==========================================
library(readxl)
library(tidyverse)
library(tidygraph)
library(ggraph)

# ==========================================
# 2. Data Loading and Preprocessing
# ==========================================
file_path <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Archaea_MGE_nt_prok.xlsx"
df <- read_excel(file_path, sheet = "co-related")

# Aggregate counts for nodes
df_count <- df %>%
  count(queryid, Species, Sample_type, Continent, name = "number")

# Define target IDs and their specific MGE-ARG colors
target_mge <- list(
  "Homologs of Tn_lunAN2_mefEN2" = "#FF6B6B",
  "Homologs of IS_copR"           = "#E377C2",
  "Homologs of plasmid_catA"      = "#9467BD"
)

# ==========================================
# 3. Define Plotting Function
# ==========================================
plot_tree <- function(data, id_name, mge_color) {
  
  # Filter data for specific ID
  sub <- data %>% filter(queryid == id_name)
  
  # Prepare node data
  continent_nodes <- sub %>%
    group_by(Continent) %>%
    summarise(number = sum(number)) %>%
    mutate(type = "Continent") %>%
    rename(name = Continent)
  
  sample_nodes <- sub %>%
    group_by(Sample_type) %>%
    summarise(number = sum(number)) %>%
    mutate(type = "Sample_type") %>%
    rename(name = Sample_type)
  
  species_nodes <- sub %>%
    group_by(Species) %>%
    summarise(number = sum(number)) %>%
    mutate(type = "Species") %>%
    rename(name = Species)
  
  mge_node <- tibble(
    name = id_name,
    number = sum(sub$number),
    type = "MGE-ARG"
  )
  
  nodes <- bind_rows(mge_node, species_nodes, sample_nodes, continent_nodes)
  
  # Construct edges (Tree hierarchy)
  edges <- bind_rows(
    sub %>% select(from = queryid, to = Species),
    sub %>% select(from = Species, to = Sample_type),
    sub %>% select(from = Sample_type, to = Continent)
  ) %>% distinct()
  
  # Create graph object
  graph <- tbl_graph(nodes = nodes, edges = edges)
  
  # Visualization
  ggraph(graph, layout = "tree") +
    geom_edge_link(colour = "grey75", alpha = 0.6) +
    geom_node_point(aes(size = number, colour = type), alpha = 0.9) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, fontface = "bold") + 
    scale_size_continuous(range = c(3, 15), name = "Count", guide = guide_legend(order = 2)) +
    scale_colour_manual(
      values = c(
        "MGE-ARG"     = mge_color, # Dynamic color based on input
        "Species"     = "#4575B4",
        "Sample_type" = "#66C2A5",
        "Continent"   = "#FDAE61"
      ),
      name = "Type",
      guide = guide_legend(order = 1)
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    ) +
    ggtitle(id_name)
}

# ==========================================
# 4. Automated Batch Processing
# ==========================================
output_dir <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/"

# Loop through list to plot and save
for (mge_name in names(target_mge)) {
  
  # Get color for current MGE
  current_color <- target_mge[[mge_name]]
  
  # Generate plot
  p <- plot_tree(df_count, mge_name, current_color)
  
  # Clean filename (replace spaces/special chars)
  file_name <- paste0(gsub(" ", "_", mge_name), ".png")
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, file_name),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    units = "in",
    bg = "white"
  )
  
  message(paste("Saved:", file_name))
}
