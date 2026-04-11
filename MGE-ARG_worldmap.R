# ==========================================
# 1. Load Required Libraries
# ==========================================
library(readxl)
library(tidyverse)
library(maps)
library(mapdata)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# ==========================================
# 2. Data Loading & Preprocessing
# ==========================================
file_path <- "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/Archaea_MGE_nt_prok.xlsx"
df <- read_excel(file_path, sheet = "co-related")

# Calculate count per Query ID for each Country
df_summary <- df %>%
  group_by(queryid, Country) %>%
  summarise(count = n(), .groups = 'drop')

# ==========================================
# 3. Geographic Data Preparation
# ==========================================
# Fetch medium-scale world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Calculate precise country centroids for point placement
country_centers <- world %>%
  st_centroid() %>%
  mutate(
    base_long = st_coordinates(geometry)[,1],
    base_lat = st_coordinates(geometry)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(name, base_long, base_lat)

# Join summarized data with geographic coordinates
df_plot <- df_summary %>%
  left_join(country_centers, by = c("Country" = "name"))

# Identify country names that failed to match geographic coordinates
unmatched <- df_plot %>% 
  filter(is.na(base_long)) %>% 
  distinct(Country) %>% 
  pull(Country)

if(length(unmatched) > 0) {
  cat("\nUnmatched country names (excluded from map):\n")
  print(unmatched)
  
  unmatched_counts <- df_summary %>%
    filter(Country %in% unmatched) %>%
    group_by(Country) %>%
    summarise(total_count = sum(count), .groups = 'drop')
  
  cat("\nSample counts from unmatched countries:\n")
  print(unmatched_counts)
}

# Remove entries without valid coordinates
df_plot <- df_plot %>% filter(!is.na(base_long))

cat("\nSuccessfully matched countries:", length(unique(df_plot$Country)), "\n")
cat("Total sample count:", sum(df_plot$count), "\n")

# ==========================================
# 4. Coordinate Offset Logic (Avoid Overlap)
# ==========================================
# Distance to offset points from the centroid
offset_dist <- 2.5

# Unique IDs to map
unique_ids <- unique(df_plot$queryid)
n_ids <- length(unique_ids)

# Generate dynamic offset directions based on the number of IDs
offset_dirs <- list()
if(n_ids == 2) {
  # Layout for 2 IDs: Left and Right
  offset_dirs[[unique_ids[1]]] <- list(dx = -offset_dist, dy = 0)
  offset_dirs[[unique_ids[2]]] <- list(dx = offset_dist, dy = 0)
} else if(n_ids == 3) {
  # Layout for 3 IDs: Triangular
  offset_dirs[[unique_ids[1]]] <- list(dx = 0, dy = offset_dist)
  offset_dirs[[unique_ids[2]]] <- list(dx = -offset_dist * cos(pi/6), 
                                       dy = -offset_dist * sin(pi/6))
  offset_dirs[[unique_ids[3]]] <- list(dx = offset_dist * cos(pi/6), 
                                       dy = -offset_dist * sin(pi/6))
} else {
  # Layout for >3 IDs: Circular distribution
  for(i in 1:n_ids) {
    angle <- 2 * pi * (i-1) / n_ids
    offset_dirs[[unique_ids[i]]] <- list(
      dx = offset_dist * cos(angle),
      dy = offset_dist * sin(angle)
    )
  }
}

# Identify countries containing multiple Query IDs
multi_id_countries <- df_plot %>%
  group_by(Country) %>%
  summarise(n_ids = n_distinct(queryid), .groups = 'drop') %>%
  filter(n_ids > 1) %>%
  pull(Country)

# Apply calculated offsets to coordinates
df_plot <- df_plot %>%
  rowwise() %>%
  mutate(
    long = ifelse(Country %in% multi_id_countries,
                  base_long + offset_dirs[[queryid]]$dx,
                  base_long),
    lat = ifelse(Country %in% multi_id_countries,
                 base_lat + offset_dirs[[queryid]]$dy,
                 base_lat)
  ) %>%
  ungroup()

# Define specific color mapping for Homologs
colors <- c(
  "Homologs of Tn_lunAN2_mefEN2" = "#FF6B6B",
  "Homologs of IS_copR"          = "#1F77B4",
  "Homologs of plasmid_catA"     = "#31A354"
)

# ==========================================
# 5. Map Visualization
# ==========================================
p_final <- ggplot() +
  # Draw base world map layers
  geom_sf(data = world, fill = "gray90", color = "white", size = 0.2) +
  # Add offset point layers
  geom_point(data = df_plot,
             aes(x = long, y = lat, 
                 size = count, 
                 color = queryid,
                 fill = queryid),
             alpha = 0.7) +
  # Legend and Color Scales
  scale_color_manual(values = colors, name = "Query ID") +
  scale_fill_manual(values = colors, name = "Query ID") +
  scale_size_continuous(range = c(2, 8), name = "Count") +
  # Layout and Theme Customization
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    fill = guide_legend(override.aes = list(size = 5)),
    size = guide_legend(override.aes = list(color = "black"))
  )

# Display result
print(p_final)

# ==========================================
# 6. Save Export
# ==========================================
ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/MGE_world_map.png",     
  plot = p_final,      
  width = 14,        
  height = 8,        
  dpi = 300,        
  units = "in",
  bg = "white" 
)
