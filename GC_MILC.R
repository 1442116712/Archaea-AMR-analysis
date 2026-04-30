library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Read data
df <- read.table("C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/1.txt",
                 header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE,
                 na.strings = c("NA", "", " "))

# 2. Long format + remove NA
df_long <- df %>%
  pivot_longer(cols = -ARG_ID, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric_Type = case_when(
      grepl("^GC3",  Metric) ~ "GC3",
      grepl("^MILC", Metric) ~ "MILC"
    ),
    Bacteria = gsub("^(GC3|MILC)\\s+", "", Metric)
  ) %>%
  filter(!is.na(Value))

# 3. Wide format + absolute difference of GC3
df_paired <- df_long %>%
  select(ARG_ID, Bacteria, Metric_Type, Value) %>%
  pivot_wider(names_from = Metric_Type, values_from = Value)

arg_gc3 <- df_paired %>%
  filter(Bacteria == "ARG") %>%
  select(ARG_ID, ARG_GC3 = GC3)

df_paired <- df_paired %>%
  left_join(arg_gc3, by = "ARG_ID") %>%
  mutate(GC3_absdiff = abs(GC3 - ARG_GC3)) %>%   # ★ Calculate absolute value
  filter(Bacteria != "ARG")

# 4. Sorting (After swapping x and y, place Archaea host on the far left)
host_first <- "Archaea host"
others <- sort(setdiff(unique(df_paired$Bacteria), host_first))   # Ascending alphabetical order

# x-axis levels from left to right: first = leftmost
df_paired$Bacteria <- factor(df_paired$Bacteria,
                             levels = c(host_first, others))

# y-axis: ARGs are plotted from bottom to top; to put the first ARG at the top, use rev()
df_paired$ARG_ID <- factor(df_paired$ARG_ID, 
                           levels = rev(unique(df_paired$ARG_ID)))

# 5. Bubble chart (x=Bacteria, y=ARG_ID)
p_bubble <- ggplot(df_paired, aes(x = Bacteria, y = ARG_ID)) +
  geom_point(aes(size = MILC, fill = GC3_absdiff),
             shape = 21, color = "grey30", stroke = 0.3) +
  scale_fill_gradient(
    name   = "|GC3 difference|%",
    low    = "#F7F4FB",
    high   = "#3C3489",
    limits = c(0, max(df_paired$GC3_absdiff, na.rm = TRUE)),
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_size_continuous(
    name   = "MILC",
    range  = c(2, 8),
    breaks = c(0.5, 1.0, 1.5, 2.0)
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                    size = 10, face = "bold"),  # Species names, italic
    axis.text.y      = element_text(size = 9, face = "bold"),   # ARG names, bold
    axis.title       = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey60", fill = NA, linewidth = 0.4),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9,  face = "bold")
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = 6),
    size = guide_legend(order = 2)
  )

print(p_bubble)

# 6. Export (Since x and y are swapped, width and height proportions are also swapped)
n_args    <- nlevels(df_paired$ARG_ID)
n_species <- nlevels(df_paired$Bacteria)

ggsave(
  "C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/ARG_align_prok/GC3absdiff_MILC_bubble_matrix.png",
  plot   = p_bubble,
  width  = 16,   # ★ Now width is determined by the number of species
  height = 7,    # ★ Height is determined by the number of ARGs
  dpi    = 1000, units = "in", bg = "white"
)
