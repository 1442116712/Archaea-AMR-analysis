# =====================================================================
# Minimal demo — reproduce a mini version of Figure 6H
# (MILC + |ΔGC3| bubble plot for 4 prioritised archaeal ARGs against
#  their archaeal host and selected bacterial comparators).
#
# Data source: subset of Supplementary Table 18 (see manuscript).
#
# Expected runtime: < 30 seconds on a standard laptop.
# Expected output: expected_output.png (created in this folder).
# =====================================================================

# ---- 1. Dependencies ------------------------------------------------
# Install once (skip if already installed):
#   install.packages(c("ggplot2", "dplyr", "readr"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---- 2. Load small demo table ---------------------------------------
input_path <- "demo_input.tsv"                # if running from demo/ folder
if (!file.exists(input_path)) input_path <- "demo/demo_input.tsv"  # if running from repo root
df <- read_tsv(input_path, show_col_types = FALSE)

# Ensure factor ordering (rows = targets, columns = references)
df$target    <- factor(df$target,
                       levels = c("Msmithii_tetC","Mmazei_TEM","Mmazei_catA","Sislandicus_copR"))
df$reference <- factor(df$reference,
                       levels = unique(df$reference))

cat("Loaded", nrow(df), "target-reference comparisons.\n")

# ---- 3. Plot mini Figure 6H bubble ----------------------------------
p <- ggplot(df, aes(x = reference, y = target)) +
  geom_point(aes(size = MILC, fill = delta_GC3_pp),
             shape = 21, colour = "grey20", stroke = 0.4) +
  scale_size_continuous(name  = "MILC",
                        range = c(2, 12),
                        limits = c(0, max(df$MILC) * 1.05)) +
  scale_fill_gradient(name  = expression("|" * Delta * "GC3| (pp)"),
                      low   = "#F7F7F7",
                      high  = "#08519C") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Demo — mini Figure 6H",
       subtitle = "MILC (bubble size) + |ΔGC3| (bubble colour)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 0, face = "italic"),
        axis.text.y  = element_text(face = "italic"),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title     = element_text(face = "bold"),
        plot.subtitle  = element_text(color = "grey40"))

# ---- 4. Save PNG ----------------------------------------------------
out_png <- "expected_output.png"
if (dirname(input_path) == "demo") out_png <- "demo/expected_output.png"
ggsave(out_png, p, width = 6.5, height = 3.6, dpi = 200, bg = "white")
cat("Saved:", normalizePath(out_png, mustWork = FALSE), "\n")
