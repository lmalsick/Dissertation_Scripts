############################
# Script is used to look at WT WNV and sfRNA mutant and make the gene correlation graphs
# Helpful for looking at differential expressed genes in each condition that are opposititely egulated
# High expression in WNV, low expression in sfRNA mut and vice versa
# Author: Lauren Malsick
# Date: Final 11/30/2025
# Need all gene lists for this script, run only on one timepoint per comparison
# Looks at p.adjust < 0.005
#
############################
library(readr)
library(dplyr)
library(ggplot2)
library(MASS)       
library(ggrepel)      
library(viridis)     
library(scales)       


# ---- Load data ----
# Replace these with your actual file paths
sfrna1 <- read_csv("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/Gene_Lists/virus_timepoint_816.72_vs_Mock.72_allGenes.csv")
west <- read_csv("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/Gene_Lists/virus_timepoint_744.72_vs_Mock.72_allGenes.csv")

# ---- Merge datasets ----
merged <- inner_join(sfrna1, west, by = "Gene", suffix = c("_sfrna1", "_west")) %>%
  filter(!is.na(log2FoldChange_sfrna1), !is.na(log2FoldChange_west))

merged_clean <- merged %>%
  filter(is.finite(log2FoldChange_sfrna1), is.finite(log2FoldChange_west))

# Annotate significance status
merged_clean <- merged_clean %>%
  mutate(sig_status = case_when(
    padj_sfrna1 < 0.005 & padj_west < 0.005 ~ "Both significant",
    padj_sfrna1 < 0.005 ~ "sfRNA1 Mutant only",
    padj_west < 0.005 ~ "WNV only",
    TRUE ~ "Not significant"
  ))

# Focus on genes that are close to 0 in one condition but strongly regulated in the other
differential_genes <- merged_clean %>%
  mutate(
    # Calculate differential impact score
    wt_near_zero = abs(log2FoldChange_west) < 1.5,  # WT close to 0
    mutant_near_zero = abs(log2FoldChange_sfrna1) < 1.5,  # Mutant close to 0
    wt_strong = abs(log2FoldChange_west) > 3,  # WT strongly regulated
    mutant_strong = abs(log2FoldChange_sfrna1) > 3,  # Mutant strongly regulated
    wt_significant = padj_west < 0.05,  # WT significant
    mutant_significant = padj_sfrna1 < 0.05,  # Mutant significant
    
    # Define differential categories (now with significance requirement)
    diff_category = case_when(
      wt_near_zero & mutant_strong & mutant_significant ~ "Mutant-specific regulation",
      mutant_near_zero & wt_strong & wt_significant ~ "WT-specific regulation", 
      TRUE ~ "Other"
    )
  ) %>%
  filter(diff_category != "Other")

# Select top differential genes for labeling (with additional significance check)
n_labels_per_category <- 6
top_mutant_specific <- differential_genes %>%
  filter(diff_category == "Mutant-specific regulation") %>%
  # Additional check to ensure overall significance for coloring
  filter(padj_sfrna1 < 0.005 | padj_west < 0.005) %>%
  arrange(desc(abs(log2FoldChange_sfrna1))) %>%
  head(n_labels_per_category)

top_wt_specific <- differential_genes %>%
  filter(diff_category == "WT-specific regulation") %>%
  # Additional check to ensure overall significance for coloring  
  filter(padj_sfrna1 < 0.005 | padj_west < 0.005) %>%
  arrange(desc(abs(log2FoldChange_west))) %>%
  head(n_labels_per_category)

highlight_genes <- bind_rows(top_mutant_specific, top_wt_specific)

cat("Highlighting", nrow(highlight_genes), "genes with differential regulation patterns (significant and colored only)\n")

# Print summary of what we're highlighting
if(nrow(highlight_genes) > 0) {
  color_check <- highlight_genes %>%
    mutate(
      will_be_colored = case_when(
        padj_sfrna1 < 0.005 & padj_west < 0.005 ~ "Both significant (RED)",
        padj_sfrna1 < 0.005 ~ "sfRNA1 Mutant only (BLUE)",
        padj_west < 0.005 ~ "WNV only (GREEN)",
        TRUE ~ "Not significant (GREY)"
      )
    )
  cat("Color breakdown of highlighted genes:\n")
  print(table(color_check$will_be_colored))
}


# Custom theme for publication-ready plots
theme_publication <- function() {
  theme_minimal() +
    theme(
      # Text elements
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30", margin = margin(b = 15)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      
      # Panel and grid
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", size = 0.5, linetype = "dotted"),
      panel.grid.minor = element_line(color = "grey95", size = 0.3, linetype = "dotted"),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.8),
      
      # Legend
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey80", size = 0.5),
      legend.margin = margin(10, 10, 10, 10),
      legend.key.size = unit(0.8, "cm"),
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Enhanced color palette
colors_enhanced <- c(
  "Both significant" = "#E31A1C",      # Bright red
  "sfRNA1 Mutant only" = "#1F78B4",   # Blue  
  "WNV only" = "#33A02C",             # Green
  "Not significant" = "grey75"         # Light grey
)

########################### ---- PLOT 1:  ----
p1_enhanced <- ggplot(merged_clean, aes(x = log2FoldChange_west, y = log2FoldChange_sfrna1)) +
  # Add subtle density contours first (background layer)
  geom_contour(data = dens_df, aes(x = x, y = y, z = z), 
               color = "grey80", alpha = 0.6, bins = 8, size = 0.3) +
  
  # Main points - moved to front, no background interference
  geom_point(aes(color = sig_status, alpha = sig_status), 
             size = 1.2, stroke = 0) +
  
  # Custom alpha values for significance levels
  scale_alpha_manual(values = c(
    "Both significant" = 0.9,
    "sfRNA1 Mutant only" = 0.8,
    "WNV only" = 0.8,
    "Not significant" = 0.4
  ), guide = "none") +
  
  # Enhanced color scheme
  scale_color_manual(values = colors_enhanced, name = "Significance\n(padj < 0.005)") +
  
  # Highlight differential genes with subtle rings - made more transparent and thinner
  geom_point(data = highlight_genes, 
             shape = 21, color = "black", fill = NA, 
             size = 4, stroke = 0.8, alpha = 0.6) +
  
  # Gene labels for differential genes
  geom_text_repel(data = highlight_genes,
                  aes(label = Gene),
                  size = 3.5, 
                  fontface = "bold",
                  color = "black",
                  bg.color = "white",
                  bg.r = 0.1,
                  max.overlaps = 100,
                  min.segment.length = 0.1,
                  segment.size = 0.3,
                  segment.alpha = 0.6) +
  
  # Diagonal reference line
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "grey40", size = 0.8) +
  
  # Vertical and horizontal lines at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.6) +
  
  
  # Enhanced axis formatting
  scale_x_continuous(breaks = pretty_breaks(n = 8)) +
  scale_y_continuous(breaks = pretty_breaks(n = 8)) +
  
  # Labels
  labs(
    title = "Gene Expression Correlation Analysis",
    subtitle = "WT WNV vs sfRNA1 Mutant",
    x = "log2 Fold Change (WT WNV)",
    y = "log2 Fold Change (WNV sfRNA1 Mutant)"
  ) +
  
  # Apply custom theme
  theme_publication() +
  
  # Fine-tune legend
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 0.9),
    title.position = "top"
  ))

write_csv(merged_clean, "WNVandsfrna1_correlation_plot_data.csv")

print(p1_enhanced)

################################### ---- PLOT 2:  ----
diagonal_crossers <- merged_clean %>%
  filter(sign(log2FoldChange_sfrna1) != sign(log2FoldChange_west)) %>%
  filter(padj_sfrna1 < 0.05 | padj_west < 0.05) %>%
  mutate(cross_direction = case_when(
    log2FoldChange_sfrna1 > 0.1 & log2FoldChange_west < 0.1 ~ "Up in sfRNA1, Down in WNV",
    log2FoldChange_sfrna1 < 0.1 & log2FoldChange_west > 0.1 ~ "Down in sfRNA1, Up in WNV"
  )) %>%
  filter(!is.na(cross_direction))

cat("Diagonal crossers found:", nrow(diagonal_crossers), "\n")

# Select top genes to label based on differential regulation 
n_labels <- 5
top_sfrna1 <- diagonal_crossers %>%
  top_n(n_labels, wt = abs(log2FoldChange_sfrna1))
top_west <- diagonal_crossers %>%
  top_n(n_labels, wt = abs(log2FoldChange_west))
label_genes <- unique(c(top_sfrna1$Gene, top_west$Gene))

p2_enhanced <- ggplot(diagonal_crossers, aes(x = log2FoldChange_west, y = log2FoldChange_sfrna1)) +
  # Add quadrant backgrounds
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           fill = "#FF6B35", alpha = 0.05) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, 
           fill = "#4ECDC4", alpha = 0.05) +
  
  # Main points with enhanced styling
  geom_point(aes(color = cross_direction), 
             size = 1.8, alpha = 0.8, stroke = 0) +
  
  # Color scheme
  scale_color_manual(values = c(
    "Up in sfRNA1, Down in WNV" = "#1F78B4",
    "Down in sfRNA1, Up in WNV" = "#33A02C"
  ), name = "Regulation Pattern") +
  
  # Highlight labeled genes with subtle rings - made more transparent  
  geom_point(data = filter(diagonal_crossers, Gene %in% label_genes), 
             shape = 21, color = "black", fill = NA, 
             size = 4, stroke = 0.8, alpha = 0.6) +
  
  # Labels for selected genes
  geom_text_repel(
    data = filter(diagonal_crossers, Gene %in% label_genes),
    aes(label = Gene),
    size = 3.5,
    fontface = "bold",
    color = "black",
    bg.color = "white",
    bg.r = 0.1,
    max.overlaps = 25,
    min.segment.length = 0.05,
    segment.size = 0.4,
    segment.alpha = 0.7,
    force = 2,
    force_pull = 0.1
  ) +
  
  # Reference lines
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", size = 0.6) +
  
  # Enhanced formatting with extended axis limits
  scale_x_continuous(breaks = pretty_breaks(n = 8), 
                     expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_continuous(breaks = pretty_breaks(n = 8), 
                     expand = expansion(mult = c(0.15, 0.15))) +
  
  # Labels
  labs(
    title = "Oppositely Regulated Genes",
    subtitle = "Divergent Expression Patterns Between sfRNA1 Mutant and WT",
    x = "log2 Fold Change (WT WNV)",
    y = "log2 Fold Change (WNV sfRNA1 Mutant)"
  ) +
  
  # Apply theme
  theme_publication() +
  
  # Legend adjustments
  guides(color = guide_legend(
    override.aes = list(size = 3),
    title.position = "top"
  ))
write_csv(diagonal_crossers, "WNVandsfrna2_opposite_regulation_genes.csv")

print(p2_enhanced)

# ---- Summary of differential genes ----
cat("\nSummary of highlighted significant differential genes:\n")
if(nrow(highlight_genes) > 0) {
  highlight_summary <- highlight_genes %>%
    mutate(
      regulation_type = case_when(
        abs(log2FoldChange_west) < 1.5 & abs(log2FoldChange_sfrna1) > 3 ~ "Mutant-specific",
        abs(log2FoldChange_sfrna1) < 1.5 & abs(log2FoldChange_west) > 3 ~ "WT-specific",
        TRUE ~ "Other"
      )
    ) %>%
    select(Gene, log2FoldChange_west, log2FoldChange_sfrna1, padj_west, padj_sfrna1, regulation_type) %>%
    arrange(regulation_type, desc(abs(log2FoldChange_sfrna1) + abs(log2FoldChange_west)))
  
  print(highlight_summary)
}

# ---- Save high-quality plots ----
ggsave("wt_mutant_differential_correlation.png", p1_enhanced, 
       width = 12, height = 9, dpi = 300, bg = "white")

ggsave("wt_mutant_diagonal_crossers.png", p2_enhanced, 
       width = 12, height = 9, dpi = 300, bg = "white")
# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()
