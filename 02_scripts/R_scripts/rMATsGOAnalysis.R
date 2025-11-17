##########################################
# Name: rMATS GO Analysis
# Description: 
# This script goes through the rMATs data previously organized and cleaned up by junction numbers
# As a reminder, PSI values without junction counts are not super useful because there could be less than 5-10 reads per junction
# with less reads, more possibility for genes to pop up significant but not biologically relevant
# 
# You must first run rmatsFilterCounts.R, then that output file is used in this script
# Author: Lauren Malsick
# Date: 11/17/2025
#########################################

# Load necessary libraries
library(clusterProfiler)
library(org.Ajamaicensis.eg.db)
library(ggplot2)
library(tidyverse)
library(enrichplot)
library(forcats)
library(stats)
library(dplyr)

setwd("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing//rMATs_analysis/rMATS_output_744_mock_72hpi")
output_dir <- dir.create("GSEA", recursive = TRUE)
rna_data <- read.csv("rMATS_filtered_events.csv", header = TRUE)

df <- rna_data

df$PValue <- suppressWarnings(as.numeric(df$PValue))
df$IncLevelDifference <- suppressWarnings(as.numeric(df$IncLevelDifference))

# Remove missing values
df <- df[!is.na(df$PValue) & df$PValue > 0, ]
df <- df[!is.na(df$IncLevelDifference), ]

#################################
# Biological significance threshold
# Filter for |deltaPSI| > 0.05 (5% change in splicing)
# Adjust threshold as needed: 
#################################

delta_psi_threshold <- 0.05

cat("Total events before filtering:", nrow(df), "\n")

df_filtered <- df %>%
  filter(abs(IncLevelDifference) >= delta_psi_threshold)

cat("Events after |deltaPSI| >=", delta_psi_threshold, "filter:", nrow(df_filtered), "\n")
cat("Genes with multiple events:", 
    df_filtered %>% count(GeneID) %>% filter(n > 1) %>% nrow(), "\n")

#################################
#This part of the script selects for the genes with the LARGEST |deltaPSI|
#################################

df_collapsed <- df_filtered %>%
  group_by(GeneID) %>%
  slice_max(order_by = abs(IncLevelDifference), n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Unique genes after collapsing:", nrow(df_collapsed), "\n\n")

# Calculate ranking score: deltaPSI * -log10(p-value)
df_collapsed$rankScore <- df_collapsed$IncLevelDifference * -log10(df_collapsed$PValue)

# Remove non-finite values
df_collapsed <- df_collapsed[is.finite(df_collapsed$rankScore), ]

# Create ranked gene list
ranked_genes <- df_collapsed$rankScore
names(ranked_genes) <- df_collapsed$GeneID

# Sort high to low
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Save ranked gene list
write.csv(
  data.frame(
    GeneID = names(ranked_genes), 
    RankScore = ranked_genes,
    IncLevelDifference = df_collapsed$IncLevelDifference[match(names(ranked_genes), df_collapsed$GeneID)],
    PValue = df_collapsed$PValue[match(names(ranked_genes), df_collapsed$GeneID)]
  ),
  paste0(output_dir, "/ranked_genes_maxAbsEffect.csv"), 
  row.names = FALSE
)

pdf(paste0(output_dir, "/Quality_Control.pdf"), width=12, height=8)

par(mfrow=c(2,2))

# Rank score distribution - This looks for a bell curve which will contain the deltaPSI*-log10(pvalues) -- shows spread
# Frequency of the # of events in the data
hist(ranked_genes, breaks=50, main="Rank Score Distribution", 
     xlab="IncLevelDiff * -log10(p-value)", col="skyblue")
abline(v=0, col="red", lwd=2, lty=2)

# IncLevelDifference distribution - This shows the deltaPSI value distribution in the sample with the filter (0.05 or 5%) as the threshold
hist(df_collapsed$IncLevelDifference, breaks=50, 
     main="IncLevelDifference Distribution", 
     xlab="IncLevelDifference", col="lightgreen")
abline(v=c(-delta_psi_threshold, delta_psi_threshold), col="red", lwd=2, lty=2)

# -log10(p-value) distribution -- shows the spread of p-values
hist(-log10(df_collapsed$PValue), breaks=50, 
     main="-log10(p-value) Distribution", 
     xlab="-log10(p-value)", col="lightcoral")

# Scatter/volcano: effect size vs significance -- Shows deltaPSI graphed against the -log10(pvalue)
#Basically shows the distribution in a different way
plot(df_collapsed$IncLevelDifference, -log10(df_collapsed$PValue),
     xlab="IncLevelDifference", ylab="-log10(p-value)",
     main="Effect Size vs Significance",
     pch=16, col=rgb(0,0,0,0.3))
abline(v=c(-delta_psi_threshold, delta_psi_threshold), col="red", lwd=2, lty=2)

dev.off()

# -------------------------------
# Run GSEA for GO categories
# -------------------------------
goCats <- c("BP", "MF", "CC")
gsea_results <- list()

for (cat in goCats) {
  message("Running GSEA GO:", cat)
  
  g <- gseGO(
    geneList = ranked_genes,
    OrgDb = org.Ajamaicensis.eg.db,
    keyType = "GID",
    ont = cat,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  if (is.null(g) || nrow(g@result) == 0) {
    message("No enriched terms for: ", cat, " — skipping plots.")
    next
  }
  
  gsea_results[[cat]] <- g
  
  # Save results table
  write.csv(g@result,
            paste0(output_dir, "/GSEA_", cat, "_results.csv"),
            row.names = FALSE)
  
  # Dotplot
  pdf(paste0(output_dir, "/GSEA_", cat, "_dotplot.pdf"), width=10, height=8)
  print(dotplot(g, showCategory=10) + ggtitle(paste("GSEA:", cat)))
  dev.off()
  
  # Enrichment plot for top pathways -- downward curves for negative NES aka decreasing, then the lines show where genes cluster in the list
  # bottom again shows positive/negative NES and the curves will correspond to the genes and where they fall in the matrix
  if (nrow(g@result) > 0) {
    pdf(paste0(output_dir, "/GSEA_", cat, "_enrichplot.pdf"), width=12, height=10)
    print(gseaplot2(g, geneSetID = 1:min(5, nrow(g@result)), pvalue_table = TRUE))
    dev.off()
  }
}

# -------------------------------
# Summary report
# -------------------------------
cat("\n=== GSEA Results Summary ===\n")
for (cat in goCats) {
  if (!is.null(gsea_results[[cat]])) {
    sig_terms <- sum(gsea_results[[cat]]@result$p.adjust < 0.05)
    cat(cat, ":", sig_terms, "significant terms (p.adj < 0.05)\n")
  }
}

message("Files saved to: ", output_dir)
message("\nFiltering parameters:")
message("|deltaPSI| threshold: ", delta_psi_threshold)
message("Ranking strategy: Maximum absolute effect size per gene")
message("Ranking metric: deltaPSI × -log10(p-value)")

# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()

