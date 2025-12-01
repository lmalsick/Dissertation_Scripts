##########################################
# Name: rMATS GO Analysis
# Description: 
# This script goes through the rMATs data previously organized and cleaned up by junction numbers
# As a reminder, PSI values without junction counts are not super useful because there could be less than 5-10 reads per junction
# with less reads, more possibility for genes to pop up significant but not biologically relevant
# 
# If GSEA has low pathways, try ORA below
#
# !!!!!!!!!!!! You must first run rmatsFilterCounts.R, then that output file is used in this script
#
# Author: Lauren Malsick
# Date: 11/17/2025
#########################################

# Load libraries
library(clusterProfiler)
library(org.Ajamaicensis.eg.db)
library(ggplot2)
library(tidyverse)
library(enrichplot)
library(forcats)
library(stats)
library(dplyr)


setwd("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/rMATs_analysis/rMATS_output_744_mock_72hpi")
output_dir <- "GSEA"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
rna_data <- read.csv("rMATS_filtered_event_744_Mock_72hr.csv", header = TRUE)

# Convert all columns to appropriate types
df <- rna_data

# Force conversion of key columns to numeric/character
df$GeneID <- as.character(df$GeneID)
df$PValue <- suppressWarnings(as.numeric(as.character(df$PValue)))
df$IncLevelDifference <- suppressWarnings(as.numeric(as.character(df$IncLevelDifference)))
df$FDR <- suppressWarnings(as.numeric(as.character(df$FDR)))

# Select only columns needed for GSEA (base R way)
df_clean <- df[, c("GeneID", "PValue", "IncLevelDifference", "FDR")]
df_clean <- as.data.frame(df_clean)

# Remove missing values
df_clean <- df_clean[!is.na(df_clean$PValue) & df_clean$PValue > 0, ]
df_clean <- df_clean[!is.na(df_clean$IncLevelDifference), ]

#################################
# Biological significance threshold
# Filter for |deltaPSI| > 0.05 (5% change in splicing)
# Adjust threshold as needed: 
#################################

delta_psi_threshold <- 0.05

cat("Total events before filtering:", nrow(df_clean), "\n")

# Filter by deltaPSI threshold (base R)
df_filtered <- df_clean[abs(df_clean$IncLevelDifference) >= delta_psi_threshold, ]

cat("Events after |deltaPSI| >=", delta_psi_threshold, "filter:", nrow(df_filtered), "\n")

# Count genes with multiple events (base R)
gene_counts <- table(df_filtered$GeneID)
cat("Genes with multiple events:", sum(gene_counts > 1), "\n")

#################################
# This part selects genes with the LARGEST |deltaPSI|
# IF THERE ARE TIES: takes the largest absolute |psi| with the smallest p-value
#################################

# Collapse to one event per gene using base R
df_collapsed <- do.call(rbind, lapply(split(df_filtered, df_filtered$GeneID), function(x) {
  # Sort by absolute IncLevelDifference (descending), then PValue (ascending)
  x <- x[order(-abs(x$IncLevelDifference), x$PValue), ]
  # Return first row
  x[1, , drop = FALSE]
}))

# Reset row names
rownames(df_collapsed) <- NULL

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

###############################
# Gene Set Enrichment Analysis (GSEA)
###############################

goCats <- c("BP", "MF", "CC")
gsea_results <- list()

for (cat in goCats) {
  message("Running GSEA GO:", cat)
  
  set.seed(2025)
  g <- gseGO(
    geneList = ranked_genes,
    OrgDb = org.Ajamaicensis.eg.db,
    keyType = "SYMBOL",
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

for (cat in goCats) {
  if (!is.null(gsea_results[[cat]])) {
    sig_terms <- sum(gsea_results[[cat]]@result$p.adjust < 0.05)
    cat(cat, ":", sig_terms, "significant terms (p.adj < 0.05)\n")
  }
}

message("Files saved to: ", output_dir)
message("\nFiltering parameters:")
message("|deltaPSI| threshold: ", delta_psi_threshold)
message("Ranking metric: deltaPSI × -log10(p-value)")

#IF GSEA has limited results, try ORA
# GSEA is usually for larger datasets, at least ORA can show us trends in patterns

################################
# Over-Representation Analysis (ORA)
###############################

# Output directory for ORA
output_dir_ora <- "ORA"
dir.create(output_dir_ora, recursive = TRUE, showWarnings = FALSE)

sig_events <- df_collapsed %>%
  filter(FDR < 0.05)

sig_genes <- sig_events$GeneID

cat("Significant genes for ORA:", length(sig_genes), "\n")

# -------------------------------
# 1) Build background gene universe from rMATS files
# This is just if I wanted to use rMATs background and not universal backgraound
# This background would answer the question of out of all the spliced genes, what pathways are the most affected by alternative splicing
# -------------------------------
# # Find all MATS.JCEC.txt files
# files <- list.files(".", pattern = "MATS.JCEC.txt$", full.names = TRUE)
# cat("Found", length(files), "rMATS event files\n")
# 
# background_list <- lapply(files, function(f) {
#   df <- read.table(f, header = TRUE, sep = "\t", 
#                    quote = "\"",  # ✅ Change from "" to "\""
#                    comment.char = "", 
#                    stringsAsFactors = FALSE)
#   
#   if ("GeneID" %in% colnames(df)) {
#     return(df$GeneID)
#   } else if ("geneSymbol" %in% colnames(df)) {
#     return(df$geneSymbol)
#   } else {
#     return(NULL)
#   }
# })
# 
# background_genes <- unique(unlist(background_list))
# background_genes <- background_genes[!is.na(background_genes) & background_genes != ""]
# 
# # Remove any lingering quotes just in case
# background_genes <- gsub('^"|"$', '', background_genes)

# cat("Background genes:", length(background_genes), "\n")
# cat("Overlap with sig_genes:", sum(sig_genes %in% background_genes), "\n")
# 
# # If overlap is 0, use database as background
# if (sum(sig_genes %in% background_genes) == 0) {
#   cat("WARNING: No overlap between sig_genes and background!\n")
#   cat("Using all genes in database as background instead.\n")
#   background_genes <- keys(org.Ajamaicensis.eg.db, keytype = "SYMBOL")
#   cat("New background size:", length(background_genes), "\n")
#   cat("New overlap:", sum(sig_genes %in% background_genes), "\n")
# }
# 
# # Test if genes are in database
# test_genes <- head(sig_genes, 5)
# cat("\nTesting first 5 sig genes in database:\n")
# for (g in test_genes) {
#   exists <- g %in% keys(org.Ajamaicensis.eg.db, keytype = "SYMBOL")
#   cat("  ", g, ":", exists, "\n")
# }

############### 1) This background will be the entire AJ database, so looking for biological processes affected by splicing

background_genes <- keys(org.Ajamaicensis.eg.db, keytype = "SYMBOL")
cat("Background universe genes:", length(background_genes), "\n")
cat("Overlap with sig_genes:", sum(sig_genes %in% background_genes), "\n")

################ 2) Run ORA

ora_results <- list()

for (cat in goCats) {
  message("\nRunning ORA GO:", cat)
  
  g <- enrichGO(
    gene          = sig_genes,
    universe      = background_genes,
    OrgDb         = org.Ajamaicensis.eg.db,
    keyType       = "SYMBOL",
    ont           = cat,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 5,
    maxGSSize     = 500,
    readable      = TRUE
  )
  
  if (is.null(g) || nrow(g@result) == 0) {
    message("  No enriched terms for: ", cat)
    next
  }
  
  sig_count <- sum(g@result$p.adjust < 0.05)
  message("  Found ", sig_count, " significant terms (FDR < 0.05)")
  
  ora_results[[cat]] <- g
  
  # Save results
  write.csv(
    g@result,
    paste0(output_dir_ora, "/ORA_", cat, "_results.csv"),
    row.names = FALSE
  )
  
  # Plots
  if (sig_count > 0) {
    pdf(paste0(output_dir_ora, "/ORA_", cat, "_dotplot.pdf"), width=10, height=14)
    print(dotplot(g, showCategory = 20) + ggtitle(paste("ORA:", cat)))
    dev.off()
    
    pdf(paste0(output_dir_ora, "/ORA_", cat, "_barplot.pdf"), width=10, height=14)
    print(barplot(g, showCategory = 20) + ggtitle(paste("ORA:", cat)))
    dev.off()
  }
}

# ORA Summary
cat("\n=== ORA Results Summary ===\n")
for (cat in goCats) {
  if (!is.null(ora_results[[cat]])) {
    sig_terms <- sum(ora_results[[cat]]@result$p.adjust < 0.05)
    cat(cat, ":", sig_terms, "significant terms (FDR < 0.05)\n")
  } else {
    cat(cat, ": No enrichment\n")
  }
}

cat("GSEA results saved to:", output_dir, "\n")
cat("ORA results saved to:", output_dir_ora, "\n")


# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()

