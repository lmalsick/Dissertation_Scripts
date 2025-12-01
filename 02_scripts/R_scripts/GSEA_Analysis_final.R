##################################################
#The following library is a custom package which must be manually installed from source. 
#This can be performed with the following command, where "org.Aaegypti.eg.db" is the path to the org.Aaegypti.eg.db folder
#install.packages("org.Aaegypti.eg.db", repos=NULL,type="source")
#This folder is found in the Rosenberg NAS Informatics Resources -> Aedes_aegypti
# For this script we will use the Artibeus jamaicensis genome:
# install.packages("org.Ajamaicensis.eg.db", repos=NULL,type="source")
library(org.Ajamaicensis.eg.db)
#
##FUNCTIONS IN THIS FILE
#
#rnaGSEA: perform GSEA analysis on single file, which should be differential expression output (generally DESeq2)
#Should be used for GO annotations
#Requires an OrgDB library
#
# Script Author: Hunter Ogg with modifications by Lauren Malsick
# Script: GSEA_Analysis_final.R
# Date: Final version 11/30/2025
#
####################################################

library(tidyverse)
library(patchwork)
library(enrichplot)
library(clusterProfiler)
library(org.Ajamaicensis.eg.db)

getwd()
setwd("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/GSEA")

# Setup output directory
output_dir <- "GSEA_Publication_Figures_bothRank"
gsea_classic_dir <- file.path(output_dir, "GSEA_Classic_Plots_both")
cnet_dir <- file.path(output_dir, "Gene_Concept_Networks_both")

for (dir in c(output_dir, gsea_classic_dir, cnet_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

########################## START GSEA Function

#rnaFile is a csv which contains the results of differential gene expression. Generally the output of DEseq analysis
#rankMetric is the the choice of ranking (pval, log2fc or both). Currently only supports a single, but more will be added in future
#goCat represents which of the three subontologies (CC,MF,BP) to include in the analysis. All three are included by default

rnaGSEA = function(rnaFile, orgData, rankMetric = "both", goCat = "All", 
                   minSize = 15, maxSize = 500, keycol = "SYMBOL") {
  myDEresults = read.csv(rnaFile)
  
  if (rankMetric == "pval") {
    myDEresults = myDEresults[!is.na(myDEresults$pvalue), ]
    newRank_pvalueAndFC = -log10(myDEresults$pvalue) * sign(myDEresults$log2FoldChange)
    names(newRank_pvalueAndFC) = myDEresults$Gene
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC, decreasing = TRUE)]
  }
  if (rankMetric == "log2fc") {
    myDEresults = myDEresults[!is.na(myDEresults$pvalue), ]
    newRank_pvalueAndFC = myDEresults$log2FoldChange
    names(newRank_pvalueAndFC) = myDEresults$Gene
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC, decreasing = TRUE)]
  }
  if (rankMetric == "both") {
    myDEresults = myDEresults[!is.na(myDEresults$pvalue), ]
    newRank_pvalueAndFC = myDEresults$log2FoldChange * -log10(myDEresults$pvalue)
    names(newRank_pvalueAndFC) = myDEresults$Gene
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC, decreasing = TRUE)]
  }
  
  newRank_pvalueAndFC = newRank_pvalueAndFC[is.finite(newRank_pvalueAndFC)]
  
  set.seed(2025)
  egoCC <- gseGO(geneList     = newRank_pvalueAndFC,
                 OrgDb        = orgData,
                 keyType      = keycol,
                 ont          = goCat,
                 minGSSize    = minSize,
                 maxGSSize    = maxSize,
                 pvalueCutoff = 0.5,
                 verbose      = FALSE)
  
  return(egoCC)
}

###################################
#
# Start graphs
#
######################################
# Create dotplot
create_gsea_dotplot <- function(gsea_obj, title, top_n = 20) {
  
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    return(NULL)
  }
  
  # Filter significant results
  sig_results <- gsea_obj@result %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = top_n)
  
  if (nrow(sig_results) == 0) {
    cat("    ⚠ No significant results for", title, "\n")
    return(NULL)
  }
  
  # Calculate GeneRatio properly from the core_enrichment genes
  # GeneRatio = (# core enriched genes) / (total genes in pathway)
  sig_results <- sig_results %>%
    mutate(
      core_count = sapply(strsplit(as.character(core_enrichment), "/"), length),
      GeneRatio = core_count / setSize
    )
  
  # Order by GeneRatio for plotting (increasing size)
  sig_results <- sig_results %>%
    arrange(GeneRatio) %>%
    mutate(Description = factor(Description, levels = Description))
  
  p <- ggplot(sig_results, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = setSize, color = pvalue)) +
    scale_color_gradient(
      low = "red", 
      high = "blue",
      name = "P-value",
      trans = "log10",
      labels = scales::scientific
    ) +
    scale_size_continuous(
      name = "Gene Set\nSize",
      range = c(3, 10)
    ) +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
  
  return(p)
}

# Create GSEA running score plots
create_classic_gsea_plots <- function(gsea_obj, comparison_name, go_cat, top_n = 5) {
  
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    return(NULL)
  }
  
  sig_results <- gsea_obj@result %>%
    filter(pvalue < 0.05) %>%
    arrange(pvalue) %>%
    slice_head(n = top_n)
  
  if (nrow(sig_results) == 0) {
    cat("No significant results for classic GSEA plots:", comparison_name, "\n")
    return(NULL)
  }
  
  # Create gseaplot2 for top pathways
  tryCatch({
    p <- gseaplot2(gsea_obj, 
                   geneSetID = 1:min(3, nrow(sig_results)), 
                   pvalue_table = TRUE,
                   title = paste(comparison_name, "-", go_cat))
    
    filename <- paste0("GSEA_Classic_", 
                       gsub("[^A-Za-z0-9]", "_", comparison_name), 
                       "_", go_cat)
    
    ggsave(file.path(gsea_classic_dir, paste0(filename, ".pdf")), 
           p, width = 10, height = 8, dpi = 300)
    ggsave(file.path(gsea_classic_dir, paste0(filename, ".png")), 
           p, width = 10, height = 8, dpi = 300)
    
    cat("GSEA plot saved:", filename, "\n")
    
  }, error = function(e) {
    cat("Error creating classic GSEA plot:", e$message, "\n")
  })
}

# Create gene-concept network
create_gene_concept_network <- function(gsea_obj, comparison_name, go_cat, 
                                        show_category = 5, gene_list = NULL) {
  
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    return(NULL)
  }
  
  sig_results <- gsea_obj@result %>%
    filter(pvalue < 0.05) %>%
    arrange(pvalue)
  
  if (nrow(sig_results) < 2) {
    cat("Not enough significant results for cnetplot:", comparison_name, "\n")
    return(NULL)
  }
  
  tryCatch({
    # Create cnetplot
    if (!is.null(gene_list)) {
      p <- cnetplot(gsea_obj, 
                    showCategory = min(show_category, nrow(sig_results)),
                    foldChange = gene_list,
                    colorEdge = TRUE,
                    cex_label_gene = 0.7,
                    cex_label_category = 1.0) +
        ggtitle(paste("Gene-Concept Network:", comparison_name, "-", go_cat))
    } else {
      p <- cnetplot(gsea_obj, 
                    showCategory = min(show_category, nrow(sig_results)),
                    colorEdge = TRUE,
                    cex_label_gene = 0.7,
                    cex_label_category = 1.0) +
        ggtitle(paste("Gene-Concept Network:", comparison_name, "-", go_cat))
    }
    
    filename <- paste0("CNet_", 
                       gsub("[^A-Za-z0-9]", "_", comparison_name), 
                       "_", go_cat)
    
    ggsave(file.path(cnet_dir, paste0(filename, ".pdf")), 
           p, width = 12, height = 10, dpi = 300)
    ggsave(file.path(cnet_dir, paste0(filename, ".png")), 
           p, width = 12, height = 10, dpi = 300)
    
    cat("Gene-concept network saved:", filename, "\n")
    
  }, error = function(e) {
    cat("Error creating gene-concept network:", e$message, "\n")
  })
}

#####################
# 
# Define the paths for the all genes locations
#
#####################

base_path <- "~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/Gene_Lists"

# Define all your comparisons with their file paths
comparisons <- list(
  list(
    name = "WT 12h vs Mock",
    file = file.path(base_path, "virus_timepoint_744.12_vs_Mock.12_allgenes.csv"),
    label = "WT WNV vs Mock (12h)"
  ),
  list(
    name = "WT 24h vs Mock",
    file = file.path(base_path, "virus_timepoint_744.24_vs_Mock.24_allgenes.csv"),
    label = "WT WNV vs Mock (24h)"
  ),
  list(
    name = "WT 48h vs Mock",
    file = file.path(base_path, "virus_timepoint_744.48_vs_Mock.48_allgenes.csv"),
    label = "WT WNV vs Mock (48h)"
  ),
  list(
    name = "WT 72h vs Mock",
    file = file.path(base_path, "virus_timepoint_744.72_vs_Mock.72_allgenes.csv"),
    label = "WT WNV vs Mock (72h)"
  ),
  list(
    name = "sfRNA1 12h vs Mock",
    file = file.path(base_path, "virus_timepoint_816.12_vs_Mock.12_allgenes.csv"),
    label = "sfRNA1 vs Mock (12h)"
  ),
  list(
    name = "sfRNA1 24h vs Mock",
    file = file.path(base_path, "virus_timepoint_816.24_vs_Mock.24_allgenes.csv"),
    label = "sfRNA1 vs Mock (24h)"
  ),
  list(
    name = "sfRNA1 48h vs Mock",
    file = file.path(base_path, "virus_timepoint_816.48_vs_Mock.48_allgenes.csv"),
    label = "sfRNA1 vs Mock (48h)"
  ),
  list(
    name = "sfRNA1 72h vs Mock",
    file = file.path(base_path, "virus_timepoint_816.72_vs_Mock.72_allgenes.csv"),
    label = "sfRNA1 vs Mock (72h)"
  ),
  list(
    name = "sfRNA1and2 12h vs Mock",
    file = file.path(base_path, "virus_timepoint_853.12_vs_Mock.12_allgenes.csv"),
    label = "sfRNA1and2 vs Mock (12h)"
  ),
  list(
    name = "sfRNA1and2 24h vs Mock",
    file = file.path(base_path, "virus_timepoint_853.24_vs_Mock.24_allgenes.csv"),
    label = "sfRNA1and2 vs Mock (24h)"
  ),
  list(
    name = "sfRNA1and2 48h vs Mock",
    file = file.path(base_path, "virus_timepoint_853.48_vs_Mock.48_allgenes.csv"),
    label = "sfRNA1and2 vs Mock (48h)"
  ),
  list(
    name = "sfRNA1and2 72h vs Mock",
    file = file.path(base_path, "virus_timepoint_853.72_vs_Mock.72_allgenes.csv"),
    label = "sfRNA1and2 vs Mock (72h)"
  ),
  list(
    name = "ENTV 12h vs Mock",
    file = file.path(base_path, "virus_timepoint_ENTV.12_vs_Mock.12_allgenes.csv"),
    label = "ENTV vs Mock (12h)"
  ),
  list(
    name = "ENTV 24h vs Mock",
    file = file.path(base_path, "virus_timepoint_ENTV.24_vs_Mock.24_allgenes.csv"),
    label = "ENTV vs Mock (24h)"
  ),
  list(
    name = "ENTV 48h vs Mock",
    file = file.path(base_path, "virus_timepoint_ENTV.48_vs_Mock.48_allgenes.csv"),
    label = "ENTV vs Mock (48h)"
  ),
  list(
    name = "ENTV 72h vs Mock",
    file = file.path(base_path, "virus_timepoint_ENTV.72_vs_Mock.72_allgenes.csv"),
    label = "ENTV vs Mock (72h)"
  )
)

################
# RUN FUNCTION

# Store all GSEA results
all_gsea_results <- list()

for (comp in comparisons) {
  cat("\n=== Processing:", comp$name, "===\n")
  
  # Check if file exists
  if (!file.exists(comp$file)) {
    cat("File not found:", comp$file, "\n")
    next
  }
  
  # Read DESeq results to get gene list for cnetplot
  de_results <- read.csv(comp$file)
  gene_list <- setNames(de_results$log2FoldChange, de_results$Gene)
  gene_list <- gene_list[is.finite(gene_list)]
  
  # Run GSEA for each GO category
  for (go_cat in c("BP", "MF", "CC")) {
    cat("  Processing GO category:", go_cat, "\n")
    
    # Run GSEA
    gsea_result <- rnaGSEA(
      rnaFile = comp$file,
      orgData = org.Ajamaicensis.eg.db,
      rankMetric = "both",
      goCat = go_cat,
      keycol = "SYMBOL"
    )
    
    # Store result
    result_name <- paste0(comp$name, "_", go_cat)
    all_gsea_results[[result_name]] <- gsea_result
    
    if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
      cat("No results for", comp$name, go_cat, "\n")
      next
    }
    
    # Create dotplot
    cat("Creating dotplot...\n")
    dotplot_obj <- create_gsea_dotplot(gsea_result, 
                                       paste(comp$label, "-", go_cat), 
                                       top_n = 20)
    
    if (!is.null(dotplot_obj)) {
      filename <- paste0("Dotplot_", gsub("[^A-Za-z0-9]", "_", comp$name), "_", go_cat)
      ggsave(file.path(output_dir, paste0(filename, ".pdf")), 
             dotplot_obj, width = 10, height = 8, dpi = 300)
      ggsave(file.path(output_dir, paste0(filename, ".png")), 
             dotplot_obj, width = 10, height = 8, dpi = 300)
      cat("Dotplot saved\n")
    }
    
    # Create GSEA plots
    cat("    Creating GSEA plots...\n")
    create_classic_gsea_plots(gsea_result, comp$label, go_cat, top_n = 5)
    
    # Create gene-concept network
    cat("    Creating gene-concept network...\n")
    create_gene_concept_network(gsea_result, comp$label, go_cat, 
                                show_category = 5, gene_list = gene_list)
  }
}


cat("Output directories:\n")
cat("  Main:", output_dir, "\n")
cat("  Classic GSEA plots:", gsea_classic_dir, "\n")
cat("  Gene-concept networks:", cnet_dir, "\n\n")

# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()
