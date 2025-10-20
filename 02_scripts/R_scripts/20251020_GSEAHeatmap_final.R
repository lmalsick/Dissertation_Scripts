#############################################
# Heatmap creation of GSEA data for WNV and WNV mutants
# Created: 10/20/2025
# Author: Lauren Malsick and Hunter Ogg
# Notes: Adapted from Hunter Ogg's Aaegypti script with a database created by EggnogMapper
# Known bug: PDF file crashes with the greek delta symbol but PNG works fine
# Files needed: .csv file outputs from DESeq2 
# Confirm your organism database has the correct keycol (aka "SYMBOL")
# Note 2: Not all libraries may be used but I don't have time to determine which ones I decided not to
##############################################

getwd()

setwd("C:/Users/lmals/Documents/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/GSEA_Artibeus")

getwd()

#The following library is a custom package which must be manually installed from source. 
#This can be performed with the following command, where "org.Aaegypti.eg.db" is the path to the org.Aaegypti.eg.db folder
#install.packages("org.Ajamaicensis.eg.db", repos=NULL,type="source")

library(org.Ajamaicensis.eg.db)

library(readr)
library(clusterProfiler)
library(tidyverse)
library(fgsea)
library(topGO)
library(Rgraphviz)
library(GenomicFeatures)
library(BiocParallel)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(purrr)
library(stringr)

#rnaFile is a csv which contains the results of differential gene expression. Generally the output of DEseq analysis
#rankMetric is the the choice of ranking for . Currently only supports a single, but more will be added in future
#goCat represents which of the three subontologies (CC,MF,BP) to include in the analysis. All three are included by default
rnaGSEA = function(rnaFile,orgData,rankMetric = "pval",goCat="All",minSize=15,maxSize=500,rankedList=c(),keycol="GID"){
  myDEresults=read.csv(rnaFile)
  if(rankMetric == "pval"){
    myDEresults=read.csv(rnaFile)
    myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
    newRank_pvalueAndFC = -log10(myDEresults$pvalue) *sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
    names(newRank_pvalueAndFC) = myDEresults$Gene
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  }
  if(rankMetric == "log2fc"){
    myDEresults=read.csv(rnaFile)
    myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
    newRank_pvalueAndFC = myDEresults$log2FoldChange #*sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
    names(newRank_pvalueAndFC) = myDEresults$Gene
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  }
  if(rankMetric == "both"){
    myDEresults=read.csv(rnaFile)
    myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
    newRank_pvalueAndFC = myDEresults$log2FoldChange * -log10(myDEresults$pvalue) #*sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
    names(newRank_pvalueAndFC) = myDEresults$Gene
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  }
  ## ----> ADD THIS BLOCK here to clean non-finite values:
  newRank_pvalueAndFC = newRank_pvalueAndFC[is.finite(newRank_pvalueAndFC)]
  
  set.seed(2025)
  egoCC <- gseGO(geneList     = newRank_pvalueAndFC,
                 OrgDb        = orgData,
                 keyType=keycol,
                 ont          = goCat,
                 minGSSize    = minSize,
                 maxGSSize    = maxSize,
                 pvalueCutoff = 0.5,
                 verbose      = TRUE)
  
  
  #head(egoCC)
  #goplot(egoCC)
  #dotplot(egoCC)+theme(axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"))#+theme(text = element_text(size = ))
  return(egoCC)
}

####################################
#All timepoints Heatmap

WT_72hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_744.72_vs_Mock.72_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")
WT_48hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_744.48_vs_Mock.48_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")
WT_24hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_744.24_vs_Mock.24_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")

sfRNA1_72hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_816.72_vs_Mock.72_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")
sfRNA1_48hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_816.48_vs_Mock.48_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")
sfRNA1_24hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_816.24_vs_Mock.24_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")

sfRNA2_72hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_853.72_vs_Mock.72_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")
sfRNA2_48hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_853.48_vs_Mock.48_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")
sfRNA2_24hpi=rnaGSEA("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/03_output/Virus_timepoint_to_mock_logfoldchange_0.5/virus_timepoint_853.24_vs_Mock.24_allgenes.csv",goCat="BP",orgData = org.Ajamaicensis.eg.db, keycol="SYMBOL")

# 1) Define your conditions and timepoints
conds   <- c("WT WNV", "sfRNA\u0394\u0031", "sfRNA\u0394\u0031&2")
timepts <- c("24hpi", "48hpi", "72hpi")

# 2) Pull together all your GSEA results into a list‑column
all_results <- tibble(
  cond = c(rep("WT WNV", 3), rep("sfRNA\u0394\u0031", 3), rep("sfRNA\u0394\u0031&2", 3)),
  tp = rep(c("24hpi", "48hpi", "72hpi"), 3),
  obj_name = c("WT_24hpi", "WT_48hpi", "WT_72hpi",
               "sfRNA1_24hpi", "sfRNA1_48hpi", "sfRNA1_72hpi",
               "sfRNA2_24hpi", "sfRNA2_48hpi", "sfRNA2_72hpi")
) %>%
  mutate(
    df = map(obj_name, ~ as.data.frame(get(.x)@result))
  ) %>%
  dplyr::select(cond, tp, df)

# 3) Compute signed –log10(FDR) per (cond, tp, pathway)
long_scores <- all_results %>%
  mutate(
    df2 = map(df, ~ .x %>%
                transmute(
                  ID,
                  Description,
                  score = sign(NES) * -log10(p.adjust)
                )
    )
  ) %>%
  dplyr::select(cond, tp, df2) %>%
  unnest(df2)

# 4) Pivot to a wide matrix: rows = pathways, cols = cond_tp
wide_scores <- long_scores %>%
  mutate(col = paste(cond, tp, sep = "_")) %>%
  dplyr::select(Description, col, score) %>%
  pivot_wider(names_from = col, values_from = score) %>%
  column_to_rownames("Description") %>%
  as.matrix()

# 5) Pick top 15 up and down pathways by absolute score
###### Here's where you could modify the # of pathways on the heatmap
summary_df <- long_scores %>%
  group_by(Description) %>%
  summarize(
    max_score = max(score, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE)
  ) %>%
  ungroup()

top_up   <- summary_df %>% arrange(desc(max_score)) %>% slice_head(n = 15) %>% pull(Description)
top_down <- summary_df %>% arrange(min_score)          %>% slice_head(n = 15) %>% pull(Description)
selected <- unique(c(top_up, top_down))

# subset to those rows
score_mat_sel <- wide_scores[selected, , drop = FALSE]

# ─────────────────────────────────────────────────────────────────────────────
# 6) Reorder the columns so WT WNV is first, then sfRNA1 and then sfRNA1&2
ordered_cols <- c(
  paste0("WT WNV_",    timepts),
  paste0("sfRNA\u0394\u0031_", timepts),
  paste0("sfRNA\u0394\u0031&2_", timepts)
)

score_mat_sel <- score_mat_sel[, ordered_cols]

score_mat_sel <- score_mat_sel[, ordered_cols]

# 7) This just adds a space to the column names so there's no longer an "_" underscore
# (Optional)
colnames(score_mat_sel) <- gsub("_", "\n", colnames(score_mat_sel))

# 8) (Optional) wrap long pathway names
rownames(score_mat_sel) <- str_wrap(rownames(score_mat_sel), width = 50)

# 9) Define colors of the heatmap
col_fun <- colorRamp2(c(-6, 0, 6), c("steelblue", "white", "firebrick"))

# 10a) Draw the heatmap WITH numbers (for personal use)
ht_with_numbers <- Heatmap(
  score_mat_sel,
  name                = "signed –log10(FDR)",
  col                 = col_fun,
  show_row_names      = TRUE,
  row_names_gp        = gpar(fontsize = 6, fontface = "bold", just = "center"),
  row_names_centered  = FALSE,
  row_names_max_width = unit(10, "cm"),
  show_column_names   = TRUE,
  column_names_gp     = gpar(fontsize = 8),
  column_names_rot    = 45,
  column_names_side   = "top",
  cluster_rows        = TRUE,
  cluster_columns     = FALSE,
  row_gap             = unit(1, "mm"),
  height              = unit(15, "cm"),
  width               = unit(8, "cm"),
  na_col              = "grey90",
  heatmap_legend_param = list(
    title          = "signed –log10(FDR)",
    at             = c(-6, 0, 6),
    direction      = "horizontal",
    title_position = "leftcenter"
  ),
  layer_fun = function(j, i, x, y, w, h, fill, ...) {
    vals <- score_mat_sel[cbind(i, j)]
    not_na <- !is.na(vals)
    if(any(not_na)) {
      grid.text(sprintf("%.2f", vals[not_na]), 
                x[not_na], y[not_na], 
                gp = gpar(fontsize = 7))
    }
  }
)

# 10b) Heatmap WITHOUT numbers
ht_without_numbers <- Heatmap(
  score_mat_sel,
  name                = "signed –log10(FDR)",
  col                 = col_fun,
  show_row_names      = TRUE,
  row_names_gp        = gpar(fontsize = 6, fontface = "bold", just = "center"),
  row_names_centered  = FALSE,
  row_names_max_width = unit(10, "cm"),
  show_column_names   = TRUE,
  column_names_gp     = gpar(fontsize = 8),
  column_names_rot    = 45,
  column_names_side   = "top",
  cluster_rows        = TRUE,
  cluster_columns     = FALSE,
  row_gap             = unit(1, "mm"),
  height              = unit(15, "cm"),
  width               = unit(8, "cm"),
  na_col              = "grey90",
  heatmap_legend_param = list(
    title          = "signed –log10(FDR)",
    at             = c(-6, 0, 6),
    direction      = "horizontal",
    title_position = "leftcenter"
  )
)

# 11a) Save WITH numbers - PDF
pdf("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/GSEA_Artibeus/gsea_heatmap_with_numbers.pdf", 
    width = 10, height = 12)
draw(ht_with_numbers, heatmap_legend_side = "bottom")
grid::grid.text(
  "Gray cells = not significant (FDR > 0.05)", 
  x = unit(0.98, "npc"), 
  y = unit(0.02, "npc"),
  just = c("right", "bottom"),
  gp = gpar(fontsize = 8)
)
dev.off()

# 11a) Save WITH numbers - PNG
png("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/GSEA_Artibeus/gsea_heatmap_with_numbers.png", 
    width = 2400, height = 2400, res = 300)
draw(ht_with_numbers, heatmap_legend_side = "bottom")
grid::grid.text(
  "Gray cells = not significant (FDR > 0.05)", 
  x = unit(0.98, "npc"), 
  y = unit(0.02, "npc"),
  just = c("right", "bottom"),
  gp = gpar(fontsize = 8)
)
dev.off()

# 11b) Save WITHOUT numbers - PDF
pdf("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/GSEA_Artibeus/gsea_heatmap_without_numbers.pdf", 
    width = 10, height = 12)
draw(ht_without_numbers, heatmap_legend_side = "bottom")
grid::grid.text(
  "Gray cells = not significant (FDR > 0.05)", 
  x = unit(0.98, "npc"), 
  y = unit(0.02, "npc"),
  just = c("right", "bottom"),
  gp = gpar(fontsize = 8)
)
dev.off()

# 11b) Save WITHOUT numbers - PNG
png("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/GSEA_Artibeus/gsea_heatmap_without_numbers.png", 
    width = 2400, height = 2800, res = 300)
draw(ht_without_numbers, heatmap_legend_side = "bottom")
grid::grid.text(
  "Gray cells = not significant (FDR > 0.05)", 
  x = unit(0.98, "npc"), 
  y = unit(0.02, "npc"),
  just = c("right", "bottom"),
  gp = gpar(fontsize = 8)
)
dev.off()
write.csv(score_mat_sel, file = "~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/GSEA_Artibeus/heatmap_scores.csv", row.names = TRUE)

###################################### END

# Will now have saved a PDF/PNG of the heatmap with the values printed on, and another heatmap
# with no numbers saved to PDF/PNG
# Reminder: This heatmap uses -log10(FDR) with the sign of NES
# SMALLER FDR = larger -log10(FDR) so more highly significant, this is then multiplied by the sign of NES to indicate direction of regulation

# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()

