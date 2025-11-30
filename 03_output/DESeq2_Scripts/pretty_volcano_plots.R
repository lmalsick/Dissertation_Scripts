#######################################

#__Date:__ December 2nd 2024
#__Author:__ Erin Osborne Nishimura (edited by Lauren Malsick)
#__Script:__ 20241202_DEseq2_loopdraft.R
#__Project:__ To analyze RNA-seq data for Genewiz sequencing of AJ6 cells infected with 744/816/853/ENTV or mock 
#second project with IFN alpha pre-treatment in comparison to mock
#__Requires:__ 
# 
# + R (4.4.1)
# + DESeq2 (1.44.0)   
# + corrplot (0.95)
# + RColorBrewer (1.1-3)
# + pheatmap (1.0.12)
# + apeglm (1.20.0)
# Need: metadata_infected.txt file and counts_infected.txt (zipped in github folder)
######################################


######### FOR FIRST TIME USE ONLY ##############
######### After use, comment this section ##############

# If you don't have bioconductor, install it. This version works with R version 4.1. Use "3.11" for R version 4.0. Use "3.18" for R version 4.3. Use "3.20" for R 4.4
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")
# 
# 
# # Install required packages/libraries:
# 
# # Install DESeq2:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# 
# 
# # Install 'apeglm'
# BiocManager::install("apeglm")
# 
# # install corrplot:
# install.packages("corrplot")
# 
# #Install pretty heatmap - pheatmap: https://cran.r-project.org/web/packages/pheatmap/index.html
# install.packages("pheatmap")
# 
# #Install RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
# install.packages("RColorBrewer")
# 
# install.packages("shiny")
# 
# install.packages("plotly")
# 
# if (!requireNamespace("rentrez", quietly = TRUE)) {
#   install.packages("rentrez")
# }
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

###########  LOAD PACKAGE  #####################

# Do this every time.
# Load packages:
library(DESeq2)
library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(apeglm)
library(ggplot2)
library(shiny)
library(plotly)
library(VennDiagram)
library(ggvenn)
library(trelliscopejs)
library(plotly)
library(htmlwidgets)
library(venn)
library(clusterProfiler)
library(BiocManager)

################################################

getwd()

setwd()

getwd()

countsData <- read.table(file = "../01_input/counts_infected.txt", header = FALSE, row.names = 1, skip = 2) # 

head(countsData)
dim(countsData)
class(countsData)

# Read in the metadata
metadata1 <- read.table(file = "../01_input/metadata_infected.txt", header = FALSE) # import the data

metadata1


# Organize the metadata file by adding some column headers
colnames(metadata1) <- c("fasta1", "fasta2", "names1", "names2", "timepoint", "virus", "rep")
metadata1

# The last names will be names for each sample. We can pull those names from metadata1:
as.vector(metadata1$names2)

# Name countsData columns headers:
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata1$names2))

# OK, our task will be to generate a table called "cts" out of the countsData table.
# Subset the countsData 
head(countsData)
dim(countsData)
head(countsData[,6:65])
dim(countsData[,6:65])


# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:65])
head(cts)

metadata1
rownames(metadata1)<- metadata1$names2
metadata1

coldata <- metadata1[,c("timepoint", "virus", "rep")]
coldata$virus <- as.factor(coldata$virus)
coldata$rep <- as.factor(coldata$rep)
coldata

rownames(coldata)
colnames(cts)
class(coldata)
class(cts)

#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))
head(colnames(cts)) 

################################Pairwise Timepoint Comparisons #########################

coldata$virus_timepoint <- interaction(coldata$virus, coldata$timepoint, sep=".")

ddsPairwise <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ virus_timepoint)

ddsPairwise <- DESeq(ddsPairwise)

vsd <- vst(ddsPairwise, blind = FALSE)

plotPCA(vsd, intgroup=c("virus", "timepoint"))

resPairwise <- results(ddsPairwise)

resultsNames(ddsPairwise)
summary(resPairwise)

############################# For loop begins here ###########################
#This will go through each of the listed combinations, relisting the target conditions and performing deSeq and shrinkage
# A labeled volcano plot, and plotly volcano will be produced for each target condition

# Define the target conditions as pairs of "virus_timepoint" and "Mock" levels

target_conditions <- list(
  c("744.12", "Mock.12"),
  c("744.24", "Mock.24"),
  c("744.48", "Mock.48"),
  c("744.72", "Mock.72"),
  c("816.12", "Mock.12"),
  c("816.24", "Mock.24"),
  c("816.48", "Mock.48"),
  c("816.72", "Mock.72"),
  c("853.12", "Mock.12"),
  c("853.24", "Mock.24"),
  c("853.48", "Mock.48"),
  c("853.72", "Mock.72"),
  c("ENTV.12", "Mock.12"),
  c("ENTV.24", "Mock.24"),
  c("ENTV.48", "Mock.48"),
  c("ENTV.72", "Mock.72")
)

upregulated_genes_by_virus_timepoint <- list()
downregulated_genes_by_virus_timepoint <- list()

# Iterate over each target condition pair
for (condition_pair in target_conditions) {
  # Extract the mock and virus timepoint for the pair
  mock_level <- condition_pair[1]
  timepoint_level <- condition_pair[2]
  
  # Relevel the 'virus_timepoint' factor to use the specific timepoint as the reference
  ddsPairwise$virus_timepoint <- relevel(ddsPairwise$virus_timepoint, ref = timepoint_level)
  
  # Re-run DESeq2 with the new reference level
  ddsPairwise <- DESeq(ddsPairwise)
  
  # Removes samples with less than 10 reads
  keep <- rowSums(counts(ddsPairwise)) >= 10
  ddsPairwise <- ddsPairwise[keep, ]
  
  # Get the result for the comparison between the mock level and the timepoint level
  result_name <- paste0("virus_timepoint_", mock_level, "_vs_", timepoint_level)
  res <- results(ddsPairwise, name = result_name)
  
  # Perform LFC shrinkage via apeglm
  res_shrunk <- lfcShrink(ddsPairwise, coef = result_name, type = "apeglm")
  
  # Convert res_shrunk to a data frame to use dplyr functions
  res_shrunk_df <- as.data.frame(res_shrunk)
  
  # Add a column for labeling significance
  res_shrunk_df$significance <- "Not Significant"
  res_shrunk_df$significance[res_shrunk_df$padj < 0.01 & res_shrunk_df$log2FoldChange > 2] <- "Upregulated"
  res_shrunk_df$significance[res_shrunk_df$padj < 0.01 & res_shrunk_df$log2FoldChange < -2] <- "Downregulated"
  
  # Calculate -log10(padj)
  res_shrunk_df$neg_log10_padj <- -log10(res_shrunk_df$padj)
  
  # Add labels for significant points
  res_shrunk_df$label <- ifelse(
    res_shrunk_df$padj < 0.01 & abs(res_shrunk_df$log2FoldChange) > 2, 
    rownames(res_shrunk_df), 
    NA
  )
  
  # Filter out NAs in your data for label assignments (only for labels)
  volcano_data_clean <- res_shrunk_df %>%
    filter(!is.na(label))  # Only keep rows where the label is not NA for labeling
  
  # Volcano Plot
  volcano_plot <- ggplot(res_shrunk_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    # Plot all points, color-coded by significance
    geom_point(aes(color = significance), alpha = 0.6, size = 1.2) +
    
    # Add text labels only for significant points
    geom_text_repel(
      data = volcano_data_clean,  # Only significant points will have labels
      aes(label = label), size = 3, max.overlaps = 5
    ) +
    
    # Custom colors for significance
    scale_color_manual(
      values = c("Upregulated" = "forestgreen", "Downregulated" = "darkorchid4", "Not Significant" = "gray")
    ) +
    
    # Add threshold lines (dashed)
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred") +  # p-value threshold
    
    # Title, axis labels, and limits
    labs(
      title = paste("Volcano Plot:", result_name),  # Add a title
      x = "Log2 Fold Change",
      y = "-log10(padj)"
    ) +
    
    # Adjust the plot range
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 100)) +
    
    # Apply a clean theme with a white background
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),  # White panel background
      plot.background = element_rect(fill = "white", color = NA),       # White overall background
      panel.grid.major = element_blank(),                               # Remove major gridlines
      panel.grid.minor = element_blank(),                               # Remove minor gridlines
      axis.line = element_line(color = "black"),                        # Add black axis lines
      axis.text = element_text(color = "black"),                     
      axis.title = element_text(color = "black")                        # Black axis titles
    )
  
  # Save the volcano plot as a PDF
  ggsave(paste0("../03_output/", result_name, "_volcano_plot.pdf"), plot = volcano_plot, width = 8, height = 6)
  
  # Convert ggplot to plotly for interactivity
  interactive_volcano <- ggplotly(volcano_plot, tooltip = "text")
  
  # Save the interactive plot as an HTML file
  htmlwidgets::saveWidget(interactive_volcano, file = paste0("../03_output/", result_name, "_volcano_plot.html"))
}


# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()
