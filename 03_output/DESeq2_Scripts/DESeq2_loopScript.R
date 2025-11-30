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
# Loops through all pairwise comparisons and makes an MA plot, histogram, volcano plot, and csv file
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

countsDataGSEA_csv <- "../03_output/countsDataGSEA.csv"

# Write the countsData data frame to a CSV file
write.csv(countsData, file = countsDataGSEA_csv, row.names = TRUE)

################################Pairwise Timepoint Comparisons #########################

coldata$virus_timepoint <- interaction(coldata$virus, coldata$timepoint, sep=".")

ddsPairwise <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ virus_timepoint)

ddsPairwise <- DESeq(ddsPairwise)

vsd <- vst(ddsPairwise, blind = FALSE)

plotPCA(vsd, intgroup=c("virus", "timepoint"))


############################# For loop begins here ###########################
#This will go through each of the listed combinations, relisting the target conditions and performing deSeq and shrinkage
# An MA plot, histogram, volcano plot, and csv file (with up and down regulated genes) will be produced for each target condition

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
  
  #removes samples with less than 10 reads
  keep <- rowSums(counts(ddsPairwise)) >= 10
  
  ddsPairwise <- ddsPairwise[keep,]
  
  # Get the result for the comparison between the mock level and the timepoint level
  result_name <- paste0("virus_timepoint_", mock_level, "_vs_", timepoint_level)
  res <- results(ddsPairwise, name = result_name)
  
  # Perform LFC shrinkage for apeglm
  res_shrunk <- lfcShrink(ddsPairwise, coef = result_name, type = "apeglm")
  
  histogram_output_file <- paste0("../03_output/", result_name, "_baseMean_histogram.pdf")
  pdf(histogram_output_file)
  
  # Apply log10 transformation to baseMean values
  log_baseMean <- log10(res_shrunk$baseMean)
  
  # Plot the histogram of log10-transformed baseMean values
  hist(log_baseMean, breaks = 50, main = paste("Distribution of log10(baseMean) for", result_name),
       xlab = "log10(baseMean)", col ="lightblue", border = "black")
  
  abline(v = 10, col = "red", lwd = 0.5, lty = 0.5) # Example threshold
  dev.off()
  #We could remove low base mean counts here with this code: filtered_res <- res_shrunk[res_shrunk$baseMean > 20, ]
  #but some online argue that baseMean is a bit uninterpretable in this way and encourage just filtering out samples 
  #that have less than 10 reads which I did before the histogram
  
  # Filter results to include only those with LFC > 0.5 (upregulated) and LFC < -0.5 (downregulated)
  ############### You can change the filtering here to be more stringent for log2FC
  res_up <- subset(res_shrunk, padj < 0.01 & log2FoldChange > 0.5)
  res_down <- subset(res_shrunk, padj < 0.01 & log2FoldChange < -0.5)
  
  # Store upregulated genes in the list, using timepoint as the key
  upregulated_genes_by_virus_timepoint[[paste(mock_level, timepoint_level, sep = "_vs_")]] <- rownames(res_up)
  
  # Store downregulated genes in the list, using timepoint as the key
  downregulated_genes_by_virus_timepoint[[paste(mock_level, timepoint_level, sep = "_vs_")]] <- rownames(res_down)
  
  # Create a dataframe for upregulated genes
  upregulated_genes_df <- data.frame(
    Gene = rownames(res_up),
    baseMean = res_up$baseMean,
    log2FoldChange = res_up$log2FoldChange,
    lfcSE = res_up$lfcSE,
    pvalue = res_up$pvalue,
    padj = res_up$padj
  )
  
  # Create a dataframe for downregulated genes
  downregulated_genes_df <- data.frame(
    Gene = rownames(res_down),
    baseMean = res_down$baseMean,
    log2FoldChange = res_down$log2FoldChange,
    lfcSE = res_down$lfcSE,
    pvalue = res_down$pvalue,
    padj = res_down$padj
  )
  
  # Define output filenames for CSVs
  upregulated_csv <- paste0("../03_output/Virus_timepoint_to_mock_logfoldchange_0.5/", result_name, "_upregulated_genes.csv")
  downregulated_csv <- paste0("../03_output/Virus_timepoint_to_mock_logfoldchange_0.5/", result_name, "_downregulated_genes.csv")
  allGenes_csv <- paste0("../03_output/Virus_timepoint_to_mock_logfoldchange_0.5/", result_name, "_allGenes.csv")
  
  # Save the upregulated and downregulated genes to CSV files
  write.csv(upregulated_genes_df, file = upregulated_csv)
  write.csv(downregulated_genes_df, file = downregulated_csv)
  write.csv(res_shrunk, file = allGenes_csv)
  
  ######################### Create and save MA plots
  
  pdf(paste0("../03_output/", result_name, "_MAplots.pdf"), width = 14, height = 7)
  
  par(mfrow = c(1, 2)) # Set up a side-by-side plotting layout
  
  # MA plot for unshrunken results
  plotMA(res, main = paste(result_name, "\nunshrunken"),
         ylim = c(-7, 10), ylab = "log2 fold change",
         xlab = "mean of normalized counts")
  
  # MA plot for shrunken results
  plotMA(res_shrunk, main = paste(result_name, "\nshrunken"),
         ylim = c(-7, 10), ylab = "log2 fold change",
         xlab = "mean of normalized counts")
  
  dev.off() 
  
  
  ##### Volcano Plot creation
  
  significantLFC <- rbind(
    subset(res_shrunk, padj < 0.01 & log2FoldChange > 0.5 ),
    subset(res_shrunk, padj < 0.01 & log2FoldChange < -0.5 )
  )
  
  # Get the significant points to plot
  significant_points_to_plot <- res_shrunk[which(rownames(res_shrunk) %in% rownames(significantLFC)), ]
  
  # Maxed out points (p-adj very small)
  maxedout <- subset(res_shrunk, padj < 10e-100)
  
  # Create and save the volcano plot as a PDF
  pdf(paste0("../03_output/", result_name, "_volcano_plot.pdf"))
  
  # Draw the volcano plot
  plot(
    res_shrunk$log2FoldChange, -log10(res_shrunk$padj),
    main=paste("Volcano plot", result_name), 
    xlab="Effect size: log2(fold-change)", 
    ylab="-log10(adjusted p-value)", 
    pch=20, cex = 0.4, ylim = c(0, 50), xlim = c(-6,6), 
    col = "#00000030"
  )
  
  # Add points for maxed out values (those with extremely small p-values)
  points(
    maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)),
    pch=17, cex = 0.4, ylim = c(0, 100), col = "darkorchid4"
  )
  #darkorchid4 for log2fc 0.5
  # Highlight significant points (up and downregulated)
  points(
    significant_points_to_plot$log2FoldChange, -log10(significant_points_to_plot$padj),
    pch=20, cex = 0.4, col = "darkorchid4"
  )
  
  # Add threshold lines for visualization
  abline(v=0, col = "darkred")  # Fold change = 1
  abline(v=0.5, col = "darkred", lty = "dashed") 
  abline(v=-0.5, col = "darkred", lty = "dashed")  
  abline(h=-log10(0.01), col = "darkred", lty = "dashed")  # p-value threshold
  
  # Close the PDF
  dev.off()
  
  interactive_volcano_file <- paste0("../03_output/", result_name, "_volcano_plot.html")
  
  # Prepare data for Plotly volcano plot
  plotly_data <- data.frame(
    Gene = rownames(res_shrunk),
    log2FoldChange = res_shrunk$log2FoldChange,
    padj = res_shrunk$padj,
    negLog10Padj = -log10(res_shrunk$padj)
  )
  
  # Add a significance column
  plotly_data$Significance <- ifelse(
    plotly_data$padj < 0.01 & abs(plotly_data$log2FoldChange) > 0.5, 
    "Significant", 
    "Not Significant"
  )
  
  # Create the Plotly volcano plot
  volcano_plot <- plot_ly(
    data = plotly_data,
    x = ~log2FoldChange,
    y = ~negLog10Padj,
    text = ~paste("Gene:", Gene, "<br>log2FC:", round(log2FoldChange, 0.5), "<br>p-adj:", signif(padj, 3)),
    color = ~Significance,
    colors = c("Not Significant" = "gray", "Significant" = "forestgreen"),
    type = "scatter",
    mode = "markers",
    marker = list(size = 5, opacity = 0.7)
  ) %>%
    layout(
      title = paste("Interactive Volcano Plot for", result_name),
      xaxis = list(title = "Effect size: log2(Fold Change)", range = c(-6, 6)),
      yaxis = list(title = "-log10(Adjusted P-Value)", range = c(0, 50)),
      shapes = list(
        list(type = "line", x0 = 2, x1 = 2, y0 = 0, y1 = 50, line = list(dash = "dash", color = "red")),
        list(type = "line", x0 = -2, x1 = -2, y0 = 0, y1 = 50, line = list(dash = "dash", color = "red")),
        list(type = "line", x0 = -6, x1 = 6, y0 = -log10(0.01), y1 = -log10(0.01), line = list(dash = "dash", color = "red"))
      )
    )
  
  # Save the volcano plot as an HTML file
  htmlwidgets::saveWidget(volcano_plot, file = interactive_volcano_file)
  #Function to create and save Venn diagrams
  # Function to create and save Venn diagrams and CSVs for each intersection
  save_venn_and_genes <- function(upregulated_genes, downregulated_genes, regulation_type, virus) {
    # Define the correct order of timepoints
    timepoints <- c("12", "24", "48", "72")
    
    upregulated_genes <- upregulated_genes[timepoints]
    downregulated_genes <- downregulated_genes[timepoints]
    
    # Filter out empty lists for upregulated and downregulated genes
    upregulated_genes <- upregulated_genes[sapply(upregulated_genes, length) > 0]
    downregulated_genes <- downregulated_genes[sapply(downregulated_genes, length) > 0]
    
    # Check if there are any valid gene sets
    if (length(upregulated_genes) == 0 && length(downregulated_genes) == 0) {
      print(paste0("No genes for ", regulation_type, " regulation in this comparison"))
      return(NULL)  # Skip if both gene lists are empty
    }
    
    # If there are upregulated genes, create the Venn diagram
    if (length(upregulated_genes) > 0) {
      venn.fill.colors <- c("purple", "orange", "cyan", "brown")[1:length(upregulated_genes)]  # Adjust color length based on the sets
      venn.diagram(upregulated_genes, 
                   filename = paste0("../03_output/venn_diagrams/venn_diagram_upregulated_", virus, ".png"),
                   main = paste("Venn diagram for upregulated genes for", virus),
                   fill = venn.fill.colors,  # Set fill colors based on the number of sets
                   alpha = 0.5, cex = 1, fontface = "bold", 
                   cat.fontface = "bold", cat.cex = 1.2)  # Customize further as needed
      
      # Save each intersection of upregulated genes as CSV
      save_venn_intersections(upregulated_genes, virus, "upregulated")
    }
    
    # If there are downregulated genes, create the Venn diagram
    if (length(downregulated_genes) > 0) {
      venn.fill.colors <- c("purple", "orange", "cyan", "brown")[1:length(downregulated_genes)]  # Adjust color length based on the sets
      venn.diagram(downregulated_genes, 
                   filename = paste0("../03_output/venn_diagrams/venn_diagram_downregulated_", virus, ".png"),
                   main = paste("Venn diagram for downregulated genes for", virus),
                   fill = venn.fill.colors,  # Set fill colors based on the number of sets
                   alpha = 0.5, cex = 1, fontface = "bold", 
                   cat.fontface = "bold", cat.cex = 1.2)  # Customize further as needed
      
      # Save each intersection of downregulated genes as CSV
      save_venn_intersections(downregulated_genes, virus, "downregulated")
    }
    
    # Print the Venn data for debugging purposes
    cat("Venn diagram for", virus, "created with the following gene counts:\n")
    print(sapply(upregulated_genes, length))
    print(sapply(downregulated_genes, length))
  }
  
  # Function to save the intersections of the Venn diagram as CSVs
  save_venn_intersections <- function(gene_sets, virus, regulation_type) {
    # Define the output directory
    output_dir2 <- "../03_output/intersections/"
    
    # Create the directory if it doesn't exist
    if (!dir.exists(output_dir2)) {
      dir.create(output_dir2, recursive = TRUE)  # Create all necessary parent directories
    }
    # Get the names of the sets (e.g., "12", "24", "48", "72")
    set_names <- names(gene_sets)
    
    # Create a list to store the intersections
    intersections <- list()
    
    # Loop over all possible combinations of the gene sets
    for (i in 1:length(set_names)) {
      for (j in i:length(set_names)) {
        # Get the intersection of the two sets
        intersection_name <- paste(set_names[i], set_names[j], sep = "_and_")
        intersection_genes <- intersect(gene_sets[[set_names[i]]], gene_sets[[set_names[j]]])
        intersections[[intersection_name]] <- intersection_genes
      }
    }
    # Compute the intersection of all sets (the center of the Venn diagram)
    all_intersection_name <- paste(set_names, collapse = "_and_")
    all_intersection_genes <- Reduce(intersect, gene_sets)
    intersections[[all_intersection_name]] <- all_intersection_genes
    
    # Save the intersections as CSV files
    for (intersection_name in names(intersections)) {
      intersection_genes <- intersections[[intersection_name]]
      
      if (length(intersection_genes) > 0) {
        intersection_df <- data.frame(Gene = intersection_genes)
        # Define the file path
        file_path <- file.path(output_dir2, paste0("venn_", regulation_type, "_genes_", virus, "_", intersection_name, ".csv"))
        
        # Write the genes to the CSV file
        write.csv(intersection_df, file_path, row.names = FALSE)
      } else {
        cat(paste("No genes to save for intersection:", intersection_name, "\n"))
      }
    }
  }
  
  # Loop over each virus and generate Venn diagrams and CSVs for upregulated and downregulated genes across timepoints
  for (virus in c("744", "816", "853", "ENTV")) {
    
    # Prepare lists for upregulated and downregulated genes across the timepoints
    upregulated_genes <- list()
    downregulated_genes <- list()
    
    # Get upregulated and downregulated genes for each timepoint
    for (timepoint in c("12", "24", "48", "72")) {
      result_name <- paste0(virus, ".", timepoint, "_vs_", "Mock.", timepoint)
      
      # Retrieve the upregulated and downregulated genes from the pre-existing results
      genes_up <- upregulated_genes_by_virus_timepoint[[result_name]]
      genes_down <- downregulated_genes_by_virus_timepoint[[result_name]]
      
      # Store the genes in the respective lists (only if the lists are not empty)
      if (length(genes_up) > 0) upregulated_genes[[timepoint]] <- genes_up
      if (length(genes_down) > 0) downregulated_genes[[timepoint]] <- genes_down
    }
    
    # Create the Venn diagram and save the intersections for upregulated and downregulated genes across all timepoints for the current virus
    save_venn_and_genes(upregulated_genes, downregulated_genes, "regulated", virus)
  }
} 


# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()
