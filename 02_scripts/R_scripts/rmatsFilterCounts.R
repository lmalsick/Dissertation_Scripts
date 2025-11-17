#############################################
# This script is used to filter rMATs outputs
# Needs to reference the txt files that come from the rMATs analysis
# I use the JCEC statistical model, though you can use JC instead just change the references to JCEC
# This current script filters by: 
#       FDR > 0.05
#       deltaPSI > 0.1 (or 10%)
#       MinCounts (set!)
# The top sections are to make the functions, please cite GeneStructureTools (DOI: 10.18129/B9.bioc.GeneStructureTools)
# Date: 11/17/2025
# Lauren Malsick
##############################################

library(data.table)
#Install then comment out
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GeneStructureTools")

library(GeneStructureTools)

#' Import RMATS Turbo results files as a rmatsDataSet
#' @param filePath path to RMATS differential splicing output files
#' @param type type of counts to use. "JC" or "JCEC"
#' @return rmatsDataSet
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
#' rmats_filePath <- system.file("extdata", "rmats_small/", package = "GeneStructureTools")
#' rds <- readRmatsDataSet(rmats_filePath)
readRmatsDataSet <- function(filePath, type = "JC") {
  rds <- new("rmatsDataSet", filePath = "~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/rMATs_analysis/rMATS_output_744_mock_72hpi")
  rmatsEventTypes <- c("SE", "MXE", "RI", "A3SS", "A5SS")
  rmatsFileList <- paste0(rmatsEventTypes, ".MATS.", type, ".txt")
  allFiles <- list.files(filePath, full.names = TRUE)
  diffSpliceFiles <- allFiles[basename(allFiles) %in% rmatsFileList]
  
  if (length(diffSpliceFiles) == 0) {
    stop("no rmats files in the specified filePath")
  } else if (!all(rmatsFileList %in% basename(diffSpliceFiles))) {
    for (f in rmatsFileList[which(!(rmatsFileList %in% basename(diffSpliceFiles)))]) {
      message(paste0("Can't find file: ", f, " please check if this file should exist in the filePath"))
    }
  }
  
  for (eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")) {
    if (any(basename(diffSpliceFiles) == paste0(eventType, ".MATS.", type, ".txt"))) {
      slot(rds, eventType) <- fread(diffSpliceFiles[basename(diffSpliceFiles) == rmatsFileList[rmatsEventTypes == eventType]], header = TRUE, data.table = FALSE)
    }
  }
  return(rds)
}



#' Filter out significant events from a RMATS dataset
#' @param rmatsDataSet rmatsDataSet generated from \code{readRmatsDataSet()}
#' @param FDR maximum FDR required to call event as significant
#' @param psiDelta minimum change in psi required to call an event as significant
#' @param idList (optional) list of gene ids to filter for
#' @param minCounts minumum number of counts for all replicates
#' in at least one condition to call an event as significant
#' @param medianCounts median count for all replicates
#' in at least one condition to call an event as significant
#' @param sampleTable data.frame with sample names and conditions.
#' Only needed if filtering with counts.
#' @return filtered rmatsDataSet
#' @export
#' @importFrom stats median
#' @import stringr
#' @family rmats data processing
#' @author Beth Signal
#' @examples
#' rmats_directory <- system.file("extdata", "rmats_small/", package = "GeneStructureTools")
#' rds <- readRmatsDataSet(rmats_directory)
#' rds.filtered <- filterRmatsEvents(rds, FDR = 0.01, psiDelta = 0.1)
#' # filter by gene name/id
#' rds.Tmem208 <- filterRmatsEvents(rds, idList = "Tmem208", FDR = 1, psiDelta = 0)
filterRmatsEvents <- function(rmatsDataSet,
                              FDR = 0.05,
                              psiDelta = 0.1,
                              idList = NA,
                              minCounts = NA,
                              medianCounts = NA,
                              sampleTable) {
  
  # set FDR/psiDelta if NA
  if (is.na(FDR[1])) {
    FDR <- 1
  }
  if (is.na(psiDelta[1])) {
    psiDelta <- 0
  }
  
  for (eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")) {
    tmp <- slot(rmatsDataSet, eventType)
    if (nrow(tmp) > 0) {
      significantEventsIndex <- which(tmp$FDR <= FDR & abs(tmp$IncLevelDifference) > psiDelta)
      slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex, ]
    }
  }
  
  if (!is.na(idList[1])) {
    for (eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")) {
      tmp <- slot(rmatsDataSet, eventType)
      if (nrow(tmp) > 0) {
        significantEventsIndex <- which(tmp$GeneID %in% idList | tmp$geneSymbol %in% idList)
        slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex, ]
      }
    }
  }
  
  if (!is.na(minCounts[1])) {
    for (eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")) {
      tmp <- slot(rmatsDataSet, eventType)
      if (nrow(tmp) > 0) {
        cond1 <- mapply(function(x, y) as.numeric(x) + as.numeric(y), x = stringr::str_split(tmp$IJC_SAMPLE_1, ","), y = stringr::str_split(tmp$SJC_SAMPLE_1, ","))
        cond2 <- mapply(function(x, y) as.numeric(x) + as.numeric(y), x = stringr::str_split(tmp$IJC_SAMPLE_2, ","), y = stringr::str_split(tmp$SJC_SAMPLE_2, ","))
        significantEventsIndex <- which(apply(rbind(cond1, cond2), 2, function(x) all(x >= minCounts)))
        slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex, ]
      }
    }
  }
  if (!is.na(medianCounts[1])) {
    for (eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")) {
      tmp <- slot(rmatsDataSet, eventType)
      if (nrow(tmp) > 0) {
        cond1 <- mapply(function(x, y) as.numeric(x) + as.numeric(y), x = stringr::str_split(tmp$IJC_SAMPLE_1, ","), y = stringr::str_split(tmp$SJC_SAMPLE_1, ","))
        cond2 <- mapply(function(x, y) as.numeric(x) + as.numeric(y), x = stringr::str_split(tmp$IJC_SAMPLE_2, ","), y = stringr::str_split(tmp$SJC_SAMPLE_2, ","))
        significantEventsIndex <- which(apply(cond1, 2, median) >= medianCounts | apply(cond2, 2, median) >= medianCounts)
        slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex, ]
      }
    }
  }
  
  return(rmatsDataSet)
}


#' Compare open reading frames for RMATS differentially spliced events
#' @param rmatsDataSet rmatsDataSet generated from \code{readRmatsDataSet()}
#' @param eventTypes which event type to filter for? default="all"
#' @param exons GRanges gtf annotation of exons
#' @param BSgenome BSGenome object containing the genome for the species analysed
#' @param NMD Use NMD predictions? (Note: notNMD must be installed to use this feature)
#' @param exportGTF file name to export alternative isoform GTFs (default=NULL)
#' @return data.frame containing significant RMATS differential splicing data and ORF change summaries
#' @export
#' @importFrom rtracklayer export.gff
#' @import GenomicRanges
#' @family rmats data processing
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata", "gencode.vM25.small.gtf", package = "GeneStructureTools"))
#' exons <- gtf[gtf$type == "exon"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' rmats_directory <- system.file("extdata", "rmats_small/", package = "GeneStructureTools")
#' rds <- readRmatsDataSet(rmats_directory)
#' rds.filtered <- filterRmatsEvents(rds, FDR = 0.01, psiDelta = 0.1)
#' rmats_summary <- rmatsTranscriptChangeSummary(rmatsDataSet = rds.filtered, exons, BSgenome = g)
rmatsTranscriptChangeSummary <- function(rmatsDataSet,
                                         exons = NULL,
                                         eventTypes = "all",
                                         BSgenome,
                                         NMD = TRUE,
                                         exportGTF = NULL) {
  if (eventTypes[1] == "all") {
    eventTypes <- c("SE", "MXE", "RI", "A3SS", "A5SS")
  }
  
  allTranscripts <- GRanges()
  orfChanges <- NULL
  
  diffSplice.SE.signif <- extractEvent(rmatsDataSet, "SE")
  if (nrow(diffSplice.SE.signif) > 0 & "SE" %in% eventTypes) {
    message(paste0("Creating isoforms for ", nrow(diffSplice.SE.signif), " SE events"))
    isoforms.SE <- skipExonByJunction(diffSplice.SE.signif, eventType = "SE", exons = exons)
    orfChanges.SE <- transcriptChangeSummary(isoforms.SE[isoforms.SE$set == "included_exon"],
                                             isoforms.SE[isoforms.SE$set == "skipped_exon"],
                                             BSgenome = BSgenome, NMD = NMD, exportGTF = exportGTF, dataSet = rmatsDataSet
    )
    m <- match(unlist(lapply(stringr::str_split(orfChanges.SE$id, "[-]"), "[[", 1)), diffSplice.SE.signif$ID)
    orfChanges <- rbind(orfChanges, cbind(diffSplice.SE.signif[m, c("ID", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference")], type = "SE", orfChanges.SE))
    allTranscripts <- c(allTranscripts, isoforms.SE)
  }
  
  diffSplice.MXE.signif <- extractEvent(rmatsDataSet, "MXE")
  if (nrow(diffSplice.MXE.signif) > 0 & "MXE" %in% eventTypes) {
    message(paste0("Creating isoforms for ", nrow(diffSplice.MXE.signif), " MXE events"))
    isoforms.MXE <- skipExonByJunction(diffSplice.MXE.signif, eventType = "MXE", exons = exons)
    orfChanges.MXE <- transcriptChangeSummary(isoforms.MXE[isoforms.MXE$set == "included_exon1"],
                                              isoforms.MXE[isoforms.MXE$set == "included_exon2"],
                                              BSgenome = BSgenome, NMD = NMD, exportGTF = exportGTF, dataSet = rmatsDataSet
    )
    m <- match(unlist(lapply(stringr::str_split(orfChanges.MXE$id, "[-]"), "[[", 1)), diffSplice.MXE.signif$ID)
    orfChanges <- rbind(orfChanges, cbind(diffSplice.MXE.signif[m, c("ID", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference")], type = "MXE", orfChanges.MXE))
    allTranscripts <- c(allTranscripts, isoforms.MXE)
  }
  
  diffSplice.RI.signif <- extractEvent(rmatsDataSet, "RI")
  if (nrow(diffSplice.RI.signif) > 0 & "RI" %in% eventTypes) {
    message(paste0("Creating isoforms for ", nrow(diffSplice.RI.signif), " RI events"))
    isoforms.RI <- altIntronRmats(diffSplice.RI.signif, exons = exons)
    orfChanges.RI <- transcriptChangeSummary(isoforms.RI[isoforms.RI$set == "spliced_intron"],
                                             isoforms.RI[isoforms.RI$set == "retained_intron"],
                                             BSgenome = BSgenome, NMD = NMD, exportGTF = exportGTF, dataSet = rmatsDataSet
    )
    m <- match(unlist(lapply(stringr::str_split(orfChanges.RI$id, "[-]"), "[[", 1)), diffSplice.RI.signif$ID)
    orfChanges <- rbind(orfChanges, cbind(diffSplice.RI.signif[m, c("ID", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference")], type = "RI", orfChanges.RI))
    allTranscripts <- c(allTranscripts, isoforms.RI)
  }
  
  diffSplice.A3SS.signif <- extractEvent(rmatsDataSet, "A3SS")
  if (nrow(diffSplice.A3SS.signif) > 0 & "A3SS" %in% eventTypes) {
    message(paste0("Creating isoforms for ", nrow(diffSplice.A3SS.signif), " A3SS events"))
    isoforms.A3SS <- altSpliceSiteRmats(diffSplice.A3SS.signif, eventType = "A3SS", exons = exons)
    orfChanges.A3SS <- transcriptChangeSummary(isoforms.A3SS[isoforms.A3SS$set == "alt3_splicesite_long"],
                                               isoforms.A3SS[isoforms.A3SS$set == "alt3_splicesite_short"],
                                               BSgenome = BSgenome, NMD = NMD, exportGTF = exportGTF, dataSet = rmatsDataSet
    )
    m <- match(unlist(lapply(stringr::str_split(orfChanges.A3SS$id, "[-]"), "[[", 1)), diffSplice.A3SS.signif$ID)
    orfChanges <- rbind(orfChanges, cbind(diffSplice.A3SS.signif[m, c("ID", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference")], type = "A3SS", orfChanges.A3SS))
    allTranscripts <- c(allTranscripts, isoforms.A3SS)
  }
  
  diffSplice.A5SS.signif <- extractEvent(rmatsDataSet, "A5SS")
  if (nrow(diffSplice.A5SS.signif) > 0 & "A5SS" %in% eventTypes) {
    message(paste0("Creating isoforms for ", nrow(diffSplice.A5SS.signif), " A5SS events"))
    isoforms.A5SS <- altSpliceSiteRmats(diffSplice.A5SS.signif, eventType = "A5SS", exons = exons)
    orfChanges.A5SS <- transcriptChangeSummary(isoforms.A5SS[isoforms.A5SS$set == "alt5_splicesite_long"],
                                               isoforms.A5SS[isoforms.A5SS$set == "alt5_splicesite_short"],
                                               BSgenome = BSgenome, NMD = NMD, exportGTF = exportGTF, dataSet = rmatsDataSet
    )
    m <- match(unlist(lapply(stringr::str_split(orfChanges.A5SS$id, "[-]"), "[[", 1)), diffSplice.A5SS.signif$ID)
    orfChanges <- rbind(orfChanges, cbind(diffSplice.A5SS.signif[m, c("ID", "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference")], type = "A5SS", orfChanges.A5SS))
    allTranscripts <- c(allTranscripts, isoforms.A5SS)
  }
  
  
  if (!is.null(exportGTF)) {
    rtracklayer::export.gff(allTranscripts,
                            con = exportGTF,
                            format = "gtf"
    )
  }
  
  return(orfChanges)
}

######################## After functions defined run this:
setwd("~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/rMATs_analysis/rMATS_output_744_mock_72hpi")

setClass(
  "rmatsDataSet",
  slots = list(
    filePath = "character",
    SE = "data.frame",
    MXE = "data.frame",
    RI = "data.frame",
    A3SS = "data.frame",
    A5SS = "data.frame"
  )
)

#here's where you edit PATH to .txt files
rds <- readRmatsDataSet(
       filePath = "~/Lauren_PhD/Geiss_Lab_PhD_Project/Bat_cell_project/Genewiz_sequencing/rMATs_analysis/rMATS_output_744_mock_72hpi",
       type = "JCEC")
before_counts <- sapply(slotNames(rds)[-1], function(x) {
  nrow(slot(rds, x))
})

#filter -- change to preferences, as a default I used FDR = 0.05, deltaPSI = 0.05, mincounts = 10
rds.filtered <- filterRmatsEvents(
  rmatsDataSet = rds,
  FDR = 0.05,
  psiDelta = 0.05,
  minCounts = 10
)

after_counts <- sapply(slotNames(rds.filtered)[-1], function(x) {
  nrow(slot(rds.filtered, x))
})


filter_summary <- data.frame(
  EventType = names(before_counts),
  Before = before_counts,
  After = after_counts,
  Removed = before_counts - after_counts,
  Percent_Removed = round(100 * (before_counts - after_counts) / before_counts, 2)
)

### Show results for filtering process to see how many reads were removed
filter_summary

######################
# Now you can make a CSV file of your reads so that you can then run rMATsGOAnalysis.R
# This will bind all splicing events into one csv files, so needs to make some extra columns for the information needed for other splicing
# Specifically MXE will make more columns for junction locations-- if not applicable like for SE will just have NA
######################

event_list <- lapply(slotNames(rds.filtered)[-1], function(eventType) {
  df <- slot(rds.filtered, eventType)
  if (nrow(df) > 0) df$EventType <- eventType
  df
})
names(event_list) <- slotNames(rds.filtered)[-1]

all_cols <- unique(unlist(lapply(event_list, names)))

event_list_std <- lapply(event_list, function(df) {
  missing_cols <- setdiff(all_cols, names(df))
  df[missing_cols] <- NA
  df <- df[all_cols]  # reorder for consistent rbind
  df
})

df_out <- do.call(rbind, event_list_std)

write.csv(df_out, "rMATS_filtered_events.csv", row.names = FALSE)
