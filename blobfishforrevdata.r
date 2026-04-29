library(netZooR)
library(fgsea)
library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
library(grid)
library(igraph)
library(org.Hs.eg.db)

# Set source directory
sourceDir <- NULL
monsterDir <- NULL
pathwayFile <- NULL
nullDir <- NULL

# Read the PANDAs
LGGFemale <- read.table(paste0(sourceDir, 'LGG_Female_Panda_output.txt'), sep = " ")
GBMFemale <- read.table(paste0(sourceDir, 'GBM_Female_Panda_output.txt'), sep = " ")
LGGMale <- read.table(paste0(sourceDir, 'LGG_Male_Panda_output.txt'), sep = " ")
GBMMale <- read.table(paste0(sourceDir, 'GBM_Male_Panda_output.txt'), sep = " ")

# Read the TFs
tfsDF <- read.csv(paste0(monsterDir, "gbmFemaleSpecificMONSTER.csv"), row.names = 1)
tfs <- tfsDF[which(tfsDF[,2] != "Remove"), "x"]

# Filter the PANDAs
filterPanda <- function(panda){
  pandaFilt <- panda[which(panda$V1 %in% tfs),]
  pandaFilt <- pandaFilt[,c(1,2,4)]
  colnames(pandaFilt) <- c("tf", "gene", "score")
  return(pandaFilt)
}

lggFemaleFilt <- filterPanda(LGGFemale)
gbmFemaleFilt <- filterPanda(GBMFemale)
lggMaleFilt <- filterPanda(LGGMale)
gbmMaleFilt <- filterPanda(GBMMale)

# Convert ENSEMBL to gene symbols
ensemblToSymbol <- function(data){
  # Map ENSEMBL to SYMBOL
  results <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                                   keys = data$gene, 
                                   column = "SYMBOL", 
                                   keytype = "ENSEMBL",
                                   multiVals = "first")
  
  # Create a new data frame
  newData <- data.frame(tf = data$tf, 
                       gene = results, 
                       score = data$score,
                       stringsAsFactors = FALSE)
  
  # Remove NAs
  newData <- newData[!is.na(newData$gene), ]
  
  pairNames <- paste(newData$tf, newData$gene, sep = "__")

  # Remove duplicate pairs in the data frame
  pairCounts <- table(pairNames)
  dupPairs <- names(pairCounts)[which(pairCounts > 1)]
  
  # Only process duplicates if they exist
  if(length(dupPairs) > 0){
    toRemove <- unlist(lapply(1:length(dupPairs), function(i){
      pair <- dupPairs[i]
      whichPairNames <- which(pairNames == pair)
      whichMaxScore <- which.max(newData[whichPairNames, "score"])
      whichToRemove <- whichPairNames[setdiff(1:length(whichPairNames), whichMaxScore)]
      if(i %% 100 == 0) print(paste(i, "out of", length(dupPairs)))
      return(whichToRemove)
    }))
    newDataDedup <- newData[setdiff(1:nrow(newData), toRemove),]
  } else {
    newDataDedup <- newData
  }

  rownames(newDataDedup) <- paste(newDataDedup$tf, newDataDedup$gene, sep = "__")
  return(newDataDedup)
}

lggFemaleTarget <- ensemblToSymbol(lggFemaleFilt)
write.csv(lggFemaleTarget, paste0(sourceDir, "lggFemaleTargetFull.csv"))

gbmFemaleTarget <- ensemblToSymbol(gbmFemaleFilt)
write.csv(gbmFemaleTarget, paste0(sourceDir, "gbmFemaleTargetFull.csv"))

lggMaleTarget <- ensemblToSymbol(lggMaleFilt)
write.csv(lggMaleTarget, paste0(sourceDir, "lggMaleTargetFull.csv"))

gbmMaleTarget <- ensemblToSymbol(gbmMaleFilt)
write.csv(gbmMaleTarget, paste0(sourceDir, "gbmMaleTargetFull.csv"))

# Read pathway result files
GBMM_GBMF <- read.csv(paste0(sourceDir, "GBMMale_Female_PathwayFile.csv"))
GBMF_GBMM <- read.csv(paste0(sourceDir, "GBMFemale_Male_PathwayFile.csv"))

# Get all genes in the pathways of interest
gmtFile = fgsea::gmtPathways(pathwayFile)

mrnaPathways <- c("REACTOME_METABOLISM_OF_RNA", "REACTOME_SPLICING", "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA")
arPathways <- "PID_AR_NONGENOMIC_PATHWAY"
immunePathways <- c("REACTOME_NEUTROPHIL_DEGRANULATION", "REACTOME_INNATE_IMMUNE_SYSTEM", "KEGG_LYSOSOME")
carbPathways <- c("WP_METABOLIC_PATHWAYS_OF_FIBROBLASTS", "WP_AEROBIC_GLYCOLYSIS", "WP_N_GLYCAN_BIOSYNTHESIS")
ecmPathways <- c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_COLLAGEN_FORMATION")
cancerPathways <- c("WP_TYPE_2_PAPILLARY_RENAL_CELL_CARCINOMA", "WP_CLEAR_CELL_RENAL_CELL_CARCINOMA_PATHWAYS")
hypoxiaPathways <- "PID_HIF1_TFPATHWAY"

getGenesInPathways <- function(pathwayNames, pathwayFile){
  geneStrings <- pathwayFile[which(pathwayFile$pathway %in% pathwayNames), "leadingGenesDrivingEnrichment"]
  geneLists <- lapply(geneStrings, function(string){return(strsplit(string, "; ")[[1]])})
  geneSet <- unique(unlist(geneLists))
  return(geneSet)
}

mrnaGenes <- getGenesInPathways(mrnaPathways, GBMM_GBMF)
arGenes <- getGenesInPathways(arPathways, GBMM_GBMF)
immuneGenes <- getGenesInPathways(immunePathways, GBMF_GBMM)
carbGenes <- getGenesInPathways(carbPathways, GBMF_GBMM)
ecmGenes <- getGenesInPathways(ecmPathways, GBMF_GBMM)
cancerGenes <- getGenesInPathways(cancerPathways, GBMF_GBMM)
hypoxiaGenes <- getGenesInPathways(hypoxiaPathways, GBMF_GBMM)
genesOfInterest <- unique(c(mrnaGenes, arGenes, immuneGenes, carbGenes, ecmGenes, cancerGenes, hypoxiaGenes))

# Run BLOBFISH on the PANDAs
null <- readRDS(paste0(nullDir, "nullPANDASubset.RDS"))

lggFemaleBlobfish <- netZooR::RunBLOBFISH(networks = list(lggFemaleTarget), 
                                          geneSet = genesOfInterest, 
                                          hopConstraint = 2, 
                                          alpha = 0.05, 
                                          nullDistribution = null,
                                          pValueFile = paste0(sourceDir, "lggFemaleBlobfishPvals"))
write.csv(lggFemaleBlobfish, paste0(sourceDir, 'LGG_Female_Panda_BLOBFISH_Full.csv'))

gbmFemaleBlobfish <- netZooR::RunBLOBFISH(networks = list(gbmFemaleTarget), 
                                          geneSet = genesOfInterest, 
                                          hopConstraint = 2, 
                                          alpha = 0.05, 
                                          nullDistribution = null,
                                          pValueFile = paste0(sourceDir, "gbmFemaleBlobfishPvals"))
write.csv(gbmFemaleBlobfish, paste0(sourceDir, 'GBM_Female_Panda_BLOBFISH_Full.csv'))

lggMaleBlobfish <- netZooR::RunBLOBFISH(networks = list(lggMaleTarget), 
                                        geneSet = genesOfInterest, 
                                        hopConstraint = 2, 
                                        alpha = 0.05, 
                                        nullDistribution = null,
                                        pValueFile = paste0(sourceDir, "lggMaleBlobfishPvals"))
write.csv(lggMaleBlobfish, paste0(sourceDir, 'LGG_Male_Panda_BLOBFISH_Full.csv'))

gbmMaleBlobfish <- netZooR::RunBLOBFISH(networks = list(gbmMaleTarget), 
                                        geneSet = genesOfInterest, 
                                        hopConstraint = 2, 
                                        alpha = 0.05, 
                                        nullDistribution = null,
                                        pValueFile = paste0(sourceDir, "gbmMaleBlobfishPvals"))
write.csv(gbmMaleBlobfish, paste0(sourceDir, 'GBM_Male_Panda_BLOBFISH_Full.csv'))

# Find the GBM Female and Male specific edges
gbmFemaleSpecificBlobfish <- gbmFemaleBlobfish[setdiff(rownames(gbmFemaleBlobfish),
                                                       c(rownames(lggFemaleBlobfish),
                                                         rownames(gbmMaleBlobfish),
                                                         rownames(lggMaleBlobfish))),]

gbmMaleSpecificBlobfish <- gbmMaleBlobfish[setdiff(rownames(gbmMaleBlobfish),
                                                   c(rownames(lggFemaleBlobfish),
                                                     rownames(gbmFemaleBlobfish),
                                                     rownames(lggMaleBlobfish))),]

# Define pathway colors
mrnaColor <- rgb(red = 1, green = 0, blue = 0, alpha = 0.5)
arColor <- rgb(red = 0, green = 0, blue = 1, alpha = 0.5)
immuneColor <- rgb(red = 0, green = 1, blue = 0, alpha = 0.5)
carbColor <- rgb(red = 1, green = 1, blue = 0, alpha = 0.5)
ecmColor <- rgb(red = 1, green = 0, blue = 1, alpha = 0.5)
cancerColor <- rgb(red = 0, green = 1, blue = 1, alpha = 0.5)
hypoxiaColor <- rgb(red = 0, green = 0, blue = 0, alpha = 0.5)

# Assign pathways to Female edges
gbmFemaleToPathway <- gbmFemaleSpecificBlobfish
gbmFemaleToPathway$pathway <- "placeholder"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% mrnaGenes), "pathway"] <- "Splicing"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% arGenes), "pathway"] <- "Androgen Receptor"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% immuneGenes), "pathway"] <- "Immune"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% carbGenes), "pathway"] <- "Carbohydrate Metabolism"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% ecmGenes), "pathway"] <- "Extracellular Matrix"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% cancerGenes), "pathway"] <- "Targets of HIF1A"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% hypoxiaGenes), "pathway"] <- "Hypoxia"

# Assign pathways to Male edges
gbmMaleToPathway <- gbmMaleSpecificBlobfish
gbmMaleToPathway$pathway <- "placeholder"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% mrnaGenes), "pathway"] <- "Splicing"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% arGenes), "pathway"] <- "Androgen Receptor"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% immuneGenes), "pathway"] <- "Immune"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% carbGenes), "pathway"] <- "Carbohydrate Metabolism"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% ecmGenes), "pathway"] <- "Extracellular Matrix"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% cancerGenes), "pathway"] <- "Targets of HIF1A"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% hypoxiaGenes), "pathway"] <- "Hypoxia"

# Plot female network
uniqueGenesFemale <- unique(gbmFemaleToPathway$pathway)
geneColorMappingFemale <- data.frame(gene = uniqueGenesFemale, color = rep("gray", length(uniqueGenesFemale)))
geneColorMappingFemale[which(uniqueGenesFemale == "Splicing"), "color"] <- mrnaColor
geneColorMappingFemale[which(uniqueGenesFemale == "Androgen Receptor"), "color"] <- arColor
geneColorMappingFemale[which(uniqueGenesFemale == "Immune"), "color"] <- immuneColor
geneColorMappingFemale[which(uniqueGenesFemale == "Carbohydrate Metabolism"), "color"] <- carbColor
geneColorMappingFemale[which(uniqueGenesFemale == "Extracellular Matrix"), "color"] <- ecmColor
geneColorMappingFemale[which(uniqueGenesFemale == "Targets of HIF1A"), "color"] <- cancerColor
geneColorMappingFemale[which(uniqueGenesFemale == "Hypoxia"), "color"] <- hypoxiaColor

gbmFemaleEdgeNames <- paste(gbmFemaleToPathway$tf, gbmFemaleToPathway$pathway, sep = "_")
gbmFemaleSimplified <- do.call(rbind, lapply(unique(gbmFemaleEdgeNames), function(edge){
  firstInstance <- which(gbmFemaleEdgeNames == edge)[1]
  return(data.frame(tf = gbmFemaleToPathway[firstInstance, "tf"], 
                   pathway = gbmFemaleToPathway[firstInstance, "pathway"]))
}))
colnames(gbmFemaleSimplified)[2] <- "gene"
write.csv(gbmFemaleSimplified, paste0(sourceDir, "gbmFemaleSimplifiedBLOBFISH.csv"))

PlotNetwork(gbmFemaleSimplified, 
            geneColorMapping = geneColorMappingFemale,
            layoutBipartite = TRUE, 
            nodeSize = 6, 
            tfColor = "gray", 
            vertexLabels = gbmFemaleSimplified[,2])

# Plot male network
uniqueGenesMale <- unique(gbmMaleToPathway$pathway)
geneColorMappingMale <- data.frame(gene = uniqueGenesMale, color = rep("gray", length(uniqueGenesMale)))
geneColorMappingMale[which(uniqueGenesMale == "Splicing"), "color"] <- mrnaColor
geneColorMappingMale[which(uniqueGenesMale == "Androgen Receptor"), "color"] <- arColor
geneColorMappingMale[which(uniqueGenesMale == "Immune"), "color"] <- immuneColor
geneColorMappingMale[which(uniqueGenesMale == "Carbohydrate Metabolism"), "color"] <- carbColor
geneColorMappingMale[which(uniqueGenesMale == "Extracellular Matrix"), "color"] <- ecmColor
geneColorMappingMale[which(uniqueGenesMale == "Targets of HIF1A"), "color"] <- cancerColor
geneColorMappingMale[which(uniqueGenesMale == "Hypoxia"), "color"] <- hypoxiaColor

gbmMaleEdgeNames <- paste(gbmMaleToPathway$tf, gbmMaleToPathway$pathway, sep = "_")
gbmMaleSimplified <- do.call(rbind, lapply(unique(gbmMaleEdgeNames), function(edge){
  firstInstance <- which(gbmMaleEdgeNames == edge)[1]
  return(data.frame(tf = gbmMaleToPathway[firstInstance, "tf"], 
                   pathway = gbmMaleToPathway[firstInstance, "pathway"]))
}))
colnames(gbmMaleSimplified)[2] <- "gene"
write.csv(gbmMaleSimplified, paste0(sourceDir, "gbmMaleSimplifiedBLOBFISH.csv"))

PlotNetwork(gbmMaleSimplified, 
            geneColorMapping = geneColorMappingMale,
            layoutBipartite = FALSE, 
            nodeSize = 6, 
            tfColor = "gray",
            vertexLabels = gbmMaleSimplified[,2])

# Obtain distributions
femaleDistrib <- table(gbmFemaleSimplified$gene) / length(unique(gbmFemaleBlobfish$tf))
maleDistrib <- table(gbmMaleSimplified$gene) / length(unique(gbmMaleBlobfish$tf))
distribDF <- data.frame(
  pathwayCategory = rep(names(femaleDistrib), 2),
  Sex = rep(c("female", "male"), each = length(femaleDistrib)),
  percentOfTFs = c(femaleDistrib, maleDistrib)
)

ggplot(distribDF, aes(x = pathwayCategory, y = percentOfTFs, fill = Sex)) +
  geom_col(position = "dodge") +
  labs(x = "Pathway Category", 
       y = "Percent of Significant TFs Targeting Pathway Category", 
       title = "Sex-Specific Pathway Targeting") +
  theme_minimal() + 
  coord_flip() + 
  scale_fill_manual(
    values = c("female" = rgb(red = 252 / 255, green = 182 / 255, blue = 195 / 255),
               "male" = rgb(red = 189 / 255, green = 190 / 255, blue = 255 / 255))
  )

# Make UpSet plots
gbmFemaleSimplifiedGraph <- igraph::graph_from_data_frame(gbmFemaleSimplified)
gbmFemaleSimplifiedAdj <- igraph::as_adjacency_matrix(gbmFemaleSimplifiedGraph, 
                                                      sparse = FALSE,
                                                      type = "upper")
gbmFemaleSimplifiedAdjSub <- as.data.frame(gbmFemaleSimplifiedAdj[unique(gbmFemaleSimplified$tf),
                                                                  unique(gbmFemaleSimplified$gene)])

gbmMaleSimplifiedGraph <- igraph::graph_from_data_frame(gbmMaleSimplified)
gbmMaleSimplifiedAdj <- igraph::as_adjacency_matrix(gbmMaleSimplifiedGraph, 
                                                    sparse = FALSE,
                                                    type = "upper")
gbmMaleSimplifiedAdjSub <- as.data.frame(gbmMaleSimplifiedAdj[unique(gbmMaleSimplified$tf),
                                                              unique(gbmMaleSimplified$gene)])

# Set the UpSet plot intersections
setIntersectionsBinary <- c("1111110", "1111111", "1110110", "1110111", "1101111",
                            "1100110", "1011111", "1000111", "0110000", "0111110", 
                            "0111111", "0111101", "0110100", "0110110", "0110111",
                            "0110010", "0101110", "0101111", "0100110", "0100111",
                            "0100101", "0100010", "0011100", "0011110", "0011111", 
                            "0011101", "0010100", "0010110", "0010111", "0010101",
                            "0001100", "0001110", "0001111", "0001001", "0000100",
                            "0000110", "0000111", "0000101", "0000010")

# Set up the matrices
gbmFemaleSimplifiedAdjSubMat <- as.matrix(gbmFemaleSimplifiedAdjSub)
gbmMaleSimplifiedAdjSubMat <- as.matrix(gbmMaleSimplifiedAdjSub)
mode(gbmFemaleSimplifiedAdjSubMat) <- "logical"
mode(gbmMaleSimplifiedAdjSubMat) <- "logical"

gbmFemaleSimplifiedAdjComb <- make_comb_mat(gbmFemaleSimplifiedAdjSubMat, mode = "distinct")
gbmMaleSimplifiedAdjComb <- make_comb_mat(gbmMaleSimplifiedAdjSubMat, mode = "distinct")
gbmFemaleSimplifiedAdjComb <- gbmFemaleSimplifiedAdjComb[sort(rownames(gbmFemaleSimplifiedAdjComb)),]
gbmMaleSimplifiedAdjComb <- gbmMaleSimplifiedAdjComb[sort(rownames(gbmMaleSimplifiedAdjComb)),]
gbmFemaleSimplifiedAdjComb <- gbmFemaleSimplifiedAdjComb[,setIntersectionsBinary]
gbmMaleSimplifiedAdjComb <- gbmMaleSimplifiedAdjComb[,setIntersectionsBinary]

# Print differences
sexDiffs <- attr(gbmFemaleSimplifiedAdjComb, "comb_size") - attr(gbmMaleSimplifiedAdjComb, "comb_size")
print(attr(gbmFemaleSimplifiedAdjComb, "dimnames")[[1]][which(sexDiffs > 100)])

# Define colors for UpSet plots
ylim_range <- c(0, 150)
maleColor <- rgb(red = 189 / 255, green = 190 / 255, blue = 255 / 255)
femaleColor <- rgb(red = 252 / 255, green = 182 / 255, blue = 195 / 255)
maleColorSat <- rgb(red = 130 / 255, green = 141 / 255, blue = 255 / 255)
femaleColorSat <- rgb(red = 255 / 255, green = 115 / 255, blue = 147 / 255)

barsToHighlightFemales <- c("1111111", "0111111", "0110111", "0010111")
barsToHighlightMales <- c("0111110", "0110110", "0011110", "0100110", "0010110", "0000110", "0000010")

combColors <- rep("black", length(setIntersectionsBinary))
combColors[which(setIntersectionsBinary %in% barsToHighlightFemales)] <- femaleColorSat
combColors[which(setIntersectionsBinary %in% barsToHighlightMales)] <- maleColorSat

# Create annotations for Male
taMale <- HeatmapAnnotation(
  "Co-Regulator Count" = anno_barplot(
    comb_size(gbmMaleSimplifiedAdjComb),
    ylim = ylim_range,
    gp = gpar(fill = maleColor, col = NA),
    border = FALSE
  ),
  annotation_name_side = "left",
  annotation_name_rot = 0,
  annotation_height = unit(4, "cm"),
  annotation_name_gp = gpar(fontface = "bold")
)

raMale <- HeatmapAnnotation(
  "Regulator Count" = anno_barplot(
    set_size(gbmMaleSimplifiedAdjComb),
    border = FALSE,
    gp = gpar(fill = maleColor, col = NA)
  ),
  which = "row",
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontface = "bold"),
  annotation_width = unit(3, "cm")
)

# Create annotations for Female
taFemale <- HeatmapAnnotation(
  "Co-Regulator Count" = anno_barplot(
    comb_size(gbmFemaleSimplifiedAdjComb),
    ylim = ylim_range,
    gp = gpar(fill = femaleColor, col = NA),
    border = FALSE
  ),
  annotation_name_side = "left",
  annotation_name_rot = 0,
  annotation_height = unit(4, "cm"),
  annotation_name_gp = gpar(fontface = "bold")
)

raFemale <- HeatmapAnnotation(
  "Regulator Count" = anno_barplot(
    set_size(gbmFemaleSimplifiedAdjComb),
    border = FALSE,
    gp = gpar(fill = femaleColor, col = NA)
  ),
  which = "row",
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontface = "bold"),
  annotation_width = unit(3, "cm")
)

# Create UpSet plots
grid.newpage()
UpSet(gbmFemaleSimplifiedAdjComb, 
      top_annotation = taFemale, 
      right_annotation = raFemale,
      set_order = order(rownames(gbmFemaleSimplifiedAdjComb)),
      comb_col = combColors)
gridFemale <- grid.grab()

grid.newpage()
UpSet(gbmMaleSimplifiedAdjComb, 
      top_annotation = taMale, 
      right_annotation = raMale,
      set_order = order(rownames(gbmFemaleSimplifiedAdjComb)),
      comb_col = combColors)
gridMale <- grid.grab()

# Display both plots
grid.arrange(grobs = list(gridFemale, gridMale), nrow = 2)

cat("\n=== ANALYSIS COMPLETE ===\n")