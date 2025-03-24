
library("igraph") 
library("fgsea")

# Modify these paths.
localDir <-NULL
pathwayDBFile <- NULL

LGGFemaleFilePath <- paste0(localDir, 'LGG_Female_Panda_output.txt')
LGGFemaleFile <- read.table(LGGFemaleFilePath, sep="", header=FALSE)
LGGMaleFilePath <- paste0(localDir, 'LGG_Male_Panda_output.txt')
LGGMaleFile <- read.table(LGGMaleFilePath, sep="", header=FALSE)
GBMFemaleFilePath <- paste0(localDir, 'GBM_Female_Panda_output.txt')
GBMFemaleFile <- read.table(GBMFemaleFilePath, sep="", header=FALSE)
GBMMaleFilePath <- paste0(localDir, 'GBM_Male_Panda_output.txt')
GBMMaleFile <- read.table(GBMMaleFilePath, sep="", header=FALSE)
GBMMaleFemalePathwayFile <- paste0(localDir, 'GBMMale_Female_PathwayFile.csv')
LGGMaleFemalePathwayFile <- paste0(localDir, 'LGGMale_Female_PathwayFile.csv')
GBMFemaleMalePathwayFile <- paste0(localDir, 'GBMFemale_Male_PathwayFile.csv')
GBMMaleLGGMalePathwayFile <- paste0(localDir, 'GBMMale_LGGMale_PathwayFile.csv')
LGGFemaleMalePathwayFile <- paste0(localDir, 'LGGFemale_Male_PathwayFile.csv')
GBMFemaleLGGFemalePathwayFile <- paste0(localDir, 'GBMFemale_LGGFemale_PathwayFile.csv')
LGGMaleGBMMalePathwayFile <- paste0(localDir, 'LGGMale_GBMMale_PathwayFile.csv')
LGGFemaleGBMFemalePathwayFile <- paste0(localDir, 'LGGFemale_GBMFemale_PathwayFile.csv')

# Function to compute targeting of each gene.
ComputeTargeting <- function(network){

  # Compute targeting for each gene.
  colnames(network) <- c("source", "target", "motif", "weight")
  str(igraph::graph_from_data_frame(network, directed = TRUE))
  str(igraph::as_adjacency_matrix(igraph::graph_from_data_frame(network, directed = TRUE), attr = "weight"))
  str(as.matrix(igraph::as_adjacency_matrix(igraph::graph_from_data_frame(network, directed = TRUE), attr = "weight")))
  netAdj <- as.matrix(igraph::as_adjacency_matrix(igraph::graph_from_data_frame(network,
                                                                                directed = TRUE),
                                                  attr = "weight"))
  uniqueTF <- unique(network$source)
  uniqueGene <- unique(network$target)
  uniqueTFInRows <- uniqueTF[which(uniqueTF %in% rownames(netAdj))]
  uniqueGeneInCols <- uniqueGene[which(uniqueGene %in% colnames(netAdj))]
  netAdjTfGene <- netAdj[uniqueTFInRows, uniqueGeneInCols]
  sumByGene <- colSums(netAdjTfGene)
  names(sumByGene) <- uniqueGeneInCols
  sumByGeneDF <- data.frame(x = sumByGene)
  return(sumByGeneDF)
}

# Compute targeting of each group.
lggFemaleTargeting <- ComputeTargeting(LGGFemaleFile)
lggMaleTargeting <- ComputeTargeting(LGGMaleFile)
gbmFemaleTargeting <- ComputeTargeting(GBMFemaleFile)
gbmMaleTargeting <- ComputeTargeting(GBMMaleFile)

#Save compute targeting results 
write.csv(lggFemaleTargeting,paste0(localDir, 'lggFemaleTargetingFile.csv'))
write.csv(lggMaleTargeting, paste0(localDir, 'lggMaleTargetingFile.csv'))
write.csv(gbmFemaleTargeting, paste0(localDir, 'gbmFemaleTargetingFile.csv'))
write.csv(gbmMaleTargeting, paste0(localDir, 'gbmMaleTargetingFile.csv'))

#Ranking targeting scores
lggFemaleTargeting$percentile <- findInterval(lggFemaleTargeting$x, quantile(lggFemaleTargeting$x, seq(0,1, by=.001)))/1000
lggMaleTargeting$percentile <- findInterval(lggMaleTargeting$x, quantile(lggMaleTargeting$x, seq(0,1, by=.001)))/1000
gbmFemaleTargeting$percentile <- findInterval(gbmFemaleTargeting$x, quantile(gbmFemaleTargeting$x, seq(0,1, by=.001)))/1000
gbmMaleTargeting$percentile <- findInterval(gbmMaleTargeting$x, quantile(gbmMaleTargeting$x, seq(0,1, by=.001)))/1000

write.csv(lggFemaleTargeting, paste0(localDir, 'lggFemaleTargeting.csv'))
write.csv(lggMaleTargeting, paste0(localDir, 'lggMaleTargeting.csv'))
write.csv(gbmFemaleTargeting, paste0(localDir, 'gbmFemaleTargeting.csv'))
write.csv(gbmMaleTargeting, paste0(localDir, 'gbmMaleTargeting.csv'))

#taking differences in percentiles
GBMMaleFemalediff <- gbmMaleTargeting$percentile - gbmFemaleTargeting$percentile
names(GBMMaleFemalediff) <- rownames(gbmMaleTargeting)
LGGMaleFemalediff <- lggMaleTargeting$percentile - lggFemaleTargeting$percentile
names(LGGMaleFemalediff) <- rownames(gbmMaleTargeting)
GBMFemaleMalediff <- gbmFemaleTargeting$percentile - gbmMaleTargeting$percentile
names(GBMFemaleMalediff) <- rownames(gbmMaleTargeting)
LGGFemaleMalediff <- lggFemaleTargeting$percentile - lggMaleTargeting$percentile
names(LGGFemaleMalediff) <- rownames(gbmMaleTargeting)
GBMMaleLGGMalediff <- gbmMaleTargeting$percentile - lggMaleTargeting$percentile
names(GBMMaleLGGMalediff) <- rownames(gbmMaleTargeting)
GBMFemaleLGGFemale <- gbmFemaleTargeting$percentile - lggFemaleTargeting$percentile
names(GBMFemaleLGGFemale) <- rownames(gbmMaleTargeting)
LGGMaleGBMMalediff <- lggMaleTargeting$percentile - gbmMaleTargeting$percentile
names(LGGMaleGBMMalediff) <- rownames(gbmMaleTargeting)
LGGFemaleGBMFemalediff <- lggFemaleTargeting$percentile - lggMaleTargeting$percentile
names(LGGFemaleGBMFemalediff) <- rownames(gbmMaleTargeting)
# GBMMaleFemalediff <- gbmMaleTargeting$x - gbmFemaleTargeting$x
# names(GBMMaleFemalediff) <- rownames(gbmMaleTargeting)
# LGGMaleFemalediff <- lggMaleTargeting$x - lggFemaleTargeting$x
# names(LGGMaleFemalediff) <- rownames(gbmMaleTargeting)
# GBMFemaleMalediff <- gbmFemaleTargeting$x - gbmMaleTargeting$x
# names(GBMFemaleMalediff) <- rownames(gbmMaleTargeting)
# LGGFemaleMalediff <- lggFemaleTargeting$x - lggMaleTargeting$x
# names(LGGFemaleMalediff) <- rownames(gbmMaleTargeting)
# GBMMaleLGGMalediff <- gbmMaleTargeting$x - lggMaleTargeting$x
# names(GBMMaleLGGMalediff) <- rownames(gbmMaleTargeting)
# GBMFemaleLGGFemale <- gbmFemaleTargeting$x - lggFemaleTargeting$x
# names(GBMFemaleLGGFemale) <- rownames(gbmMaleTargeting)
# LGGMaleGBMMalediff <- lggMaleTargeting$x - gbmMaleTargeting$x
# names(LGGMaleGBMMalediff) <- rownames(gbmMaleTargeting)
# LGGFemaleGBMFemalediff <- lggFemaleTargeting$x - lggMaleTargeting$x
# names(LGGFemaleGBMFemalediff) <- rownames(gbmMaleTargeting)

#Add SYMBOLS.
library(org.Hs.eg.db)
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=rownames(gbmMaleTargeting), 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
symbols <- ens2symbol[which(!is.na(ens2symbol$SYMBOL)),]
whichInEnsembl <- which(rownames(gbmMaleTargeting) %in% symbols$ENSEMBL)
GBMMaleFemalediffSymbol <- GBMMaleFemalediff[whichInEnsembl]
LGGMaleFemalediffSymbol <- LGGMaleFemalediff[whichInEnsembl]
GBMFemaleMalediffSymbol <- GBMFemaleMalediff[whichInEnsembl]
LGGFemaleMalediffSymbol <- LGGFemaleMalediff[whichInEnsembl]
GBMMaleLGGMalediffSymbol <- GBMMaleLGGMalediff[whichInEnsembl]
GBMFemaleLGGFemaleSymbol <- GBMFemaleLGGFemale[whichInEnsembl]
LGGMaleGBMMalediffSymbol <- LGGMaleGBMMalediff[whichInEnsembl]
LGGFemaleGBMFemalediffSymbol <- LGGFemaleGBMFemalediff[whichInEnsembl]


mappedSymbols <- unlist(lapply(names(GBMMaleFemalediffSymbol), function(ensembl){
  return(symbols[which(symbols$ENSEMBL == ensembl)[[1]], "SYMBOL"])
}))
names(GBMMaleFemalediffSymbol) <- mappedSymbols
names(LGGMaleFemalediffSymbol) <- mappedSymbols
names(GBMFemaleMalediffSymbol) <- mappedSymbols
names(LGGFemaleMalediffSymbol) <- mappedSymbols
names(GBMMaleLGGMalediffSymbol) <- mappedSymbols
names(GBMFemaleLGGFemaleSymbol) <- mappedSymbols 
names(LGGMaleGBMMalediffSymbol) <- mappedSymbols
names(LGGFemaleGBMFemalediffSymbol) <- mappedSymbols

# Remove duplicates.
symbolCounts <- table(mappedSymbols)
dupSymbols <- names(symbolCounts)[which(symbolCounts > 1)]
removeDups <- function(diffData){
  
  # Find indices to remove.
  toRemove <- unlist(lapply(dupSymbols, function(symbol){
    whichSymbol <- which(mappedSymbols == symbol)
    absVals <- abs(diffData[whichSymbol])
    whichToKeep <- which.max(absVals)
    whichToRemove <- whichSymbol[setdiff(1:length(absVals), whichToKeep)]
    return(whichToRemove)
  }))
  
  # Update the data set.
  str(setdiff(1:length(diffData), toRemove))
  return(diffData[setdiff(1:length(diffData), toRemove)])
}
GBMMaleFemalediffSymbolDedup <- removeDups(GBMMaleFemalediffSymbol)
LGGMaleFemalediffSymbolDedup <- removeDups(LGGMaleFemalediffSymbol)
GBMFemaleMalediffSymbolDedup <- removeDups(GBMFemaleMalediffSymbol)
LGGFemaleMalediffSymbolDedup <- removeDups(LGGFemaleMalediffSymbol)
GBMMaleLGGMalediffSymbolDedup <- removeDups(GBMMaleLGGMalediffSymbol)
GBMFemaleLGGFemaleSymbolDedup <- removeDups(GBMFemaleLGGFemaleSymbol)
LGGMaleGBMMalediffSymbolDedup <- removeDups(LGGMaleGBMMalediffSymbol)
LGGFemaleGBMFemalediffSymbolDedup <- removeDups(LGGFemaleGBMFemalediffSymbol)

# Run the pathway analysis with a "greater" and "less" group.
pathways = fgsea::gmtPathways(pathwayDBFile)
RunPathwayAnalysis <- function(networkScores, pathwayFile){

  # Run pathway analysis.
  pathwayResult <- fgsea::fgsea(pathways = pathways, stats = networkScores,
                                scoreType = "pos")
  
  # Compile the "leading edge" vector into a list.
  leadingEdge <- unlist(lapply(1:length(pathwayResult$leadingEdge), function(i){
    return(paste(pathwayResult$leadingEdge[i][[1]], collapse = "; "))
  }))
  pathwayResultDf <- data.frame(pathway = pathwayResult$pathway,
                                pval = pathwayResult$pval,
                                padj = pathwayResult$padj,
                                pvalErrorSD = pathwayResult$log2err,
                                enrichmentScore = pathwayResult$ES,
                                normalizedEnrichmentScore = pathwayResult$NES,
                                remainingGeneCount = pathwayResult$size,
                                leadingGenesDrivingEnrichment = leadingEdge)
  print(GBMMaleFemalePathwayFile)
  write.csv(pathwayResultDf, pathwayFile)

  # Return result.
  return(pathwayResultDf)
}

# Run for each comparison.
GBMMaleFemalePathways <- RunPathwayAnalysis(GBMMaleFemalediffSymbolDedup, GBMMaleFemalePathwayFile)
LGGMaleFemaledPathways <- RunPathwayAnalysis(LGGMaleFemalediffSymbolDedup, LGGMaleFemalePathwayFile)
GBMFemaleMalePathways <- RunPathwayAnalysis(GBMFemaleMalediffSymbolDedup, GBMFemaleMalePathwayFile)
LGGFemaleMalePathways <- RunPathwayAnalysis(LGGFemaleMalediffSymbolDedup, LGGFemaleMalePathwayFile)
GBMMaleLGGMalePathways <- RunPathwayAnalysis(GBMMaleLGGMalediffSymbolDedup, GBMMaleLGGMalePathwayFile)
GBMFemaleLGGFemalePathways <- RunPathwayAnalysis(GBMFemaleLGGFemaleSymbolDedup, GBMFemaleLGGFemalePathwayFile)
LGGMaleGBMMalePathways <- RunPathwayAnalysis(LGGMaleGBMMalediffSymbolDedup, LGGMaleGBMMalePathwayFile)
LGGFemaleGBMFemalePathways <- RunPathwayAnalysis(LGGFemaleGBMFemalediffSymbolDedup, LGGFemaleGBMFemalePathwayFile)

# Paste together.
PastePathways <- function(up, down, pathwayFile){
  # Reformat for Cytoscape.
  pathwayResultDfCytoscapeUp <- data.frame(id = up$pathway,
                                         name = up$pathway,
                                         pval = up$pval,
                                         fdr = up$padj,
                                         phenotype = rep("+1", nrow(up)),
                                         geneList = up$leadingGenesDrivingEnrichment)
  pathwayResultDfCytoscapeDown <- data.frame(id = down$pathway,
                                           name = down$pathway,
                                           pval = down$pval,
                                           fdr = down$padj,
                                           phenotype = rep("-1", nrow(down)),
                                           geneList = down$leadingGenesDrivingEnrichment)
  pathwayResultDf <- rbind(pathwayResultDfCytoscapeUp, pathwayResultDfCytoscapeDown)
  write.table(pathwayResultDf, paste0(pathwayFile), quote = FALSE, row.names = FALSE,
              sep = "\t")
}
PastePathways(GBMMaleFemalePathways, GBMFemaleMalePathways,
              paste0(localDir, "GBMMaleFemale.txt"))
PastePathways(LGGMaleFemaledPathways, LGGFemaleMalePathways,
              paste0(localDir, "LGGMaleFemale.txt"))
PastePathways(GBMMaleLGGMalePathways, LGGMaleGBMMalePathways,
              paste0(localDir, "GBMLGGMale.txt"))
PastePathways(GBMFemaleLGGFemalePathways, LGGFemaleGBMFemalePathways,
              paste0(localDir, "GBMLGGFemale.txt"))