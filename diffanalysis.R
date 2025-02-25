
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
  return(sumByGene)
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
LGGMaleFemalediff <- lggMaleTargeting$percentile - lggFemaleTargeting$percentile
GBMFemaleMalediff <- gbmFemaleTargeting$percentile - gbmMaleTargeting$percentile
LGGFemaleMalediff <- lggFemaleTargeting$percentile - lggMaleTargeting$percentile
GBMMaleLGGMalediff <- gbmMaleTargeting$percentile - lggMaleTargeting$percentile
GBMFemaleLGGFemale <- gbmFemaleTargeting$percentile - lggFemaleTargeting$percentile
LGGMaleGBMMalediff <- lggMaleTargeting$percentile - gbmMaleTargeting$percentile
LGGFemaleGBMFemalediff <- lggFemaleTargeting$percentile - lggMaleTargeting$percentile

#Add SYMBOLS.
library(org.Hs.eg.db)
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=gbmMaleTargeting$X, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
symbols <- ens2symbol[which(!is.na(ens2symbol$SYMBOL)),]
GBMMaleFemalediffSymbol <- GBMMaleFemalediff[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]
LGGMaleFemalediffSymbol <- LGGMaleFemalediff[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]
GBMFemaleMalediffSymbol <- GBMFemaleMalediff[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]
LGGFemaleMalediffSymbol <- LGGFemaleMalediff[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]
GBMMaleLGGMalediffSymbol <- GBMMaleLGGMalediff[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]
GBMFemaleLGGFemaleSymbol <- GBMFemaleLGGFemale[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]
LGGMaleGBMMalediffSymbol <- LGGMaleGBMMalediff[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]
LGGFemaleGBMFemalediffSymbol <- LGGFemaleGBMFemalediff[which(gbmMaleTargeting$X %in% symbols$ENSEMBL)]


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


# Run the pathway analysis with a "greater" and "less" group.
pathways = fgsea::gmtPathways(pathwayDBFile)
RunPathwayAnalysis <- function(networkScores, pathwayFile){

  # Run pathway analysis.
  pathwayResult <- fgsea::fgsea(pathways = pathways, stats = networkScores,
                                scoreType = "pos")
  str(pathways)
  str(networkScores)
  str(pathwayResult)
  
  # Compile the "leading edge" vector into a list.
  leadingEdge <- unlist(lapply(1:length(pathwayResult$leadingEdge), function(i){
    return(paste(pathwayResult$leadingEdge[i][[1]], collapse = "; "))
  }))
  str(leadingEdge)
  pathwayResultDf <- data.frame(pathway = pathwayResult$pathway,
                                pval = pathwayResult$pval,
                                padj = pathwayResult$padj,
                                pvalErrorSD = pathwayResult$log2err,
                                enrichmentScore = pathwayResult$ES,
                                normalizedEnrichmentScore = pathwayResult$NES,
                                remainingGeneCount = pathwayResult$size,
                                leadingGenesDrivingEnrichment = leadingEdge)
  str(pathwayResultDf)
  write.csv(pathwayResultDf, pathwayFile)



  # Return result.
  return(pathwayResultDf)
}

# Run for each comparison.
GBMMaleFemalePathways <- RunPathwayAnalysis(GBMMaleFemalediffSymbol, GBMMaleFemalePathwayFile)
LGGMaleFemaledPathways <- RunPathwayAnalysis(LGGMaleFemalediffSymbol, LGGMaleFemalePathwayFile)
GBMFemaleMalePathways <- RunPathwayAnalysis(GBMFemaleMalediffSymbol, GBMFemaleMalePathwayFile)
LGGFemaleMalePathways <- RunPathwayAnalysis(LGGFemaleMalediffSymbol, LGGFemaleMalePathwayFile)
GBMMaleLGGMalePathways <- RunPathwayAnalysis(GBMMaleLGGMalediffSymbol, GBMMaleLGGMalePathwayFile)
GBMFemaleLGGFemalePathways <- RunPathwayAnalysis(GBMFemaleLGGFemaleSymbol, GBMFemaleLGGFemalePathwayFile)
LGGMaleGBMMalePathways <- RunPathwayAnalysis(LGGMaleGBMMalediffSymbol, LGGMaleGBMMalePathwayFile)
LGGFemaleGBMFemalePathways <- RunPathwayAnalysis(LGGFemaleGBMFemalediffSymbol, LGGFemaleGBMFemalePathwayFile)