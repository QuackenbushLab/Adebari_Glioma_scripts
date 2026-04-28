# Read pathway result files.
pathwayDir <- NULL
tcgaDir <- NULL
rembrandtDir <- NULL
gmtFilePath <- NULL
GBMM_GBMF <- read.csv(paste0(pathwayDir, "GBMMale_Female_PathwayFile.csv"))
GBMF_GBMM <- read.csv(paste0(pathwayDir, "GBMFemale_Male_PathwayFile.csv"))

# Get all genes in the pathways of interest.
gmtFile = fgsea::gmtPathways(gmtFilePath)
cancerPathways <- c("WP_TYPE_2_PAPILLARY_RENAL_CELL_CARCINOMA", "WP_CLEAR_CELL_RENAL_CELL_CARCINOMA_PATHWAYS")
getGenesInPathways <- function(pathwayNames, pathwayFile){
  geneStrings <- pathwayFile[which(pathwayFile$pathway %in% pathwayNames), "leadingGenesDrivingEnrichment"]
  geneLists <- lapply(geneStrings, function(string){return(strsplit(string, "; ")[[1]])})
  geneSet <- unique(unlist(geneLists))
  return(geneSet)
}
cancerGenes <- getGenesInPathways(cancerPathways, GBMF_GBMM)

#Read BLOBFISH results
lggFemaleBlobfish <- read.csv(paste0(rembrandtDir, "LGG_Female_Panda_BLOBFISH_Full.csv"),
                              row.names = 1)
gbmFemaleBlobfish <- read.csv(paste0(rembrandtDir, "GBM_Female_Panda_BLOBFISH_Full.csv"),
                              row.names = 1)
lggMaleBlobfish <- read.csv(paste0(rembrandtDir, "LGG_Male_Panda_BLOBFISH_Full.csv"),
                            row.names = 1)
gbmMaleBlobfish <- read.csv(paste0(rembrandtDir, "GBM_Male_Panda_BLOBFISH_Full.csv"),
                            row.names = 1)

# Find the GBM Female and Male specific edges.
gbmFemaleSpecificBlobfish <- gbmFemaleBlobfish[setdiff(rownames(gbmFemaleBlobfish),
                                                       c(rownames(lggFemaleBlobfish),
                                                         rownames(gbmMaleBlobfish),
                                                         rownames(lggMaleBlobfish))),]
gbmMaleSpecificBlobfish <- gbmMaleBlobfish[setdiff(rownames(gbmMaleBlobfish),
                                                   c(rownames(lggFemaleBlobfish),
                                                     rownames(gbmFemaleBlobfish),
                                                     rownames(lggMaleBlobfish))),]
gbmFemaleHIF1A <- gbmFemaleSpecificBlobfish[which(gbmFemaleSpecificBlobfish$gene %in% cancerGenes),]
gbmMaleHIF1A <- gbmMaleSpecificBlobfish[which(gbmMaleSpecificBlobfish$gene %in% cancerGenes),]

# Do the same for TCGA.
lggFemaleBlobfishTCGA <- read.csv(paste0(tcgaDir, "LGG_Female_Panda_BLOBFISH_Full.csv"),
                              row.names = 1)
gbmFemaleBlobfishTCGA <- read.csv(paste0(tcgaDir, "GBM_Female_Panda_BLOBFISH_Full.csv"),
                              row.names = 1)
lggMaleBlobfishTCGA <- read.csv(paste0(tcgaDir, "LGG_Male_Panda_BLOBFISH_Full.csv"),
                            row.names = 1)
gbmMaleBlobfishTCGA <- read.csv(paste0(tcgaDir, "GBM_Male_Panda_BLOBFISH_Full.csv"),
                            row.names = 1)
gbmFemaleSpecificBlobfishTCGA <- gbmFemaleBlobfishTCGA[setdiff(rownames(gbmFemaleBlobfishTCGA),
                                                       c(rownames(lggFemaleBlobfishTCGA),
                                                         rownames(gbmMaleBlobfishTCGA),
                                                         rownames(lggMaleBlobfishTCGA))),]
gbmMaleSpecificBlobfishTCGA <- gbmMaleBlobfishTCGA[setdiff(rownames(gbmMaleBlobfishTCGA),
                                                   c(rownames(lggFemaleBlobfishTCGA),
                                                     rownames(gbmFemaleBlobfishTCGA),
                                                     rownames(lggMaleBlobfishTCGA))),]
gbmFemaleHIF1A_TCGA <- gbmFemaleSpecificBlobfishTCGA[which(gbmFemaleSpecificBlobfishTCGA$gene %in% cancerGenes),]
gbmMaleHIF1A_TCGA <- gbmMaleSpecificBlobfishTCGA[which(gbmMaleSpecificBlobfishTCGA$gene %in% cancerGenes),]

# Print the tables.
table(gbmMaleHIF1A$gene)
table(gbmFemaleHIF1A$gene)
table(gbmMaleHIF1A_TCGA$gene)
table(gbmFemaleHIF1A_TCGA$gene)