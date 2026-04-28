library(AnnotationDbi)
library(hgu133plus2.db)

# Read in the data.
sourceDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Tomi/REMBRANDT/"
pccFile <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/gencode.v49.basic.annotation.gtf"
clinical <- read.table(paste0(sourceDir, "GSE108474_REMBRANDT_clinical.data.txt"), sep = "\t",
                       header = TRUE)
expression <- read.table(paste0(sourceDir, "GSE108474_REMBRANDT_GeneExpression.txt"), 
                         sep = "\t", header = TRUE, row.names = 1)

# Convert Affy ID to ENSEMBL ID.
mapped <- AnnotationDbi::mapIds(
  hgu133plus2.db,
  keys = rownames(expression),
  column = "ENSEMBL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Limit to only protein coding genes.
mapping <- read.table(pccFile, sep = "\t")
mapping$gene_type <- unlist(lapply(mapping$V9, function(id){return(strsplit(strsplit(id, "gene_type ", fixed = TRUE)[[1]][2], ";")[[1]][1])}))
proteinCoding <- mapping[which(mapping$gene_type == "protein_coding"),]
proteinCoding$gene_name <- unlist(lapply(proteinCoding$V9, function(id){return(strsplit(strsplit(id, "gene_name ", fixed = TRUE)[[1]][2], ";")[[1]][1])}))
proteinCoding$gene_id <- unlist(lapply(proteinCoding$V9, function(id){return(strsplit(strsplit(strsplit(id, "gene_id ", fixed = TRUE)[[1]][2], ";")[[1]][1],
                                                                             ".", fixed = TRUE)[[1]][1])}))

write.csv(proteinCoding, paste0(sourceDir, "proteinCodingGenes.csv"))
mappedProteinCoding <- which(mapped %in% proteinCoding$gene_id)
expressionProteinCoding <- expression[mappedProteinCoding,]
mappingForExpression <- data.frame(affy = rownames(expressionProteinCoding),
                                   ensembl = mapped[mappedProteinCoding])
write.csv(mappingForExpression, paste0(sourceDir, "mappingForExpression.csv"))

# Change expression data to include only ENSEMBL genes.
ensemblIds <- table(mappingForExpression[,"ensembl"])
uniqueIds <- names(ensemblIds)
expressionNoDup <- do.call(rbind, lapply(uniqueIds, function(id){
  
  # Check where the Affy probe maps to the current ID.
  whichId <- which(mappingForExpression$ensembl == id)
  affyWhichId <- mappingForExpression[whichId, "affy"]

  # If only one location, return it with the Ensembl ID.
  # If multiple, average together the values and assign the Ensembl ID.
  ensembl <- NULL
  if(length(whichId) == 1){
    ensembl <- as.data.frame(expressionProteinCoding[affyWhichId,])
    colnames(ensembl) <- colnames(expressionProteinCoding)
    rownames(ensembl) <- mappingForExpression[whichId, "ensembl"]
  }else{
    ensembl <- as.data.frame(t(colSums(expressionProteinCoding[affyWhichId,]) / length(whichId)))
    colnames(ensembl) <- colnames(expressionProteinCoding)
    rownames(ensembl) <- mappingForExpression[whichId, "ensembl"][1]
  }
  return(ensembl)
}))
write.csv(expressionNoDup, paste0(sourceDir, "expressionProteinCodingEnsembl.csv"))

# Filter clinical data to separate out female GBM, male GBM, female LGG, and male LGG.
femaleGBM <- clinical[intersect(which(clinical$DISEASE_TYPE == "GBM"),
                                   which(clinical$GENDER == "FEMALE")),]
maleGBM <- clinical[intersect(which(clinical$DISEASE_TYPE == "GBM"),
                                 which(clinical$GENDER == "MALE")),]
femaleLGG <- clinical[intersect(which(clinical$WHO_GRADE == "II"),
                                   which(clinical$GENDER == "FEMALE")),]
maleLGG <- clinical[intersect(which(clinical$WHO_GRADE == "II"),
                                 which(clinical$GENDER == "MALE")),]

# Filter expression data.
expressionProteinCodingID <- expressionNoDup
colnames(expressionProteinCodingID) <- unlist(lapply(colnames(expressionNoDup),
                                                              function(name){
                                                                return(strsplit(name, "_")[[1]][1])
                                                              }))
expressionClinicalIntersectFemaleGBM <- intersect(colnames(expressionProteinCodingID), femaleGBM$SUBJECT_ID)
expressionClinicalIntersectFemaleLGG <- intersect(colnames(expressionProteinCodingID), maleGBM$SUBJECT_ID)
expressionClinicalIntersectMaleGBM <- intersect(colnames(expressionProteinCodingID), femaleLGG$SUBJECT_ID)
expressionClinicalIntersectMaleLGG <- intersect(colnames(expressionProteinCodingID), maleLGG$SUBJECT_ID)

femaleGBMExpression <- expressionProteinCodingID[,expressionClinicalIntersectFemaleGBM]
maleGBMExpression <- expressionProteinCodingID[,expressionClinicalIntersectFemaleLGG]
femaleLGGExpression <- expressionProteinCodingID[,expressionClinicalIntersectMaleGBM]
maleLGGExpression <- expressionProteinCodingID[,expressionClinicalIntersectMaleLGG]
write.csv(femaleGBMExpression, paste0(sourceDir, "femaleGBM_REMBRANDT.csv"))
write.csv(maleGBMExpression, paste0(sourceDir, "maleGBM_REMBRANDT.csv"))
write.csv(femaleLGGExpression, paste0(sourceDir, "femaleLGG_REMBRANDT.csv"))
write.csv(maleLGGExpression, paste0(sourceDir, "maleLGG_REMBRANDT.csv"))

write.csv(femaleGBM[expressionClinicalIntersectFemaleGBM,], paste0(sourceDir, "femaleGBM_clinical.csv"))
write.csv(maleGBM[expressionClinicalIntersectMaleGBM,], paste0(sourceDir, "maleGBM_clinical.csv"))
write.csv(femaleLGG[expressionClinicalIntersectFemaleLGG,], paste0(sourceDir, "femaleLGG_clinical.csv"))
write.csv(maleLGG[expressionClinicalIntersectMaleLGG,], paste0(sourceDir, "maleLGG_clinical.csv"))
