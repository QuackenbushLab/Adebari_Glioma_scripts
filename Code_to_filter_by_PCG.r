mappingFile <- NULL
matrixFile <- NULL
filteredMatrixFile <- NULL

mapping <- read.csv(mappingFile, row.names = 1)
normConcatMatrix <- read.csv(matrixFile, row.names = 1) 
mapping$gene_id_nodot <- unlist(lapply(mapping$gene_id, function(id){return(strsplit(id, ".", fixed = TRUE)[[1]][1])}))
sharedGenes <- intersect(mapping[which(mapping$gene_type == "protein_coding"), "gene_id_nodot"], rownames(normConcatMatrix)) 
proteinCodingMatrix <- normConcatMatrix[sharedGenes,]
write.csv(mapping, filteredMatrixFile)