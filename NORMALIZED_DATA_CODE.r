concatenatedFile <- NULL
normalizedFile <- NULL

data<-read.table(concatenatedFile, header = TRUE, sep = "\t")
normalized <- voom(data)
write.csv(normalized, normalizedFile)