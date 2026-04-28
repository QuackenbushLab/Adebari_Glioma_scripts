# Packages
library("org.Hs.eg.db")
library("AnnotationDbi")

# Read the motif prior.
motifDir <- NULL
femaleName <- "MotifPriorGencode_p5_female_PARonX.txt"
maleName <- "MotifPriorGencode_p5.txt"

modifyMotif <- function(motifName){
  
  # Read the file.
  motif <- read.table(paste0(motifDir, motifName))
  # Map the targets to gene symbols.
  motifTargets <- unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys=motif$V2, 
                                               column=c("SYMBOL"), keytype = "ENSEMBL"))
  
  # Construct and save the new motif.
  newMotif <- motif
  newMotif$V2 <- motifTargets
  newMotif <- newMotif[which(!is.na(newMotif$V2)),]
  write.table(newMotif, paste0(motifDir, "Symbol_", motifName), sep = "\t", col.names = FALSE,
              quote = FALSE, row.names = FALSE)
}
modifyMotif(femaleName)
modifyMotif(maleName)
