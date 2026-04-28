library(netZooR)

# Read pathway result files.
sourceDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Tomi/"
GBMM_GBMF <- read.csv(paste0(sourceDir, "GBMMale_Female_PathwayFile.csv"))
GBMF_GBMM <- read.csv(paste0(sourceDir, "GBMFemale_Male_PathwayFile.csv"))

# Get all genes in the pathways of interest.
gmtFile = fgsea::gmtPathways(paste0(sourceDir, "c2.cp.v2023.2.Hs.symbols.gmt"))
mrnaPathways <- c("REACTOME_METABOLISM_OF_RNA", "REACTOME_Splicing_SPLICING", "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_Splicing")
immunePathways <- c("REACTOME_NEUTROPHIL_DEGRANULATION", "REACTOME_INNATE_IMMUNE_SYSTEM",
                    "KEGG_LYSOSOME")
carbPathways <- c("WP_METABOLIC_PATHWAYS_OF_FIBROBLASTS", "WP_AEROBIC_GLYCOLYSIS",
                  "WP_N_GLYCAN_BIOSYNTHESIS")
ecmPathways <- c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                 "REACTOME_COLLAGEN_FORMATION")
cancerPathways <- c("WP_TYPE_2_PAPILLARY_RENAL_CELL_CARCINOMA", "WP_CLEAR_CELL_RENAL_CELL_CARCINOMA_PATHWAYS")
hypoxiaPathways <- "PID_HIF1_TFPATHWAY"
getGenesInPathways <- function(pathwayNames, pathwayFile){
  geneStrings <- pathwayFile[which(pathwayFile$pathway %in% pathwayNames), "leadingGenesDrivingEnrichment"]
  geneLists <- lapply(geneStrings, function(string){return(strsplit(string, "; ")[[1]])})
  geneSet <- unique(unlist(geneLists))
  return(geneSet)
}
mrnaGenes <- getGenesInPathways(mrnaPathways, GBMM_GBMF)
immuneGenes <- getGenesInPathways(immunePathways, GBMF_GBMM)
carbGenes <- getGenesInPathways(carbPathways, GBMF_GBMM)
ecmGenes <- getGenesInPathways(ecmPathways, GBMF_GBMM)
cancerGenes <- getGenesInPathways(cancerPathways, GBMF_GBMM)
hypoxiaGenes <- getGenesInPathways(hypoxiaPathways, GBMF_GBMM)
genesOfInterest <- unique(c(mrnaGenes, arGenes, immuneGenes, carbGenes, ecmGenes, cancerGenes,
                     hypoxiaGenes))

#Read BLOBFISH results
lggFemaleBlobfish <- read.csv(paste0(sourceDir, "REMBRANDT/Blobfish/LGG_Female_Panda_BLOBFISH_Full.csv"),
                              row.names = 1)
gbmFemaleBlobfish <- read.csv(paste0(sourceDir, "REMBRANDT/Blobfish/GBM_Female_Panda_BLOBFISH_Full.csv"),
                              row.names = 1)
lggMaleBlobfish <- read.csv(paste0(sourceDir, "REMBRANDT/Blobfish/LGG_Male_Panda_BLOBFISH_Full.csv"),
                              row.names = 1)
gbmMaleBlobfish <- read.csv(paste0(sourceDir, "REMBRANDT/Blobfish/GBM_Male_Panda_BLOBFISH_Full.csv"),
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

# Plot these networks, color-coding by pathway (use geneColorMapping for this.)
# Also, include the TF labels.
mrnaColor <- rgb(red = 1, green = 0, blue = 0, alpha = 0.5)
immuneColor <- rgb(red = 0, green = 1, blue = 0, alpha = 0.5)
carbColor <- rgb(red = 1, green = 1, blue = 0, alpha = 0.5)
ecmColor <- rgb(red = 1, green = 0, blue = 1, alpha = 0.5)
cancerColor <- rgb(red = 0, green = 1, blue = 1, alpha = 0.5)
hypoxiaColor <- rgb(red = 0, green = 0, blue = 0, alpha = 0.5)

gbmFemaleToPathway <- gbmFemaleSpecificBlobfish
gbmFemaleToPathway$pathway <- "placeholder"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% mrnaGenes), "pathway"] <- "Splicing"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% immuneGenes), "pathway"] <- "Immune"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% carbGenes), "pathway"] <- "Carbohydrate Metabolism"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% ecmGenes), "pathway"] <- "Extracellular Matrix"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% cancerGenes), "pathway"] <- "Targets of HIF1A"
gbmFemaleToPathway[which(gbmFemaleToPathway$gene %in% hypoxiaGenes), "pathway"] <- "Hypoxia"

gbmMaleToPathway <- gbmMaleSpecificBlobfish
gbmMaleToPathway$pathway <- "placeholder"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% mrnaGenes), "pathway"] <- "Splicing"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% immuneGenes), "pathway"] <- "Immune"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% carbGenes), "pathway"] <- "Carbohydrate Metabolism"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% ecmGenes), "pathway"] <- "Extracellular Matrix"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% cancerGenes), "pathway"] <- "Targets of HIF1A"
gbmMaleToPathway[which(gbmMaleToPathway$gene %in% hypoxiaGenes), "pathway"] <- "Hypoxia"

# Plot female.
uniqueGenesFemale <- unique(gbmFemaleToPathway$pathway)
geneColorMappingFemale <- data.frame(gene = uniqueGenesFemale, color = rep("gray", length(uniqueGenesFemale)))
geneColorMappingFemale[which(uniqueGenesFemale == "Splicing"), "color"] <- mrnaColor
geneColorMappingFemale[which(uniqueGenesFemale == "Immune"), "color"] <- immuneColor
geneColorMappingFemale[which(uniqueGenesFemale == "Carbohydrate Metabolism"), "color"] <- carbColor
geneColorMappingFemale[which(uniqueGenesFemale == "Extracellular Matrix"), "color"] <- ecmColor
geneColorMappingFemale[which(uniqueGenesFemale == "Targets of HIF1A"), "color"] <- cancerColor
geneColorMappingFemale[which(uniqueGenesFemale == "Hypoxia"), "color"] <- hypoxiaColor

gbmFemaleSimplified <- read.csv(paste0(sourceDir, "REMBRANDT/Blobfish/gbmFemaleSimplifiedBLOBFISH.csv"),
                              row.names = 1)
gbmFemaleSimplified <- gbmFemaleSimplified[which(gbmFemaleSimplified$gene != "Androgen Receptor"),]
PlotNetwork(gbmFemaleSimplified, geneColorMapping = geneColorMappingFemale,
            layoutBipartite = TRUE, nodeSize = 6, tfColor = "gray", 
            vertexLabels = c(gbmFemaleSimplified[,2]))

# Plot male.
uniqueGenesMale <- unique(gbmMaleToPathway$pathway)
geneColorMappingMale <- data.frame(gene = uniqueGenesMale, color = rep("gray", length(uniqueGenesMale)))
geneColorMappingMale[which(uniqueGenesMale == "Splicing"), "color"] <- mrnaColor
geneColorMappingMale[which(uniqueGenesMale == "Immune"), "color"] <- immuneColor
geneColorMappingMale[which(uniqueGenesMale == "Carbohydrate Metabolism"), "color"] <- carbColor
geneColorMappingMale[which(uniqueGenesMale == "Extracellular Matrix"), "color"] <- ecmColor
geneColorMappingMale[which(uniqueGenesMale == "Targets of HIF1A"), "color"] <- cancerColor
geneColorMappingMale[which(uniqueGenesMale == "Hypoxia"), "color"] <- hypoxiaColor
gbmMaleSimplified <- read.csv(paste0(sourceDir, "REMBRANDT/Blobfish/gbmMaleSimplifiedBLOBFISH.csv"),
                                row.names = 1)
PlotNetwork(gbmMaleSimplified, geneColorMapping = geneColorMappingMale,
            layoutBipartite = TRUE, nodeSize = 6, tfColor = "gray",
            vertexLabels = c(gbmMaleSimplified[,2]))

# Obtain distributions.
femaleDistrib <- table(gbmFemaleSimplified$gene) / length(unique(gbmFemaleBlobfish$tf))
maleDistrib <- table(gbmMaleSimplified$gene) / length(unique(gbmMaleBlobfish$tf))
maleDistrib["Androgen Receptor"] <- 0
maleDistrib <- maleDistrib[names(femaleDistrib)]
distribDF <- data.frame(
  pathwayCategory  = rep(names(femaleDistrib), 2),
  Sex = rep(c("female", "male"), each = length(femaleDistrib)),
  percentOfTFs = c(femaleDistrib, maleDistrib)
)
ggplot(distribDF, aes(x = pathwayCategory, y = percentOfTFs, fill = Sex)) +
  geom_col(position = "dodge") +
  labs(x = "Pathway Category", y = "Percent of Significant TFs Targeting Pathway Category", title = "Sex-Specific Pathway Targeting") +
  theme_minimal() + 
  coord_flip() + 
  scale_fill_manual(
    values = c("female" = rgb(red = 252 / 255, green = 182 / 255, blue = 195 / 255, alpha = 1),
               "male" = rgb(red = 189 / 255, green = 190 / 255, blue = 255 / 255, alpha = 1))
  )

# Make an UpSet plot.
library(ComplexHeatmap)
gbmFemaleSimplifiedGraph <- igraph::graph_from_data_frame(gbmFemaleSimplified)
gbmFemaleSimplifiedAdj <- igraph::as_adjacency_matrix(gbmFemaleSimplifiedGraph, sparse = FALSE,
                                                      type = "upper")
gbmFemaleSimplifiedAdjSub <- as.data.frame(gbmFemaleSimplifiedAdj[unique(gbmFemaleSimplified$tf),
                                                                  unique(gbmFemaleSimplified$gene)])
gbmFemaleSimplifiedAdjPerc <- gbmFemaleSimplifiedAdjSub / nrow(gbmFemaleSimplifiedAdjSub)

gbmMaleSimplifiedGraph <- igraph::graph_from_data_frame(gbmMaleSimplified)
gbmMaleSimplifiedAdj <- igraph::as_adjacency_matrix(gbmMaleSimplifiedGraph, sparse = FALSE,
                                                      type = "upper")
gbmMaleSimplifiedAdjSub <- as.data.frame(gbmMaleSimplifiedAdj[unique(gbmMaleSimplified$tf),
                                                              unique(gbmMaleSimplified$gene)])


library(gridExtra)
library(grid)

# Set the UpSet plot intersections.
setIntersections <- c("CarbohydrateMetabolism&ExtracellularMatrix&Hypoxia&Immune&mrnaSplicing&RCCHIF1AActivity",
                      "CarbohydrateMetabolism&ExtracellularMatrix&Hypoxia&Immune&mrnaSplicing",
                      "CarbohydrateMetabolism&ExtracellularMatrix&Hypoxia&Immune&RCCHIF1AActivity",
                      "CarbohydrateMetabolism&ExtracellularMatrix&Immune",
                      "CarbohydrateMetabolism&ExtracellularMatrix&Immune&mrnaSplicing",
                      "CarbohydrateMetabolism&ExtracellularMatrix&Immune&mrnaSplicing&RCCHIF1AActivity",
                      "CarbohydrateMetabolism&ExtracellularMatrix&mrnaSplicing",
                      "CarbohydrateMetabolism&ExtracellularMatrix",
                      "CarbohydrateMetabolism&Hypoxia&mrnaSplicing&Immune&RCCHIF1AActivity",
                      "CarbohydrateMetabolism&Hypoxia&mrnaSplicing&Immune",
                      "CarbohydrateMetabolism&Immune&mrnaSplicing",
                      "CarbohydrateMetabolism&Immune&mrnaSplicing&RCCHIF1AActivity",
                      "CarbohydrateMetabolism&Immune&RCCHIF1AActivity",
                      "CarbohydrateMetabolism&mrnaSplicing",
                      "ExtracellularMatrix&Hypoxia&Immune&mrnaSplicing&RCCHIF1AActivity",
                      "ExtracellularMatrix&Hypoxia&Immune&RCCHIF1AActivity",
                      "ExtracellularMatrix&Immune",
                      "ExtracellularMatrix&Immune&Hypoxia",
                      "ExtracellularMatrix&Immune&Hypoxia&mrnaSplicing",
                      "ExtracellularMatrix&Immune&mrnaSplicing",
                      "ExtracellularMatrix&Immune&mrnaSplicing&RCCHIF1AActivity",
                      "ExtracellularMatrix&Immune&RCCHIF1AActivity",
                      "Hypoxia&Immune",
                      "Hypoxia&Immune&mrnaSplicing&RCCHIF1AActivity",
                      "Hypoxia&Immune&RCCHIF1AActivity",
                      "Hypoxia&Immune&mrnaSplicing",
                      "Immune",
                      "Immune&ExtracellularMatrix",
                      "Immune&Hypoxia",
                      "Immune&mrnaSplicing",
                      "Immune&mrnaSplicing&RCCHIF1AActivity",
                      "Immune&RCCHIF1AActivity",
                      "mrnaSplicing")
setIntersectionsBinary <- c("111111", "111110", 
                            "111101", "110100", "110110", "110111",
                            "110010", "110000", "101110", "101111", "100110", "100111",
                            "100101", "100010", "011100", "011110", "011111", 
                            "011101", "010100", "010110", "010111", "010101",
                            "001100", "001110", "001111", "001001", "000100",
                            "000110", "000111", "000101", "000010")


# Set up the matrices accordingly.
gbmFemaleSimplifiedAdjSubMat <- as.matrix(gbmFemaleSimplifiedAdjSub)
gbmMaleSimplifiedAdjSubMat <- as.matrix(gbmMaleSimplifiedAdjSub)
mode(gbmFemaleSimplifiedAdjSubMat) <- "logical"
mode(gbmMaleSimplifiedAdjSubMat) <- "logical"
gbmFemaleSimplifiedAdjComb <- make_comb_mat(gbmFemaleSimplifiedAdjSubMat, 
                                            mode = "distinct")
gbmMaleSimplifiedAdjComb <- make_comb_mat(gbmMaleSimplifiedAdjSubMat, 
                                          mode = "distinct")
gbmFemaleSimplifiedAdjComb <- gbmFemaleSimplifiedAdjComb[sort(rownames(gbmFemaleSimplifiedAdjComb)),]
gbmMaleSimplifiedAdjComb <- gbmMaleSimplifiedAdjComb[sort(rownames(gbmMaleSimplifiedAdjComb)),]
gbmFemaleSimplifiedAdjComb <- gbmFemaleSimplifiedAdjComb[,setIntersectionsBinary]
gbmMaleSimplifiedAdjComb <- gbmMaleSimplifiedAdjComb[,setIntersectionsBinary]

# Print differences.
femalePercentages <- attr(gbmFemaleSimplifiedAdjComb, "comb_size") / sum(attr(gbmFemaleSimplifiedAdjComb, "comb_size"))
malePercentages <- attr(gbmMaleSimplifiedAdjComb, "comb_size") / sum(attr(gbmMaleSimplifiedAdjComb, "comb_size"))
sexDiffs <- attr(gbmFemaleSimplifiedAdjComb, "comb_size") - attr(gbmMaleSimplifiedAdjComb, "comb_size")
print(attr(gbmFemaleSimplifiedAdjComb, "dimnames")[[1]][which(sexDiffs > 100)])

# Make plots.
ylim_range <- c(0, 150)
maleColor <- rgb(red = 189 / 255, green = 190 / 255,
                 blue = 255 / 255)
femaleColor <- rgb(red = 252 / 255, green = 182 / 255,
                   blue = 195 / 255)
maleColorSat <- rgb(red = 130 / 255, green = 141 / 255,
                 blue = 255 / 255)
femaleColorSat <- rgb(red = 255 / 255, green = 115 / 255,
                   blue = 147 / 255)
combColors <- rep("black", length(setIntersectionsBinary))

taBoth <- HeatmapAnnotation(
  "Co-Regulator Count" = anno_barplot(
    as.matrix(data.frame(male = comb_size(gbmMaleSimplifiedAdjComb),
                         female = comb_size(gbmFemaleSimplifiedAdjComb))),          # the bar heights
    ylim = ylim_range,      # fix the y-axis range
    gp = gpar(fill = c(maleColor, femaleColor), col = NA),  # bar color + remove border
    border = FALSE
  ),
  annotation_name_side = "left",   # put label on left
  annotation_name_rot = 0,
  annotation_height = unit(4, "cm"),
  annotation_name_gp   = gpar(fontface = "bold")
)
raBoth <- HeatmapAnnotation(
  "Regulator Count" = anno_barplot(as.matrix(data.frame(
    male = set_size(gbmMaleSimplifiedAdjComb),
    female = set_size(gbmFemaleSimplifiedAdjComb)
  )),                 # bar lengths = set sizes
    border = FALSE,
    gp = gpar(fill = c(maleColor, femaleColor), col = NA)
  ),
  which = "row",                   # <- row annotation
  annotation_name_side = "bottom", # put label under x-axis
  annotation_name_rot  = 0,        # horizontal
  annotation_name_gp   = gpar(fontface = "bold"),
  annotation_width     = unit(3, "cm")
)
grid.newpage()
UpSet(gbmFemaleSimplifiedAdjComb, top_annotation = taBoth, right_annotation = raBoth,
      set_order = order(rownames(gbmFemaleSimplifiedAdjComb)),
      comb_col = combColors)
gridFemale <- grid.grab()
grid.newpage()
UpSet(gbmMaleSimplifiedAdjComb, top_annotation = taBoth, right_annotation = raBoth,
      set_order = order(rownames(gbmFemaleSimplifiedAdjComb)),
      comb_col = combColors)
gridMale <- grid.grab()
grid.arrange(grobs = list(gridFemale, gridMale), nrow = 2)   # two rows

