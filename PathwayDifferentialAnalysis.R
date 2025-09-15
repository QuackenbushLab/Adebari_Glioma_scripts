sourceDir <- NULL
GBMF_LGGF <- read.csv(paste0(sourceDir, "GBMFemale_LGGFemale_PathwayFile.csv"))
GBMM_LGGM <- read.csv(paste0(sourceDir, "GBMMale_LGGMale_PathwayFile.csv"))
GBMM_GBMF <- read.csv(paste0(sourceDir, "GBMMale_Female_PathwayFile.csv"))
GBMF_GBMM <- read.csv(paste0(sourceDir, "GBMFemale_Male_PathwayFile.csv"))
LGGM_LGGF <- read.csv(paste0(sourceDir, "LGGMale_Female_PathwayFile.csv"))
LGGF_LGGM <- read.csv(paste0(sourceDir, "LGGFemale_Male_PathwayFile.csv"))

# Get significant.
getSignificant <- function(pathways){
  return(pathways[which(pathways$padj < 0.05), "pathway"])
}
GBMF_LGGF_sig <- getSignificant(GBMF_LGGF)
GBMM_LGGM_sig <- getSignificant(GBMM_LGGM)
GBMM_GBMF_sig <- getSignificant(GBMM_GBMF)
GBMF_GBMM_sig <- getSignificant(GBMF_GBMM)
LGGM_LGGF_sig <- getSignificant(LGGM_LGGF)
LGGF_LGGM_sig <- getSignificant(LGGF_LGGM)

# Filter.
maleGBMSpecificPathways <- intersect(setdiff(GBMM_GBMF_sig, c(GBMF_GBMM_sig, LGGM_LGGF_sig)),
                                     GBMM_LGGM_sig)
femaleGBMSpecificPathways <- intersect(setdiff(GBMF_GBMM_sig, c(GBMM_GBMF_sig, LGGF_LGGM_sig)),
                                       GBMF_LGGF_sig)

# Make barplot of pathway significance.
GBMF_GBMM <- GBMF_GBMM[,1:8]
maleSigPathways <- GBMM_GBMF[which(GBMM_GBMF$pathway %in% maleGBMSpecificPathways),]
femaleSigPathways <- GBMF_GBMM[which(GBMF_GBMM$pathway %in% femaleGBMSpecificPathways),]
dfToInput <- rbind(maleSigPathways, femaleSigPathways)
dfToInput$sex <- c(rep("Male", length(maleGBMSpecificPathways)), rep("Female", length(femaleGBMSpecificPathways)))
dfToInput$logpval <- -1 * log10(dfToInput$padj)

order <- rev(c("REACTOME_METABOLISM_OF_RNA", 
           "REACTOME_MRNA_SPLICING", 
           "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA",
           "PID_AR_NONGENOMIC_PATHWAY", 
           "REACTOME_NEUTROPHIL_DEGRANULATION", 
           "REACTOME_INNATE_IMMUNE_SYSTEM",
           "KEGG_LYSOSOME", 
           "WP_METABOLIC_PATHWAYS_OF_FIBROBLASTS", 
           "WP_AEROBIC_GLYCOLYSIS",
           "WP_N_GLYCAN_BIOSYNTHESIS",
           "KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
           "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
           "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", 
           "REACTOME_COLLAGEN_DEGRADATION",
           "REACTOME_COLLAGEN_FORMATION",
           "PID_HIF1_TFPATHWAY",
           "WP_2Q37_COPY_NUMBER_VARIATION_SYNDROME",
           "WP_TYPE_2_PAPILLARY_RENAL_CELL_CARCINOMA",
           "WP_CLEAR_CELL_RENAL_CELL_CARCINOMA_PATHWAYS"))

alpha <- 0.05
thr   <- -log10(alpha)
maleColor <- rgb(red = 189 / 255, green = 190 / 255,
                 blue = 255 / 255)
femaleColor <- rgb(red = 252 / 255, green = 182 / 255,
                   blue = 195 / 255)
rownames(dfToInput) <- dfToInput$pathway
dfToInput <- dfToInput[order,]
dfToInput$pathway <- factor(dfToInput$pathway, levels = order)
ggplot(dfToInput, aes(pathway, logpval, fill = sex)) +
  geom_col(width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = thr, linetype = 2, linewidth = 0.4) +
  scale_fill_manual(values = c(Male = maleColor, Female = femaleColor)) +
  labs(x = NULL, y = expression(-log[10](p)), fill = NULL,
       title = "GBM-Specific \nEnriched Pathways") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(sourceDir, "pathwayBarplot.pdf"), plot = last_plot(),
       width = 10, height = 10, units = "in")

# Function to calculate dTFI without plotting it.
monsterdTFIEmpiricalPValues <- function(monsterObj, filename){
  num.iterations <- length(monsterObj@nullTM)
  # Calculate the off-diagonal squared mass for each transition matrix
  null.SSODM <- lapply(monsterObj@nullTM,function(x){
    apply(x,2,function(y){t(y)%*%y})
  })
  null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
  null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))
  
  ssodm <- apply(monsterObj@tm,2,function(x){t(x)%*%x})
  
  seqssdom <- seq_along(ssodm)
  names(seqssdom) <- names(ssodm)
  p.values <- 1-pnorm(vapply(seqssdom,function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE))
  
  # Set up as data frame.
  pValDF <- data.frame(pval = p.values)
  pValDF$padj <- p.adjust(pValDF$pval, method = "fdr")
  rownames(pValDF) <- names(p.values)
  write.csv(pValDF, filename)
  return(pValDF)
}

# Read MONSTER files.
GBMM_GBMF_monster <- readRDS(paste0(sourceDir, "Latest_Percentile_GBM_Male_Female_MonsterFile.rds"))
LGGM_LGGF_monster <- readRDS(paste0(sourceDir, "Latest_Percentile_LGG_Male_Female_MonsterFile.rds"))
GBMM_LGGM_monster <- readRDS(paste0(sourceDir, "Latest_Percentile_GBM_LGG_Male_MonsterFile.rds"))
GBMF_LGGF_monster <- readRDS(paste0(sourceDir, "Latest_Percentile_GBM_LGG_Female_MonsterFile.rds"))

# Get dTFI.
GBMM_GBMF_dTFI <- monsterdTFIEmpiricalPValues(GBMM_GBMF_monster,
                                              paste0(sourceDir, "Latest_Percentile_GBM_Male_Female_MonsterFile_pvals.csv"))
LGGM_LGGF_dTFI <- monsterdTFIEmpiricalPValues(LGGM_LGGF_monster,
                                              paste0(sourceDir, "Latest_Percentile_LGG_Male_Female_MonsterFile_pvals.csv"))
GBMM_LGGM_dTFI <- monsterdTFIEmpiricalPValues(GBMM_LGGM_monster,
                                              paste0(sourceDir, "Latest_Percentile_GBM_LGG_Male_MonsterFile_pvals.csv"))
GBMF_LGGF_dTFI <- monsterdTFIEmpiricalPValues(GBMF_LGGF_monster,
                                              paste0(sourceDir, "Latest_Percentile_GBM_LGG_Female_MonsterFile_pvals.csv"))

# Filter.
GBM_only <- setdiff(rownames(GBMM_GBMF_dTFI)[which(GBMM_GBMF_dTFI$padj < 0.05)],
                    rownames(LGGM_LGGF_dTFI)[which(LGGM_LGGF_dTFI$padj < 0.05)])
GBM_only_female <- rownames(GBMF_LGGF_dTFI)[which(GBMF_LGGF_dTFI$padj < 0.05)]
GBM_only_male <- rownames(GBMM_LGGM_dTFI)[which(GBMM_LGGM_dTFI$padj < 0.05)]
GBM_overlap_female_specific <- intersect(GBM_only, GBM_only_female)
GBM_overlap_male_specific <- intersect(GBM_only, GBM_only_male)
GBM_female_specific <- setdiff(GBM_only_female, GBM_only_male)
GBM_male_specific <- setdiff(GBM_only_male, GBM_only_female)
write.csv(GBM_female_specific, paste0(sourceDir, "gbmFemaleSpecificMONSTER.csv"))