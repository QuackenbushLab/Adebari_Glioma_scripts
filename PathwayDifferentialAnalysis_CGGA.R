sourceDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Tomi/CGGA/Differential_Analysis_results/"
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
maleGBMSpecificPathwaysStringent <- intersect(setdiff(GBMM_GBMF_sig, c(GBMF_GBMM_sig, LGGM_LGGF_sig)),
                                     GBMM_LGGM_sig)
femaleGBMSpecificPathwaysStringent <- intersect(setdiff(GBMF_GBMM_sig, c(GBMM_GBMF_sig, LGGF_LGGM_sig)),
                                       GBMF_LGGF_sig)

# Filter more leniently.
maleGBMSpecificPathways <- setdiff(GBMM_GBMF_sig, c(GBMF_GBMM_sig, LGGM_LGGF_sig))
femaleGBMSpecificPathways <- setdiff(GBMF_GBMM_sig, c(GBMM_GBMF_sig, LGGF_LGGM_sig))

