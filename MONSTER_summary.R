library(netZooR)
localDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Tomi/CGGA/"
gbm_mf <- readRDS(paste0(localDir, "Percentile_GBM_Male_Female_MonsterFile.rds"))
lgg_mf <- readRDS(paste0(localDir, "Percentile_LGG_Male_Female_MonsterFile.rds"))
f_both <- readRDS(paste0(localDir, "Percentile_GBM_LGG_Female_MonsterFile.rds"))
m_both <- readRDS(paste0(localDir, "Percentile_GBM_LGG_Male_MonsterFile.rds"))

# Calculate stats.
constructStatDF <- function(monsterFile, outFile){
  stats <- monsterCalculateTmStats(monsterFile)$p.values
  statsFDR <- p.adjust(stats, method = "fdr")
  df <- data.frame(tf = names(stats), p = stats, padj = statsFDR)
  print(df[which(df$padj < 0.05),])
  write.csv(df, outFile, row.names = FALSE)
}
constructStatDF(gbm_mf, paste0(localDir, "GBM_MF_Monster_sig.csv"))
constructStatDF(lgg_mf, paste0(localDir, "LGG_MF_Monster_sig.csv"))
constructStatDF(m_both, paste0(localDir, "M_Monster_sig.csv"))
constructStatDF(f_both, paste0(localDir, "F_Monster_sig.csv"))