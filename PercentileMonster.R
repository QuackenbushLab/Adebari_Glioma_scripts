library(netZooR)
library(reshape2)

#Read in PANDA output paths
localDir
LGGFemalePath <- paste0(localDir,'/LGG_Female_Panda_output.txt')
GBMFemalePath <- paste0(localDir,'/GBM_Female_Panda_output.txt')
LGGMalePath <- paste0(localDir,'/LGG_Male_Panda_output.txt')
GBMMalePath <- paste0(localDir,'/GBM_Male_Panda_output.txt')

#Read as space delimited files
lggFemale <-read.table(LGGFemalePath, sep = " ")
lggMale <-read.table(LGGMalePath, sep = " ")
gbmFemale <-read.table(GBMFemalePath, sep = " ")
gbmMale <-read.table(GBMMalePath, sep = " ")


#Monster Results Path
GBM_LGG_Male_FilePath <- paste0(localDir,'/Percentile_GBM_LGG_Male_MonsterFile.rds')
GBM_LGG_Female_FilePath <- paste0(localDir,'/Percentile_GBM_LGG_Female_MonsterFile.rds')
GBM_Male_Female_FilePath <- paste0(localDir,'/Percentile_GBM_Male_Female_MonsterFile.rds')
LGG_Male_Female_FilePath <- paste0(localDir,'/Percentile_LGG_Male_Female_MonsterFile.rds')


#Percentile
lggFemale$percentile <- findInterval(lggFemale$V4, quantile(lggFemale$V4, seq(0,1, by=.001)))/1000
lggMale$percentile <- findInterval(lggMale$V4, quantile(lggMale$V4, seq(0,1, by=.001)))/1000
gbmFemale$percentile <- findInterval(gbmFemale$V4, quantile(gbmFemale$V4, seq(0,1, by=.001)))/1000
gbmMale$percentile <- findInterval(gbmMale$V4, quantile(gbmMale$V4, seq(0,1, by=.001)))/1000


#Monster Function              
domonster <- function(exp_graph, control_graph){
  combdf <- as.data.frame(cbind(control_graph, exp_graph))
  combdes <- c(rep(0, ncol(control_graph)),rep(1, ncol(exp_graph)))
  res <- netZooR::monster(expr = combdf,
                 design = combdes,
                 motif = NA,
                 nullPerms = 100,
                 numMaxCores = 4,
                 mode = 'regNet')
  return(res)
}

lggFemale <- lggFemale[,c(1:2,5)]
lggMale <- lggMale[,c(1:2,5)]
gbmFemale <- gbmFemale[,c(1:2,5)]
gbmMale <- gbmMale[,c(1:2,5)]


#convert to matrix
lggFemale_Matrix <- reshape2::dcast(lggFemale, V1~V2)
rownames(lggFemale_Matrix) <- lggFemale_Matrix[,"V1"]
lggFemale_Matrix <- lggFemale_Matrix[,2:ncol(lggFemale_Matrix)]


lggMale_Matrix <- reshape2::dcast(lggMale, V1~V2)
rownames(lggMale_Matrix) <- lggMale_Matrix[,"V1"]
lggMale_Matrix <- lggMale_Matrix[,2:ncol(lggMale_Matrix)]

gbmMale_Matrix <- reshape2::dcast(gbmMale, V1~V2)
rownames(gbmMale_Matrix) <- gbmMale_Matrix[,"V1"]
gbmMale_Matrix <- gbmMale_Matrix[,2:ncol(gbmMale_Matrix)]

gbmFemale_Matrix <- reshape2::dcast(gbmFemale, V1~V2)
rownames(gbmFemale_Matrix) <-gbmFemale_Matrix[,"V1"]
gbmFemale_Matrix <- gbmFemale_Matrix[,2:ncol(gbmFemale_Matrix)]

#Call Monster 
GBM_LGG_Male_Monster <- domonster(gbmMale_Matrix, lggMale_Matrix)
GBM_LGG_Female_Monster <- domonster(gbmFemale_Matrix, lggFemale_Matrix)
GBM_Male_Female_Monster <- domonster(gbmMale_Matrix, gbmFemale_Matrix)
LGG_Male_Female_Monster <- domonster(lggMale_Matrix, lggFemale_Matrix)


#Save files
saveRDS(GBM_LGG_Male_Monster, file = GBM_LGG_Male_FilePath)
saveRDS(GBM_LGG_Female_Monster, file = GBM_LGG_Female_FilePath)
saveRDS(GBM_Male_Female_Monster, file = GBM_Male_Female_FilePath)
saveRDS(LGG_Male_Female_Monster, file = LGG_Male_Female_FilePath)

print('The monsters are out!')

