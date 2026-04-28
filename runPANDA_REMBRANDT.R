library(netZooR)

# Transpose the files.
fileDir <- "../REMBRANDT/"
write.table(read.csv(paste0(fileDir, "femaleGBM_REMBRANDT.csv"), row.names = 1),
            sep = "\t",
            paste0(fileDir, "femaleGBM_REMBRANDT.txt"))
write.table(read.csv(paste0(fileDir, "maleGBM_REMBRANDT.csv"), row.names = 1),
            sep = "\t",
            paste0(fileDir, "maleGBM_REMBRANDT.txt"))
write.table(read.csv(paste0(fileDir, "femaleLGG_REMBRANDT.csv"), row.names = 1),
            sep = "\t",
            paste0(fileDir, "femaleLGG_REMBRANDT.txt"))
write.table(read.csv(paste0(fileDir, "maleLGG_REMBRANDT.csv"), row.names = 1),
            sep = "\t",
            paste0(fileDir, "maleLGG_REMBRANDT.txt"))

# Call Python PANDA for each input.
femaleGBM <- netZooR::pandaPy(expr_file = paste0(fileDir, "femaleGBM_REMBRANDT.txt"), 
                              motif_file = paste0(fileDir, "Symbol_MotifPriorGencode_p5_female_PARonX_nh.txt"), 
                              ppi_file = paste0(fileDir, "ppi_Gencode_p5_nh.txt"),
                              modeProcess = "intersection", with_header = TRUE)
saveRDS(femaleGBM, paste0(fileDir, "femaleGBM_PANDA.RDS"))
femaleLGG <- netZooR::pandaPy(expr_file = paste0(fileDir, "femaleLGG_REMBRANDT.txt"), 
                              motif_file = paste0(fileDir, "Symbol_MotifPriorGencode_p5_female_PARonX_nh.txt"), 
                              ppi_file = paste0(fileDir, "ppi_Gencode_p5_nh.txt"),
                              modeProcess = "intersection", with_header = TRUE)
saveRDS(femaleLGG, paste0(fileDir, "femaleLGG_PANDA.RDS"))
maleGBM <- netZooR::pandaPy(expr_file = paste0(fileDir, "maleGBM_REMBRANDT.txt"), 
                              motif_file = paste0(fileDir, "Symbol_MotifPriorGencode_p5_nh.txt"), 
                              ppi_file = paste0(fileDir, "ppi_Gencode_p5_nh.txt"),
                              modeProcess = "intersection", with_header = TRUE)
saveRDS(maleGBM, paste0(fileDir, "maleGBM_PANDA.RDS"))
maleLGG <- netZooR::pandaPy(expr_file = paste0(fileDir, "maleLGG_REMBRANDT.txt"), 
                              motif_file = paste0(fileDir, "Symbol_MotifPriorGencode_p5_nh.txt"), 
                              ppi_file = paste0(fileDir, "ppi_Gencode_p5_nh.txt"),
                              modeProcess = "intersection", with_header = TRUE)
saveRDS(maleLGG, paste0(fileDir, "maleLGG_PANDA.RDS"))

