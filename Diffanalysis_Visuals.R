library("forcats")

#Read files in.
localDir <- NULL
localDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Tomi/"
LGG_Male_Female <- read.csv(paste0(localDir, 'LGGMale_Female_PathwayFile.csv'))
LGG_Female_Male <- read.csv(paste0(localDir, 'LGGFemale_Male_PathwayFile.csv'))
LGG_Female_GBM_Female <- read.csv(paste0(localDir, 'LGGFemale_GBMFemale_PathwayFile.csv'))
LGG_Male_GBM_Male <- read.csv(paste0(localDir, 'LGGMale_GBMMale_PathwayFile.csv'))
GBM_Male_Female <- read.csv(paste0(localDir, 'GBMMale_Female_PathwayFile.csv'))
GBM_Female_Male <-read.csv(paste0(localDir, 'GBMFemale_Male_PathwayFile.csv'))
GBM_Female_LGG_Female <- read.csv(paste0(localDir, 'GBMFemale_LGGFemale_PathwayFile.csv'))
GBM_Male_LGG_Male <- read.csv(paste0(localDir, 'GBMMale_LGGMale_PathwayFile.csv'))

# Sort pathways by p-value for each.
GBM_Male_LGG_Male_sorted <- rbind(GBM_Male_LGG_Male, LGG_Male_GBM_Male)
GBM_Male_LGG_Male_sorted$direction <- c(rep("up", nrow(GBM_Male_LGG_Male)),
                                        rep("down", nrow(LGG_Male_GBM_Male)))
GBM_Male_LGG_Male_sorted <- GBM_Male_LGG_Male_sorted[order(GBM_Male_LGG_Male_sorted$padj),]

GBM_Female_LGG_Female_sorted <- rbind(GBM_Female_LGG_Female, LGG_Female_GBM_Female)
GBM_Female_LGG_Female_sorted$direction <- c(rep("up", nrow(GBM_Female_LGG_Female)),
                                        rep("down", nrow(LGG_Female_GBM_Female)))
GBM_Female_LGG_Female_sorted <- GBM_Female_LGG_Female_sorted[order(GBM_Female_LGG_Female_sorted$padj),]

GBM_Male_Female_sorted <- rbind(GBM_Male_Female, GBM_Female_Male)
GBM_Male_Female_sorted$direction <- c(rep("up", nrow(GBM_Male_Female)),
                                            rep("down", nrow(GBM_Female_Male)))
GBM_Male_Female_sorted <- GBM_Male_Female_sorted[order(GBM_Male_Female_sorted$padj),]

LGG_Male_Female_sorted <- rbind(LGG_Male_Female, LGG_Female_Male)
LGG_Male_Female_sorted$direction <- c(rep("up", nrow(LGG_Male_Female)),
                                      rep("down", nrow(LGG_Female_Male)))
LGG_Male_Female_sorted <- LGG_Male_Female_sorted[order(LGG_Male_Female_sorted$padj),]

#Visualize GBM_M vs LGG_M
pathwayFlat <- data.frame(padj = GBM_Male_LGG_Male_sorted[1:50, "padj"],
                          pathwayNames = GBM_Male_LGG_Male_sorted[1:50, "pathway"],
                          direction = GBM_Male_LGG_Male_sorted[1:50, "direction"])
group.colors <- c(up= "red", down = "blue")
pathwayFlat %>% ggplot() + geom_col(aes(x = padj, 
                                              y = fct_reorder(pathwayNames, padj,
                                                              .desc = TRUE),
                                              fill = direction)) +
  scale_fill_manual(values=group.colors) + 
  theme(axis.title.y = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank())


#Visualize GBM_F vs LGG_F
pathwayFlat <- data.frame(padj = GBM_Female_LGG_Female_sorted[1:50, "padj"],
                          pathwayNames = GBM_Female_LGG_Female_sorted[1:50, "pathway"],
                          direction = GBM_Female_LGG_Female_sorted[1:50, "direction"])
group.colors <- c(up= "red", down = "blue")
pathwayFlat %>% ggplot() + geom_col(aes(x = padj, 
                                        y = fct_reorder(pathwayNames, padj,
                                                        .desc = TRUE),
                                        fill = direction)) +
  scale_fill_manual(values=group.colors) + 
  theme(axis.title.y = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank())

#Visualize GBM_M vs GBM_F
pathwayFlat <- data.frame(padj = GBM_Male_Female_sorted[1:50, "padj"],
                          pathwayNames = GBM_Male_Female_sorted[1:50, "pathway"],
                          direction = GBM_Male_Female_sorted[1:50, "direction"])
group.colors <- c(up= "red", down = "blue")
pathwayFlat %>% ggplot() + geom_col(aes(x = padj, 
                                        y = fct_reorder(pathwayNames, padj,
                                                        .desc = TRUE),
                                        fill = direction)) +
  scale_fill_manual(values=group.colors) + 
  theme(axis.title.y = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank())


#Visualize LGG_M vs LGG_F
pathwayFlat <- data.frame(padj = LGG_Male_Female_sorted[1:50, "padj"],
                          pathwayNames = LGG_Male_Female_sorted[1:50, "pathway"],
                          direction = LGG_Male_Female_sorted[1:50, "direction"])
group.colors <- c(up= "red", down = "blue")
pathwayFlat %>% ggplot() + geom_col(aes(x = padj, 
                                        y = fct_reorder(pathwayNames, padj,
                                                        .desc = TRUE),
                                        fill = direction)) +
  scale_fill_manual(values=group.colors) + 
  theme(axis.title.y = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank())