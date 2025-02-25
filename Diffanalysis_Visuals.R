
#Read files in.
localDir <- NULL
LGG_Male_Female <- read.csv(paste0(localDir, 'LGGMale_Female_PathwayFile.csv'))
LGG_Female_Male <- read.csv(paste0(localDir, 'LGGFemale_Male_PathwayFile.csv'))
LGG_Female_GBM_Female <- read.csv(paste0(localDir, 'LGGFemale_GBMFemale_PathwayFile.csv'))
LGG_Male_GBM_Male <- read.csv(paste0(localDir, 'LGGMale_GBMMale_PathwayFile.csv'))
GBM_Male_Female <- read.csv(paste0(localDir, 'GBMMale_Female_PathwayFile.csv'))
GBM_Female_Male <-read.csv(paste0(localDir, 'GBMFemale_Male_PathwayFile.csv'))
GBM_Female_LGG_Female <- read.csv(paste0(localDir, 'GBMFemale_LGGFemale_PathwayFile.csv'))
GBM_Male_LGG_Male <- read.csv(paste0(localDir, 'GBMMale_LGGMale_PathwayFile.csv'))


#Visualize GBM_M vs LGG_M
pathwayFlat <- data.frame(padj = c(GBM_Male_LGG_Male[which(GBM_Male_LGG_Male$padj < 0.000005), "padj"],
                                           LGG_Male_GBM_Male[which(LGG_Male_GBM_Male$padj < 0.000005), "padj"]),
                                  pathwayNames = c(GBM_Male_LGG_Male[which(GBM_Male_LGG_Male$padj < 0.000005), "pathway"],
                                           LGG_Male_GBM_Male[which(LGG_Male_GBM_Male$padj < 0.000005), "pathway"]),
                                  direction = c(rep("up", length(which(GBM_Male_LGG_Male$padj < 0.000005))),
                                  rep("down", length(which(LGG_Male_GBM_Male$padj < 0.000005)))))
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

pathwayFlat <- data.frame(padj = c(GBM_Female_LGG_Female[which(GBM_Female_LGG_Female$padj < 0.000005), "padj"],
                                           LGG_Female_GBM_Female[which(LGG_Female_GBM_Female$padj < 0.000005), "padj"]),
                                  pathwayNames = c(GBM_Female_LGG_Female[which(GBM_Female_LGG_Female$padj < 0.000005), "pathway"],
                                           LGG_Female_GBM_Female[which(LGG_Female_GBM_Female$padj < 0.000005), "pathway"]),
                                  direction = c(rep("up", length(which(GBM_Female_LGG_Female$padj < 0.000005))),
                                  rep("down", length(which(LGG_Female_GBM_Female$padj < 0.000005)))))
pathwayFlat %>% ggplot() + geom_col(aes(x = padj,
                                              y = fct_reorder(pathwayNames, padj,
                                                              .desc = TRUE),
                                              fill = direction)) +
  scale_fill_manual(values=group.colors) + 
  theme(axis.title.y = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank())


#Visualize GBM_M vs GBM_F

pathwayFlat <- data.frame(padj = c(GBM_Male_Female[which(GBM_Male_Female$padj < 0.005), "padj"],
                                           GBM_Female_Male[which(GBM_Female_Male$padj < 0.005), "padj"]),
                                  pathwayNames = c(GBM_Male_Female[which(GBM_Male_Female$padj < 0.005), "pathway"],
                                           GBM_Female_Male[which(GBM_Female_Male$padj < 0.005), "pathway"]),
                                  direction = c(rep("up", length(which(GBM_Male_Female$padj < 0.005))),
                                  rep("down", length(which(GBM_Female_Male$padj < 0.005)))))
pathwayFlat %>% ggplot() + geom_col(aes(x = padj,
                                              y = fct_reorder(pathwayNames, padj,
                                                              .desc = TRUE),
                                              fill = direction)) +
  scale_fill_manual(values=group.colors) + 
  theme(axis.title.y = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank())


#Visualize LGG_M vs LGG_F

pathwayFlat <- data.frame(padj = c(LGG_Male_Female[which(LGG_Male_Female$padj < 0.00005), "padj"],
                                           LGG_Female_Male[which(LGG_Female_Male$padj < 0.00005), "padj"]),
                                  pathwayNames = c(LGG_Male_Female[which(LGG_Male_Female$padj < 0.00005), "pathway"],
                                           LGG_Female_Male[which(LGG_Female_Male$padj < 0.00005), "pathway"]),
                                  direction = c(rep("up", length(which(LGG_Male_Female$padj < 0.00005))),
                                  rep("down", length(which(LGG_Female_Male$padj < 0.00005)))))
pathwayFlat %>% ggplot() + geom_col(aes(x = padj,
                                              y = fct_reorder(pathwayNames, padj,
                                                              .desc = TRUE),
                                              fill = direction)) +
  scale_fill_manual(values=group.colors) + 
  theme(axis.title.y = element_blank(), panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), panel.border = element_blank())