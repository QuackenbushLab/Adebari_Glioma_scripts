library(netZooR)

# Read all MONSTER files.
localDir <- NULL
LGG_Male_Female <- readRDS(paste0(localDir,"Latest_Percentile_LGG_Male_Female_MonsterFile.rds"))
GBM_Male_Female <- readRDS(paste0(localDir,"Latest_Percentile_GBM_Male_Female_MonsterFile.rds"))
LGG_GBM_Female <- readRDS(paste0(localDir,"Latest_Percentile_GBM_LGG_Female_MonsterFile.rds"))
LGG_GBM_Male <- readRDS(paste0(localDir,"Latest_Percentile_GBM_LGG_Male_MonsterFile.rds"))

# Make plots.
tfsOfInterest <- c("ZNF418", "ZNF879", "ZNF235", "ZNF496", "ZNF225", "ARNTL", "ARNT2", "GABPA")
monsterdTFIPlotCustom <- function(monsterObj, rescale='none', plot.title=NA, highlight.tfs=tfsOfInterest,
                            nTFs=-1, filename, w, h){
  if(is.na(plot.title)){
    plot.title <- "Differential TF Involvement"
  }
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
  
  t.values <- vapply(seqssdom,function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
  
  # Process the data for ggplot2
  combined.mat <- cbind(null.ssodm.matrix, ssodm)
  colnames(combined.mat) <- c(rep('Null',num.iterations),"Observed")
  
  
  if (rescale == 'significance'){
    combined.mat <- t(apply(combined.mat,1,function(x){
      (x-mean(x[-(num.iterations+1)]))/sd(x[-(num.iterations+1)])
    }))
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-t.values)]
    x.axis.size    <- 10 # pmin(15,7-log(p.values[order(p.values)]))
  } else if (rescale == 'none'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])
    x.axis.size    <- pmin(15,7-log(p.values))
  } else if (rescale == 'magnitude'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-combined.mat[, dim(combined.mat)[2]])]
    x.axis.size    <- pmin(15,7-log(p.values))
  }
  if(nTFs==-1){
    nTFs = length(x.axis.order)
  }
  null.SSODM.melt <- melt(combined.mat)[,-1][,c(2,1)]
  null.SSODM.melt$TF<-rep(rownames(monsterObj@nullTM[[1]]),num.iterations+1)
  
  ## Plot the data
  ggplot(null.SSODM.melt, aes(x=TF, y=value))+
    geom_point(aes(color=factor(Var2), alpha = .5*as.numeric(factor(Var2))), size=2) +
    scale_color_manual(values = c("blue", "red")) +
    scale_alpha(guide = "none") +
    scale_x_discrete(limits = x.axis.order[seq_len(nTFs)] ) +
    ylim(-5, 25) +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text.x = element_text(colour = 1+x.axis.order%in%highlight.tfs, 
                                     angle = 90, hjust = 1, 
                                     size=x.axis.size,face="bold")) +
    ylab("dTFI") +
    ggtitle(plot.title)
  
  ggsave(filename, width = w, height = h)
  
}
monsterdTFIPlotCustom(LGG_Male_Female, nTFs = 25, rescale = "significance",
                      filename = paste0(localDir, "monsterLGGMaleFemale.pdf"), 
                      w = 5, h = 3)
monsterdTFIPlotCustom(GBM_Male_Female, nTFs = 25, rescale = "significance",
                      filename = paste0(localDir, "monsterGBMMaleFemale.pdf"), 
                      w = 5, h = 3)
monsterdTFIPlotCustom(LGG_GBM_Female, nTFs = 25, rescale = "significance",
                      filename = paste0(localDir, "monsterLGGGBMFemale.pdf"), 
                      w = 5, h = 3)
monsterdTFIPlotCustom(LGG_GBM_Male, nTFs = 25, rescale = "significance",
                      filename = paste0(localDir, "monsterLGGGBMMale.pdf"), 
                      w = 5, h = 3)

# Make plots with subset information.
monsterdTFIPlotCustom <- function(monsterObj, rescale='none', plot.title=NA, highlight.tfs=NA,
                                  nTFs=-1, filename, w, h){
  if(is.na(plot.title)){
    plot.title <- "Differential TF Involvement"
  }
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
  
  t.values <- vapply(seqssdom,function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
  
  # Process the data for ggplot2
  combined.mat <- cbind(null.ssodm.matrix, ssodm)
  colnames(combined.mat) <- c(rep('Null',num.iterations),"Observed")
  
  
  if (rescale == 'significance'){
    combined.mat <- t(apply(combined.mat,1,function(x){
      (x-mean(x[-(num.iterations+1)]))/sd(x[-(num.iterations+1)])
    }))
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-t.values)]
    x.axis.size    <- 10 # pmin(15,7-log(p.values[order(p.values)]))
  } else if (rescale == 'none'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])
    x.axis.size    <- pmin(15,7-log(p.values))
  } else if (rescale == 'magnitude'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-combined.mat[, dim(combined.mat)[2]])]
    x.axis.size    <- pmin(15,7-log(p.values))
  }
  if(nTFs==-1){
    nTFs = length(x.axis.order)
  }
  null.SSODM.melt <- melt(combined.mat)[,-1][,c(2,1)]
  null.SSODM.melt$TF<-rep(rownames(monsterObj@nullTM[[1]]),num.iterations+1)
  
  ## Plot the data
  ggplot(null.SSODM.melt, aes(x=TF, y=value))+
    geom_point(aes(color=factor(Var2), alpha = .5*as.numeric(factor(Var2))), size=2) +
    scale_color_manual(values = c("blue", "red")) +
    scale_alpha(guide = "none") +
    scale_x_discrete(limits = x.axis.order[seq_len(nTFs)] ) +
    ylim(-5, 25) +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text.x = element_text(colour = 1+x.axis.order%in%highlight.tfs, 
                                     angle = 90, hjust = 1, 
                                     size=x.axis.size,face="bold")) +
    ylab("dTFI") +
    ggtitle(plot.title)
  
  ggsave(filename, width = w, height = h)
  
}
monsterdTFIPlotCustom(LGG_Male_Female, nTFs = 100, rescale = "significance",
                      filename = paste0(localDir, "monsterLGGMaleFemale.pdf"), 
                      w = 15, h = 3)
monsterdTFIPlotCustom(GBM_Male_Female, nTFs = 100, rescale = "significance",
                      filename = paste0(localDir, "monsterGBMMaleFemale.pdf"), 
                      w = 15, h = 3)
monsterdTFIPlotCustom(LGG_GBM_Female, nTFs = 100, rescale = "significance",
                      filename = paste0(localDir, "monsterLGGGBMFemale.pdf"), 
                      w = 15, h = 3)
monsterdTFIPlotCustom(LGG_GBM_Male, nTFs = 100, rescale = "significance",
                      filename = paste0(localDir, "monsterLGGGBMMale.pdf"), 
                      w = 15, h = 3)

# Custom function based on netZooR code.
monsterCSV <- function(monsterObj, filename = NA){
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
  
  t.values <- vapply(seqssdom,function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
  
  # Process the data for ggplot2
  combined.mat <- cbind(null.ssodm.matrix, ssodm)
  colnames(combined.mat) <- c(rep('Null',num.iterations),"Observed")
  
  
  combined.mat <- t(apply(combined.mat,1,function(x){
    (x-mean(x[-(num.iterations+1)]))/sd(x[-(num.iterations+1)])
  }))
  x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-t.values)]
  x.axis.size    <- 10 # pmin(15,7-log(p.values[order(p.values)]))
  null.SSODM.melt <- melt(combined.mat)[,-1][,c(2,1)]
  null.SSODM.melt$TF<-rep(rownames(monsterObj@nullTM[[1]]),num.iterations+1)
  
  ## Save the data
  print(table(factor(null.SSODM.melt$Var2)))
  observed <- null.SSODM.melt[which(null.SSODM.melt$Var2 == "Observed"), c("TF", "value")]
  write.csv(observed, filename)
  
}
monsterCSV(LGG_Male_Female, filename = paste0(localDir,"Latest_Percentile_LGG_Male_Female_MonsterFile.csv"))
monsterCSV(GBM_Male_Female, paste0(localDir,"Latest_Percentile_GBM_Male_Female_MonsterFile.csv"))
monsterCSV(LGG_GBM_Female, paste0(localDir,"Latest_Percentile_LGG_GBM_Female_MonsterFile.csv"))
monsterCSV(LGG_GBM_Male, paste0(localDir,"Latest_Percentile_LGG_GBM_Male_MonsterFile.csv"))
