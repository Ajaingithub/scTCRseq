TCR_AJ_EBV <- function(object, savedir, objname){
  library(taRifx)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(Seurat)
  
  dir.create(savedir,showWarnings = FALSE)
  setwd(savedir)
  
  Idents(object) <- object@meta.data$seurat_clusters
  ### In here we will do TRA CDR3 and TRB CDR3 separately
  TRB_CDR3 <- table(object@meta.data$TCR) %>% as.data.frame()
  TRB_clonal <- TRB_CDR3[order(-TRB_CDR3$Freq),]
  TRB_clonal_gt_2 <- TRB_clonal[TRB_clonal$Freq > 1,]
  TRB_clonal_gt_2$Var1 <- factor(TRB_clonal_gt_2$Var1, levels = TRB_clonal_gt_2$Var1)
  TRB_clonal_gt_2 <- TRB_clonal_gt_2[grep("^_$",TRB_clonal_gt_2$Var1,invert = TRUE),]
  
  p <- ggplot(TRB_clonal_gt_2, aes(x=Var1,y=Freq, fill = Var1)) + geom_bar(stat="identity", color="black") +
    ggtitle(paste(objname," TCR CDR3 aa gt 1",sep = "")) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.minor = element_line(colour = "grey"),
          panel.grid.major = element_line(colour = "grey")) + NoLegend()
  
  pdf(paste(savedir,"/",objname,"_TRB_clonal_expanded_gt_1.pdf",sep = ""), width = 20, height = 6)
  print(p)
  dev.off()
  
  rm(plot_list)
  plot_list <- list()
  for(i in 1:length(TRB_clonal_gt_2$Var1)){
    clonal_expand <- rownames(object@meta.data[grep(TRB_clonal_gt_2$Var1[i],object@meta.data$TCR),])
    p <- DimPlot(object,cells.highlight = clonal_expand, reduction = "wnn.umap", label = TRUE, label.size = 3, sizes.highlight = 0.5) + 
      ggtitle(TRB_clonal_gt_2$Var1[i]) + 
      NoLegend() +
      theme(plot.title = element_text(hjust = 0.5))
    # pdf(paste(savedir,"object/object_TRB_",TRB_clonal_gt_2$Var1[i],".pdf",sep = ""))
    # print(p)
    # dev.off()
    plot_list[[i]] <- p
  }
  
  pdf(paste(savedir,"/",objname,"_TCR_all_UMAP.pdf",sep = ""), width = 4.3, height = 4)
  print(plot_list)
  dev.off()
  
  ## TRB greater than 1
  clonal_expand <- rownames(object@meta.data[grep(paste(TRB_clonal_gt_2$Var1,collapse = "|"),object@meta.data$TCR),])
  p <- DimPlot(object,cells.highlight = clonal_expand, reduction = "wnn.umap", label = TRUE, sizes.highlight = 0.3) + ggtitle("TRB clones gt 1") + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
  pdf(paste(savedir,objname,"_TCR_all_UMAP_combined.pdf",sep = ""), width = 5, height = 5)
  print(p)
  dev.off()
  
  n = levels(object@meta.data$seurat_clusters)
  df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
  colnames(df) <- paste("C",n,sep = "")
  rownames(df) <- paste("C",n,sep = "")
  
  object_TCR_cluster <- table(object@meta.data$TCR, object@meta.data$seurat_clusters)
  # head(object_TCR_cluster)
  object_TCR_cluster <- object_TCR_cluster[grep("^_$",rownames(object_TCR_cluster),invert=TRUE),]
  object_TCR_cluster[object_TCR_cluster > 0] <- 1
  # test <- head(object_TCR_cluster)
  
  for (i in 1:length(df)) {
    for (j in 1:length(df)) {
      df[i,j] <- as.data.frame(table(object_TCR_cluster[,i] + object_TCR_cluster[,j] == 2)[2])[,1]
    }
  }
  
  df$cluster <- rownames(df)
  hm <- melt(df)
  colnames(hm) <- c("Cluster1","Cluster2",'Clones')
  
  hm$Cluster1 <- factor(hm$Cluster1,levels = paste("C",c(0:12),sep = ""))
  hm$Cluster2 <- factor(hm$Cluster2,levels = paste("C",c(0:12),sep = ""))
  
  library(ggplot2)
  p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=Clones)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="RdYlBu", direction=1) +
    scale_fill_gradientn(limits = c(0,100),colours=c(ArchRPalettes$solarExtra)) +
    # guides(fill=F) + # removing legend for `fill`
    labs(title = "Value distribution") + # using a title instead
    geom_text(aes(label=Clones), color="black") # printing values
  
  dir.create(paste(savedir,"/","heatmap",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"/","heatmap/",objname,"_TCR_number.pdf",sep=""))
  print(p)
  dev.off()
  
  write.table(df,paste(savedir,"/","heatmap/",objname,"_TCR_number.txt",sep=""), 
              quote = F, row.names = T, col.names = T)
  
  
  #### Celltypes Number ####
  n = unique(object@meta.data$Celltypes)
  df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
  # colnames(df) <- n
  # rownames(df) <- n
  
  object_TCR_cluster <- table(object@meta.data$TCR, object@meta.data$Celltypes)
  object_TCR_cluster <- object_TCR_cluster[grep("^_$",rownames(object_TCR_cluster),invert=TRUE),]
  head(object_TCR_cluster)
  object_TCR_cluster[object_TCR_cluster > 0] <- 1
  # test <- head(object_TCR_cluster)
  
  colnames(df) <- colnames(object_TCR_cluster)
  rownames(df) <- colnames(object_TCR_cluster)
  
  for (i in 1:length(df)) {
    for (j in 1:length(df)) {
      df[i,j] <- as.data.frame(table(object_TCR_cluster[,i] + object_TCR_cluster[,j] == 2)[2])[,1]
    }
  }
  
  # colnames(df) <- c("CM_EM","EM","Naive","TEMRA","Cycling_2")
  # rownames(df)  <- c("CM_EM","EM","Naive","TEMRA","Cycling_2")
  
  df$cluster <- rownames(df)
  hm <- melt(df)
  colnames(hm) <- c("Cluster1","Cluster2",'Clones')
  hm$Cluster1 <- factor(hm$Cluster1,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))
  hm$Cluster2 <- factor(hm$Cluster2,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))
  
  library(ggplot2)
  p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=Clones)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="RdYlBu", direction=1) +
    # guides(fill=F) + # removing legend for `fill`
    scale_fill_gradientn(limits = c(0,100),colours=c(ArchRPalettes$solarExtra)) +
    labs(title = "TRB Value distribution") + # using a title instead
    geom_text(aes(label=Clones), color="black") # printing values
  
  
  dir.create(paste(savedir,"/","heatmap",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"/","heatmap/",objname,"_celltypes_number.pdf",sep=""))
  print(p)
  dev.off()
  
  write.table(df,paste(savedir,"/","heatmap/",objname,"_Celltypes_number.txt",sep=""), 
              quote = F, row.names = T, col.names = T)
  
  
  ### Using Bhattacharya cofficient ####
  # from this paper https://www.sciencedirect.com/science/article/pii/S0092867421013337?via%3Dihub#bib28
  library(stringr)
  
  n = levels(object@meta.data$seurat_clusters)
  df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
  colnames(df) <- paste("C",n,sep = "")
  rownames(df) <- paste("C",n,sep = "")
  
  object_TCR_cluster <- table(object@meta.data$TCR, object@meta.data$seurat_clusters)
  object_TCR_cluster <- object_TCR_cluster[grep("^_$",rownames(object_TCR_cluster),invert=TRUE),]
  head(object_TCR_cluster)
  # object_TCR_cluster[object_TCR_cluster > 0] <- 1
  # test <- head(object_TCR_cluster,20) %>% tail(10) ## this is more interesting
  column_sum <- colSums(object_TCR_cluster) %>% as.vector()
  vec <- c()
  
  for (i in 1:length(df)) {
    for (j in 1:length(df)) {
      req <- object_TCR_cluster[object_TCR_cluster[,i] + object_TCR_cluster[,j] > 1,]
      if(length(req) > 0){
        if (length(req) > length(column_sum)) {
          req2 <- req[,c(i,j)]
          vec <- c()
          for (k in 1:nrow(req2)) {
            vec[k] <- sqrt((req2[k,1]/column_sum[i])*(req2[k,2]/column_sum[j]))
          }
        }
        else if (length(req) == length(column_sum)) {
          req_df <- as.data.frame(req) 
          req2 <- req_df[c(i,j),]
          vec <- c()
          vec <- sqrt((req2[1]/column_sum[i])*(req2[2]/column_sum[j]))
        }
        df[i,j] <- sum(vec)
      }
    }
  }
  
  # colnames(df) <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")
  # rownames(df)  <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")
  
  df$cluster <- rownames(df)
  hm <- melt(df)
  colnames(hm) <- c("Cluster1","Cluster2",'TRSS')
  hm$TRSS <- round(hm$TRSS, digits = 2)
  hm$Cluster1 <- factor(hm$Cluster1,levels = paste("C",c(0:12),sep = ""))
  hm$Cluster2 <- factor(hm$Cluster2,levels = paste("C",c(0:12),sep = ""))
  
  library(ggplot2)
  p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=TRSS)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="RdYlBu", direction=1) +
    # guides(fill=F) + # removing legend for `fill`
    scale_fill_gradientn(limits = c(0,0.60),colours=c(ArchRPalettes$solarExtra)) +
    labs(title = "TRB Value distribution") + # using a title instead
    geom_text(aes(label=TRSS), color="black") # printing values
  
  
  dir.create(paste(savedir,"/","heatmap",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"/","heatmap/",objname,"_cluster_TCR_only_proportion_TRSS.pdf",sep=""))
  print(p)
  dev.off()
  
  write.table(df,paste(savedir,"/","heatmap/",objname,"_cluster_TRB_only_proportion_TRSS_table.txt",sep=""), quote = F, 
              row.names = T, col.names = T, sep = "\t")
  
  ### BC Celltypes ####
  library(stringr)
  
  n = unique(object@meta.data$Celltypes)
  df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
  # colnames(df) <- n
  # rownames(df) <- n
  
  object_TCR_cluster <- table(object@meta.data$TCR, object@meta.data$Celltypes)
  object_TCR_cluster <- object_TCR_cluster[grep("^_$",rownames(object_TCR_cluster),invert=TRUE),]
  head(object_TCR_cluster)
  # object_TCR_cluster[object_TCR_cluster > 0] <- 1
  # test <- head(object_TCR_cluster,20) %>% tail(10) ## this is more interesting
  column_sum <- colSums(object_TCR_cluster) %>% as.vector()
  vec <- c()
  
  colnames(df) <- colnames(object_TCR_cluster)
  rownames(df) <- colnames(object_TCR_cluster)
  
  
  for (i in 1:length(df)) {
    for (j in 1:length(df)) {
      req <- object_TCR_cluster[object_TCR_cluster[,i] + object_TCR_cluster[,j] > 1,]
      if(length(req) > 0){
        if (length(req) > length(column_sum)) {
          req2 <- req[,c(i,j)]
          vec <- c()
          for (k in 1:nrow(req2)) {
            vec[k] <- sqrt((req2[k,1]/column_sum[i])*(req2[k,2]/column_sum[j]))
          }
        }
        else if (length(req) == length(column_sum)) {
          req_df <- as.data.frame(req) 
          req2 <- req_df[c(i,j),]
          vec <- c()
          vec <- sqrt((req2[1]/column_sum[i])*(req2[2]/column_sum[j]))
        }
        df[i,j] <- sum(vec)
      }
    }
  }
  
  # colnames(df) <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")
  # rownames(df)  <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")
  
  df$cluster <- rownames(df)
  hm <- melt(df)
  colnames(hm) <- c("Cluster1","Cluster2",'TRSS')
  hm$TRSS <- round(hm$TRSS, digits = 2)
  hm$Cluster1 <- factor(hm$Cluster1,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))
  hm$Cluster2 <- factor(hm$Cluster2,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))
  
  library(ggplot2)
  p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=TRSS)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="RdYlBu", direction=1) +
    # guides(fill=F) + # removing legend for `fill`
    scale_fill_gradientn(limits = c(0,0.50),colours=c(ArchRPalettes$solarExtra)) +
    labs(title = "TRB Value distribution") + # using a title instead
    geom_text(aes(label=TRSS), color="black") # printing values
  
  
  dir.create(paste(savedir,"/","heatmap",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"/","heatmap/",objname,"_celltypes_TCR_only_proportion_TRSS.pdf",sep=""))
  print(p)
  dev.off()
  
  write.table(df,paste(savedir,"/","heatmap/",objname,"_celltypes_TCR_only_proportion_TRSS_table.txt",sep=""), quote = F,
              row.names = T, col.names = T, sep = "\t")
  
  
}




