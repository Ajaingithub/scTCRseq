# In this object is the seurat object
# savedir is the directory where you save
# clus_col is the seurat clusters
# TCR_col is where the TCR are present
# group_col for which you want to do clonal sharing like vaccine, condition etc
# group_val within that group_col which vaccine you want to extract like Z, S, or DMSO
# split_col which column to split
# column_name on splitting the column what will be the column name

TCR_AJ_VZV <- function(object, savedir, clus_col="seurat_clusters",
                       TCR_col, group_col, group_val, split_col="orig.ident", 
                       column_name, splitting = "_"){
  # library(taRifx)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  # library(Seurat)
  library(reshape2)
  library(ArchR)
  library(stringr)
  
  dir.create(savedir,showWarnings = FALSE)
  setwd(savedir)
  
  ### Extracting the metadata from the object for the TCR analysis
  require_df <- object@meta.data[,c(split_col,clus_col,TCR_col,group_col)]
  require_df_noNA <- require_df[!is.na(require_df[,TCR_col]),]
  
  ### Splitting the dataframe 
  split_to_dataframe <- function(x) {
    split_elements <- strsplit(x, splitting)[[1]]
    data.frame(t(split_elements))
  }
  
  require_df_split <- do.call(rbind, lapply(require_df_noNA[,split_col], split_to_dataframe))
  colnames(require_df_split) <- column_name
  
  if("condition" %in% colnames(require_df_split) & "vaccine" %in% colnames(require_df_split)){
    require_df_split[grep("DMSO",require_df_split$condition),"vaccine"] <- "DMSO"
  }
  
  for (i in 1:length(column_name)) {
    require_df_noNA[,column_name[i]] <- require_df_split[,column_name[i]]
  }
  
  ### Now our df is ready
  ## Subsetting the group_val ####
  message(paste("subsetting ",group_val," from ",group_col))
  require_df_noNA_subset <- require_df_noNA[grep(paste("^",group_val,"$",sep = ""),require_df_noNA[,group_col]),]
  
  # require_df_noNA_subset[,clus_col] <- factor(require_df_noNA_subset[,clus_col],levels = unique(require_df_noNA_subset[,clus_col])[order(unique(require_df_noNA_subset[,clus_col]))])
  total_clusters <- length(levels(require_df_noNA_subset$seurat_clusters))
  require_df_noNA_subset[,clus_col] <- factor(require_df_noNA_subset[,clus_col],levels = c(0:(total_clusters-1)))
  clusters = unique(require_df_noNA_subset[,clus_col])
  # n <- clusters[order(clusters)]
  n <- c(0:(total_clusters-1))
  df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
  colnames(df) <- paste("C",n,sep = "")
  rownames(df) <- paste("C",n,sep = "")
  
  object_TCR_cluster <- table(require_df_noNA_subset[,TCR_col], require_df_noNA_subset[,clus_col])
  object_TCR_cluster <- object_TCR_cluster[grep("^_$",rownames(object_TCR_cluster),invert=TRUE),]
  
  
  # clusters = unique(require_df_noNA_subset[,clus_col])
  # n <- clusters[order(clusters)]
  # df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
  # colnames(df) <- paste("C",n,sep = "")
  # rownames(df) <- paste("C",n,sep = "")
  # 
  # object_TCR_cluster <- table(require_df_noNA_subset[,TCR_col], require_df_noNA_subset[,clus_col])
  # object_TCR_cluster <- object_TCR_cluster[grep("^_$",rownames(object_TCR_cluster),invert=TRUE),]
  
  object_TCR_cluster[object_TCR_cluster > 0] <- 1
  
  for (i in 1:length(df)) {
    for (j in 1:length(df)) {
      df[i,j] <- as.data.frame(table(object_TCR_cluster[,i] + object_TCR_cluster[,j] == 2)[2])[,1]
    }
  }
  
  # df[is.na(df)] <- 0
  
  df$cluster <- rownames(df)
  hm <- melt(df)
  colnames(hm) <- c("Cluster1","Cluster2",'Clones')
  
  hm$Cluster1 <- factor(hm$Cluster1,levels = paste("C",c(0:(length(n))),sep = ""))
  hm$Cluster2 <- factor(hm$Cluster2,levels = paste("C",c(0:(length(n))),sep = ""))
  
  library(ggplot2)
  p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=Clones)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="RdYlBu", direction=1) +
    scale_fill_gradientn(limits = c(0,100),colours=c(ArchRPalettes$solarExtra)) +
    # guides(fill=F) + # removing legend for `fill`
    labs(title = paste(group_val,"_from_",group_col,sep = "")) + # using a title instead
    geom_text(aes(label=Clones), color="black") # printing values
  
  dir.create(paste(savedir,"/","clonal_sharing",sep = ""),showWarnings = FALSE)
  dir.create(paste(savedir,"/","clonal_sharing/",group_col,"/",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"/","clonal_sharing/",group_col,"/",group_val,"_from_",group_col,"_TCR_number.pdf",sep=""), width = 15, height = 10)
  print(p)
  dev.off()
  
  write.table(df,paste(savedir,"/","clonal_sharing/",group_col,"/",group_val,"_from_",group_col,"_TCR_number.txt",sep=""), 
              quote = F, row.names = T, col.names = T)
  
  
  ### Using Bhattacharya cofficient ####
  # from this paper https://www.sciencedirect.com/science/article/pii/S0092867421013337?via%3Dihub#bib28
  
  rm(df)
  df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
  colnames(df) <- paste("C",n,sep = "")
  rownames(df) <- paste("C",n,sep = "")
  
  message(paste("Running Bhattacharya cofficient for ",group_val," of ",group_col))
  object_TCR_cluster <- table(require_df_noNA_subset[,TCR_col], require_df_noNA_subset[,clus_col])
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
  
  # df[is.na(df)] <- 0

  df$cluster <- rownames(df)
  hm <- melt(df)
  colnames(hm) <- c("Cluster1","Cluster2",'TRSS')
  hm$TRSS <- round(hm$TRSS, digits = 2)
  hm$Cluster1 <- factor(hm$Cluster1,levels = paste("C",c(0:(length(n))),sep = ""))
  hm$Cluster2 <- factor(hm$Cluster2,levels = paste("C",c(0:(length(n))),sep = ""))
  
  library(ggplot2)
  p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=TRSS)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="RdYlBu", direction=1) +
    # guides(fill=F) + # removing legend for `fill`
    scale_fill_gradientn(limits = c(0,0.08),colours=c(ArchRPalettes$solarExtra)) +
    labs(title = paste(group_val,"_from_",group_col,sep = "")) + # using a title instead
    geom_text(aes(label=TRSS), color="black") # printing values
  
  
  dir.create(paste(savedir,"/","clonal_sharing",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"/","clonal_sharing/",group_col,"/",group_val,"_from_",group_col,"_cluster_TCR_only_proportion_TRSS.pdf",sep=""), width = 15, height = 10)
  print(p)
  dev.off()
  
  write.table(df,paste(savedir,"/","clonal_sharing/",group_col,"/",group_val,"_from_",group_col,"_TRSS_table.txt",sep=""), quote = F,
              row.names = T, col.names = T, sep = "\t")
}




