# H = -Σ (P_i * log2(P_i))
# In total cluster please count from 0 so if you have 21 cluster then the length should be 22
# The Shannon diversity index measures the uncertainty of the probability distribution of species abundance. A higher Shannon diversity index
# indicates a more diverse community, because it means that there is more uncertainty about which species will be encountered. The Shannon
# diversity index is calculated as follows:
# H = -∑_i p_i log_2 p_i
# where:
# H is the Shannon diversity index
# p_i is the proportion of individuals of species i in the community
# log_2 is the logarithm to base 2

Shannon_AJ <- function(object, savedir, clus_col="seurat_clusters",
                       TCR_col, group_col, group_val, split_col="orig.ident",column_name, 
                       total_clusters, downsampled=FALSE,splitting="_", cellnumber=20){
  
  library(ggplot2)
  # library(Seurat)
  library(dplyr)
  library(Seurat)
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
  rm(sha_ent)
  message(paste("subsetting ",group_val," from ",group_col))
  require_df_noNA_subset <- require_df_noNA[grep(paste("^",group_val,"$",sep = ""),require_df_noNA[,group_col]),]
  require_df_noNA_subset[,clus_col] <- factor(require_df_noNA_subset[,clus_col],levels = unique(require_df_noNA_subset[,clus_col])[order(unique(require_df_noNA_subset[,clus_col]))])
  
  size <- nrow(require_df_noNA_subset)
  sha_ent <- rep(NA,length(total_clusters)+1)
  
  if (size > cellnumber) {
    clone_size <- as.data.frame(table(require_df_noNA_subset[,TCR_col]))[,2]
    total_sum <- sum(clone_size)
    clone_size_prop <- (clone_size/total_sum)
    sha_ent[1] = -(sum(clone_size_prop * (log2(clone_size_prop))))
    
    cluster_num <- total_clusters
    
    for (i in 1:length(cluster_num)) {
      tryCatch({
        require_df_noNA_subset_clus <- require_df_noNA_subset[grep(paste("^",cluster_num[i],"$",sep = ""), require_df_noNA_subset[,clus_col]),]
        clone_size <- as.data.frame(table(require_df_noNA_subset_clus[,TCR_col]))[,2]
        total_sum <- sum(clone_size)
        clone_size_prop <- (clone_size/total_sum)
        sha_ent[i+1]=-(sum(clone_size_prop * (log2(clone_size_prop))))
        # print(-(sum(clone_size_prop * (log2(clone_size_prop)))))
      }, error = function(e) {
        # Error handling: If an error occurs, assign NA to the result
        sha_ent[i+1] <- NA
      })
    }
  }
  return(sha_ent)
}

