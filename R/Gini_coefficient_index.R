### Gini Coefficient index is calculated mainly calculated for https://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
## We can use desctools to calculate the Gini index
# Sample clonotype counts (replace with your own data)
## In this it does not matter if you take the counts or proportion give a similar output since it is divided with its sum


Gini_coef_AJ <- function(object, savedir, clus_col="seurat_clusters",
                         TCR_col, group_col, group_val, split_col="orig.ident",column_name,
                         total_clusters, splitting = "_", cellnumber=20){
  
  library(ggplot2)
  library(Seurat)
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
    
    sorted_proportions <- sort(clone_size_prop)
    n <- length(sorted_proportions)
    denominator <- n*(sum(sorted_proportions))
    
    rm(combine)
    combine <- vector()
    
    for (j in 1:length(sorted_proportions)) {
      combine[j] <- (2*j -n -1)*sorted_proportions[j]
    }
    
    sha_ent[1] = sum(combine)/denominator
    
    cluster_num <- total_clusters
    
    for (i in 1:length(cluster_num)) {
      tryCatch({
        require_df_noNA_subset_clus <- require_df_noNA_subset[grep(paste("^",cluster_num[i],"$",sep = ""), require_df_noNA_subset[,clus_col]),]
        clone_size <- as.data.frame(table(require_df_noNA_subset_clus[,TCR_col]))[,2]
        total_sum <- sum(clone_size)
        clone_size_prop <- (clone_size/total_sum)
        sorted_proportions <- sort(clone_size_prop)
        n <- length(sorted_proportions)
        denominator <- n*(sum(sorted_proportions))
        
        rm(combine)
        combine <- vector()
        
        for (j in 1:length(sorted_proportions)) {
          combine[j] <- (2*j -n -1)*sorted_proportions[j]
        }
        sha_ent[i+1] = sum(combine)/denominator
      }, error = function(e) {
        # Error handling: If an error occurs, assign NA to the result
        sha_ent[i+1] <- NA
      })
    }
  }
  return(sha_ent)
}

