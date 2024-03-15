# D = Σ p_i^2 ---> Simpson Index
# 1-D --> Gini-Simpson Index or Simpson diversity Index
# Inverse Simpson index 1/D
# H = -Σ (P_i * log2(P_i)) ---> Shannon
# group_col which column you want to use to perform the diversity vaccine
# group_val within the group_val each sample that you wanted use like DMSO
# column_name when you split the sample id what column name that you wanted to provide.
# The Simpson diversity index measures the probability that two individuals randomly selected from the community will be of the same species. A lower Simpson diversity index indicates a more diverse community, because it means that it is less likely that two individuals will be of the same species. The Simpson diversity index is calculated as follows:
#   D = 1 - ∑_i p_i^2
#   where:
#   D is the Simpson diversity index
#   p_i is the proportion of individuals of species i in the community

Simpson_CT_AJ <- function(object, savedir, clus_col="celltypes",
                       TCR_col, group_col, group_val, split_col="orig.ident",column_name, 
                       total_clusters, splitting = "_", clus_fac){
  
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
  require_df_noNA_subset[,clus_col] <- factor(require_df_noNA_subset[,clus_col],levels = clus_fac)
  sha_ent <- rep(NA,length(clus_fac)+1)
  
  clone_size <- as.data.frame(table(require_df_noNA_subset[,TCR_col]))[,2]
  total_sum <- sum(clone_size)
  clone_size_prop <- (clone_size/total_sum)^2
  sha_ent[1] = sum(clone_size_prop)
  
  cluster_num <- clus_fac
  
  for (i in 1:length(cluster_num)) {
    tryCatch({
      require_df_noNA_subset_clus <- require_df_noNA_subset[grep(paste("^",cluster_num[i],"$",sep = ""), require_df_noNA_subset[,clus_col]),]
      clone_size <- as.data.frame(table(require_df_noNA_subset_clus[,TCR_col]))[,2]
      total_sum <- sum(clone_size)
      clone_size_prop <- (clone_size/total_sum)^2
      sha_ent[i+1]=sum(clone_size_prop)
    }, error = function(e) {
      # Error handling: If an error occurs, assign NA to the result
      sha_ent[i+1] <- NA
    })
  }
  
  return(sha_ent)
}

