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

Shannon_AJ_perm <- function(object, savedir, clus_col="seurat_clusters",
                       TCR_col, group_col, group_val, split_col="orig.ident",column_name, 
                       total_clusters, remove_sample = NULL, perm_num=10){
  
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(Seurat)
  library(reshape2)
  library(ArchR)
  library(stringr)
  
  dir.create(savedir,showWarnings = FALSE)
  setwd(savedir)
  
  df2 <- data.frame(matrix(ncol = total_clusters+1, nrow = perm_num))
  rownames(df2) <- paste("perm",1:perm_num,sep = "")
  colnames(df2) <- c("all",paste("C",0:(total_clusters-1),sep = ""))
  
  # perm_table = rep(NA,perm_num)
  sha_ent <- rep(NA,total_clusters+1)
  for (a in 1:perm_num) {
    rm(total_cells)
    rm(cells)
    total_cells = vector()
    cells = vector()
    remove_samples <- c("DMSO","Run01")
    metadata <- object@meta.data
    metadata_remove <- metadata[grep(paste(remove_samples,collapse = "|"),rownames(metadata),invert = TRUE),]
    metadata_remove_TCR <- metadata_remove[!is.na(metadata_remove[,TCR_col]),]
    downsample_number <- min(as.data.frame(table(metadata_remove_TCR$orig.ident))[,"Freq"])
    sampleids <- unique(metadata_remove_TCR$orig.ident)
    for (b in 1:length(sampleids)) {
      cells=metadata_remove_TCR[sample(grep(sampleids[b],metadata_remove_TCR$orig.ident), size=downsample_number),"rows"]
      print(head(cells,1))
      (head(total_cells,1))
      total_cells = c(total_cells,cells)
    }
    
    # metadata_remove_equal<- metadata_remove_TCR[match(total_cells, metadata_remove_TCR$rows),]
    # 
    # print(head(metadata_remove_equal$rows,3))
    print(head(total_cells,1))
    
    # cellnames=metadata_remove_equal$rows
    object_subset <- subset(object, cells = total_cells)
    
    ### Extracting the metadata from the object_subset for the TCR analysis
    require_df <- object_subset@meta.data[,c(split_col,clus_col,TCR_col,group_col)]
    require_df_noNA <- require_df[!is.na(require_df[,TCR_col]),]
    
    ### Splitting the dataframe 
    split_to_dataframe <- function(x) {
      split_elements <- strsplit(x, "_")[[1]]
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
    
    
    clone_size <- as.data.frame(table(require_df_noNA_subset[,TCR_col]))[,2]
    total_sum <- sum(clone_size)
    clone_size_prop <- (clone_size/total_sum)
    df2[a,1] = -(sum(clone_size_prop * (log2(clone_size_prop))))
    
    cluster_num <- 0:(total_clusters-1)
    
    for (i in 1:length(cluster_num)) {
      tryCatch({
        require_df_noNA_subset_clus <- require_df_noNA_subset[grep(paste("^",cluster_num[i],"$",sep = ""), require_df_noNA_subset[,clus_col]),]
        clone_size <- as.data.frame(table(require_df_noNA_subset_clus$TRB_CDR3))[,2]
        total_sum <- sum(clone_size)
        clone_size_prop <- (clone_size/total_sum)
        df2[a,i+1]=-(sum(clone_size_prop * (log2(clone_size_prop))))
        # print(-(sum(clone_size_prop * (log2(clone_size_prop)))))
      }, error = function(e) {
        # Error handling: If an error occurs, assign NA to the result
        df2[a,i+1] <- NA
      })
    }
    
  }
  
  
  
  return(sha_ent)
}

