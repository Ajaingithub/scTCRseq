# In this object is the seurat object
# savedir is the directory where you save
# clus is the seurat clusters
# TCR col is where the TCR are present
# group for which you want to do clonal sharing like vaccine, condition etc
# group_val within that group_col which vaccine you want to extract like Z, S, or DMSO
# split_col which column to split
# column_name on splitting the column what will be the column name
# breaking how to break the expansion default is interval (1 Unique), (1,4]-expanded, (4,10] hyperexpanded, (>10 Extremely hyperexpanded)
# df labels is Unique, expanded, hyperexpanded and extremely expanded 

TCR_expansion_AJ_VZV_percent <- function(object, savedir, clus_col="seurat_clusters",
                                 TCR_col, group_col, group_val, split_col="orig.ident", column_name, 
                                 df_labels = c("expanded (abs1,0.1]","hyperexpanded (0.1,1]","extremely_hyperexpanded >1"),
                                 breaking=c(0.1,1)){
  # library(taRifx)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(Seurat)
  library(reshape2)
  library(ArchR)
  library(stringr)
  
  dir.create(savedir,showWarnings = FALSE)
  setwd(savedir)
  
  ### Absolute Number ###
  # 2. >1 clone size <= 0.1% expanded
  # 3. >0.1% clone size <= 1 hyperexpanded
  # 4. >1% clone size
  
  ### Extracting the metadata from the object for the TCR analysis
  require_df <- object@meta.data[,c(split_col,clus_col,TCR_col)]
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
  message(paste("subsetting ",group_val," from ",group_col))
  require_df_noNA_subset <- require_df_noNA[grep(paste("^",group_val,"$",sep = ""),require_df_noNA[,group_col]),]
  
  require_df_noNA_subset[,clus_col] <- factor(require_df_noNA_subset[,clus_col],levels = unique(require_df_noNA_subset[,clus_col])[order(unique(require_df_noNA_subset[,clus_col]))])
  
  clone_size <- as.data.frame(table(require_df_noNA_subset$TRB_CDR3))[,2]
  total_sum <- sum(clone_size)
  clone_size_percent <- (clone_size/total_sum)*100
  first_break <- (1/total_sum)*100
  clone_size_percent_remove_1 <- as.numeric(grep(first_break,clone_size_percent,value=TRUE,invert=TRUE))
  expanded <- cut(clone_size_percent_remove_1, breaks=c(-Inf,breaking,+Inf), labels=df_labels) %>% table() %>% as.data.frame
  return(expanded$Freq)
}










