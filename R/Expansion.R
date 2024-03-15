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

TCR_expansion_AJ_VZV <- function(object, savedir, clus_col="seurat_clusters",
                       TCR_col, group_col, group_val, split_col="orig.ident", column_name, 
                       df_labels = c("Unique [1]","expanded (1,4]","hyperexpanded (4,10]","extremely_hyperexpanded >10"), 
                       breaking=c(1,4,10)){
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
  # 1. clone size 1
  # 2. >1 clone size < 4
  # 3. >4 clone size < 10
  # 4. >10 clone size
  
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
  expanded <- cut(clone_size, breaks=c(-Inf,breaking,+Inf), labels=df_labels) %>% table() %>% as.data.frame
  return(expanded$Freq)
}


## To run this 
# rm(df)
# vaccine_id <- CD4@meta.data$vaccine %>% unique()
# df_labels <- c("Unique [1]","expanded (1,4]","hyperexpanded (4,10]","extremely_hyperexpanded >10")
# breaking <- c(1,4,10)
# 
# df <- data.frame(matrix(ncol = length(df_labels), nrow = length(vaccine_id)))
# rownames(df) <- vaccine_id
# colnames(df) <- df_labels
# dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
# 
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_expansion.R")
# for (i in 1:length(vaccine_id)) {
#   df[i,] <- TCR_expansion_AJ_VZV(object = CD4, savedir = savedir, clus_col = "seurat_clusters2", TCR_col = "TRB_CDR3", 
#                                  group_col = "vaccine", group_val = vaccine_id[i], split_col = "orig.ident", column_name = c("condition","vaccine","age","id","gender","Run"),
#                                  df_labels = df_labels, breaking = breaking)
#   
# }
# 
# write.table(df, paste(savedir,"expansion/vaccine.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
# 
# df$vaccine <- row.names(df)
# df_melted <- melt(df)
# colnames(df_melted) <- c("vaccine","clone_size","clone_val")
# p <- ggplot(df_melted, aes(x=clone_size,y=clone_val, fill = vaccine)) + 
#   geom_bar(stat = "identity",color="black", position = "dodge") +
#   ggtitle(paste("Clone Size for ",group_col)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.grid.major = element_line(colour = "white")) 
# 
# pdf(paste(savedir,"expansion/vaccine_bargraph.pdf",sep = ""), width = 6, height = 6)
# p
# dev.off()







