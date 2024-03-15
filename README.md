# scTCRseq
TCR Downstream analysis clonal sharing, diversity, and clonality.

The argument information is written in the R funtion

## Clonal Sharing: 
To identify sharing between different clusters, celltypes for different samples

To Run the code
        
        source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing.R")
        
        TCR_AJ_VZV(object = CD8, 
                   savedir = "./output/", 
                   clus_col = "seurat_clusters", 
                   TCR_col = "CDR3",
                   group_col = "Age_Ag", 
                   group_val = "O_EBNA1", 
                   split_col = "orig.ident",
                   column_name = c("id","Ag","Age","Age_number","gender","Run"))

## 
                   
