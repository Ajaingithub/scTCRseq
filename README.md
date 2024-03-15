# scTCRseq
TCR Downstream analysis clonal sharing, diversity, and clonality.

The argument information is written in the R funtion

### Clonal Sharing: 
Clone haring between different clusters, celltypes for different samples

**Code**
        
        source("./R/VZV_clonal_sharing.R")
        
        TCR_AJ_VZV(object = CD8, 
                   savedir = "./output/", 
                   clus_col = "seurat_clusters", 
                   TCR_col = "CDR3",
                   group_col = "Age_Ag", 
                   group_val = "O_EBNA1", 
                   split_col = "orig.ident",
                   column_name = c("id","Ag","Age","Age_number","gender","Run"))

### Clonal Diversity
Clonal diversity between individuals within the clusters or celltypes

**Shannon Diveristy**: The Shannon diversity index measures the uncertainty of the probability distribution of species abundance. A higher Shannon diversity index indicates a more diverse community, because it means that there is more uncertainty about which species will be encountered. The Shannon diversity index is calculated as follows:

H = -∑_i p_i log_2 p_i

where:
H is the Shannon diversity index

p_i is the proportion of individuals of species i in the community

log_2 is the logarithm to base 2


**Code**

        source("./R/Shannon_Diversity.R")
        
        Shannon_AJ(object = CD8, 
                   savedir = "./output/", 
                   clus_col = "seurat_clusters", 
                   TCR_col = "CDR3",
                   group_col = "orig.ident", 
                   group_val = "BB2003", 
                   split_col = "orig.ident",
                   column_name = c("id","Ag","Age","Age_number","gender","Run"), 
                   total_clusters = length(cluster_num))

**Inverse Simpson Diversity**: The Simpson diversity index measures the probability that two individuals randomly selected from the community will be of the same species. A lower Simpson diversity index indicates a more diverse community, because it means that it is less likely that two individuals will be of the same species. The Simpson diversity index is calculated as follows:

D = Σ p_i^2 ---> Simpson Index

1-D = Gini-Simpson Index or Simpson diversity Index

1/D = Inverse Simpson index

where: D is the Simpson diversity index;  p_i is the proportion of individuals of species i in the community

**Code**

        source("./R/Simpson_Diversity.R")
        
        Simpson_AJ(object = CD8, 
                   savedir = savedir,
                   clus_col = "seurat_clusters", 
                   TCR_col = "CDR3",
                   group_col = "orig.ident", 
                   group_val = sample_id[i], 
                   split_col = "orig.ident",
                   column_name = c("id","Ag","Age","Age_number","gender","Run"), 
                   total_clusters = length(cluster_num))
                             


                   

                   
