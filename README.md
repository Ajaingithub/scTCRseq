# scTCRseq
TCR Downstream analysis clonal sharing, diversity, and clonality.

The argument information is written in the R funtion

### Clonal Sharing: 
Clone sharing between different clusters, celltypes for different samples using TCR Receptor Similarity Score or Bhattacharya Coefficient

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
Clonal diversity between individuals within the clusters or celltypes using Shannon or Simpson diversity index

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

                   
### Gini Clonality Index: 
This method calculates the Gini coefficient (G) of inequality with bootstrap confidence intervals. A Lorenz plot is produced when a single variable is specified for analysis, otherwise the summary statistics alone are displayed for a group of variables.

![Screenshot 2024-03-15 at 4 44 35 PM](https://github.com/Ajaingithub/scTCRseq/assets/37553954/289fa13b-e792-4b07-a671-23531f877df2)

where x is an observed value, n is the number of values observed and i is the rank of values in ascending order.

        source("./R/Gini_coefficient_index.R")
        
        Gini_coef_AJ(object = CD8, 
                     savedir = savedir, 
                     clus_col = "seurat_clusters", 
                     TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                     group_col = "orig.ident", 
                     group_val = sample_id[i], 
                     split_col = "orig.ident",
                     column_name = c("id","Ag","Age","Age_number","gender","Run"), 
                     total_clusters = cluster_name)



These would calculate the Clonal Sharing, Shannon Diversity, Simpson Diversity and Gini Clonality Index




                             


                   

                   
