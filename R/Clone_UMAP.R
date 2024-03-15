#### This function will highlight the different clusters for same clone
## Object: Seurat object
# reductions: Dimensionality reduction (PCA, UMAP)
# umap1: x axis of dimensionality reduction
# umap2: y axis of dimensionality reduction
# highlight: The cluster column that needs to be highlighted
# savedir: saving directory
# colour: color for different for different clusters

clone_UMAP_AJ <- function(object, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight, savedir, colour = c("red","green","blue")){
  p <- DimPlot(object, group.by = highlight, reduction = reductions)
  clone_data <- p$data
  size <- rep(1,length(unique(clone_data[,highlight]))-1)
  clone_data[,highlight] <- as.character(clone_data[,highlight])
  clone_data[is.na(clone_data[,highlight]),highlight] <- "ZZZ"
  indexing <- sort(clone_data[,highlight],decreasing=TRUE, index=TRUE)
  
  clone_data_ordered <- clone_data[indexing$ix,]
  
  p <- ggplot(clone_data_ordered, aes_string(umap1, umap2, color=highlight, size = highlight)) + geom_point() + 
    scale_size_manual(values = c(size,0.1)) +
    scale_color_manual(values = c(colour,"gray90")) +
    theme_bw()
  
  pdf(paste(savedir, "UMAP_cells_highlighted.pdf",sep = ""))
  print(p)
  dev.off()
}
