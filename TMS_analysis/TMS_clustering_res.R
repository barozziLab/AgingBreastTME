

#~~~~~~~~~~~~~~~~~~~~~
#Pre-Analysis pipeline
######################

if(run_preAnalysis | !file.exists(paste0(resultsFolder, "so_pre.rds"))) {
  
  so <- SCTransform(so, variable.features.n = n_features)
  so <- RunPCA(so, verbose = FALSE, npcs = 100)
  so <- RunUMAP(so, dims = 1:pca_dim_sel)
  so <- FindNeighbors(so, dims = 1:pca_dim_sel)
  
  saveRDS(object = so, file = paste0(resultsFolder, "so_pre.rds"))

} else {
  
  so <- readRDS(paste0(resultsFolder, "so_pre.rds"))

}

########
# assess optimal clustering resolution using silhouette score

if (run_silhouette) {
  
  stats_res <- c()
  
  # Silhouette score - identify best clustering resolution
  w <- rep(T, ncol(so))
  
  # subsample, if needed
  if(subsample_n) {
    subsample_n_real <- min(subsample_n, length(w))
    set.seed(subsample_n_real)
    w_s <- sample(w, subsample_n_real)
  } else {
    w_s <- w
  }
  
  # run clustering and compute silhouette score
  
  # Resolutions for leiden (Silhouette)
  res_grid <- seq(0.1, 2, 0.1)
  
  for (res in res_grid) {
    
    # print progress
    print(paste(format(Sys.time(), "%H:%M:%S"), 
                "Calculating Silhouette score for clustering res", res, "...", sep = " "))
    
    so <- FindClusters(so, resolution = res, algorithm = 4)
    meta_id_col <- paste("SCT_snn_res.", res, sep = "")
    
    d <- so@reductions$pca@cell.embeddings[w_s, 1:60]
    cls <- as.integer(so@meta.data[w_s,meta_id_col])
    cls <- as.integer(as.factor(cls))
    
    #keep only clusters with 2+ elements
    cls_valid <- as.integer(names(table(cls)[table(cls) >= 2]))
    w_valid <- cls %in% cls_valid
    cls <- cls[w_valid]
    d <- d[w_valid,]
    
    n_cls <- length(unique(cls))
    
    if (n_cls > 1) {
      intCriteria_res <- intCriteria(d, cls, c("Silhouette", "Dunn"))
    } else {
      intCriteria_res <- list()
      intCriteria_res$silhouette <- NA
      intCriteria_res$dunn <- NA
    }
    
    stats_res <- rbind(stats_res, c(res, n_cls, intCriteria_res$silhouette, intCriteria_res$dunn))
    
  }
  
  gc() # clean environment
  
  colnames(stats_res) <- c("resolution", "n_clusters", "silhouette", "dunn")
  
  stats_res_tbl <- as_tibble(stats_res) %>% 
    mutate(across(c(resolution, n_clusters, silhouette, dunn), as.numeric))
  
  # plotting
  silhouette <- ggplot(stats_res_tbl, aes(x = resolution, y = silhouette)) + 
    geom_bar(position = "dodge", stat = "identity", fill = "cornflowerblue") +
    ggtitle("Silhouette score") +
    ylab("Silhouette score") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = round(seq(min(stats_res_tbl$resolution), max(stats_res_tbl$resolution), by = 0.2),1)) +
    scale_y_continuous(expand = c(0, 0))
  
  dunn <- ggplot(stats_res_tbl, aes(x = resolution, y = dunn)) + 
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle("Dunn index") +
    ylab("Silhouette score") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = round(seq(min(stats_res_tbl$resolution), max(stats_res_tbl$resolution), by = 0.2),1)) +
    scale_y_continuous(expand = c(0, 0))
  
  n_clusters <- ggplot(stats_res_tbl, aes(x = resolution, y = n_clusters)) +
    geom_bar(position = "dodge", stat = "identity", fill = "coral2") +
    ggtitle("Number of clusters") +
    ylab("Number of clusters") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = round(seq(min(stats_res_tbl$resolution), max(stats_res_tbl$resolution), by = 0.2), 1)) +
    scale_y_continuous(expand = c(0, 0))
  
  silhouette_summary <- (silhouette | n_clusters)
  
  clustering_summary <- list()
  for (res in res_grid) {
    meta_id_col <- paste("SCT_snn_res.", res, sep = "")
    
    n_cells_celltype <- so@meta.data %>% group_by(!!sym(meta_id_col), cell_ontology_class) %>% summarize(n_cells = n())
    n_cells_age <- so@meta.data %>% group_by(!!sym(meta_id_col), age) %>% summarize(n_cells = n())
    
    n_cells_clusters_celltype <- ggplot(n_cells_celltype, aes(fill = cell_ontology_class, y = n_cells, x = !!sym(meta_id_col))) + 
      geom_bar(position = "stack", stat = "identity") +
      theme_classic() +
      scale_fill_brewer(palette = "Spectral") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
      ggtitle(paste0("resolution=", res))
    
    n_cells_clusters_age <- ggplot(n_cells_age, aes(fill = age, y = n_cells, x = !!sym(meta_id_col))) + 
      geom_bar(position = "stack", stat = "identity") +
      theme_classic() +
      scale_fill_brewer(palette = "Spectral") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
      ggtitle(paste0("resolution=", res))
    
    clustering_summary[[meta_id_col]] <- (silhouette | n_clusters | n_cells_clusters_celltype | n_cells_clusters_age)
  }
  
  pdf(file = paste0(resultsFolder, "silhouette_results.pdf"), width = 20, height = 4)
  for(i in 1:length(clustering_summary)) {
    plot(clustering_summary[[i]])
  }
  dev.off()
  
  write_tsv(x = stats_res_tbl, file = paste0(resultsFolder, "stats_res_tbl.tsv"))
  saveRDS(object = so, file = paste0(resultsFolder, "so_silhouette.rds"))
  
}

########
# clustree analysis

determined_res <- 1.0        # selected resolution to end clustree plots at

if (run_clustree) {
  
  so[["PC_1"]] <- so@reductions$pca@cell.embeddings[,"PC_1"]
  so[["PC_2"]] <- so@reductions$pca@cell.embeddings[,"PC_2"]
  so[["UMAP_1"]] <- so@reductions$umap@cell.embeddings[,"UMAP_1"]
  so[["UMAP_2"]] <- so@reductions$umap@cell.embeddings[,"UMAP_2"]
  
  full_tree <- clustree(so, prefix = "SCT_snn_res.")
  umap_full_tree <- clustree_overlay(so, prefix = "SCT_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2")
  
  res_sel <- paste0("SCT_snn_res.", res)
  
  so_tmp <- so
  res <- seq(determined_res + 0.1, 2.0, 0.1)
  so_tmp@meta.data <- so[[]] %>% select(-all_of(res_sel))
  
  sub_tree <- clustree(so_tmp, prefix = "SCT_snn_res.")
  
  umap_sub_tree <- clustree_overlay(so_tmp, prefix = "SCT_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2")
  
  ggsave(filename = paste0(resultsFolder, "clustree_full_tree.pdf"), full_tree, width = 22, height = 16)
  ggsave(filename = paste0(resultsFolder, "clustree_sub_tree.pdf"), sub_tree, width = 14, height = 7)
  ggsave(filename = paste0(resultsFolder, "clustree_umap_full_tree.png"), umap_full_tree + theme_classic(), width = 16, height = 16)
  ggsave(filename = paste0(resultsFolder, "clustree_umap_sub_tree.png"), umap_sub_tree + theme_classic(), width = 12, height = 10)
  
  rm(so_tmp)
  gc()
  
}

########
# multiK analysis

if (run_multiK) {
  
  emb <- Embeddings(so, reduction = "pca")
  emb <- emb[, 1:pca_dim_sel]
  
  saveRDS(object = so, file = paste0(resultsFolder, "so_tmp.rds"))
  rm(so)
  gc()
  
  multiK_res <- MultiK_par_conserveMem(input.matr = emb,
                                       reps = 120, 
                                       pSample = 0.8,
                                       numCores_first = n_cores,
                                       numCores_second = n_cores,
                                       numChunks = 8)
  
  pdf(file = paste0(resultsFolder, "multiK_plot.pdf"), width = 16, height = 6)
  plot(multiK_res$plots)
  dev.off()
  
  saveRDS(object = multiK_res, file = paste0(resultsFolder, "multiK_res.rds"))
  
  so <- readRDS(paste0(resultsFolder, "so_tmp.rds"))
  
}
