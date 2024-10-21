library(dplyr)
library(Seurat)
library(clustree)
library(patchwork)

so.split.proc <- readRDS(file = "combined_analysis/split_objects/so.split.processed.rds")

compartments <- c("epcam" = "epcam+", "immune" = "cd45+", "stromal" = "stromal")
determined_res <- c("epcam" = 0.4, "immune" = 1.0, "stromal" = 0.5)

for (comp in names(compartments)) {
  so <- so.split.proc[[compartments[[comp]]]]
  
  so[["PC_1"]] <- so@reductions$pca@cell.embeddings[,"PC_1"]
  so[["PC_2"]] <- so@reductions$pca@cell.embeddings[,"PC_2"]
  so[["UMAP_1"]] <- so@reductions$umap@cell.embeddings[,"UMAP_1"]
  so[["UMAP_2"]] <- so@reductions$umap@cell.embeddings[,"UMAP_2"]
  
  full_tree <- clustree(so, prefix = "SCT_snn_res.")
  umap_full_tree <- clustree_overlay(so, prefix = "SCT_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2")
  
  res <- seq(determined_res[[comp]] + 0.1, 2.0, 0.1)
  res_sel <- paste0("SCT_snn_res.", res)
  
  so@meta.data <- so[[]] %>% select(-all_of(res_sel))
  
  sub_tree <- clustree(so, prefix = "SCT_snn_res.")
  
  umap_sub_tree <- clustree_overlay(so, prefix = "SCT_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2")
  
  out_folder <- paste0("combined_analysis/split_objects/ds_analyses/", comp)
  if (!dir.exists(out_folder)) {dir.create(comp)}
  
  ggsave(filename = file.path(comp, "clustree_full_tree.pdf"), full_tree, width = 22, height = 16)
  ggsave(filename = file.path(comp, "clustree_sub_tree.pdf"), sub_tree, width = 14, height = 7)
  ggsave(filename = file.path(comp, "clustree_umap_full_tree.png"), umap_full_tree + theme_classic(), width = 16, height = 16)
  ggsave(filename = file.path(comp, "clustree_umap_sub_tree.png"), umap_sub_tree + theme_classic(), width = 10, height = 8)
  
}