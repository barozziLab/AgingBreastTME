# downstream analyses for split objects in combined analysis
#1. Removal of artefact (potential doublets) clusters
#2. General analyses: statistical testing of cluster composition (young/old cells, treated/untreated cells)
#3. different DEG analyses in each of the split objects
#4. DEG analyses between old and young cells on Cluster Label level (after removing contaminated and unclear clusters)

library(tidyverse)
library(Seurat)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(readxl)
source("TME.ds_functions.R")

#~~~~~~~~~~~~~~~~~~~~~
#Load and prepare data
######################
data_folder <- "combined_analysis/split_objects"
so.split.proc <- readRDS(file.path(data_folder, "so.split.processed.rds"))
annotation_file <- "combined_analysis/split_objects/results_20230925_metadata_final_v2.xlsx"

cluster_anno <- read_xlsx(annotation_file) %>%
  mutate(compartment = ifelse(compartment == "cd45-", "stromal", compartment)) %>% 
  mutate(cluster_label = paste(Cell_Type_Label, cluster, sep = "."))
color_table <- read_tsv("plotting_colors.tsv")

# update metadata slots
for (so_name in names(so.split.proc)) {
  ca <- cluster_anno %>% filter(compartment == so_name) %>%
    mutate(cluster_label = paste(Cell_Type_Label, cluster, sep = ".")) %>%
    mutate(cluster = factor(cluster))
  
  so.split.proc[[so_name]]$cluster_label <- so.split.proc[[so_name]]@meta.data %>%
    left_join(ca, by = c("seurat_clusters" = "cluster")) %>%
    pull(cluster_label)
  
  so.split.proc[[so_name]]$age_tm <- so.split.proc[[so_name]]@meta.data %>%
    mutate(age_tm = paste(age, treatment, sep = " ")) %>%
    pull(age_tm)
  
  Idents(so.split.proc[[so_name]]) <- "cluster_label"
  so.split.proc[[so_name]]@active.ident <- factor(so.split.proc[[so_name]]@active.ident, levels = ca$cluster_label)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#UMAP Plots with final colors
#############################
out_folder <- "combined_analysis/split_objects/"
for (so_name in names(so.split.proc)) {
  
  if (so_name == "epcam+") {
    Idents(so.split.proc[[so_name]]) <- "seurat_clusters"
    umap <- UMAPPlot(so.split.proc[[so_name]], label = T, label.box = TRUE)
    ggsave(umap, file = file.path(out_folder, paste0(so_name, "_umap_final.pdf")), height = 6, width = 8)
  } else {
    col <- cluster_anno %>% filter(compartment == so_name) %>% pull(Color_Scheme_1)
    names(col) <- cluster_anno %>% filter(compartment == so_name) %>% pull(cluster_label)
    
    Idents(so.split.proc[[so_name]]) <- "cluster_label"
    umap <- UMAPPlot(so.split.proc[[so_name]], group.by = "cluster_label", label = T,
                     cols  = col[so.split.proc[[so_name]][[]] %>% pull(cluster_label) %>% unique()])
    
    Idents(so.split.proc[[so_name]]) <- "seurat_clusters"
    umap_lab <- UMAPPlot(so.split.proc[[so_name]], label = T, label.box = TRUE, group.by = "seurat_clusters", 
                         cols = cluster_anno %>% filter(compartment == so_name) %>% pull(Color_Scheme_1))
      
    if(so_name == "stromal") ggsave(umap, file = file.path(out_folder, paste0(so_name, "_umap_final.pdf")), height = 6, width = 12)
    if(so_name == "cd45+") ggsave(umap, file = file.path(out_folder, paste0(so_name, "_umap_final.pdf")), height = 6, width = 12)
    if(so_name == "stromal") ggsave(umap_lab, file = file.path(out_folder, paste0(so_name, "_umap_final_labels.pdf")), height = 6, width = 8)
    if(so_name == "cd45+") ggsave(umap_lab, file = file.path(out_folder, paste0(so_name, "_umap_final_labels.pdf")), height = 6, width = 8)
  }
  
}

#~~~~~~~~~~~~~~~~~~~~
#Epcam: clusters 7-13
#####################
out_folder <- "combined_analysis/split_objects/ds_analyses/epcam"
if(!dir.exists(out_folder)) dir.create(out_folder)

# Assess whether clusters 7-13 are artefacts in Epcam object
# Plot coverage of their marker genes across all the cells

markers_epcam <- read_tsv(file.path(data_folder, "markers_epcam_soupx.tsv"))
markers_epcam_filtered <- markers_epcam %>%
  filter(p_val_adj <= 0.05)

# Earlier exploratory approaches

# for (cl in 7:13) {
#   mg <- markers_epcam_filtered %>% 
#     filter(cluster == cl) %>% 
#     pull(gene)
#   
#   plt <- VlnPlot(so.split.proc[["epcam+"]], features = mg)
#   
#   plt <- (plt) + plot_annotation(title = paste0("Expression of Markers for cluster ", cl),
#                                  theme = theme(plot.title = element_text(size = 19,
#                                                                          face = "bold",
#                                                                          hjust = 0.5)))
#   
#   ggsave(filename = file.path(out_folder, paste0("epcam_marker_expr_cl_", cl, ".pdf")), plot = plt)
# }
# 
# dp_small <- DotPlot(so.split.proc[["epcam+"]], features = mg_small) + 
#   theme(axis.text.x = element_text(angle = -45, hjust = 0))
# 
# ggsave(filename = file.path(out_folder, "markers_small_dotplot.pdf"), plot = dp_small, width = 10, height = 8)
# 
# mg_big_log2FC <- markers_epcam_filtered %>% 
#   filter(cluster %in% c(1, 2)) %>%
#   group_by(cluster) %>%
#   slice_max(avg_log2FC, n = 10) %>%
#   pull(gene)
# 
# mg_big_pct.ratio <- markers_epcam_filtered %>% 
#   filter(cluster %in% c(1, 2)) %>%
#   group_by(cluster) %>%
#   slice_max(pct.ratio, n = 10) %>%
#   pull(gene)
# 
# dp_big_FC <- DotPlot(so.split.proc[["epcam+"]], features = mg_big_log2FC) + 
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle("sorted by log2FC")
# 
# dp_big_perc <- DotPlot(so.split.proc[["epcam+"]], features = mg_big_pct.ratio) + 
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle("sorted by perc.ratio")
# 
# ggsave(filename = file.path(out_folder, "markers_cl1_2_FC_dotplot.pdf"), plot = dp_big_FC, width = 10, height = 8)
# ggsave(filename = file.path(out_folder, "markers_cl1_2_perc_dotplot.pdf"), plot = dp_big_perc, width = 10, height = 8)
# 
# 
# dp_big_FC_filt <- DotPlot(so.split.proc[["epcam+"]], features = mg_big_logFC_filt) + 
#   theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
#   ggtitle("sorted by log2FC, filtered")

# potential figures for paper
mg_small <- markers_epcam_filtered %>% 
  filter(cluster %in% c(7:13)) %>% 
  filter(!grepl("^Rpl|^Rps|^mt-|Malat1|Neat1", gene)) %>%
  pull(gene)

mg_big_logFC_filt <-
  markers_epcam_filtered %>%
  filter(cluster %in% c(1, 2)) %>%
  filter(!grepl("^Rpl|^Rps|^mt-|Malat1|Neat1", gene)) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  pull(gene)

Idents(so.split.proc[["epcam+"]]) <- "cluster_label"
dp_combined_filtered <- custom_dotplot(so = so.split.proc[["epcam+"]], 
                                       features = c(mg_small, mg_big_logFC_filt),
                                       color_pal = "Greys")

ggsave(filename = file.path(out_folder, "markers_full_filtered_dotplot.pdf"), plot = dp_combined_filtered, width = 12, height = 4)

dp_small_cl_markers_adapted <- custom_dotplot(so = so.split.proc[["epcam+"]], 
                                              features = mg_small, 
                                              color_pal = "Greys")

ggsave(filename = file.path(out_folder, "markers_small_dotplot.pdf"), plot = dp_small_cl_markers_adapted, width = 8, height = 4)

#~~~~~~~~~~~~~~~~~~~~~~~
#Stromal DEG comparisons
########################
out_folder <- "combined_analysis/split_objects/ds_analyses/stromal"
if(!dir.exists(out_folder)) dir.create(out_folder)
  
#Figure 4A (top2 non-Haemoglobin markers in each stromal cluster) and S4A
stromal_markers <- read_tsv("combined_analysis/split_objects/markers_stromal_soupx.tsv")
top2_markers <- stromal_markers %>%
  filter(!grepl("^Hba|^Hbb|Malat1|^mt", gene)) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 2) %>%
  pull(gene) %>%
  unique()

Idents(so.split.proc[["stromal"]]) <- "cluster_label"
fig4a <- custom_dotplot(so = so.split.proc[["stromal"]],
                        features = top2_markers,
                        color_pal = "Reds",
                        cluster.idents = TRUE,
                        assay = "SCT")

ggsave(filename = file.path(out_folder, "fig4a_top2markers.pdf"), plot = fig4a, width = 10, height = 6)

stromal_minimal <- read_xlsx("markers/Stromal_Cell_Markers_Minimal.xlsx")

Idents(so.split.proc[["stromal"]]) <- "cluster_label"
fig4a_minimal <- custom_dotplot(so = so.split.proc[["stromal"]],
                                features = stromal_minimal$Gene,
                                color_pal = "Reds",
                                cluster.idents = TRUE,
                                assay = "SCT")

ggsave(filename = file.path(out_folder, "figs4_minimal_markers.pdf"), plot = fig4a_minimal, width = 7, height = 6)

#dot plot for all collagen genes that occur as markers for any cluster
stromal_comm_markers <- read_tsv("combined_analysis/split_objects/markers_stromal_soupx.tsv")

coll_markers <- stromal_comm_markers %>% 
  filter(p_val_adj <= 0.05, pct.difference >= 0.1) %>%
  filter(grepl("^Col[0-9]", gene)) %>%
  pull(gene) %>%
  unique() %>%
  sort()

coll_markers_dp <- custom_dotplot(so = so.split.proc[["stromal"]], 
                                  features = coll_markers, 
                                  color_pal = "Reds",
                                  assay = "SCT",
                                  cluster.idents = TRUE)

ggsave(filename = file.path(out_folder, "stromal_collagen_markers_dotplot.pdf"), plot = coll_markers_dp, width = 8, height = 5)

# Figure 4B (collagen gene expression across populations except Erythrocytes, Basophils, unclear and contaminated)
clus <- so.split.proc[["stromal"]]@meta.data$cluster_label %>% unique()
clus_keep <- clus[!grepl(pattern = "Granulocytes|Contaminated|Erythrocytes|Basophils|Pro-inflammatory", x = clus)]
so.stromal.clear <- subset(x = so.split.proc[["stromal"]],
                           subset = cluster_label %in% clus_keep)

coll_markers_filtered <- custom_dotplot(so = so.stromal.clear, 
               features = coll_markers, 
               color_pal = "Reds",
               assay = "SCT",
               cluster.idents = TRUE)

ggsave(filename = file.path(out_folder, "fig4b_stromal_collagen_markers_selected_dotplot.pdf"), plot = coll_markers_filtered, width = 7, height = 4)

# Figure 4C final panels: exclude basophils, erythrocytes, contaminated, unclear
# vln plot for Col1a1, Col3a1 and module score of all collagens
so.stromal.clear <- AddModuleScore(object = so.stromal.clear,
                                             features = list(coll_markers),
                                             ctrl = 5,
                                             name = "total_collagens")
so.stromal.clear$total_collagens <- so.stromal.clear$total_collagens1
so.stromal.clear$total_collagens1 <- NULL

Idents(so.stromal.clear) <- "cluster_label"
clus_order <- c("CAFs_Matrix.1", "CAFs_Matrix.2", "CAFs_Matrix.5", "CAFs_Matrix.6", 
                "CAFs_Matrix.12", "CAFs_Matrix.14", "CAFs_Matrix.18", "CAFs_Vascular.20", "Myofibroblasts.9",
                "Smooth_Muscle.10", "Endothelial.19")
so.stromal.clear@active.ident <- factor(x = so.stromal.clear@active.ident, levels = clus_order)

for (feat in c("Col1a1", "Col3a1", "total_collagens")) {
  pl <- VlnPlot(object = so.stromal.clear, 
          features = feat, 
          split.by = "age") #+
    #stat_summary(fun.y = median, shape = 95, size = 3, aes(color = split))
  
  ggsave(filename = file.path(out_folder, paste0("fig4c_", feat, ".pdf")), plot = pl, width = 7, height = 5)
}

# vln plots not split by clusters
clus_keep <- c("CAFs_Matrix.2", "CAFs_Matrix.6", "CAFs_Matrix.12", "CAFs_Matrix.14", "CAFs_Matrix.18", "CAFs_Vascular.20")
so.stromal.good_cafs <- subset(x = so.stromal.clear,
                           subset = cluster_label %in% clus_keep)
Idents(so.stromal.good_cafs) <- "age"
so.stromal.good_cafs@active.ident <- factor(so.stromal.good_cafs@active.ident, levels = c("old", "young"))

vln_col_age_minimal <- VlnPlot(so.stromal.good_cafs, features = c("Col1a1", "Col3a1", "total_collagens"))
ggsave(filename = file.path(out_folder, "fig4c_onlyCAF_wo_1_5_by_age.pdf"), width = 7, height = 5)

# Figure 4C exploratory: violin and dotplot for all collagens from stromal markers per cluster and per cluster split by age
# also: as a module (module score across all collagen genes and violin thereof, split by old and young)
vln_col <- VlnPlot(object = so.split.proc[["stromal"]], features = coll_markers)
dot_col <- custom_dotplot(so = so.split.proc[["stromal"]],
                          features = coll_markers,
                          assay = "SCT",
                          color_pal = "Reds")
vln_col_age <- VlnPlot(object = so.split.proc[["stromal"]], features = coll_markers, split.by = "age", split.plot = TRUE) +
  plot_annotation(title = 'Collagen expression in stroma, split by age')
dot_col_age <- custom_dotplot(so = so.split.proc[["stromal"]],
                              features = coll_markers,
                              assay = "SCT",
                              split.by = "age",
                              color_pal = "Reds") +
  ggtitle("Collagen expression in stromal clusters, split by age")

vln_col_tm <- VlnPlot(object = so.split.proc[["stromal"]], features = coll_markers, split.by = "treatment", split.plot = TRUE) +
  plot_annotation(title = "Collagen expression in stromal clusters, split by treatment")
dot_col_tm <- custom_dotplot(so = so.split.proc[["stromal"]],
                              features = coll_markers,
                              assay = "SCT",
                              split.by = "treatment",
                              color_pal = "Reds") +
  ggtitle("Collagen expression in stromal clusters, split by treatment")


ggsave(filename = file.path(out_folder, "vln_collagens.png"), plot = vln_col, width = 20, height = 17)
ggsave(filename = file.path(out_folder, "vln_collagens_age.png"), plot = vln_col_age, width = 20, height = 17)
ggsave(filename = file.path(out_folder, "vln_collagens_tm.png"), plot = vln_col_tm, width = 20, height = 17)

ggsave(filename = file.path(out_folder, "dot_collagens.pdf"), plot = dot_col, width = 7, height = 5)
ggsave(filename = file.path(out_folder, "dot_collagens_age.pdf"), plot = dot_col_age, width = 8, height = 9)
ggsave(filename = file.path(out_folder, "dot_collagens_tm.pdf"), plot = dot_col_tm, width = 8, height = 9)

# module score
so.split.proc[["stromal"]] <- AddModuleScore(object = so.split.proc[["stromal"]],
                                             features = list(coll_markers),
                                             ctrl = 5,
                                             name = "total_collagens")

vln_col_module <- VlnPlot(so.split.proc[["stromal"]], features = "total_collagens1", split.by = "age", split.plot = T)
ggsave(filename = file.path(out_folder, "vln_collagens_module.png"), plot = vln_col_module, width = 7, height = 5)

# Figure S4B - Feature Plots for collagen markers
for (col in coll_markers) {
  fp <- FeaturePlot(object = so.split.proc[["stromal"]], features = col)
  ggsave(filename = file.path(out_folder, paste0("FeaturePlot_", col, ".pdf")), plot = fp, width = 7, height = 6)
}

# exclude clusters 15 and 21
so.stromal.filt <- subset(x = so.split.proc[["stromal"]], 
                          subset = seurat_clusters %in% c(15, 21), 
                          invert = TRUE)

# all vascular cells vs all others
Idents(so.stromal.filt) <- "seurat_clusters"
alltbi_vs_rest.markers <- FindMarkers(so.stromal.filt, 
                                       ident.1 = c(3, 4, 5, 7, 8, 9, 10, 16, 17, 20), 
                                       ident.2 = c(1, 2, 6, 11, 12, 13, 14, 18, 19, 22), 
                                       min.pct = 0.1)

alltbi_vs_rest.markers <- alltbi_vs_rest.markers %>% 
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")
write_tsv(alltbi_vs_rest.markers, file = file.path(out_folder, "all_tbi_vs_rest_markers.tsv"), col_names = T)


# each of vascular vs all others
count = 0
for (cl in c(3, 4, 5, 7, 8, 9, 10, 16, 17, 20)) {
  tmp_markers <- FindMarkers(so.stromal.filt,
                         ident.1 = cl,
                         ident.2 = c(1, 2, 6, 11, 12, 13, 14, 18, 19, 22), 
                         min.pct = 0.1) %>%
    mutate(cluster = cl)
  
  if (count == 0) {
    tbi_single.markers <- tmp_markers
  } else {
    tbi_single.markers <- rbind(tbi_single.markers, tmp_markers)
  }
  count = count + 1
}    

tbi_single.markers <- tbi_single.markers %>% 
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")
write_tsv(tbi_single.markers, file = file.path(out_folder, "tbi_single_vs_rest_markers.tsv"), col_names = T)

# 1 vs 2
markers.1vs2 <- FindMarkers(so.stromal.filt, 
                            ident.1 = 1, 
                            ident.2 = 2, 
                            min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")
# 1 vs 6
markers.1vs6 <- FindMarkers(so.stromal.filt, 
                            ident.1 = 1, 
                            ident.2 = 6, 
                            min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")
# 1 vs 2+6
markers.1vs2_6 <- FindMarkers(so.stromal.filt, 
                            ident.1 = 1, 
                            ident.2 = c(2, 6), 
                            min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")
# 2 vs 6
markers.2vs6 <- FindMarkers(so.stromal.filt, 
                            ident.1 = 2, 
                            ident.2 = 6, 
                            min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")
# 4 vs 8
markers.4vs8 <- FindMarkers(so.stromal.filt, 
                            ident.1 = 4, 
                            ident.2 = 8, 
                            min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

write_tsv(markers.1vs2, file = file.path(out_folder, "markers.1vs2.tsv"), col_names = T)
write_tsv(markers.1vs6, file = file.path(out_folder, "markers.1vs6.tsv"), col_names = T)
write_tsv(markers.1vs2_6, file = file.path(out_folder, "markers.1vs2_6.tsv"), col_names = T)
write_tsv(markers.2vs6, file = file.path(out_folder, "markers.2vs6.tsv"), col_names = T)
write_tsv(markers.4vs8, file = file.path(out_folder, "markers.4vs8.tsv"), col_names = T)

#~~~
# DEG analyses to finalize validations: young vs. old
####
annotation <- read_xlsx(annotation_file) %>%
  mutate(cluster = as.factor(cluster))
anno_stroma <- annotation %>% filter(compartment == "cd45-")
all_CAFs <- anno_stroma %>% filter(sub_compartment == "CAF")

# in all CAFs
so.stromal.CAFs <- subset(x = so.stromal.filt, 
                          subset = seurat_clusters %in% all_CAFs$cluster)

Idents(so.stromal.CAFs) <- "age"
markers_old_vs_young_allCAF <- FindMarkers(so.stromal.CAFs, 
                                           ident.1 = "old", 
                                           ident.2 = "young", 
                                           min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

# in all CAFs, except cluster 5
so.stromal.CAFs_wo_5 <- subset(x = so.stromal.CAFs, 
                          subset = seurat_clusters == 5,
                          invert = TRUE)
markers_old_vs_young_allCAF_wo_5 <- FindMarkers(so.stromal.CAFs_wo_5, 
                                           ident.1 = "old", 
                                           ident.2 = "young", 
                                           min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

# in all CAFs, except clusters 5 and 1
so.stromal.CAFs_wo_5_1 <- subset(x = so.stromal.CAFs, 
                               subset = seurat_clusters %in% c(5, 1),
                               invert = TRUE)
markers_old_vs_young_allCAF_wo_5_1 <- FindMarkers(so.stromal.CAFs_wo_5_1, 
                                                ident.1 = "old", 
                                                ident.2 = "young", 
                                                min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

write_tsv(markers_old_vs_young_allCAF, file = file.path(out_folder, "markers.allCAFs.old_vs_young.tsv"), col_names = T)
write_tsv(markers_old_vs_young_allCAF_wo_5, file = file.path(out_folder, "markers.allCAFs_wo_5.old_vs_young.tsv"), col_names = T)
write_tsv(markers_old_vs_young_allCAF_wo_5_1, file = file.path(out_folder, "markers.allCAFs_wo_5_1.old_vs_young.tsv"), col_names = T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DEG analyses old versus young on cell type label level
#######################################################
so.stromal <- so.split.proc[["stromal"]]
annotation <- read_xlsx(annotation_file) %>%
  mutate(cluster = as.factor(cluster)) %>%
  mutate(Cluster_Label = paste0(Cell_Type_Label, "_c", cluster))

annotation_stromal <- annotation %>% filter(compartment == "cd45-")
obj_meta <- so.stromal@meta.data %>% rownames_to_column("barcodes")
joined_meta <- obj_meta %>% 
  left_join(annotation_stromal, by = c("seurat_clusters" = "cluster")) %>%
  column_to_rownames("barcodes")
so.stromal <- AddMetaData(object = so.stromal, metadata = joined_meta)

# filter (remove all unclear and contaminated clusters)
so.stromal.filt.2 <- subset(x = so.stromal, 
                          subset = sub_compartment %in% c("Contamination", "Unclear"), 
                          invert = TRUE)

# DEGs
markers_stromal_oldyoung_by_label <- list()

for (cell_type in unique(so.stromal.filt.2@meta.data$Cell_Type_Label)) {
  so.subs <- subset(x = so.stromal.filt.2, 
                    subset = Cell_Type_Label == cell_type)
  
  if(length(unique(so.subs@meta.data$age)) < 2) {
    next
  }
  
  Idents(so.subs) <- "age"
  markers_stromal_oldyoung_by_label[[cell_type]] <- FindMarkers(so.subs, 
                                                                ident.1 = "old", 
                                                                ident.2 = "young", 
                                                                min.pct = 0.1) %>%
    filter(p_val_adj <= 0.05) %>%
    mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
    rownames_to_column("gene")
  
  write_tsv(markers_stromal_oldyoung_by_label[[cell_type]], file = file.path(out_folder, paste0("markers_stromal_old_vs_young_byCTlabel_", cell_type, ".tsv")), col_names = T)
}

#~~~~~~~~~~~~~~~~~~
#Immune compartment
###################

# Investigate Tcf7 expression across T-cells.
# 1. Pooling all T-cells, DEGs between young and old conditions?
# 2. For each cluster, DEGs between young and old? (Might be tricky as some clusters contain mainly old or young)
# 3. DEGs for comparisons: cluster 17 vs 8, 1 vs 9, 1 vs 10
# 4. number of cells expressing Tcf7

out_folder <- "combined_analysis/split_objects/ds_analyses/immune"
if(!dir.exists(out_folder)) dir.create(out_folder)

# 1. Pooling all T-cells
so.tcells <- subset(so.split.proc[["cd45+"]], subset = seurat_clusters %in% c(2, 4, 5, 7, 12, 13, 14, 16, 19, 20))

Idents(so.tcells) <- "age"

markers_tcells_y_vs_o <- FindMarkers(so.tcells,
                                     ident.1 = "young",
                                     ident.2 = "old",
                                     min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

write_tsv(markers_tcells_y_vs_o, file = file.path(out_folder, "markers_tcells_y_vs_o.tsv"), col_names = T)

# 2. DEGs for each cluster individually
markers_tcells_byCluster_y_vs_o <- tibble()

for (cl in c(2, 4, 5, 7, 12, 13, 14, 16, 19, 20)) {
  c_id.1 <- rownames(so.tcells@meta.data[so.tcells@meta.data$age == "young" & so.tcells@meta.data$seurat_clusters == cl,])
  c_id.2 <- rownames(so.tcells@meta.data[so.tcells@meta.data$age == "old" & so.tcells@meta.data$seurat_clusters == cl,])
  
  if(length(c_id.1) > 5 & length(c_id.2) > 5) {
    tmp_markers <- FindMarkers(so.tcells, 
                               ident.1 = c_id.1,
                               ident.2 = c_id.2,
                               min.pct = 0.1) %>%
      filter(p_val_adj <= 0.05) %>%
      mutate(cluster = cl) %>%
      mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
      rownames_to_column("gene")
  } else {
    tmp_markers <- tibble(gene = NA, p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA, cluster = cl, pct.ratio = NA, pct.difference = NA)
  }
  markers_tcells_byCluster_y_vs_o <- rbind(markers_tcells_byCluster_y_vs_o, tmp_markers)
}

write_tsv(markers_tcells_byCluster_y_vs_o, file = file.path(out_folder, "markers_tcells_byCluster_y_vs_o.tsv"), col_names = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DEG analyses old versus young on cell type label level
#######################################################
so.immune <- so.split.proc[["cd45+"]]
annotation <- read_xlsx(annotation_file) %>%
  mutate(cluster = as.factor(cluster)) %>%
  mutate(Cluster_Label = paste0(Cell_Type_Label, "_c", cluster))

annotation_immune <- annotation %>% filter(compartment == "cd45+")
obj_meta <- so.immune@meta.data %>% rownames_to_column("barcodes")
joined_meta <- obj_meta %>% 
  left_join(annotation_immune, by = c("seurat_clusters" = "cluster")) %>%
  column_to_rownames("barcodes")
so.immune <- AddMetaData(object = so.immune, metadata = joined_meta)

# DEGs
markers_immune_oldyoung_by_label <- list()

for (cell_type in unique(so.immune@meta.data$Cell_Type_Label)) {
  so.subs <- subset(x = so.immune, 
                    subset = Cell_Type_Label == cell_type)
  
  if(length(unique(so.subs@meta.data$age)) < 2) { # skip iteration, if there is only one age group in this cell type
    next
  }
  
  Idents(so.subs) <- "age"
  
  # Avoid error for fewer than 3 cells per group
  min_cells <- so.subs@meta.data %>% group_by(age) %>% summarize(count = n()) %>% pull(count) %>% min()
  if (min_cells >= 3) {
    markers_immune_oldyoung_by_label[[cell_type]] <- FindMarkers(so.subs, 
                                                                 ident.1 = "old", 
                                                                 ident.2 = "young", 
                                                                 min.pct = 0.1) %>%
      filter(p_val_adj <= 0.05) %>%
      mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
      rownames_to_column("gene")
    
    write_tsv(markers_immune_oldyoung_by_label[[cell_type]], file = file.path(out_folder, paste0("markers_immune_old_vs_young_byCTlabel_", cell_type, ".tsv")), col_names = T)
  }
  
}


# 3. DEGs for specific comparisons
# cluster 17 vs. 8
Idents(so.immune) <- "seurat_clusters"
immune_17_vs_8 <- FindMarkers(object = so.immune, ident.1 = 17, ident.2 = 8, min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

immune_1_vs_9 <- FindMarkers(object = so.immune, ident.1 = 1, ident.2 = 9, min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

immune_1_vs_10 <- FindMarkers(object = so.immune, ident.1 = 1, ident.2 = 10, min.pct = 0.1) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

write_tsv(immune_17_vs_8, file = file.path(out_folder, "markers_immune_17vs8.tsv"))
write_tsv(immune_1_vs_9, file = file.path(out_folder, "markers_immune_1vs9.tsv"))
write_tsv(immune_1_vs_10, file = file.path(out_folder, "markers_immune_1vs10.tsv"))

# 4. Number of cells expressing Tcf7
custom_legend_text <- theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
                            legend.key.height = unit(0.4, 'cm'), #change legend key height
                            legend.key.width = unit(0.4, 'cm'), #change legend key width
                            legend.title = element_text(size = 11), #change legend title font size
                            legend.text = element_text(size = 9)) #change legend text font size


tcf7_expr <- as_tibble(as.list(so.tcells@assays$SCT@scale.data["Tcf7",])) %>% 
  pivot_longer(cols = everything(), 
               names_to = "cell_name", 
               values_to = "Tcf7_scale.data")

tcf7_data <- so.tcells@meta.data %>% 
  rownames_to_column("cell_name") %>%
  select(cell_name, treatment, age, cluster_label) %>%
  mutate(condition = paste(age, treatment, sep = "_")) %>%
  as_tibble() %>%
  inner_join(tcf7_expr, by = c("cell_name" = "cell_name"))

expr_threshold <- 0.5
Tcf7_stats <- tcf7_data %>% 
  mutate(Tcf7_expr = Tcf7_scale.data >= 0.5) %>%
  group_by(cluster_label, age, Tcf7_expr) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count))

Tcf7_count <- ggplot(Tcf7_stats %>% filter(Tcf7_expr), aes(x = age, y = count, fill = cluster_label)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Tcf7 expressing cells") +
  custom_legend_text

Tcf7_stats_plus_neg <- Tcf7_stats %>% filter(!Tcf7_expr) %>%
  group_by(age) %>%
  summarize(cluster_label = as_factor('negative'), count = sum(count), Tcf7_expr = FALSE, freq = NA) %>%
  relocate(cluster_label, age, Tcf7_expr, count, freq) %>%
  rbind(Tcf7_stats %>% filter(Tcf7_expr))

Tcf7_count_plus_neg <- ggplot(Tcf7_stats_plus_neg, aes(x = age, y = count, fill = cluster_label)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  scale_fill_manual(values = c("#a9a9a9", brewer.pal(n = 10, name = "Spectral"))) +
  ggtitle("Tcf7 positive and negative cells") +
  guides(fill = guide_legend(title = "group")) +
  custom_legend_text

freq_age <- Tcf7_stats_plus_neg %>% 
  group_by(age) %>%
  mutate(freq = count / sum(count))

Tcf7_freq_plus_neg <- ggplot(freq_age, aes(x = age, y = freq, fill = cluster_label)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  scale_fill_manual(values = c("#a9a9a9", brewer.pal(n = 10, name = "Spectral"))) +
  ggtitle("Tcf7 freq. of positive and negative cells") +
  guides(fill = guide_legend(title = "group")) +
  custom_legend_text

Tcf7_stats_filteredCounts <- Tcf7_stats %>% 
  mutate(total_count = sum(count))

Tcf7_freq <- ggplot(Tcf7_stats_filteredCounts, aes(x = cluster_label, y = freq, fill = Tcf7_expr)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~age) +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Tcf7 expression freq.") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  custom_legend_text

Tcf7_stats_adv_filteredCounts <- Tcf7_stats %>% 
  mutate(total_count = sum(count)) %>% 
  ungroup() %>%
  group_by(age) %>%
  mutate(abs_freq = count / sum(count))

Tcf7_abs_freq <- ggplot(Tcf7_stats_adv_filteredCounts, aes(x = cluster_label, y = abs_freq, fill = Tcf7_expr)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~age) +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Tcf7 expression absolute freq. (within age group)") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  custom_legend_text

Tcf7_pw <- (Tcf7_count + plot_spacer()) / (Tcf7_count_plus_neg + Tcf7_freq_plus_neg) / (Tcf7_freq) / (Tcf7_abs_freq) +
  plot_annotation(tag_levels = 'A')
Tcf7_pw
ggsave(filename = file.path(out_folder, "Tcf7_expression.pdf"), plot = Tcf7_pw, width = 8, height = 14)

#5. Number of all Tcf7 expressing cells
so.immune <- so.split.proc[["cd45+"]]

tcf7_expr_total <- as_tibble(as.list(so.immune@assays$SCT@scale.data["Tcf7",])) %>% 
  pivot_longer(cols = everything(), 
               names_to = "cell_name", 
               values_to = "Tcf7_scale.data")

tcf7_data_total <- so.immune@meta.data %>% 
  rownames_to_column("cell_name") %>%
  select(cell_name, treatment, age, cluster_label) %>%
  mutate(condition = paste(age, treatment, sep = "_")) %>%
  as_tibble() %>%
  inner_join(tcf7_expr_total, by = c("cell_name" = "cell_name"))

expr_threshold <- 0.5
Tcf7_stats_total <- tcf7_data_total %>% 
  mutate(Tcf7_expr = Tcf7_scale.data >= 0.5) %>%
  group_by(cluster_label, age, Tcf7_expr) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Supplementary: Marker gene stats
#################################
markers <- list()
for (compartment in c("cd45", "epcam", "stromal")) {
  markers[[compartment]] <- read_tsv(paste0("combined_analysis/split_objects/markers_", compartment, "_soupx.tsv")) %>%
    mutate(compartment = compartment)
}
all_markers <- rbind(markers[["cd45"]], markers[["epcam"]], markers[["stromal"]]) %>%
  mutate(direction = factor(ifelse(avg_log2FC >= 0, "up", "down")))

#Panel A: Number of markers genes (up/down-regulated) for each cluster, by compartment; stacked based on pct.difference groups

# filtering
all_significant_markers <- all_markers %>%
  filter(p_val_adj <= 0.05)

library(ggallin)
n_pct_markers <- all_significant_markers %>% 
  mutate(pct.diff.group = factor(ifelse(abs(pct.difference) <= 0.1, "< 0.1", 
                                 ifelse(abs(pct.difference) <= 0.2, "0.1-0.2", "> 0.2")), 
                                 levels = c("< 0.1", "0.1-0.2", "> 0.2"), ordered = TRUE)) %>%
  group_by(compartment, cluster, direction, pct.diff.group) %>%
  summarize(count = n()) %>%
  mutate(count = ifelse(direction == "down", -count, count)) %>%
  mutate(compartment = ifelse(compartment == "cd45", "cd45+", ifelse(compartment == "epcam", "epcam+", "stromal")))

n_pct_genes <- list()
for (comp in c("stromal", "cd45+", "epcam+")) {
  if(comp == "epcam+") {
    n_pct_genes[[comp]] <- ggplot(n_pct_markers %>% filter(compartment == comp), aes(x = cluster, y = count, fill = pct.diff.group)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~compartment) +
      theme_bw() +
      geom_hline(yintercept = 0) +
      ylim(c(-800, 1300))
  } else {
    n_pct_genes[[comp]] <- ggplot(n_pct_markers %>% filter(compartment == comp), aes(x = cluster, y = count, fill = pct.diff.group)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~compartment) +
      theme_bw() +
      theme(legend.position = "none") +
      geom_hline(yintercept = 0) +
      ylim(c(-800, 1300))
  }
}
n_pct_markers_plt <- wrap_plots(n_pct_genes, ncol = 3)
ggsave(filename = file.path(out_folder, "number_pct_diff_markers_per_cluster.pdf"), 
       n_pct_markers_plt,
       width = 10, height = 3)

# alternative way to plot it
# ggplot(n_pct_markers, aes(x = cluster, y = count, fill = pct.diff.group)) +
#   geom_bar(stat = "identity", position = "stack") +
#   facet_wrap(~compartment) +
#   theme_bw() +
#   geom_hline(yintercept = 0)

#Panel B: Number of markers genes (up/down-regulated) across clusters, by compartment; boxplot
n_markers_boxplot <- all_significant_markers %>% 
  group_by(compartment, cluster, direction) %>%
  summarize(count = n()) %>%
  mutate(compartment = ifelse(compartment == "cd45", "cd45+", ifelse(compartment == "epcam", "epcam+", "stromal"))) %>%
  mutate(compartment = factor(compartment, levels = c("stromal", "cd45+", "epcam+"))) %>%
  mutate(direction = factor(direction, levels = c("up", "down")))

n_markers_boxplot_larger0.2 <- all_significant_markers %>% 
  group_by(compartment, cluster, direction) %>%
  filter(pct.difference > 0.2) %>%
  summarize(count = n()) %>%
  mutate(compartment = ifelse(compartment == "cd45", "cd45+", ifelse(compartment == "epcam", "epcam+", "stromal"))) %>%
  mutate(compartment = factor(compartment, levels = c("stromal", "cd45+", "epcam+"))) %>%
  mutate(direction = factor(direction, levels = c("up", "down")))

n_markers_boxplot_plt <- ggplot(n_markers_boxplot, aes(x = direction, y = count)) + 
  geom_boxplot() + 
  facet_wrap(~compartment) +
  theme_bw() +
  ggtitle("Number of marker genes across clusters")
n_markers_boxplot_plt_larger02 <- ggplot(n_markers_boxplot_larger0.2, aes(x = direction, y = count)) + 
  geom_boxplot() + 
  facet_wrap(~compartment) +
  theme_bw() +
  ggtitle("Number of marker genes across clusters with pct.difference > 0.2")

ggsave(filename = file.path(out_folder, "n_markers_across_clusters.pdf"), 
       n_markers_boxplot_plt,
       width = 6, height = 3)
ggsave(filename = file.path(out_folder, "n_markers_across_clusters_pct_diff_larger0.2.pdf"), 
       n_markers_boxplot_plt_larger02,
       width = 6, height = 3)

#~~~~~~~~~~ other versions of these plots
# n_markers <- all_significant_markers %>% 
#   group_by(compartment, cluster, direction) %>%
#   summarize(count = n()) %>%
#   mutate(count = ifelse(direction == "down", -count, count)) %>%
#   mutate(compartment = ifelse(compartment == "cd45", "cd45+", ifelse(compartment == "epcam", "epcam+", "stromal")))
# 
# ct <- color_table %>% 
#   filter(condition %in% c("old untreated", "young triple_treated")) %>%
#   mutate(condition = ifelse(condition == "old untreated", "down", "up"))
# 
# n_markers <- n_markers %>% left_join(ct, by = c("compartment" = "compartment", "direction" = "condition"))
# n_genes_facet <- ggplot(n_markers, aes(x = cluster, y = count)) +
#   geom_bar(stat = "identity", position = "identity", fill = n_markers$color) +
#   scale_color_identity() +
#   facet_wrap(~compartment) +
#   theme_bw() #+
#   # scale_y_continuous(trans = pseudolog10_trans)
# 
# n_genes <- list()
# for (comp in c("cd45+", "epcam+", "stromal")) {
#   cols <- color_table %>% filter(compartment == comp)
#   n_genes[[comp]] <- ggplot(n_markers %>% filter(compartment == comp), aes(x = cluster, y = count, fill = direction)) +
#     geom_bar(stat = "identity", position = "identity") +
#     # facet_wrap(~compartment) +
#     theme_bw() +
#     ylim(c(-1000, 2000)) +
#  #   scale_y_continuous(trans = pseudolog10_trans) +
#     scale_fill_manual(values = cols$color[2:3]) +
#     ggtitle(comp)
# }
# 
# wrap_plots(n_genes, ncol = 3)

