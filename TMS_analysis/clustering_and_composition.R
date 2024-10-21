source("TME.ds_functions.R")

subFolder <- paste0(resultsFolder, "/clustering_res_", res, "/")

if(!file.exists(subFolder)) {
  dir.create(subFolder)
}

#~~~~~~~~~~~
# Clustering
############

so <- FindClusters(so, resolution = res, algorithm = 4)
saveRDS(object = so, file = paste0(subFolder, "so_clusters.rds"))

#~~~~~~~~~~~~~~~~~~~~~~
# Cell-cycle/Cell types
#######################

# plots for cell-cycle phase and cell type per cluster

# number of cells per cluster and cell cycle phase
n_cc <- so@meta.data %>% group_by(seurat_clusters, Phase) %>% summarize(n_cells = n())

n_cells_phase <- ggplot(n_cc, aes(fill = Phase, y = n_cells, x = seurat_clusters)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))

perc_cells_phase <- ggplot(n_cc, aes(fill = Phase, y = n_cells, x = seurat_clusters)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))

phase_plots <- (n_cells_phase + perc_cells_phase) +
  plot_annotation(title = "Cell cycle phases in clusters",
                  theme = theme(plot.title = element_text(size = 18,
                                                          face = "bold",
                                                          hjust = 0.5)))

ggsave(plot = phase_plots, filename = paste0(subFolder, "cell_cycle_phase_clusters.pdf"), width = 10, height = 4)

top_phases <- so@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  summarize(n_cells = n()) %>%
  filter(n_cells == max(n_cells)) %>%
  dplyr::select(-n_cells) %>%
  dplyr::rename(top_phase = Phase)

freq_phases <- so@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  summarize(n_cells = n()) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  left_join(top_phases, by = c("seurat_clusters" = "seurat_clusters")) %>%
  dplyr::select(seurat_clusters, Phase, top_phase, freq, n_cells)

write_tsv(freq_phases, file = paste0(subFolder, "cell_cycle_phase_frequencies.tsv"))

# number of cells per cluster and cell ontology class
n_cells <- so@meta.data %>% 
  group_by(seurat_clusters, cell_ontology_class) %>% 
  summarize(n_cells = n()) %>%
  mutate(frac_cells = n_cells / sum(n_cells))

n_cells_clusters <- ggplot(n_cells, aes(fill = cell_ontology_class, y = n_cells, x = seurat_clusters)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  ylab("Number of cells") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))

perc_cells_clusters <- ggplot(n_cells, aes(fill = cell_ontology_class, y = frac_cells, x = seurat_clusters)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  ylab("Fraction of cells") +
  scale_fill_brewer(palette = "Spectral") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))

celltype_plots <- cowplot::plot_grid(n_cells_clusters, perc_cells_clusters, 
                                     rel_widths = c(1, 1.6))

ggsave(plot = celltype_plots, filename = paste0(subFolder, "cell_type_clusters.pdf"), width = 12, height = 4)
write_tsv(x = n_cells, file = paste0(subFolder, "cell_type_clusters.tsv"))

# number of cells per cluster and cell ontology class, labeled

## annotate each cluster with majority cell ontology class
major_cell_ont <- so@meta.data %>% 
  group_by(seurat_clusters, cell_ontology_class) %>% 
  summarize(n_cells = n()) %>%
  slice_max(n = 1, order_by = n_cells) %>% select(-n_cells) %>% dplyr::rename(majority_cell_ont = cell_ontology_class)

so[["majority_cell_ont"]] <- so@meta.data %>% 
  left_join(major_cell_ont, by = c("seurat_clusters" = "seurat_clusters")) %>%
  pull(majority_cell_ont)

so[["cell_ont.labels"]] <- so@meta.data %>% 
  mutate(majority_cell_ont_short = str_replace_all(majority_cell_ont, " ", "_")) %>%
  mutate(majority_cell_ont_short = ifelse(majority_cell_ont_short == "luminal_epithelial_cell_of_mammary_gland", "luminal_epithelial", majority_cell_ont_short)) %>%
  mutate(cell_ont.labels = paste(majority_cell_ont_short, seurat_clusters, sep = ".")) %>%
  pull(cell_ont.labels)

lvl <- so@meta.data %>% select(cell_ont.labels, seurat_clusters) %>% as_tibble() %>% distinct() %>% arrange(seurat_clusters)

n_cells <- so@meta.data %>% 
  group_by(cell_ont.labels, cell_ontology_class) %>% 
  summarize(n_cells = n()) %>%
  mutate(frac_cells = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(desc(n_cells))

n_cells <- n_cells %>%
  mutate(cell_ont.labels = factor(cell_ont.labels, levels = lvl$cell_ont.labels, ordered = TRUE))

n_cells_clusters <- ggplot(n_cells, aes(fill = cell_ontology_class, y = n_cells, x = cell_ont.labels)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  ylab("Number of cells") +
  xlab("cluster labels") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

perc_cells_clusters <- ggplot(n_cells, aes(fill = cell_ontology_class, y = frac_cells, x = cell_ont.labels)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  ylab("Fraction of cells") +
  xlab("cluster labels") +
  scale_fill_brewer(palette = "Spectral") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

celltype_plots_labeled <- cowplot::plot_grid(n_cells_clusters, perc_cells_clusters, 
                                     rel_widths = c(1, 1.6))

ggsave(plot = celltype_plots_labeled, filename = paste0(subFolder, "cell_type_clusters_labeled.pdf"), width = 13, height = 5)
write_tsv(x = n_cells, file = paste0(subFolder, "cell_type_clusters_labeled.tsv"))

#~~~~~~~~~~~~~
# Split object
##############

# Split stroma and immune compartments according to cell ontology class

immune_types <- c("T cell", "B cell", "macrophage")

so[["age_group"]] <- so@meta.data %>% 
  mutate(age_group = ifelse(age == "3m", "young", "old")) %>% 
  pull(age_group)

so@meta.data[["compartment"]]  <- so@meta.data %>% 
  mutate(compartment = ifelse(cell_ontology_class %in% immune_types, "immune", "stromal")) %>% 
  pull(compartment)

## identify top class for each cluster
enrichment_cluster_majority <- so@meta.data %>% 
  group_by(seurat_clusters, compartment) %>% 
  summarize(n_cells = n())

enrichment_conversion_table <- enrichment_cluster_majority %>% 
  slice_max(n = 1, order_by = n_cells) %>%
  select(-n_cells) %>%
  dplyr::rename(split_group = compartment)

## add to metadata
so[["split_group"]] <- so@meta.data %>% left_join(enrichment_conversion_table, 
                                                  by = c("seurat_clusters" = "seurat_clusters")) %>%
  pull(split_group)

## split object
so.list <- SplitObject(so, split.by = "split_group")

#~~~~~~~~~~~~~
# Marker genes
##############

# Markers by cluster (full object)
markers_clusters <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1) %>% 
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2)

write_tsv(x = markers_clusters, file = paste0(subFolder, "markers_clusters_combined_object.tsv"), col_names = T)

# Markers by cluster (split object)
for (comp in c("immune", "stromal")) {
  markers <- FindAllMarkers(so.list[[comp]], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1) %>% 
    filter(p_val_adj <= 0.05) %>%
    mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2)
  
  write_tsv(x = markers, file = paste0(subFolder, "markers_clusters_split_", comp, ".tsv"), col_names = T)
  
  # top 2 markers per cluster
  top2_markers <- markers %>%
    filter(!grepl("^Hba|^Hbb|Malat1|^mt", gene)) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 2) %>%
    pull(gene) %>%
    unique()
  
  if (comp == "immune") { # define color palette
    palette = "Greens"
  } else {
    palette = "YlOrBr"
  }
  
  top2_plt <- custom_dotplot(so = so.list[[comp]],
                             features = top2_markers,
                             color_pal = palette,
                             cluster.idents = TRUE,
                             assay = "SCT")
  
  ggsave(filename = paste0(subFolder, "top2markers_", comp, ".pdf"), plot = top2_plt, width = 10, height = 6)
  
}

# DEGs old vs. young by cluster
markers_split_old_young <- list()
markers_split_old_young_full <- tibble()
for (comp in c("immune", "stromal")) {
  
  for (cluster in unique(so.list[[comp]]@meta.data$seurat_clusters) ) {
    so.subs <- subset(x = so.list[[comp]], 
                      subset = seurat_clusters == cluster)
    
    # if no cells in one or both age_groups, skip iteration
    if(length(unique(so.subs@meta.data$age_group)) < 2) {
      next
    }
    
    # compute DEGs old vs. young
    Idents(so.subs) <- "age_group"
    id <- paste(comp, cluster, sep = "_")
    markers_split_old_young[[id]] <- FindMarkers(so.subs, 
                                                 ident.1 = "old", 
                                                 ident.2 = "young", 
                                                 min.pct = 0.1) %>%
      filter(p_val_adj <= 0.05) %>%
      mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2)
    
    if (nrow(markers_split_old_young[[id]]) > 0) {
          markers_tmp <- markers_split_old_young[[id]] %>% 
            mutate(comparison = paste(comp, cluster, "old.young", sep = "_")) %>%
            rownames_to_column(var = "genes")
    }
    
    markers_split_old_young_full <- rbind(markers_split_old_young_full, markers_tmp)
    
  }
  
}

write_tsv(markers_split_old_young_full, file = paste0(subFolder, "markers_split_old_young.tsv"), col_names = T)

#~~~~~~~~~~~~~~~~~~
# T-cell exhaustion
###################

# Feature plots
feature_plot <- FeaturePlot(so.list[["immune"]], features = c("Tcf7", "Havcr2", "Entpd1"), ncol = 3)
ggsave(filename = paste0(subFolder, "t_cell_exhaustion_featureplots.pdf"), width = 15, height = 4)

# Dot plots
dotplot_clusters <- custom_dotplot(so = so.list[["immune"]], 
                                   features = c("Tcf7", "Havcr2", "Entpd1"),
                                   assay = "SCT",
                                   color_pal = "Blues",
                                   cluster.idents = TRUE)

ggsave(filename = paste0(subFolder, "exhaustion_markers_clusters.pdf"), plot = dotplot_clusters, width = 4, height = 4)

dot_split_age <- custom_dotplot(so = so.list[["immune"]], 
                                features = c("Tcf7", "Havcr2", "Entpd1"),
                                assay = "SCT",
                                color_pal = "Blues",
                                cluster.idents = FALSE,
                                split.by = "age_group")

ggsave(filename = paste0(subFolder, "exhaustion_markers_age.pdf"), plot = dot_split_age, width = 5, height = 4)


#~~~~~~~~~~~~~~~~~~~~
# Cluster composition
#####################

for (comp in c("immune", "stromal")) {
  
  ct <- colors_tme %>% filter(compartment == paste0("tms_", comp))
  col_vec <- ct$color
  names(col_vec) <- ct$condition
  
  proptest_age <- sc_utils(so.list[[comp]])
  
  proptest_age <- permutation_test(
    proptest_age, cluster_identity = "cell_ont.labels",
    sample_1 = "young", sample_2 = "old",
    sample_identity = "age_group"
  )
  
  plot_data <- return_permutation_plotdata(proptest_age)
  cl_order_age <- plot_data %>% arrange(desc(obs_log2FD)) %>% pull(clusters)
  
  prop_table_age <- plot_data %>% arrange(desc(obs_log2FD))
  write_tsv(x = prop_table_age, file = paste0(subFolder, "cluster_composition_age_split_", comp, ".tsv"))
  
  perm_plot_age <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) + 
    theme_bw() + 
    geom_hline(yintercept = log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(values = c(ct %>% 
                                    filter(condition == "young") %>% 
                                    pull(color), 
                                  "grey")) + 
    coord_flip() +
    ylim(-10, 10)
  
  n_cells <- so.list[[comp]]@meta.data %>%
    group_by(cell_ont.labels, age_group) %>% summarize(n_cells = n()) %>%
    left_join(ct, by = c("age_group" = "condition"))
  
  n_cells_age <- ggplot(n_cells, aes(fill = age_group, y = n_cells, x = factor(cell_ont.labels, levels = cl_order_age))) + 
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    scale_fill_manual(values = col_vec) +
    ylab("number of cells") +
    labs(fill = "") +
    coord_flip()
  
  perc_age <- ggplot(n_cells, aes(fill = age_group, y = n_cells, x = factor(cell_ont.labels, levels = cl_order_age))) + 
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    scale_fill_manual(values = col_vec) +
    ylab("percentage of cells") +
    coord_flip()
  
  summary_age <- (perm_plot_age | perc_age | n_cells_age) +
    plot_annotation(title = paste0("TMS", " cluster composition: old vs. young"),
                    theme = theme(plot.title = element_text(size = 16, 
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  ggsave(filename = paste0(subFolder, "TMS_cluster_composition_age_split_", comp, ".pdf"),
         summary_age,
         height = 4,
         width = 11)
  
  # minimal version for stroma
  if(comp == "stromal") {
    
    clusters_keep <- so.list[[comp]]@meta.data %>% select(majority_cell_ont, cell_ont.labels) %>% filter(majority_cell_ont %in% c("stromal cell", "endothelial cell")) %>% pull(cell_ont.labels) %>% unique()
    
    plot_data <- plot_data %>% filter(clusters %in% clusters_keep)
    perm_plot_age <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
      geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) + 
      theme_bw() + 
      geom_hline(yintercept = log2FD_threshold, lty = 2) + 
      geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
      geom_hline(yintercept = 0) + 
      scale_color_manual(values = c(ct %>% 
                                      filter(condition == "young") %>% 
                                      pull(color), 
                                    "grey")) + 
      coord_flip() +
      ylim(-10, 10)
    
    n_cells <- so.list[[comp]]@meta.data %>%
      group_by(cell_ont.labels, age_group) %>% summarize(n_cells = n()) %>%
      left_join(ct, by = c("age_group" = "condition")) %>%
      filter(cell_ont.labels %in% clusters_keep)
    
    n_cells_age <- ggplot(n_cells, aes(fill = age_group, y = n_cells, x = factor(cell_ont.labels, levels = cl_order_age))) + 
      geom_bar(position = "stack", stat = "identity") +
      theme_classic() +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
      scale_fill_manual(values = col_vec) +
      ylab("number of cells") +
      labs(fill = "") +
      coord_flip()
    
    perc_age <- ggplot(n_cells, aes(fill = age_group, y = n_cells, x = factor(cell_ont.labels, levels = cl_order_age))) + 
      geom_bar(position = "fill", stat = "identity") +
      theme_classic() +
      theme(legend.position = "none", axis.title.y = element_blank()) +
      scale_fill_manual(values = col_vec) +
      ylab("percentage of cells") +
      coord_flip()
    
    summary_age <- (perm_plot_age | perc_age | n_cells_age) +
      plot_annotation(title = paste0("TMS", " cluster composition: old vs. young"),
                      theme = theme(plot.title = element_text(size = 16, 
                                                              face = "bold",
                                                              hjust = 0.5)))
    
    ggsave(filename = paste0(subFolder, "TMS_cluster_composition_age_split_min_", comp, ".pdf"),
           summary_age,
           height = 3,
           width = 11)
    
  }
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~
# UMAPs and visualization
#########################

# combined UMAPs
umap_age <- UMAPPlot(so, group.by = "age")
umap_clusters <- UMAPPlot(so, group.by = "seurat_clusters")
umap_cell_ont <- UMAPPlot(so, group.by = "cell_ontology_class") +
  theme(legend.position = "bottom")

ggsave(filename = paste0(subFolder, "umap_age.pdf"), plot = umap_age, height = 5, width = 6)
ggsave(filename = paste0(subFolder, "umap_clusters.pdf"), plot = umap_clusters, height = 5, width = 6)
ggsave(filename = paste0(subFolder, "umap_cell_ont.pdf"), plot = umap_cell_ont, height = 6, width = 6)

#~~~~~~~~~~~~~~~~~~~~~
# Signature expression
######################

marker_tbl <- read_xlsx("Stromal_Cell_Markers_202303.xlsx")
marker_tbl <- marker_tbl %>% 
  mutate(Cell_Type = str_replace_all(Cell_Type, "\\+", "_pos")) %>%
  mutate(Ct_Study_Id = paste(Cell_Type, Study_Id, sep = "_")) %>% 
  mutate(Ct_Study_Id = str_replace_all(Ct_Study_Id, '-', ''))

# AddModuleScore to aggregate the gene expression of markers for one celltype in a cell

# mouse markers
marker_tbl_mouse <- marker_tbl %>% filter(Species == "Mouse")

vln.list.mouse <- list()
Idents(so.list[["stromal"]]) <- "seurat_clusters"
for (ct in unique(marker_tbl_mouse$Ct_Study_Id)) {
  
  feat <- list(marker_tbl_mouse %>% filter(Ct_Study_Id == ct) %>% pull(Gene))
  
  so.list[["stromal"]] <- AddModuleScore(so.list[["stromal"]],
                                         features = feat,
                                         name = ct,
                                         assay = "SCT")
  
  vln.list.mouse[[ct]] <- VlnPlot(so.list[["stromal"]], features = paste0(ct, "1")) +
    theme(legend.position = "none")
}

pw <- patchwork::wrap_plots(vln.list.mouse, ncol = 4, nrow = 8)

ggsave(filename = paste0(subFolder, "stromal_marker_expr_mouse.png"), plot = pw, height = 32, width = 20)

# human markers
marker_tbl_human <- marker_tbl %>% filter(Species == "Human")

#human->mouse orthologs map (keep only mouse gene with a single unanbigous human orthologs)
h2m <- read_tsv("human.txt", col_names = c("human", "mouse"))
cs <- table(h2m$human)
genes_single <- names(cs)[cs == 1]
h2m_single_match <- h2m %>% filter(human %in% genes_single)

marker_tbl_human <- marker_tbl_human %>%
  left_join(h2m_single_match, by = c("Gene" = "human"))


vln.list.human <- list()
Idents(so.list[["stromal"]]) <- "seurat_clusters"
for (ct in unique(marker_tbl_human$Ct_Study_Id)) {
  
  feat <- list(marker_tbl_human %>% filter(Ct_Study_Id == ct) %>% pull(mouse))
  
  so.list[["stromal"]] <- AddModuleScore(so.list[["stromal"]],
                                         features = feat,
                                         name = ct,
                                         assay = "SCT")
  
  vln.list.human[[ct]] <- VlnPlot(so.list[["stromal"]], features = paste0(ct, "1")) +
    theme(legend.position = "none")
}

pw2 <- patchwork::wrap_plots(vln.list.human, ncol = 4, nrow = 3)

ggsave(filename = paste0(subFolder, "stromal_marker_expr_human.png"), plot = pw2, height = 10, width = 20)

# heatmap with median scores per cluster and program
marker_tbl_full <- rbind(marker_tbl_human, marker_tbl_mouse %>% mutate(mouse = Gene))
Idents(so.list[["stromal"]]) <- "seurat_clusters"
for (ct in unique(marker_tbl_full$Ct_Study_Id)) {
  
  feat <- list(marker_tbl_full %>% filter(Ct_Study_Id == ct) %>% pull(mouse))
  
  so.list[["stromal"]] <- AddModuleScore(so.list[["stromal"]],
                                         features = feat,
                                         name = ct,
                                         assay = "SCT")
}

median_tbl <- so.list[["stromal"]]@meta.data %>% 
  select(seurat_clusters, ends_with("1") & !contains("SCT_snn")) %>%
  group_by(seurat_clusters) %>%
  summarize(across(2:last_col(), median))

median_matr <- as.matrix(median_tbl[,-1])
rownames(median_matr) <- median_tbl$seurat_clusters
colnames(median_matr) <- str_sub(colnames(median_matr), 1, -2)
saveRDS(median_matr, file = paste0(subFolder, "median_matrix_stromal_signatures.rds"))

pdf(paste0(subFolder, "heatmap_modules_clusters.pdf"), width = 10, height = 10)
plot(Heatmap(median_matr, 
             column_names_max_height = unit(10, "cm"), 
             heatmap_legend_param = list(title = "module score")))
dev.off()

#~~~~~~~~~~~~~~~~~~~~
# Collagen expression
#####################

# dot plot for all collagen genes that occur as markers for any cluster (in tumor data)
stromal_comm_markers <- read_tsv("markers_stromal_soupx.tsv")

coll_markers <- stromal_comm_markers %>% 
  filter(p_val_adj <= 0.05, pct.difference >= 0.1) %>%
  filter(grepl("^Col[0-9]", gene)) %>%
  pull(gene) %>%
  unique() %>%
  sort()

coll_markers_dp <- custom_dotplot(so = so.list[["stromal"]], 
                                  features = coll_markers, 
                                  color_pal = "YlOrBr",
                                  assay = "SCT",
                                  cluster.idents = TRUE)

clusters_plot <- so.list[["stromal"]]@meta.data %>%
  filter(majority_cell_ont == "stromal cell") %>% 
  pull(seurat_clusters) %>%
  unique()

coll_markers_dp_str <- custom_dotplot(so = so.list[["stromal"]], 
                                  features = coll_markers, 
                                  idents = clusters_plot,
                                  color_pal = "YlOrBr",
                                  assay = "SCT",
                                  cluster.idents = TRUE)

ggsave(filename = paste0(subFolder, "stromal_collagen_markers_dotplot.pdf"), plot = coll_markers_dp, width = 8, height = 5)
ggsave(filename = paste0(subFolder, "stromal_minimal_collagen_markers_dotplot.pdf"), plot = coll_markers_dp_str, width = 8, height = 2.5)

# similar to Figure 4C: vln plot for Col1a1, Col3a1, Col4a1, Col4a2 and module score of all collagens
so.list[["stromal"]] <- AddModuleScore(object = so.list[["stromal"]],
                                       features = list(coll_markers),
                                       assay = "SCT",
                                       name = "total_collagens")
so.list[["stromal"]]$total_collagens <- so.list[["stromal"]]$total_collagens1
so.list[["stromal"]]$total_collagens1 <- NULL

# set idents and subset object for plotting
Idents(so.list[["stromal"]]) <- "cell_ont.labels"
clusters_keep <- so.list[[comp]]@meta.data %>% select(majority_cell_ont, cell_ont.labels) %>% filter(majority_cell_ont %in% c("stromal cell", "endothelial cell")) %>% pull(cell_ont.labels) %>% unique()
so.tmp <- subset(so.list[["stromal"]], idents = clusters_keep)

for (feat in c(coll_markers, "total_collagens")) {
  pl <- VlnPlot(object = so.tmp, 
                features = feat, 
                split.by = "age_group") #+
  #stat_summary(fun.y = median, shape = 95, size = 3, aes(color = split))
  
  ggsave(filename = paste0(subFolder, paste0("tms_stromal_", feat, ".pdf")), plot = pl, width = 7, height = 5)
}

# vln plot with all stromal cells combined and all stromal and endothelial cells combined
Idents(so.tmp) <- "majority_cell_ont"

for (feat in c(coll_markers, "total_collagens")) {
  p_str <- VlnPlot(object = so.tmp,
                   idents = c("stromal cell"),
                   features = feat, 
                   split.by = "age_group")
  
  p_endo <- VlnPlot(object = so.tmp,
                    idents = c("endothelial cell"),
                    features = feat, 
                    split.by = "age_group")
  
  ggsave(filename = paste0(subFolder, paste0("tms_stromal_combined_", feat, ".pdf")), plot = p_str, width = 3, height = 4)
  ggsave(filename = paste0(subFolder, paste0("tms_endo_combined_", feat, ".pdf")), plot = p_endo, width = 3, height = 4)
}

# delete tmp object
rm(so.tmp)
gc()
# set back idents
Idents(so.list[["stromal"]]) <- "seurat_clusters"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# stromal/CAF signature expression
##################################
library(DOSE)

tme_metadata <- read_xlsx("results_20230925_metadata_final_v2.xlsx")

# load data from our own CAF/stroma populations
stroma_markers <- read_tsv("markers_stromal_soupx.tsv") %>%
  left_join(tme_metadata %>% filter(compartment == "cd45-"), by = c("cluster" = "cluster")) %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(cluster_id = paste(Cell_Type_Label, cluster, sep = "."))

tms_stromal_markers <- read_tsv(file = paste0(subFolder, "markers_clusters_split_stromal.tsv"), col_names = T) %>%
  mutate(cluster_id = paste("stroma", cluster, sep = "."))

markers_use <- stroma_markers %>% 
  select(cluster_id, gene)

enricher_res <- list()
gsea_res <- list()
for (tms_clus in unique(tms_stromal_markers$cluster_id)) {
  
  ## Marker genes
  gene_set <- tms_stromal_markers %>% filter(cluster_id == tms_clus) %>% pull(gene)
  enricher_res[[tms_clus]] <- enricher(gene = gene_set, TERM2GENE = markers_use)
  
  gs <- tms_stromal_markers %>% filter(cluster_id == tms_clus) %>% arrange(desc(avg_log2FC))
  gs_sorted <- gs$avg_log2FC
  names(gs_sorted) <- gs$gene
  gsea_res[[tms_clus]] <- GSEA(gs_sorted, TERM2GENE = markers_use)
  
  # plotting
  p1 <- enricher_res[[tms_clus]] %>% filter(p.adjust < .05, qvalue < 0.2) %>%
    mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
    ggplot(showCategory = 20, 
           aes(richFactor, fct_reorder(Description, richFactor))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(range=c(2, 10)) +
    theme_minimal() + 
    xlab("rich factor") +
    ylab(NULL) + 
    ggtitle(tms_clus)
  
  if (nrow(gsea_res[[tms_clus]]@result) > 0) {

    p2 <- dotplot(gsea_res[[tms_clus]], showCategory = 30) + ggtitle(tms_clus)
  
    grid <- cowplot::plot_grid(p1, p2)

    ggsave(filename = paste0(subFolder, "enrichment_analysis_", tms_clus, ".pdf"),
         plot = grid, width = 10, height = 6)

  } else {

    ggsave(filename = paste0(subFolder, "enrichment_analysis_", tms_clus, ".pdf"),
         plot = p1, width = 6, height = 6)

  }
  
}

## GSEA for old vs. young comparisons
markers_old_young_stroma <- markers_split_old_young_full %>%
  mutate(compartment = str_split_i(comparison, pattern = "_", i = 1)) %>%
  mutate(cluster = str_split_i(comparison, pattern = "_", i = 2))

markers_old_young_stroma <- markers_old_young_stroma %>%
  filter(compartment == "stromal") %>%
  mutate(cluster_id = paste(compartment, cluster, sep = "."))

enricher_res_y.o <- list()
gsea_res_y.o <- list()
for (clus in unique(markers_old_young_stroma$cluster_id)) {
  
  gene_set <- markers_old_young_stroma %>% filter(cluster_id == clus) %>% pull(genes)
  enricher_res_y.o[[clus]] <- enricher(gene = gene_set, TERM2GENE = markers_use)
  
  gs <- markers_old_young_stroma %>% filter(cluster_id == clus) %>% arrange(desc(avg_log2FC))
  gs_sorted <- gs$avg_log2FC
  names(gs_sorted) <- gs$genes
  gsea_res_y.o[[clus]] <- GSEA(gs_sorted, TERM2GENE = markers_use)
  
  # plotting
  if (!is.null(enricher_res_y.o[[clus]])) {
    
    p1 <- enricher_res_y.o[[clus]] %>% filter(p.adjust < .05, qvalue < 0.2) %>%
      mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
      ggplot(showCategory = 20, 
             aes(richFactor, fct_reorder(Description, richFactor))) + 
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=p.adjust, size = Count)) +
      scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
      scale_size_continuous(range=c(2, 10)) +
      theme_minimal() + 
      xlab("rich factor") +
      ylab(NULL) + 
      ggtitle(clus)
    
  } else {

    p1 <- ggplot() + theme_void()

  }
  
  if (nrow(gsea_res_y.o[[clus]]@result) > 0) {
    
    p2 <- dotplot(gsea_res_y.o[[clus]], showCategory = 30) + ggtitle(clus)
    
    grid <- cowplot::plot_grid(p1, p2)
    
    ggsave(filename = paste0(subFolder, "enrichment_analysis_young_old_", clus, ".pdf"),
           plot = grid, width = 10, height = 6)
    
  } else {
    
    ggsave(filename = paste0(subFolder, "enrichment_analysis_young_old_", clus, ".pdf"),
           plot = p1, width = 6, height = 6)
    
  }
  
}

saveRDS(object = so.list, file = paste0(subFolder, "so_split.rds"))
