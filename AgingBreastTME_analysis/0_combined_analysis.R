#######################
## Analysis Pipeline ##
#######################

# Pre-clustering and re-assignment of cells to compartments, as described in Figure S2A

#~~~~~~~~~~~~~~~~~~~
#Parameter selection
####################

regress_haemo_genes <- FALSE
n_features <- 3000
pca_dim_sel <- 50
stromal_soupx <- TRUE # whether to use the soupx corrected counts for the stromal_4t1 dataset or not
subsample_n <- 0 # number of cells to subsample (for Silhouette score), set to 0 if no subsampling should be performed

#~~~~~~~~~~~
#Preparation
############
if (regress_haemo_genes) {
  data_folder <- file.path("../data/objects_after_QC")
  out_folder_combined <- "combined_analysis_regHaem/combined_object"
  out_folder_split <- "combined_analysis_regHaem/split_objects"
  
  if(!dir.exists("combined_analysis_regHaem")) dir.create("combined_analysis_regHaem")
  if(!dir.exists(out_folder_combined)) dir.create(out_folder_combined)
  if(!dir.exists(out_folder_split)) dir.create(out_folder_split)
  
} else {
  data_folder <- "../data/objects_after_QC"
  out_folder_combined <- "combined_object"
  out_folder_split <- "split_objects"
  
  if(!dir.exists(out_folder_combined)) dir.create(out_folder_combined)
  if(!dir.exists(out_folder_split)) dir.create(out_folder_split)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load and combine seurat objects
################################

if (stromal_soupx) {
  stromal_4t1 <- readRDS(file.path(data_folder, "stromal.rds"))
} else {
  stromal_4t1 <- readRDS(file.path(data_folder, "stromal_4t1.rds")) # This file is not provided, since it wasn't used.
}

tumor_immune_4t1 <- readRDS(file.path(data_folder, "tumor_and_immune.rds"))                       
pub_immune_4t1 <- readRDS(file.path(data_folder, "published_immune.rds"))

if (stromal_soupx) {
  stromal_4t1[["project_label"]] <- "stromal_4t1_soupx"
} else {
  stromal_4t1[["project_label"]] <- "stromal_4t1"
}

tumor_immune_4t1[["project_label"]] <- "tumor_immune_4t1"
pub_immune_4t1[["project_label"]] <- "pub_immune_4t1"

so_combined <- merge(stromal_4t1, y = c(tumor_immune_4t1, pub_immune_4t1), project = "combined_project")
rm(stromal_4t1, tumor_immune_4t1, pub_immune_4t1)
gc()

#~~~~~~~~~~~~~~~~~~~~~
#Pre-analysis pipeline
######################
genes <- rownames(so_combined@assays$RNA@counts)
haem_genes <- genes[grepl(pattern = "^Hba-", x = genes) | grepl(pattern = "^Hbb-", x = genes)]
haem_prog <- list(haem_prog = haem_genes)

so_combined <- AddModuleScore(so_combined, features = haem_prog, nbin = 20) # not enough data for 24 bins
so_combined$haem_prog <- so_combined$Cluster1
so_combined$Cluster1 <- NULL

if(regress_haemo_genes) {
  print("regressing additional variables ...")
  so_combined <- SCTransform(so_combined, variable.features.n = n_features, vars.to.regress = "haem_prog")
} else {
  print("no additional regression ...")
  so_combined <- SCTransform(so_combined, variable.features.n = n_features)
}

so_combined <- RunPCA(so_combined, verbose = FALSE, npcs = 100)

so_combined <- RunUMAP(so_combined, dims = 1:pca_dim_sel)

so_combined <- FindNeighbors(so_combined, dims = 1:pca_dim_sel)
so_combined <- FindClusters(so_combined, algorithm = 4)

# Cell-cycle scoring

# The conversion of human to mouse cell cycle marker genes is performed as described in this github post: https://github.com/satijalab/seurat/issues/2493 (second approach)
m.s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", 
                  target_organism = "mmusculus")$ortholog_name

m.g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", 
                    target_organism = "mmusculus")$ortholog_name

so_combined <- CellCycleScoring(so_combined, s.features = m.s.genes, g2m.features = m.g2m.genes)

# number of cells per cluster and cell cycle phase
n_cc <- so_combined@meta.data %>% group_by(seurat_clusters, Phase) %>% summarize(n_cells = n())

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

if (stromal_soupx) {
  ggsave(file.path(out_folder_combined, "cc_phase_clusters_plots_soupx.pdf"), phase_plots, height = 5, width = 14)
} else {
  ggsave(file.path(out_folder_combined, "cc_phase_clusters_plots.pdf"), phase_plots, height = 5, width = 14)
}

top_phases <- so_combined@meta.data %>% 
  group_by(seurat_clusters, Phase) %>% 
  summarize(n_cells = n()) %>% 
  filter(n_cells == max(n_cells)) %>%
  dplyr::select(-n_cells) %>%
  dplyr::rename(top_phase = Phase)

freq_phases <- so_combined@meta.data %>% 
  group_by(seurat_clusters, Phase) %>% 
  summarize(n_cells = n()) %>% 
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  left_join(top_phases, by = c("seurat_clusters" = "seurat_clusters")) %>%
  dplyr::select(seurat_clusters, Phase, top_phase, freq, n_cells)
  
if (stromal_soupx) {
  write_tsv(top_phases, file = file.path(out_folder_combined, "cc_top_phases_soupx.tsv"))
  write_tsv(freq_phases, file = file.path(out_folder_combined, "cc_freq_phases_soupx.tsv"))
} else {
  write_tsv(top_phases, file = file.path(out_folder_combined, "cc_top_phases.tsv"))
  write_tsv(freq_phases, file = file.path(out_folder_combined, "cc_freq_phases.tsv"))
}

# UMAP plots
umap_lbl <- UMAPPlot(so_combined, group.by = 'label', label = T)
umap_phase <- UMAPPlot(so_combined, group.by = 'Phase')
umap_clust <- UMAPPlot(so_combined, label = T)
feature_dbl <- FeaturePlot(so_combined, features = 'scDblFinder.score')
feature_mt <- FeaturePlot(so_combined, features = 'percent.mt')
umap_dbl <- UMAPPlot(so_combined, group.by = 'scDblFinder.class')

if (stromal_soupx) {
  ggsave(umap_lbl, file = file.path(out_folder_combined, "umap_lbl_soupx.pdf"), height = 6, width = 8)
  ggsave(umap_phase, file = file.path(out_folder_combined, "umap_phase_soupx.pdf"), height = 6, width = 8)
  ggsave(umap_clust, file = file.path(out_folder_combined, "umap_clust_soupx.pdf"), height = 6, width = 8)
  ggsave(feature_dbl, file = file.path(out_folder_combined, "feature_dbl_soupx.pdf"), height = 6, width = 8)
  ggsave(feature_mt, file = file.path(out_folder_combined, "feature_mt_soupx.pdf"), height = 6, width = 8)
  ggsave(umap_dbl, file = file.path(out_folder_combined, "umap_dbl_soupx.pdf"), height = 6, width = 8)
} else {
  ggsave(umap_lbl, file = file.path(out_folder_combined, "umap_lbl.pdf"), height = 6, width = 8)
  ggsave(umap_phase, file = file.path(out_folder_combined, "umap_phase.pdf"), height = 6, width = 8)
  ggsave(umap_clust, file = file.path(out_folder_combined, "umap_clust.pdf"), height = 6, width = 8)
  ggsave(feature_dbl, file = file.path(out_folder_combined, "feature_dbl.pdf"), height = 6, width = 8)
  ggsave(feature_mt, file = file.path(out_folder_combined, "feature_mt.pdf"), height = 6, width = 8)
  ggsave(umap_dbl, file = file.path(out_folder_combined, "umap_dbl.pdf"), height = 6, width = 8)
}

# number of cells per cluster and sample (label)
n_cells <- so_combined@meta.data %>% group_by(seurat_clusters, label) %>% summarize(n_cells = n())

n_cells_clusters <- ggplot(n_cells, aes(fill = label, y = n_cells, x = seurat_clusters)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Spectral")

# Marker genes
so_combined.markers <- FindAllMarkers(so_combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

if (stromal_soupx) {
  write_tsv(so_combined.markers, file = file.path(out_folder_combined, "markers_str_soupx.tsv"))
  write_tsv(n_cells, file = file.path(out_folder_combined, "ncells_clusters_str_soupx.tsv"))
} else {
  write_tsv(so_combined.markers, file = file.path(out_folder_combined, "markers.tsv"))
  write_tsv(n_cells, file = file.path(out_folder_combined, "ncells_clusters.tsv"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#scType cell-type annotation (combined_object)
##############################################

# run sctype for each of the marker gene sets
count = 0
cl_type_full <- tibble()
cl_type_best_full <- tibble()
for (run in names(db_list)) {
  gs_list = db_list[[run]]
  
  # get cell-type by cell matrix
  scores <- sctype_score(scRNAseqData = so_combined@assays$SCT@scale.data, 
                         scaled = TRUE, 
                         gs = gs_list$gs_positive, 
                         gs2 = gs_list$gs_negative,
                         gene_names_to_uppercase = F)
  
  # merge by cluster
  cl_cell_ont <- so_combined@meta.data$seurat_clusters
  
  cl_type <- do.call("rbind", lapply(unique(cl_cell_ont), cl_to_type_top10, scores, so_combined))
  
  # build summary tibble of cl_type
  if(count == 0) {
    cl_type_full <- cl_type %>% mutate(ann_strategy = run)
  } else {
    tmp <- cl_type %>% mutate(ann_strategy = run)
    cl_type_full <- rbind(cl_type_full, tmp)
  }
  
  cl_type_best <- cl_type %>% 
    group_by(cluster) %>% 
    top_n(n = 1, wt = score) %>%
    filter(duplicated(cluster) == FALSE) %>% # take care of ties in the score
    ungroup()
  
  cl_type_best <- cl_type_best %>%
    mutate(type = ifelse(score <= ncells / 4 | score < 0, "Unknown", type))
  
  # build summary tibble of cl_type_best
  if(count == 0) {
    cl_type_best_full <- cl_type_best %>% mutate(ann_strategy = run)
  } else {
    tmp <- cl_type_best %>% mutate(ann_strategy = run)
    cl_type_best_full <- rbind(cl_type_best_full, tmp)
  }
  
  tmp <- cl_type_best %>% select(cluster, type) %>% rename(!!run := type)
  so_combined@meta.data[run] <- so_combined@meta.data %>% left_join(tmp, by = c("seurat_clusters" = "cluster")) %>% pull(!!run)
  
  # increment counter
  count = count + 1
}

# save results of ct annotation
if (stromal_soupx) {
  write_tsv(cl_type_best_full, file = file.path(out_folder_combined, "cluster_sctype_annotation_full_str_soupx.txt"))
  write_tsv(cl_type_full, file = file.path(out_folder_combined, "cluster_sctype_extended_full_str_soupx.txt"))
} else {
  write_tsv(cl_type_best_full, file = file.path(out_folder_combined, "cluster_sctype_annotation_full.txt"))
  write_tsv(cl_type_full, file = file.path(out_folder_combined, "cluster_sctype_extended_full.txt"))
}

#~~~~~~~~~~~~~~~~~
#split by clusters
##################

# We are splitting the dataset by majority of cells in each cluster, meaning that
# if a cluster contains most cells from one of the captures, all cells will be assigned
# to that seurat object (the argument behind that is that probably in the sorting/enrichment
# process some contamination occured that would be "corrected" according to this
# initial clustering on the whole dataset)

# annotate new set of category (enrichment_group): cd45+, epcam+, stromal
so_combined[["enrichment_group"]] <- so_combined@meta.data %>% 
  mutate(enrichment_group = ifelse(str_detect(label, ".*cd45"), "cd45+",  # immune cell
                             ifelse(str_detect(label, ".*epcam"), "epcam+",  # cancer cell
                                    "stromal")))  %>% # stromal cell
  pull(enrichment_group)

# identify top class for each cluster
enrichment_cluster_majority <- so_combined@meta.data %>% 
  group_by(seurat_clusters, enrichment_group) %>% 
  summarize(n_cells = n())

# save table for inspection (e.g. if there are any close calls, due to bad initial clustering)
if (stromal_soupx) {
  write_tsv(enrichment_cluster_majority, 
            file = file.path(file.path(out_folder_split, "split_enrichment_cluster_majority_soupx.tsv")))
} else {
  write_tsv(enrichment_cluster_majority, 
            file = file.path(file.path(out_folder_split, "split_enrichment_cluster_majority_soupx.tsv")))
}

# add annotation to metadata
enrichment_conversion_table <- enrichment_cluster_majority %>% 
  slice_max(n = 1, order_by = n_cells) %>%
  select(-n_cells) %>%
  dplyr::rename(split_group = enrichment_group)

so_combined[["split_group"]] <- so_combined@meta.data %>% left_join(enrichment_conversion_table, 
                                    by = c("seurat_clusters" = "seurat_clusters")) %>%
  pull(split_group)

# save combined seurat object
if (stromal_soupx) {
  saveRDS(so_combined, file = file.path(out_folder_combined, "so_combined_soupx.rds"))
} else {
  saveRDS(so_combined, file = file.path(out_folder_combined, "so_combined.rds"))
}

# plot cell numbers before and after reassignment
cell_nums_original <- so_combined@meta.data %>% 
  group_by(enrichment_group, treatment, age) %>% 
  summarize(count = n()) %>% 
  dplyr::rename(compartment = enrichment_group) %>% 
  mutate(reassignment = factor("before reassignment", levels = c("before reassignment", "after reassignment"), ordered = T))
cell_nums_reassigned <- so_combined@meta.data %>% 
  group_by(split_group, treatment, age) %>% 
  summarize(count = n()) %>% 
  dplyr::rename(compartment = split_group) %>% 
  mutate(reassignment = factor("after reassignment", levels = c("before reassignment", "after reassignment"), ordered = T))

cell_nums_full <- rbind(cell_nums_original, cell_nums_reassigned) %>% 
  mutate(treatment = ifelse(treatment == "triple_treated", "treatment", "ctrl")) %>%
  mutate(condition = factor(paste(age, treatment, sep = " "), levels = c("old ctrl", "old treatment", "young ctrl", "young treatment"), ordered = T))

cell_nums_plt <- ggplot(cell_nums_full, aes(x = condition, y = count, fill = reassignment)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_classic() +
  ylab("Number of cells") +
  xlab("Condition") +
  labs(fill = "") +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0)) +
  facet_wrap(~compartment)

ggsave(cell_nums_plt, file = file.path(out_folder_split, "cell_numbers_reassignment.pdf"), device = "pdf", width = 6, height = 4)

# split seurat object
so.split <- SplitObject(so_combined, split.by = "split_group")
rm(so_combined)
gc()

# Resolutions for leiden (Silhouette)
res_grid <- seq(0.2, 2, 0.1)

# analysis pipeline
stats_res <- c()
for (so_name in names(so.split)) {
  so <- so.split[[so_name]]
  so[["seurat_clusters"]] <- NULL
  so <- RunPCA(so, verbose = FALSE, npcs = 100)
  
  # Silhouette score - identify best clustering resolution
  # parts of the code in this section were originally written by Iros Barozzi
  w <- rep(T, ncol(so))
  
  # subsample, if needed
  if(subsample_n) {
    subsample_n_real <- min(subsample_n, length(w))
    set.seed(subsample_n_real)
    w_s <- sample(w, subsample_n_real)
  } else {
    w_s <- w
  }
  
  # run clustering and calculate silhouette score
  so <- FindNeighbors(so, dims = 1:pca_dim_sel)
  
  for (res in res_grid) {
    
    # print progress
    print(paste(format(Sys.time(), "%H:%M:%S"), 
                "Calculating Silhouette score for object", 
                so_name, "and clustering res", res, "...", sep = " "))
    
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
    
    stats_res <- rbind(stats_res, c(so_name, res, n_cls, intCriteria_res$silhouette, intCriteria_res$dunn))
    
  }
  
  # update seurat object
  so.split[[so_name]] <- so
  
}

rm(so)
gc()

colnames(stats_res) <- c("enrichment_group", "resolution", "n_clusters", "silhouette", "dunn")

stats_res_tbl <- as_tibble(stats_res) %>% 
  mutate(across(c(resolution, n_clusters, silhouette, dunn), as.numeric)) %>% 
  mutate(across(c(enrichment_group), factor))

# plotting
silhouette <- ggplot(stats_res_tbl, aes(x = resolution, y = silhouette, fill = enrichment_group)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~enrichment_group) + 
  ggtitle("Silhouette score") +
  ylab("Silhouette score") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = round(seq(min(stats_res_tbl$resolution), max(stats_res_tbl$resolution), by = 0.2),1)) +
  scale_y_continuous(expand = c(0, 0))

dunn <- ggplot(stats_res_tbl, aes(x = resolution, y = dunn, fill = enrichment_group)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~enrichment_group) + 
  ggtitle("Dunn index") +
  ylab("Silhouette score") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = round(seq(min(stats_res_tbl$resolution), max(stats_res_tbl$resolution), by = 0.2),1)) +
  scale_y_continuous(expand = c(0, 0))

n_clusters <- ggplot(stats_res_tbl, aes(x = resolution, y = n_clusters, fill = enrichment_group)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~enrichment_group) +
  ggtitle("Number of clusters") +
  ylab("Number of clusters") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = round(seq(min(stats_res_tbl$resolution), max(stats_res_tbl$resolution), by = 0.2), 1)) +
  scale_y_continuous(expand = c(0, 0))

silhouette_summary <- (silhouette | n_clusters)

if (stromal_soupx) {
  ggsave(silhouette_summary, 
         file = file.path(out_folder_split, "silhouette_summary_soupx.pdf"),
         width = 16, 
         height = 6)
} else {
  ggsave(silhouette_summary, 
         file = file.path(out_folder_split, "silhouette_summary.pdf"),
         width = 16, 
         height = 6)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Analysis with inferred clustering resolution
#############################################

if(regress_haemo_genes) {
  final_res <- list("cd45+" = 1, "epcam+" = 0.4, "stromal" = 0.4)
} else {
  final_res <- list("cd45+" = 1, "epcam+" = 0.4, "stromal" = 0.5)
}


sav_n <- c("cd45+" = "cd45", "epcam+" = "epcam", "stromal" = "stromal")

n_cells_clusters =
  perc_cells_clusters =
  cells_clusters = 
  n_cells_phase = 
  perc_cells_phase =
  phase_plots =
  n_cells_dbl =
  perc_cells_dbl =
  dbl_plots =
  n_feat =
  n_count =
  mt =
  dbl_score =
  quality_stats <- list()
so.split.processed <- list()

for (so_name in names(so.split)) {
  so <- so.split[[so_name]]
  res <- final_res[[so_name]]
  n <- sav_n[so_name]
  
  # print progress
  print(paste0("Analysing ", so_name, " with resolution: ", res, " ..."))
  
  # clustering and cc scoring (PCA was already run in Silhouette workflow)
  if(regress_haemo_genes) {
    so <- FindClusters(so, resolution = res, algorithm = 4) %>% # clustering, FindNeighbors was already run for determining silhouette scores above
      CellCycleScoring(s.features = m.s.genes, g2m.features = m.g2m.genes, nbin = 20) %>% # cc scoring, nbin = 20 because the default 24 does not work for cd45+ cells
      RunUMAP(dims = 1:pca_dim_sel)
  } else {
    so <- FindClusters(so, resolution = res, algorithm = 4) %>% # clustering, FindNeighbors was already run for determining silhouette scores above
      CellCycleScoring(s.features = m.s.genes, g2m.features = m.g2m.genes, nbin = 22) %>% # cc scoring, nbin = 22 because the default 24 does not work for cd45+ cells
      RunUMAP(dims = 1:pca_dim_sel)
  }
  
  # marker genes
  markers <- FindAllMarkers(so, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1) 
  markers <- markers %>% mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2)
  
  if (stromal_soupx) {
    write_tsv(markers, file = file.path(out_folder_split, 
                                                    paste0("markers_", n, "_soupx.tsv")))
  } else {
    write_tsv(markers, file = file.path(out_folder_split,
                                                    paste0("markers_", n, ".tsv")))
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #scType cell-type annotation
  ############################
  
  # run sctype for each of the marker gene sets
  count = 0
  cl_type_full <- tibble()
  cl_type_best_full <- tibble()
  for (run in names(db_list)) {
    gs_list = db_list[[run]]
    
    # get cell-type by cell matrix
    scores <- sctype_score(scRNAseqData = so@assays$SCT@scale.data, 
                           scaled = TRUE, 
                           gs = gs_list$gs_positive, 
                           gs2 = gs_list$gs_negative,
                           gene_names_to_uppercase = F)
    
    # merge by cluster
    cl_cell_ont <- so@meta.data$seurat_clusters
    
    cl_type <- do.call("rbind", lapply(unique(cl_cell_ont), cl_to_type_top10, scores, so))
    
    # build summary tibble of cl_type
    if(count == 0) {
      cl_type_full <- cl_type %>% mutate(ann_strategy = run)
    } else {
      tmp <- cl_type %>% mutate(ann_strategy = run)
      cl_type_full <- rbind(cl_type_full, tmp)
    }
    
    cl_type_best <- cl_type %>% 
      group_by(cluster) %>% 
      top_n(n = 1, wt = score) %>%
      filter(duplicated(cluster) == FALSE) %>% # take care of ties in the score
      ungroup()
    
    cl_type_best <- cl_type_best %>%
      mutate(type = ifelse(score <= ncells / 4 | score < 0, "Unknown", type))
    
    # build summary tibble of cl_type_best
    if(count == 0) {
      cl_type_best_full <- cl_type_best %>% mutate(ann_strategy = run)
    } else {
      tmp <- cl_type_best %>% mutate(ann_strategy = run)
      cl_type_best_full <- rbind(cl_type_best_full, tmp)
    }
    
    col <- paste0(run, "_split")
    tmp <- cl_type_best %>% select(cluster, type) %>% dplyr::rename(!!col := type)
    so@meta.data[col] <- so@meta.data %>% left_join(tmp, by = c("seurat_clusters" = "cluster")) %>% pull(!!col)
    
    # increment counter
    count = count + 1
  }
  
  # save sctype results
  if(stromal_soupx) {
    write_tsv(cl_type_best_full, 
              file = file.path(out_folder_split, 
                               paste0("cluster_sctype_annotation_", n, "_soupx.txt")))
    write_tsv(cl_type_full, 
              file = file.path(out_folder_split, 
                               paste0("cluster_sctype_annotation_extended_", n, "_soupx.txt")))
  } else {
    write_tsv(cl_type_best_full, 
              file = file.path(out_folder_split, 
                               paste0("cluster_sctype_annotation_", n, ".txt")))
    write_tsv(cl_type_full, 
              file = file.path(out_folder_split, 
                               paste0("cluster_sctype_annotation_extended_", n, ".txt")))
  }
  
  #~~~~~~~~~~~~~~~
  #stats and plots
  ################
  
  # number of cells per cluster and sample (label)
  n_cells <- so@meta.data %>% group_by(seurat_clusters, label) %>% summarize(n_cells = n())
  
  n_cells_clusters[[so_name]] <- ggplot(n_cells, aes(fill = label, y = n_cells, x = seurat_clusters)) + 
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral") +
    theme(legend.position = "none")
  
  perc_cells_clusters[[so_name]] <- ggplot(n_cells, aes(fill = label, y = n_cells, x = seurat_clusters)) + 
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral")
  
  cells_clusters[[so_name]] <- (n_cells_clusters[[so_name]] + perc_cells_clusters[[so_name]]) +
    plot_annotation(title = paste0(so_name, ": Cells per sample in clusters"),
                    theme = theme(plot.title = element_text(size = 18, 
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  if (stromal_soupx) {
    ggsave(file.path(out_folder_split, paste0(n, "_labels_clusters_plots_soupx.pdf")), cells_clusters[[so_name]], height = 5, width = 14)
  } else {
    ggsave(file.path(out_folder_split, paste0(n, "_labels_clusters_plots.pdf")), cells_clusters[[so_name]], height = 5, width = 14)
  }
  
  # number of cells per cc phase
  n_cc <- so@meta.data %>% group_by(seurat_clusters, Phase) %>% summarize(n_cells = n())
  
  n_cells_phase[[so_name]] <- ggplot(n_cc, aes(fill = Phase, y = n_cells, x = seurat_clusters)) + 
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
  
  perc_cells_phase[[so_name]] <- ggplot(n_cc, aes(fill = Phase, y = n_cells, x = seurat_clusters)) + 
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral") +
    ylab("fraction of cells") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
  
  phase_plots[[so_name]] <- (n_cells_phase[[so_name]] + perc_cells_phase[[so_name]]) +
    plot_annotation(title = "Cell cycle phases in clusters",
                    theme = theme(plot.title = element_text(size = 18, 
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  # barplots doublet detection
  n_dbl <- so@meta.data %>% group_by(seurat_clusters, scDblFinder.class) %>% summarize(n_cells = n())
  
  n_cells_dbl[[so_name]] <- ggplot(n_dbl, aes(fill = scDblFinder.class, y = n_cells, x = seurat_clusters)) + 
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
  
  perc_cells_dbl[[so_name]] <- ggplot(n_dbl, aes(fill = scDblFinder.class, y = n_cells, x = seurat_clusters)) + 
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral") +
    ylab("fraction of cells") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
  
  dbl_plots[[so_name]] <- (n_cells_dbl[[so_name]] + perc_cells_dbl[[so_name]]) +
    plot_annotation(title = "Predicted doublets in clusters",
                    theme = theme(plot.title = element_text(size = 18, 
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  if (stromal_soupx) {
    ggsave(file.path(out_folder_split, paste0(n, "_cc_phase_clusters_plots_soupx.pdf")), phase_plots[[so_name]], height = 5, width = 14)
    ggsave(file.path(out_folder_split, paste0(n, "_dbl_clusters_plots_soupx.pdf")), dbl_plots[[so_name]], height = 5, width = 14)
  } else {
    ggsave(file.path(out_folder_split, paste0(n, "_cc_phase_clusters_plots.pdf")), phase_plots[[so_name]], height = 5, width = 14)
    ggsave(file.path(out_folder_split, paste0(n, "_dbl_clusters_plots.pdf")), dbl_plots[[so_name]], height = 5, width = 14)
  }
  
  # quality stats
  mt[[so_name]] <- VlnPlot(so, features = "percent.mt", group.by = "seurat_clusters", pt.size = 0.1) +
    theme(legend.position = "none") +
    ggtitle("Percentage mitochondrial reads")
  n_feat[[so_name]] <- VlnPlot(so, features = "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0.1) +
    theme(legend.position = "none") +
    ggtitle("Number of detected genes")
  n_count[[so_name]] <- VlnPlot(so, features = "nCount_RNA", group.by = "seurat_clusters", pt.size = 0.1) +
    theme(legend.position = "none") +
    ggtitle("Number of counts")
  dbl_score[[so_name]] <- VlnPlot(so, features = "scDblFinder.score", group.by = "seurat_clusters", pt.size = 0.1) +
    theme(legend.position = "none") +
    ggtitle("scDblFinder score")
  
  quality_stats[[so_name]] <- (n_count[[so_name]] + n_feat[[so_name]]) | (mt[[so_name]] + dbl_score[[so_name]])
  
  if(stromal_soupx) {
    ggsave(file.path(out_folder_split, paste0(n, "_quality_stats_soupx.png")), quality_stats[[so_name]], height = 6, width = 12)
  } else {
    ggsave(file.path(out_folder_split, paste0(n, "_quality_stats.png")), quality_stats[[so_name]], height = 6, width = 12)
  }
  
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
  
  if (stromal_soupx) {
    write_tsv(top_phases, file = file.path(out_folder_split, paste0(n, "_cc_top_phases_soupx.tsv")))
    write_tsv(freq_phases, file = file.path(out_folder_split, paste0(n, "_cc_freq_phases_soupx.tsv")))
  } else {
    write_tsv(top_phases, file = file.path(out_folder_split, paste0(n, "_cc_top_phases.tsv")))
    write_tsv(freq_phases, file = file.path(out_folder_split, paste0(n, "_cc_freq_phases.tsv")))
  }
  
  # umap plots
  umap_lbl <- UMAPPlot(so, group.by = 'label', label = T)
  umap_phase <- UMAPPlot(so, group.by = 'Phase')
  umap_clust <- UMAPPlot(so, label = T)
  feature_dbl <- FeaturePlot(so, features = 'scDblFinder.score')
  feature_mt <- FeaturePlot(so, features = 'percent.mt')
  umap_dbl <- UMAPPlot(so, group.by = 'scDblFinder.class')
  
  if (stromal_soupx) {
    ggsave(umap_lbl, file = file.path(out_folder_split, paste0(n, "_umap_lbl_soupx.pdf")), height = 6, width = 8)
    ggsave(umap_phase, file = file.path(out_folder_split, paste0(n, "_umap_phase_soupx.pdf")), height = 6, width = 8)
    ggsave(umap_clust, file = file.path(out_folder_split, paste0(n, "_umap_clust_soupx.pdf")), height = 6, width = 8)
    ggsave(feature_dbl, file = file.path(out_folder_split, paste0(n, "_feature_dbl_soupx.pdf")), height = 6, width = 8)
    ggsave(feature_mt, file = file.path(out_folder_split, paste0(n, "_feature_mt_soupx.pdf")), height = 6, width = 8)
    ggsave(umap_dbl, file = file.path(out_folder_split, paste0(n, "_umap_dbl_soupx.pdf")), height = 6, width = 8)
  } else {
    ggsave(umap_lbl, file = file.path(out_folder_split, paste0(n, "_umap_lbl.pdf")), height = 6, width = 8)
    ggsave(umap_phase, file = file.path(out_folder_split, paste0(n, "_umap_phase.pdf")), height = 6, width = 8)
    ggsave(umap_clust, file = file.path(out_folder_split, paste0(n, "_umap_clust.pdf")), height = 6, width = 8)
    ggsave(feature_dbl, file = file.path(out_folder_split, paste0(n, "_feature_dbl.pdf")), height = 6, width = 8)
    ggsave(feature_mt, file = file.path(out_folder_split, paste0(n, "_feature_mt.pdf")), height = 6, width = 8)
    ggsave(umap_dbl, file = file.path(out_folder_split, paste0(n, "_umap_dbl.pdf")), height = 6, width = 8)
  }
  
  so.split.processed[[so_name]] <- so
  
}

# save so objects (intermediate and processed)
saveRDS(so.split.processed, file = file.path(out_folder_split, "so.split.processed.rds"))
saveRDS(so.split, file = file.path(out_folder_split, "so.split.rds"))

#~~~~~~~~~~~~~~~~~~~~
# overview UMAP plots
#####################
# using the combined object and UMAP representation with clusters identified in the separate compartments

meta.data_stroma <- so.split.processed$stromal@meta.data %>% mutate(comp_cluster = paste0("stroma_", seurat_clusters))
meta.data_cd45 <- so.split.processed$`cd45+`@meta.data %>% mutate(comp_cluster = paste0("cd45_", seurat_clusters))
meta.data_epcam <- so.split.processed$`epcam+`@meta.data %>% mutate(comp_cluster = paste0("epcam_", seurat_clusters))
meta.data <- rbind(meta.data_stroma, meta.data_cd45, meta.data_epcam) %>% select(cell_id, comp_cluster)

# load combined seurat object
if (stromal_soupx) {
  so_combined <- readRDS(file = file.path(out_folder_combined, "so_combined_soupx.rds"))
} else {
  so_combined <- readRDS(file = file.path(out_folder_combined, "so_combined.rds"))
}
rn <- so_combined@meta.data %>% rownames()
so_combined[["comp_cluster"]] <- factor(meta.data[rn,]$comp_cluster, 
                                        levels = c(paste("cd45", 1:21, sep = "_"), 
                                                   paste("epcam", 1:13, sep = "_"), 
                                                   paste("stroma", 1:22, sep = "_")))
UMAP_all <- UMAPPlot(so_combined, group.by = "split_group")
UMAP_clusters <- UMAPPlot(so_combined, group.by = "comp_cluster") 
UMAP_split <- UMAPPlot(so_combined, split.by = "split_group", group.by = "split_group")

ggsave(filename = file.path(out_folder_combined, "UMAP_all_cells.pdf"), width = 7, height = 6, plot = UMAP_all)
ggsave(filename = file.path(out_folder_combined, "UMAP_all_clusters.pdf"), width = 10, height = 6, plot = UMAP_clusters)
ggsave(filename = file.path(out_folder_combined, "UMAP_split_compartments.pdf"), width = 15, height = 6, plot = UMAP_split)

# UMAPs with each stromal and immune cluster highlighted once
stromal_clusters <- paste("stroma", 1:22, sep = "_")
immune_clusters <- paste("cd45", 1:21, sep = "_")

strom <- list()
Idents(so_combined) <- "comp_cluster"
for (clus in stromal_clusters) {
  cells <- WhichCells(so_combined, idents = clus)
  cells <- rownames(so_combined@meta.data[so_combined@meta.data$comp_cluster == clus,]) # get cell names belonging to respective cluster
  strom[[clus]] <- Cell_Highlight_Plot(seurat_object = so_combined, cells_highlight = list(cluster = cells), highlight_color = "blue") +
    theme(legend.position = "none") +
    ggtitle(clus)
  ggsave(filename = file.path(out_folder_combined, paste0("UMAP_", clus, ".pdf")), width = 6, height = 6, plot = strom[[clus]])
}

immune <- list()
for (clus in immune_clusters) {
  cells <- WhichCells(so_combined, idents = clus)
  cells <- rownames(so_combined@meta.data[so_combined@meta.data$comp_cluster == clus,]) # get cell names belonging to respective cluster
  immune[[clus]] <- Cell_Highlight_Plot(seurat_object = so_combined, cells_highlight = list(cluster = cells), highlight_color = "red") +
    theme(legend.position = "none") +
    ggtitle(clus)
  ggsave(filename = file.path(out_folder_combined, paste0("UMAP_", clus, ".pdf")), width = 6, height = 6, plot = immune[[clus]])
}
