# Testing whether the difference in proportions (young/old, treated/untreated) of cells in each cluster
# are significant.
# For immune compartment also testing for differences in the composition of clusters with regard to treatment, separately for old and young mice
# Dotplot for expression of Tcf7 and exhaustion markers Havcr2 and Entpd1 in young and old mice separately

library(scProportionTest)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(Seurat)
library(readxl)
library(data.table)
source("TME.ds_functions.R")

log2FD_threshold = log2(1.5)

out.dir <- file.path("combined_analysis", "split_objects")
if(!dir.exists(out.dir)) dir.create(out.dir)
annotation_file <- "combined_analysis/split_objects/results_20230925_metadata_final_v2.xlsx"

so.split.processed <- readRDS("combined_analysis/split_objects/so.split.processed.rds")
conv <- read_tsv("label_condition_conversion.txt")
cluster_anno <- read_xlsx(annotation_file) %>%
  mutate(compartment = ifelse(compartment == "cd45-", "stromal", compartment))
color_table <- read_tsv("plotting_colors.tsv")

# update metadata slots
for (so_name in names(so.split.processed)) {
  ca <- cluster_anno %>% filter(compartment == so_name) %>%
    mutate(cluster_label = paste(Cell_Type_Label, cluster, sep = ".")) %>%
    mutate(cluster = factor(cluster))
  
  so.split.processed[[so_name]]$cluster_label <- so.split.processed[[so_name]]@meta.data %>%
    left_join(ca, by = c("seurat_clusters" = "cluster")) %>%
    pull(cluster_label)
  
  so.split.processed[[so_name]]$age_tm <- so.split.processed[[so_name]]@meta.data %>%
    mutate(age_tm = paste(age, treatment, sep = " ")) %>%
    pull(age_tm)
  
}

# loop over cell compartments, two plots each (y/o, t/ut) and for each of them
# a barplot orderd according to the cluster order in the result
for (so_name in names(so.split.processed)) {
  so <- so.split.processed[[so_name]]
  
  ct <- color_table %>% filter(compartment == so_name)
  col_vec <- ct$color
  names(col_vec) <- ct$condition
  
  proptest_age <- sc_utils(so)
  
  proptest_age <- permutation_test(
    proptest_age, cluster_identity = "cluster_label",
    sample_1 = "young", sample_2 = "old",
    sample_identity = "age"
  )
  
  plot_data <- return_permutation_plotdata(proptest_age)
  cl_order_age <- plot_data %>% arrange(desc(obs_log2FD)) %>% pull(clusters)
  
  prop_table_age <- plot_data %>% arrange(desc(obs_log2FD))
  write_tsv(x = prop_table_age, file = file.path(out.dir, paste(so_name, "cluster_composition_age.tsv", sep = "_")))
  
  perm_plot_age <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) + 
    theme_bw() + 
    geom_hline(yintercept = log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(values = c(ct %>% 
                                    filter(condition == "young untreated") %>% 
                                    pull(color), 
                                  "grey")) + 
    coord_flip() +
    ylim(-10, 10)
  
  n_cells <- so@meta.data %>% 
    left_join(conv, by = c("label" = "label")) %>%
    group_by(cluster_label, age_tm) %>% summarize(n_cells = n()) %>%
    left_join(ct, by = c("age_tm" = "condition"))
  
  n_cells_age <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_age))) + 
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    scale_fill_manual(values = col_vec) +
    ylab("number of cells") +
    labs(fill = "") +
    coord_flip()
  
  perc_age <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_age))) + 
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    scale_fill_manual(values = col_vec) +
    ylab("percentage of cells") +
    coord_flip()
  
  summary_age <- (perm_plot_age | perc_age | n_cells_age) +
    plot_annotation(title = paste0(so_name, " cells cluster composition: old vs. young"),
                    theme = theme(plot.title = element_text(size = 16, 
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  ggsave(filename = file.path(out.dir, paste(so_name, "cluster_composition_age.pdf", sep = "_")),
         summary_age,
         height = 6,
         width = 12)
  
  if(so_name == "cd45+") {
    pd <- as_tibble(plot_data)
    
    ca <- cluster_anno %>% filter(compartment == so_name) %>%
      mutate(cluster_label = paste(Cell_Type_Label, cluster, sep = ".")) %>%
      mutate(cluster = factor(cluster))
    
    pd <- pd %>% left_join(ca, by = c("clusters" = "cluster_label"))
    
    pd_filt <- pd %>% filter(obs_log2FD != -Inf & obs_log2FD != Inf)
    pd <- pd %>% mutate(obs_log2FD = ifelse(obs_log2FD == Inf, 
                                            max(pd_filt$obs_log2FD), 
                                            ifelse(obs_log2FD == -Inf, 
                                                   min(pd_filt$obs_log2FD), obs_log2FD)))
    
    pd <- pd %>% filter(sub_compartment != "Low-quality")
    subcomp_plot <- ggplot(pd, aes(x = reorder(sub_compartment, obs_log2FD, median, decreasing = T), y = obs_log2FD)) + 
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 1, width = 0.2) +
      coord_flip() +
      theme_classic() +
      xlab("Sub-compartment") +
      geom_hline(yintercept = 0, linetype = 2)
    
    ggsave(filename = file.path(out.dir, paste(so_name, "cluster_composition_by_subcom_age.pdf", sep = "_")),
           subcomp_plot,
           height = 2,
           width = 6)
  }
  
  
  # treated/untreated
  proptest_tm <- sc_utils(so)
  
  proptest_tm <- permutation_test(
    proptest_tm, cluster_identity = "cluster_label",
    sample_1 = "untreated", sample_2 = "triple_treated",
    sample_identity = "treatment"
  )
  
  plot_data <- return_permutation_plotdata(proptest_tm)
  cl_order_tm <- plot_data %>% arrange(desc(obs_log2FD)) %>% pull(clusters)
  
  prop_table_tm <- plot_data %>% arrange(desc(obs_log2FD))
  write_tsv(x = prop_table_tm, file = file.path(out.dir, paste(so_name, "cluster_composition_tm.tsv", sep = "_")))
  
  perm_plot_tm <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) + 
    theme_bw() + 
    geom_hline(yintercept = log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(values = c(ct %>% 
                                    filter(condition == "young untreated") %>% 
                                    pull(color), 
                                  "grey")) + 
    coord_flip() +
    ylim(-10, 10)
  
  n_cells <- so@meta.data %>% 
    left_join(conv, by = c("label" = "label")) %>%
    group_by(cluster_label, age_tm) %>% summarize(n_cells = n()) %>%
    left_join(ct, by = c("age_tm" = "condition")) %>%
    filter(cluster_label %in% cl_order_tm)
  
  n_cells_tm <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_tm))) + 
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
    scale_fill_manual(values = col_vec) +
    ylab("number of cells") +
    labs(fill = "") +
    coord_flip()
  
  perc_tm <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_tm))) + 
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    scale_fill_manual(values = col_vec) +
    ylab("percentage of cells") +
    coord_flip()
  
  summary_tm <- (perm_plot_tm | perc_tm | n_cells_tm) +
    plot_annotation(title = paste0(so_name, " cells cluster composition: treated vs. untreated"),
                    theme = theme(plot.title = element_text(size = 16, 
                                                            face = "bold",
                                                            hjust = 0.5)))
  
  ggsave(filename = file.path(out.dir, paste(so_name, "cluster_composition_treatment.pdf", sep = "_")),
         summary_tm,
         height = 6,
         width = 12)
}

# subsetted permutation plot for stromal cells
so_name <- "stromal"
so <- so.split.processed[[so_name]]

ct <- color_table %>% filter(compartment == so_name)
col_vec <- ct$color
names(col_vec) <- ct$condition

proptest_age <- sc_utils(so)

proptest_age <- permutation_test(
  proptest_age, cluster_identity = "cluster_label",
  sample_1 = "young", sample_2 = "old",
  sample_identity = "age"
)

plot_data <- return_permutation_plotdata(proptest_age)

prop_table_age <- plot_data %>% arrange(desc(obs_log2FD))
write_tsv(x = prop_table_age, file = file.path(out.dir, paste(so_name, "subset_cluster_composition_age.tsv", sep = "_")))

stromal_clus_sel <- 
  plot_data$clusters %>% 
  str_subset(pattern = "^Unclear.*", negate = T) %>% 
  str_subset(pattern = "^Contaminated.*", negate = T)

perm_plot_age <- ggplot(plot_data[clusters %in% stromal_clus_sel], aes(x = clusters, y = obs_log2FD)) +
  geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) + 
  theme_bw() + 
  geom_hline(yintercept = log2FD_threshold, lty = 2) + 
  geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = c(ct %>% 
                                  filter(condition == "young untreated") %>% 
                                  pull(color), 
                                "grey")) + 
  coord_flip() +
  ylim(-10, 10)

cl_order_age <- plot_data %>% arrange(desc(obs_log2FD)) %>% pull(clusters)

n_cells <- so@meta.data %>% 
  left_join(conv, by = c("label" = "label")) %>%
  group_by(cluster_label, age_tm) %>% summarize(n_cells = n()) %>%
  left_join(ct, by = c("age_tm" = "condition")) %>%
  filter(cluster_label %in% stromal_clus_sel)

n_cells_age <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_age))) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  scale_fill_manual(values = col_vec) +
  ylab("number of cells") +
  labs(fill = "") +
  coord_flip()

perc_age <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_age))) + 
  geom_bar(position = "fill", stat = "identity") +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  scale_fill_manual(values = col_vec) +
  ylab("percentage of cells") +
  coord_flip()

summary_age <- (perm_plot_age | perc_age | n_cells_age) +
  plot_annotation(title = paste0(so_name, " cells cluster composition: old vs. young"),
                  theme = theme(plot.title = element_text(size = 16, 
                                                          face = "bold",
                                                          hjust = 0.5)))

ggsave(filename = file.path(out.dir, paste(so_name, "subset_cluster_composition_age.pdf", sep = "_")),
       summary_age,
       height = 5,
       width = 12)

# treated/untreated
proptest_tm <- sc_utils(so)

proptest_tm <- permutation_test(
  proptest_tm, cluster_identity = "cluster_label",
  sample_1 = "untreated", sample_2 = "triple_treated",
  sample_identity = "treatment"
)

plot_data <- return_permutation_plotdata(proptest_tm)

prop_table_tm <- plot_data %>% arrange(desc(obs_log2FD))
write_tsv(x = prop_table_tm, file = file.path(out.dir, paste(so_name, "subset_cluster_composition_tm.tsv", sep = "_")))

perm_plot_tm <- ggplot(plot_data[clusters %in% stromal_clus_sel], aes(x = clusters, y = obs_log2FD)) +
  geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) + 
  theme_bw() + 
  geom_hline(yintercept = log2FD_threshold, lty = 2) + 
  geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = c(ct %>% 
                                  filter(condition == "young untreated") %>% 
                                  pull(color), 
                                "grey")) + 
  coord_flip() +
  ylim(-10, 10)

cl_order_tm <- perm_plot_tm$data %>% arrange(desc(obs_log2FD)) %>% pull(clusters)

n_cells <- so@meta.data %>% 
  left_join(conv, by = c("label" = "label")) %>%
  group_by(cluster_label, age_tm) %>% summarize(n_cells = n()) %>%
  left_join(ct, by = c("age_tm" = "condition")) %>%
  filter(cluster_label %in% stromal_clus_sel)

n_cells_tm <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_tm))) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  scale_fill_manual(values = col_vec) +
  ylab("number of cells") +
  labs(fill = "") +
  coord_flip()

perc_tm <- ggplot(n_cells, aes(fill = age_tm, y = n_cells, x = factor(cluster_label, levels = cl_order_tm))) + 
  geom_bar(position = "fill", stat = "identity") +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  scale_fill_manual(values = col_vec) +
  ylab("percentage of cells") +
  coord_flip()

summary_tm <- (perm_plot_tm | perc_tm | n_cells_tm) +
  plot_annotation(title = paste0(so_name, " cells cluster composition: treated vs. untreated"),
                  theme = theme(plot.title = element_text(size = 16, 
                                                          face = "bold",
                                                          hjust = 0.5)))

ggsave(filename = file.path(out.dir, paste(so_name, "subset_cluster_composition_treatment.pdf", sep = "_")),
       summary_tm,
       height = 5,
       width = 12)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tcf7 and exhaustion marker expression
#######################################
annotation <- read_xlsx(annotation_file) %>%
  mutate(cluster_label = paste(Cell_Type_Label, cluster, sep = "."))
annotation_cd45 <- annotation %>% filter(compartment == "cd45+")
annotation_cd45_tcells <- annotation %>% filter(compartment == "cd45+", sub_compartment == "Lymphoid_T")

immune <- so.split.processed[["cd45+"]]
immune$sub_compartment <- immune[[]] %>% left_join(annotation_cd45, by = c("cluster_label" = "cluster_label")) %>% pull(sub_compartment)
immune.ages <- immune %>% SplitObject(split.by = "age")

Idents(immune.ages$young) <- "cluster_label"
Idents(immune.ages$old) <- "cluster_label"

young_exh <- DotPlot(immune.ages$young, features = c("Tcf7", "Havcr2", "Entpd1"), 
                     assay = "SCT", cluster.idents = TRUE, idents = annotation_cd45_tcells$cluster_label) +
  ggtitle("Young T-cells")

old_exh <- DotPlot(immune.ages$old, features = c("Tcf7", "Havcr2", "Entpd1"), 
                   assay = "SCT", cluster.idents = TRUE, idents = annotation_cd45_tcells$cluster_label) +
  ggtitle("Old T-cells")

exh_pw <- (young_exh | old_exh) +
  plot_annotation(title = "CD45+ T-cell marker expression")
ggsave(filename = file.path(out.dir, "ds_analyses/immune/", "exhaustion_marker_expr.pdf"), plot = exh_pw, width = 14, height = 6)

# Feature plots
exh_fp <- FeaturePlot(immune, features = c("Tcf7", "Havcr2", "Entpd1"), order = TRUE, ncol = 3)
ggsave(filename = file.path(out.dir, "ds_analyses/immune/", "exhaustion_marker_featureplot.pdf"), plot = exh_fp, width = 18, height = 5)

# Violin plots
old_exh_vln <- immune.ages$old %>% VlnPlot(features = c("Tcf7", "Havcr2", "Entpd1")) +
  plot_annotation(title = "Old T-cells")

young_exh_vln <- immune.ages$young %>% VlnPlot(features = c("Tcf7", "Havcr2", "Entpd1")) +
  plot_annotation(title = "Young T-cells")

vln_pw <- old_exh_vln / young_exh_vln

Idents(immune) <- "cluster_label"
tcf7_expr <- annotation_cd45[which(annotation_cd45$cluster %in% c(2, 5 ,7 , 12, 13, 16, 19)),] %>% arrange(cluster_label)
vln_tcf7_tcells <- VlnPlot(immune, features = "Tcf7", idents = tcf7_expr$cluster_label, 
                           split.plot = TRUE, split.by = "age", group.by = "cluster_label") +
  stat_summary(fun = median, geom = 'point', size = 3, aes(color = split))
ggsave(filename = file.path(out.dir, "ds_analyses/immune/", "tcf7_tcells_age.pdf"), plot = vln_tcf7_tcells, width = 7, height = 4)

# Dot plots
dotplot_full <- DotPlot(immune, features = c("Tcf7", "Havcr2", "Entpd1"),
                        assay = "SCT", cluster.idents = TRUE) +
  ggtitle("Exhaustion marker expression in immune compartment")
ggsave(filename = file.path(out.dir, "ds_analyses/immune/", "exhaustion_dotplot_all_clusters.pdf"), plot = dotplot_full, width = 9, height = 6)

Idents(immune) <- "sub_compartment"
dotplot_subcomp <- custom_dotplot(so = immune, 
               features = c("Tcf7", "Havcr2", "Entpd1"),
               assay = "SCT",
               color_pal = "Blues",
               cluster.idents = TRUE)

ggsave(filename = file.path(out.dir, "ds_analyses/immune/", "exhaustion_markers_subcompartments.pdf"), plot = dotplot_subcomp, width = 4, height = 4)

dot_split_age <- custom_dotplot(so = immune, 
                                features = c("Tcf7", "Havcr2", "Entpd1"),
                                assay = "SCT",
                                color_pal = "Blues",
                                cluster.idents = TRUE,
                                split.by = "age")
dot_split_tm <- custom_dotplot(so = immune, 
                               features = c("Tcf7", "Havcr2", "Entpd1"),
                               assay = "SCT",
                               color_pal = "Blues",
                               cluster.idents = FALSE,
                               split.by = "treatment")

ggsave(filename = file.path(out.dir, "ds_analyses/immune/", "exhaustion_markers_subcompartments_age.pdf"), plot = dot_split_age, width = 5, height = 4)
ggsave(filename = file.path(out.dir, "ds_analyses/immune/", "exhaustion_markers_subcompartments_treatment.pdf"), plot = dot_split_tm, width = 5, height = 4)

subcomp_colors <- annotation_cd45 %>% group_by(sub_compartment) %>% summarize(Color_Scheme_1_subcomp = first(Color_Scheme_1))
subcomp_color_vec <- subcomp_colors$Color_Scheme_1_subcomp
names(subcomp_color_vec) <- subcomp_colors$sub_compartment

# Violin plots of exhaustion marker expression in subcompartments
for (gene in c("Tcf7", "Havcr2", "Entpd1")) {
  # Violin plots by subcompartment/age
  vln1 <- VlnPlot(immune, features = gene , assay = "SCT") +
    scale_fill_manual(values = subcomp_color_vec) +
    stat_summary(fun = median, geom = 'point', size = 3, color = "black")
  
  vln2 <- VlnPlot(immune, features = gene, assay = "SCT", split.plot = TRUE, split.by = "age", group.by = "sub_compartment") +
    stat_summary(fun = median, geom = 'point', size = 3, aes(color = split))
  
  vln3 <- ggplot(vln2$data, aes(x = split, y = .data[[gene]], fill = split)) +
    geom_violin() +
    geom_jitter(width = 0.3, size = 0.5) +
    theme_classic() +
    facet_wrap(~ident) +
    ggtitle(gene) +
    labs(fill = "Age") +
    xlab("")

  ggsave(filename = file.path(out.dir, "ds_analyses/immune", paste0(gene, "_by_subcompartment.pdf")), plot = vln1, width = 6, height = 4)
  ggsave(filename = file.path(out.dir, "ds_analyses/immune", paste0(gene, "_by_subcompartment_and_age.pdf")), plot = vln2, width = 6, height = 4)
  ggsave(filename = file.path(out.dir, "ds_analyses/immune", paste0(gene, "_by_subcompartment_and_age_facet.pdf")), plot = vln3, width = 6, height = 4)
  
}

# Violin plots by subcompartment and statistical test by age/treatment using wilcoxon test
VlnPlot(immune, features = "Tcf7", split.by = "age", group.by = "sub_compartment", split.plot = T)

# table of stats (p-values, log-fold-changes, and fraction of cells expressing the gene)
count <- 0
for (gene_test in c("Tcf7", "Havcr2", "Entpd1")) {
  compartments <- immune[[]]$sub_compartment %>% unique()
  
  for (compartment in compartments) {
    subs <- subset(x = immune, 
                   subset = sub_compartment == compartment)
    
    Idents(subs) <- "age"
    dim_a <- subs[[]] %>% group_by(age) %>% summarize(count = n()) %>% dim()
    if(dim_a[1] == 2) {
      tmp_age <- FindMarkers(subs, ident.1 = "old", ident.2 = "young", logfc.threshold = 0, min.pct = 0) %>% 
        mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
        rownames_to_column("gene") %>%
        filter(gene == gene_test) %>%
        mutate(comparison = "old_vs_young", sub_compartment = compartment)
    }
    
    Idents(subs) <- "treatment"
    dim_t <- subs[[]] %>% group_by(treatment) %>% summarize(count = n()) %>% dim()
    if(dim_t[1] == 2) {
      tmp_tm <- FindMarkers(subs, ident.1 = "triple_treated", ident.2 = "untreated", logfc.threshold = 0, min.pct = 0) %>%  
        mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2) %>%
        rownames_to_column("gene") %>%
        filter(gene == gene_test) %>% 
        mutate(comparison = "treated_vs_untreated", sub_compartment = compartment)
    }
    
    if(count == 0) {
      tbl <- rbind(tmp_age, tmp_tm)
    } else {
      tbl <- rbind(tbl, tmp_age, tmp_tm)
    }
    count = count + 1
  }
}

write_tsv(tbl, file = "combined_analysis/split_objects/ds_analyses/immune/exhaustion_markers_expression_by_sub_compartment.tsv")






