# LIANA cell-cell communication analysis

library(tidyverse)
library(magrittr)
library(liana)
library(cowplot)
library(ggridges)
library(gdata)

#~~~~~~~~~~~~~~~~
# explore results
#################

so.split.processed <- readRDS("../data/so.split.processed.rds")
labels <- read.xls("../data/Table_S2_plus_colors.xlsx") %>%
  mutate(compartment = ifelse(compartment == "cd45-", "stromal", compartment)) %>%
  mutate(Cluster_number = factor(Cluster_number))

metadata <- rbind(so.split.processed[["stromal"]][[]],
      so.split.processed[["epcam+"]][[]],
      so.split.processed[["cd45+"]][[]]) %>%
  left_join(labels, by = c("split_group" = "compartment", "seurat_clusters" = "Cluster_number")) %>%
  mutate(sub_compartment = ifelse(Sub_compartment == "Cancer", "tumor", Sub_compartment))

liana_all_res <- readRDS("results/liana_res_orthologs.rds")

res <- rbind(liana_all_res$subcomp$young %>% mutate(age = "young"),
             liana_all_res$subcomp$old %>% mutate(age = "old")) %>%
  select(age, everything()) %>%
  filter(!source %in% c("Low-quality", "Contamination"),
         !target %in% c("Low-quality", "Contamination"))  %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  mutate(source = factor(source)) %>%
  mutate(target = factor(target)) %>%
  mutate(interaction_score = natmi.edge_specificity * sca.LRscore) # compute interaction score (specificity * magnitude)

write_tsv(res, file = "results/liana_res_combined_y_o.tsv")

#~~~~~~~~~~~~~~~
# overview plots
################

## Mean interaction score

stats <- res %>% 
  group_by(age, source, target) %>% 
  summarize(Mean = mean(interaction_score),
            Median = median(interaction_score),
            n_interactions = n())

cell_n <- metadata %>% 
  group_by(sub_compartment, age) %>% 
  summarize(n_cells = n()) %>%
  rbind(tibble(sub_compartment = "Lymphoid_NK", age = "old", n_cells = 0))

# right barplot
target_bar_mean_o <- cell_n %>%
  filter(sub_compartment %in% levels(stats$target),
         age == "old") %>%
  mutate(sub_compartment = factor(sub_compartment, levels = levels(stats$target))) %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  ggplot(aes(x = n_cells, y = sub_compartment)) +
  geom_bar(stat = "identity", fill = "#006D2C") +
  theme_classic() +
  facet_wrap(~age) +
  scale_x_continuous(breaks = c(0, 4000)) +
  scale_y_discrete(position = "right") +
  xlab("# cells") +
  ylab("") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(23, -10, -24, 5.5), "pt"))

# left barplot
target_bar_mean_y <- cell_n %>%
  filter(sub_compartment %in% levels(stats$target),
         age == "young") %>%
  mutate(sub_compartment = factor(sub_compartment, levels = levels(stats$target))) %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  ggplot(aes(x = n_cells, y = sub_compartment)) +
  geom_bar(stat = "identity", fill = "#006D2C") +
  theme_classic() +
  facet_wrap(~age) +
  xlab("# cells") +
  ylab("") +
  scale_x_reverse(breaks = c(0, 4000)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(23, 5.5, -24, -10), "pt"))

# bottom barplot
source_bar_mean <- cell_n %>%
  filter(sub_compartment %in% levels(stats$source)) %>%
  mutate(sub_compartment = factor(sub_compartment, levels = levels(stats$source))) %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  ggplot(aes(x = sub_compartment, y = n_cells)) +
  geom_bar(stat = "identity", fill = "#006D2C") +
  theme_classic() +
  facet_wrap(~age) +
  scale_y_reverse(breaks = c(0, 4000)) +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(5.5, -1, 5, -35), "pt"))

# main plot
level = "subcompartment"
color_palette = "Greens"
interaction_score = "Mean"
mean_plt <- ggplot(stats, aes(x = source, y = target)) +
  geom_point(aes(size = n_interactions, fill = !!sym(interaction_score)), shape = 21) +
  theme_classic() +
  xlab(paste0("Sender (", level, ")")) +
  ylab(paste0("Receiver (", level, ")")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_distiller(palette = color_palette, direction = 1) +
  scale_size(range = c(0.5, 5), name = "Number of interactions") +
  guides(fill = guide_colorbar(title = paste0(interaction_score, " interaction score"))) +
  facet_wrap(~age)

# main plot adapted
mean_plt_stripped <- mean_plt +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 0, 0, 0), "pt"))

# legend
legend <- cowplot::get_legend(mean_plt +
                      theme(legend.spacing.y = unit(-5, "pt"),
                            legend.box.margin = margin(5, 5, 5, 5),
                            legend.box.background = element_rect(colour = "black"),
                            legend.title = element_text(size = 10), 
                            legend.text = element_text(size = 10)))

mean_grid <- cowplot::plot_grid(target_bar_mean_y, mean_plt_stripped, target_bar_mean_o, 
                   NULL, source_bar_mean, NULL, 
                   rel_widths = c(0.68, 1, 0.68), rel_heights = c(0.82, 1), ncol = 3)

ggsave(filename = "results/liana_stats_subcomp_mean.pdf", plot = mean_plt, width = 7.5, height = 4)
ggsave(filename = "results/liana_stats_subcomp_mean_plotgrid.pdf", plot = mean_grid, width = 8.2, height = 4.2)
ggsave(filename = "results/liana_stats_subcomp_mean_plotgrid_legend.pdf", plot = legend, width = 4, height = 6)

## Median interaction score

stats <- res %>% 
  group_by(age, source, target) %>% 
  summarize(Mean = mean(interaction_score),
            Median = median(interaction_score),
            n_interactions = n())

cell_n <- metadata %>% 
  group_by(sub_compartment, age) %>% 
  summarize(n_cells = n()) %>%
  rbind(tibble(sub_compartment = "Lymphoid_NK", age = "old", n_cells = 0))

# right barplot
target_bar_median_o <- cell_n %>%
  filter(sub_compartment %in% levels(stats$target),
         age == "old") %>%
  mutate(sub_compartment = factor(sub_compartment, levels = levels(stats$target))) %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  ggplot(aes(x = n_cells, y = sub_compartment)) +
  geom_bar(stat = "identity", fill = "#54278F") +
  theme_classic() +
  facet_wrap(~age) +
  scale_x_continuous(breaks = c(0, 4000)) +
  scale_y_discrete(position = "right") +
  xlab("# cells") +
  ylab("") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(23, -10, -24, 5.5), "pt"))

# left barplot
target_bar_median_y <- cell_n %>%
  filter(sub_compartment %in% levels(stats$target),
         age == "young") %>%
  mutate(sub_compartment = factor(sub_compartment, levels = levels(stats$target))) %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  ggplot(aes(x = n_cells, y = sub_compartment)) +
  geom_bar(stat = "identity", fill = "#54278F") +
  theme_classic() +
  facet_wrap(~age) +
  xlab("# cells") +
  ylab("") +
  scale_x_reverse(breaks = c(0, 4000)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(23, 5.5, -24, -10), "pt"))

# bottom barplot
source_bar_median <- cell_n %>%
  filter(sub_compartment %in% levels(stats$source)) %>%
  mutate(sub_compartment = factor(sub_compartment, levels = levels(stats$source))) %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  ggplot(aes(x = sub_compartment, y = n_cells)) +
  geom_bar(stat = "identity", fill = "#54278F") +
  theme_classic() +
  facet_wrap(~age) +
  scale_y_reverse(breaks = c(0, 4000)) +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(5.5, -1, 5, -35), "pt"))

# main plot
level = "subcompartment"
color_palette = "Purples"
interaction_score = "Median"
median_plt <- ggplot(stats, aes(x = source, y = target)) +
  geom_point(aes(size = n_interactions, fill = !!sym(interaction_score)), shape = 21) +
  theme_classic() +
  xlab(paste0("Sender (", level, ")")) +
  ylab(paste0("Receiver (", level, ")")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_distiller(palette = color_palette, direction = 1) +
  scale_size(range = c(0.5, 5), name = "Number of interactions") +
  guides(fill = guide_colorbar(title = paste0(interaction_score, " interaction score"))) +
  facet_wrap(~age)

# main plot adapted
median_plt_stripped <- median_plt +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 0, 0, 0), "pt"))

# legend
legend <- cowplot::get_legend(median_plt +
                                theme(legend.spacing.y = unit(-5, "pt"),
                                      legend.box.margin = margin(5, 5, 5, 5),
                                      legend.box.background = element_rect(colour = "black"),
                                      legend.title = element_text(size = 10), 
                                      legend.text = element_text(size = 10)))

median_grid <- cowplot::plot_grid(target_bar_median_y, median_plt_stripped, target_bar_median_o, 
                                NULL, source_bar_median, NULL, 
                                rel_widths = c(0.68, 1, 0.68), rel_heights = c(0.82, 1), ncol = 3)

ggsave(filename = "results/liana_stats_subcomp_median.pdf", plot = median_plt, width = 7.5, height = 4)
ggsave(filename = "results/liana_stats_subcomp_median_plotgrid.pdf", plot = median_grid, width = 8.2, height = 4.2)
ggsave(filename = "results/liana_stats_subcomp_median_plotgrid_legend.pdf", plot = legend, width = 4, height = 6)

## distributions of magnitude and specificity

distr_plots <- list()
stat <- "sca.LRscore"
col = "#006D2C"
distr_plots[[stat]] <- ggplot(res, aes(x = !!sym(stat), y = age, fill = age)) +
  geom_density_ridges() +
  theme_bw() +
  theme(legend.position = "none")

stat <- "natmi.edge_specificity"
col = "#54278F"
distr_plots[[stat]] <- ggplot(res, aes(x = !!sym(stat), y = age, fill = age)) +
  geom_density_ridges() +
  theme_bw() +
  theme(legend.position = "none")

distr <- cowplot::plot_grid(plotlist = distr_plots, rel_widths = c(0.8, 1))
ggsave(filename = "results/liana_distr_subcomp.pdf", plot = distr, width = 6, height = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Jaccard index old vs. young
#############################
# results subfolder
if (!file.exists("results/similarity_o_y")) { dir.create("results/similarity_o_y") }

# select relevant columns and only keep source/target cell type pairs that are available for both, young and old
res_min <- res %>% 
  select(age, source, target, ligand.complex, receptor.complex, aggregate_rank, interaction_score) %>%
  mutate(source.target = paste(source, target, sep = ".")) %>%
  mutate(ligand.receptor = paste(ligand.complex, receptor.complex, sep = "."))

pairs_old <- res_min %>% filter(age == "old") %>% group_by(source.target) %>% summarize(count = n()) %>% pull(source.target)
pairs_young <- res_min %>% filter(age == "young") %>% group_by(source.target) %>% summarize(count = n()) %>% pull(source.target)

pairs <- intersect(pairs_old, pairs_young)

# compute JI across a range of rank cutoffs
agg_rank_grid <- seq(0.05, 1, 0.05)
ji_tibble <- tibble()
for (agg_rank_cut in agg_rank_grid) {
  
  for (pair in pairs) {
    
    r <- res_min %>% filter(source.target == pair) %>% filter(aggregate_rank <= agg_rank_cut)
    
    n_interactions <- nrow(r)
    
    young <- r %>% filter(age == "young")
    old <- r %>% filter(age == "old")
    
    ji <- jaccard_index(young$ligand.receptor, old$ligand.receptor)
    
    row <- tibble("agg_rank_cut" = agg_rank_cut, "pair" = pair, "jaccard_index" = ji, "n_interactions" = n_interactions)
    
    ji_tibble <- rbind(ji_tibble, row)
    
  }
  
}
ji_tibble <- ji_tibble %>% mutate(agg_rank_cut = factor(agg_rank_cut))

# plot results as dotplots
for (cutoff in agg_rank_grid) {
  
  plt <- ji_tibble %>%
    filter(agg_rank_cut == cutoff) %>%
    mutate(source = str_split_i(pair, "\\.", 1)) %>%
    mutate(target = str_split_i(pair, "\\.", 2)) %>%
    ggplot(aes(x = source, y = target)) +
      geom_point(aes(size = n_interactions, fill = jaccard_index), shape = 21) +
      theme_classic() +
      xlab("Sender") +
      ylab("Receiver") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.box.margin = margin(116, 6, 6, 6)) +
      scale_fill_distiller(palette = "Spectral", direction = -1,
                           limits = c(0, 1)) +
      scale_size(range = c(1, 8), name = "Number of interactions") +
      guides(fill = guide_colorbar(title = "Jaccard Index"))
  
  ggsave(plot = plt, filename = paste0("results/similarity_o_y/jaccard_old_young_subcomp_aggrank_cutoff_", cutoff, ".pdf"),
         width = 7, height = 4)
  
}

# plot 
ji_l <- ji_tibble %>% 
  drop_na(jaccard_index) %>% 
  mutate(ligand = str_split_i(pair, pattern = "\\.", i = 1),
         receptor = str_split_i(pair, pattern = "\\.", i = 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  group_by(agg_rank_cut, ligand) %>%
  summarize(mean_ji = mean(jaccard_index), median_ji = median(jaccard_index)) %>%
  pivot_longer(cols = mean_ji:median_ji, names_to = "stat", values_to = "value") %>%
  mutate(stat_base = str_split_i(stat, pattern = "_", 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  ggplot(aes(x = agg_rank_cut, y = value, color = ligand)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess", se = F) +
  ylab("Jaccard index") +
  xlab("") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~stat, ncol = 1) +
  ylim(c(0, 0.7)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.margin = unit(c(5, 4, 5, 5), "pt"))

n_l <- ji_tibble %>% 
  drop_na(jaccard_index) %>% 
  mutate(ligand = str_split_i(pair, pattern = "\\.", i = 1),
         receptor = str_split_i(pair, pattern = "\\.", i = 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  group_by(agg_rank_cut, ligand) %>%
  summarize(mean_n = mean(n_interactions), median_n = median(n_interactions)) %>%
  pivot_longer(cols = mean_n:median_n, names_to = "stat", values_to = "value") %>%
  mutate(stat_base = str_split_i(stat, pattern = "_", 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  ggplot(aes(x = agg_rank_cut, y = value, color = ligand)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess", se = F) +
  ylab("") +
  xlab("Aggregated rank cutoff") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~stat, ncol = 1) +
  ylim(c(0, 150)) +
  theme_bw() +
  theme(plot.margin = unit(c(5, 5, 5, 0), "pt"))

ji_r <- ji_tibble %>% 
  drop_na(jaccard_index) %>% 
  mutate(ligand = str_split_i(pair, pattern = "\\.", i = 1),
         receptor = str_split_i(pair, pattern = "\\.", i = 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  group_by(agg_rank_cut, receptor) %>%
  summarize(mean_ji = mean(jaccard_index), median_ji = median(jaccard_index)) %>%
  pivot_longer(cols = mean_ji:median_ji, names_to = "stat", values_to = "value") %>%
  mutate(stat_base = str_split_i(stat, pattern = "_", 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  ggplot(aes(x = agg_rank_cut, y = value, color = receptor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess", se = F) +
  ylab("Jaccard index") +
  xlab("") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~stat, ncol = 1) +
  ylim(c(0, 0.7)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.margin = unit(c(5, 4, 5, 5), "pt"))

n_r <- ji_tibble %>% 
  drop_na(jaccard_index) %>% 
  mutate(ligand = str_split_i(pair, pattern = "\\.", i = 1),
         receptor = str_split_i(pair, pattern = "\\.", i = 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  group_by(agg_rank_cut, receptor) %>%
  summarize(mean_n = mean(n_interactions), median_n = median(n_interactions)) %>%
  pivot_longer(cols = mean_n:median_n, names_to = "stat", values_to = "value") %>%
  mutate(stat_base = str_split_i(stat, pattern = "_", 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  ggplot(aes(x = agg_rank_cut, y = value, color = receptor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess", se = F) +
  ylab("") +
  xlab("Aggregated rank cutoff") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~stat, ncol = 1) +
  ylim(c(0, 150)) +
  theme_bw() +
  theme(plot.margin = unit(c(5, 5, 5, 0), "pt"))

ji_n_grid <- cowplot::plot_grid(ji_l, n_l, ji_r, n_r, nrow = 1, rel_widths = c(0.51, 1))

ggsave(plot = ji_n_grid, filename = "results/similarity_o_y/ji_rank_ligands.receptors.pdf",
       width = 13, height = 4)


# overview of jaccard index statistics across rank cutoffs
ji_by_cut_ligands <- ji_tibble %>% 
  drop_na(jaccard_index) %>% 
  mutate(ligand = str_split_i(pair, pattern = "\\.", i = 1),
         receptor = str_split_i(pair, pattern = "\\.", i = 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  group_by(agg_rank_cut, ligand) %>%
  summarize(mean_ji = mean(jaccard_index), median_ji = median(jaccard_index), 
            mean_n = mean(n_interactions), median_n = median(n_interactions)) %>%
  pivot_longer(cols = mean_ji:median_n, names_to = "stat", values_to = "value") %>%
  mutate(stat_base = str_split_i(stat, pattern = "_", 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  ggplot(aes(x = agg_rank_cut, y = value, color = ligand)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess", se = F) +
  ylab("Jaccard index") +
  xlab("Aggregated rank cutoff") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~stat, scales = "free") +
  theme_bw()

ji_by_cut_receptors <- ji_tibble %>% 
  drop_na(jaccard_index) %>% 
  mutate(ligand = str_split_i(pair, pattern = "\\.", i = 1),
         receptor = str_split_i(pair, pattern = "\\.", i = 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  group_by(agg_rank_cut, receptor) %>%
  summarize(mean_ji = mean(jaccard_index), median_ji = median(jaccard_index), 
            mean_n = mean(n_interactions), median_n = median(n_interactions)) %>%
  pivot_longer(cols = mean_ji:median_n, names_to = "stat", values_to = "value") %>%
  mutate(stat_base = str_split_i(stat, pattern = "_", 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut))) %>%
  ggplot(aes(x = agg_rank_cut, y = value, color = receptor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess", se = F) +
  ylab("Jaccard index") +
  xlab("Aggregated rank cutoff") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~stat, scales = "free") +
  theme_bw()

ji_rank_ligands_receptors <- cowplot::plot_grid(ji_by_cut_ligands, ji_by_cut_receptors, ncol = 2)

ggsave(plot = ji_rank_ligands_receptors, filename = "results/similarity_o_y/ji_rank_ligands.receptors.pdf",
       width = 12, height = 4)

##

ji_summary <- ji_tibble %>% 
  drop_na(jaccard_index) %>% 
  group_by(agg_rank_cut) %>% 
  summarize(mean_ji = mean(jaccard_index), median_ji = median(jaccard_index), 
            mean_n = mean(n_interactions), median_n = median(n_interactions)) %>%
  pivot_longer(cols = mean_ji:median_n, names_to = "stat", values_to = "value") %>%
  mutate(stat_base = str_split_i(stat, pattern = "_", 2)) %>%
  mutate(agg_rank_cut = as.double(as.character(agg_rank_cut)))

ji_by_cut <- ggplot(ji_summary, aes(x = agg_rank_cut, y = value, color = stat)) +
  geom_line(stat = "identity", position = "identity") +
  facet_wrap(~stat_base, scales = "free") +
  theme_bw() +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 12), "pt"))

agg_rank_by_age_ridge <- ggplot(res, aes(x = aggregate_rank, y = age, fill = age)) +
  geom_density_ridges() +
  theme_bw() +
  theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"))

agg_rank_by_age_hist <- ggplot(res, aes(x = aggregate_rank)) +
  geom_histogram(color = "white", fill = "skyblue3") +
  theme_bw() +
  theme(plot.margin = unit(c(5.5, 90, 5.5, 10), "pt"))

ji_stat_plt <- cowplot::plot_grid(ji_by_cut,
                   agg_rank_by_age_ridge,
                   agg_rank_by_age_hist, ncol = 1)

ggsave(plot = ji_stat_plt, filename = "results/similarity_o_y/ji_stats.pdf", width = 6, height = 5)

#~~~~~~~~~~~~~~
# Dot plots o/y
###############
# results subfolder
if (!file.exists("results/old_young_dotplots")) {dir.create("results/old_young_dotplots")}

dp_res <- res %>% mutate(ligand.receptor = paste(ligand.complex, receptor.complex, sep = "."))

n_top <- 20

for (source_filter in unique(dp_res$source)) {
  
  subs <- dp_res %>% filter(source == source_filter) %>% 
    slice_min(n = 30, order_by = aggregate_rank)
  pairs <- subs %>% pull(ligand.receptor) %>% unique()
  pairs <- pairs[1:n_top]
  
  plt <- liana_custom_dotplot(subs %>% 
                         filter(ligand.receptor %in% pairs), 
                       source_groups = !!source_filter, 
                       plot.title = paste0("Source: ", source_filter))
  
  ggsave(plot = plt, filename = paste0("results/old_young_dotplots/pairs_", source_filter, ".pdf"),
         width = 10, height = 12)
  
}

########################