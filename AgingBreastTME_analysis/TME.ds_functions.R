###################################################################################
## Function and library collection for TME downstream analysis and visualization ##
###################################################################################

#~~~~~~~~~~
# libraries
###########

library(tidyverse)
library(Seurat)
library(clustree)
library(patchwork)
library(RColorBrewer)
library(scProportionTest)
library(viridis)
library(readxl)
library(data.table)
library(gprofiler2)
library(clusterCrit)
library(scCustomize)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#scType cell-type annotation - Preparation
##########################################

# initialize db_list
db_list <- list()

#1. Full scType-DB plus Bach (any tissue)
#########################################

source("../data/scType/gene_sets_prepare_multi_tissues.mouse.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
markers_combined_path <- "../data/scType/ScTypeDB_full_plus_bach_mm.txt"

# prepare genesets, allowing to map to any of the tissues
gs_list = gene_sets_prepare_multi_bach(path_to_db_file = markers_combined_path, cell_type = NULL)

# add to db_list
db_list[["full_sctype_bach_any"]] <- gs_list

#2. Only Immune system
######################

# prepare genesets, allowing to map to any of the tissues
gs_list = gene_sets_prepare_multi_bach(path_to_db_file = markers_combined_path, cell_type = "Immune system")

# add to db_list
db_list[["Immune_system"]] <- gs_list

#3. Immune system + Mammary Gland
#################################

# prepare genesets, allowing to map to any of the tissues
gs_list = gene_sets_prepare_multi_bach(path_to_db_file = markers_combined_path, cell_type = c("Immune system", "Mammary Gland"))

# add to db_list
db_list[["Immune_system_Mammary"]] <- gs_list

#4. Function for cluster - celltype mapping
###########################################

cl_to_type_top10 <- function(cl, scores, so) {
  cl_scores <- scores[,rownames(so@meta.data[so@meta.data$seurat_clusters == cl,])]
  cl_scores_per_type <- rowSums(cl_scores) %>% sort(decreasing = TRUE)
  summary <- tibble(cluster = cl, 
                    type = names(cl_scores_per_type), 
                    score = cl_scores_per_type, 
                    ncells = sum(so@meta.data$seurat_clusters == cl))
  return(head(summary, 10))
}

#~~~~~~~~~~~~~~~
# custom dotplot
################
# function to create custom dotplot using the plot data from seurat DotPlot output

custom_dotplot <- function(so, features, color_pal, cluster.idents = FALSE, assay = "SCT", split.by = NULL) {
  
  dp <- DotPlot(so, features = features, cluster.idents = cluster.idents, split.by = split.by)
  
  plt <- ggplot(dp$data, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
    scale_fill_distiller(palette = color_pal, direction = 1) +
    scale_size(range = c(0.5, 5), name = "Percent Expressed") +
    guides(fill = guide_colorbar(title = "Average Expression"))
  
  return(plt)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# custom permutation plot fct
#############################
# adaptation of permutation_plot function from scProportionTest package to return plot_data

return_permutation_plotdata <- function (sc_utils_obj, FDR_threshold = 0.05, log2FD_threshold = log2(1.5), 
                                         order_clusters = TRUE) 
{
  plot_data <- copy(sc_utils_obj@results$permutation)
  plot_data[, `:=`(significance, ifelse(FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold, paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s."))]
  plot_data[, `:=`(significance, factor(significance, levels = c(paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s.")))]
  if (order_clusters) {
    plot_data[, `:=`(clusters, fct_reorder(factor(clusters), desc(obs_log2FD)))]
  }
  # p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) + 
  #   geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, 
  #                       color = significance)) + theme_bw() + geom_hline(yintercept = log2FD_threshold, 
  #                                                                        lty = 2) + geom_hline(yintercept = -log2FD_threshold, 
  #                                                                                              lty = 2) + geom_hline(yintercept = 0) + scale_color_manual(values = c("salmon", 
  #                                                                                                                                                                    "grey")) + coord_flip()
  return(plot_data)
}