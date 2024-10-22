library(tidyverse)
library(Seurat)
library(liana)
library(magrittr)
library(readxl)
library(OmnipathR)


dataFolder <- "../data/"
resultsFolder <- "results/"

min_cells <- 10
n_cores <- 8

ortholog_resource <- readRDS(file = paste0(resultsFolder, "ortholog_resource.rds"))
so.list <- readRDS(paste0(dataFolder, "so.split.processed.rds"))

annotation_file <- paste0(dataFolder, "Table_S2_plus_colors.xlsx")

cluster_anno <- read_xlsx(annotation_file) %>%
  mutate(compartment = ifelse(compartment == "cd45-", "stromal", compartment)) %>% 
  mutate(cluster_label = paste(`Cell-type_label`, Cluster_number, sep = "."))

# update metadata slots
for (so_name in names(so.list)) {
  ca <- cluster_anno %>% filter(compartment == so_name) %>%
    mutate(cluster_label = paste(`Cell-type_label`, Cluster_number, sep = ".")) %>%
    mutate(Cluster_number = factor(Cluster_number))
  
  so.list[[so_name]]$cluster_label <- so.list[[so_name]]@meta.data %>%
    left_join(ca, by = c("seurat_clusters" = "Cluster_number")) %>%
    pull(cluster_label)
  
  so.list[[so_name]]$age_tm <- so.list[[so_name]]@meta.data %>%
    mutate(age_tm = paste(age, treatment, sep = " ")) %>%
    pull(age_tm)
  
  so.list[[so_name]]$Sub_compartment <- so.list[[so_name]]@meta.data %>%
    left_join(ca, by = c("seurat_clusters" = "Cluster_number")) %>%
    pull(Sub_compartment)
  
  so.list[[so_name]]$`Cell-type_label` <- so.list[[so_name]]@meta.data %>%
    left_join(ca, by = c("seurat_clusters" = "Cluster_number")) %>%
    pull(`Cell-type_label`)
  
  Idents(so.list[[so_name]]) <- "cluster_label"
  so.list[[so_name]]@active.ident <- factor(so.list[[so_name]]@active.ident, levels = ca$cluster_label)
  
}

# create combined seurat object for analysis
# untreated, run old and young separately: tumor as one cluster/group 

## combine all objects
so_full <- merge(x = so.list[["stromal"]], y = so.list[2:3])

## subset untreated
so <- subset(x = so_full, subset = treatment == "untreated")
so_young <- subset(x = so, subset = age == "young")
so_old <- subset(x = so, subset = age == "old")

so.list <- list("young" = so_young, "old" = so_old)

rm(so_full, so_young, so_old)
gc()

#~~~~~~~~~~~~~~
# cluster level
###############

# create interaction groups
## all cancer cells in one group
for (age in c("young", "old")) {
  
  # define groups
  so.list[[age]][["interaction_group"]] <- so.list[[age]]@meta.data %>% 
    mutate(interaction_group = ifelse(split_group == "epcam+", "tumor", cluster_label)) %>%
    pull(interaction_group)
  
}

#~~~~~~
# LIANA
#######
liana_res_clusters <- list()
liana_res_clusters_sub <- list()
for (age in c("young", "old")) {
  
  liana_res_clusters[[age]] <- liana_wrap(so.list[[age]],
                                          resource = 'custom', # resource has to be set to 'custom' to work with external resources
                                          external_resource = ortholog_resource, # provide orthologous resource
                                          idents_col = "interaction_group", 
                                          min_cells = min_cells, 
                                          parallelize = TRUE, 
                                          workers = n_cores) %>% 
    liana_aggregate()
  
}

#~~~~~~~~~~~~~~~~~~~~~~
# sub-compartment level
#######################

# create interaction groups
## all cancer cells in one group
for (age in c("young", "old")) {
  
  # define groups
  so.list[[age]][["interaction_group"]] <- so.list[[age]]@meta.data %>% 
    mutate(interaction_group = ifelse(split_group == "epcam+", "tumor", sub_compartment)) %>%
    pull(interaction_group)
  
}

#~~~~~~
# LIANA
#######
liana_res_subcomp <- list()

for (age in c("young", "old")) {
  
  liana_res_subcomp[[age]] <- liana_wrap(so.list[[age]],
                                         resource = 'custom', # resource has to be set to 'custom' to work with external resources
                                         external_resource = ortholog_resource, # provide orthologous resource
                                         idents_col = "interaction_group", 
                                         min_cells = min_cells, 
                                         parallelize = TRUE, 
                                         workers = n_cores) %>% 
    liana_aggregate()
  
}

#~~~~~~~~~~~~~~~~
# cell type level
#################

# create interaction groups
## all cancer cells in one group
for (age in c("young", "old")) {
  
  # define groups
  so.list[[age]][["interaction_group"]] <- so.list[[age]]@meta.data %>% 
    mutate(interaction_group = ifelse(split_group == "epcam+", "tumor", Cell_Type_Label)) %>%
    pull(interaction_group)
  
}

#~~~~~~
# LIANA
#######
liana_res_celltype <- list()

for (age in c("young", "old")) {
  
  liana_res_celltype[[age]] <- liana_wrap(so.list[[age]],
                                          resource = 'custom', # resource has to be set to 'custom' to work with external resources
                                          external_resource = ortholog_resource, # provide orthologous resource
                                          idents_col = "interaction_group", 
                                          min_cells = min_cells, 
                                          parallelize = TRUE, 
                                          workers = n_cores) %>% 
    liana_aggregate()
  
}

# aggregate and save results
liana_all_res <- list("cluster" = liana_res_clusters, "subcomp" = liana_res_subcomp, "celltype" = liana_res_celltype)
saveRDS(liana_all_res, file = paste0(resultsFolder, "liana_res_orthologs.rds"))

writeLines(capture.output(sessionInfo()), paste0(resultsFolder, "sessionInfo_hpc.txt"))

########################