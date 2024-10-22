############################################
## Run script for AgingBreastTME_analysis ##
############################################

if (!dir.exists("results")) { dir.create("results") }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load custom function library
##############################

source("TME.ds_functions.R")

#~~~~~~~~~~~~~~~~~~
# 0. Pre-processing - initial clustering, compartment re-assignment and identification of optimal clustering resolution
###################

# This script does not need to be executed, since the relevant output ("so.split.processed.rds") for downstream analyses is available on zenodo (see link in README.txt)
# The script is available mainly for documentation.
# If you decide to rerun it, be aware that the cluster assignments might not be 100% identical to our result
# and that the subsequent scripts might therefore not fully work, as the code partially uses the manual labels from the publication.

#source("0_combined_analysis.R")
so.split.processed <- readRDS(file = "../data/so.split.processed.rds")

#~~~~~~~~~~~~~~~
# Prepare object
################

annotation_file <- "../data/Table_S2_plus_colors.xlsx"
conv <- read_tsv("../data/label_condition_conversion.txt")
cluster_anno <- read_xlsx(annotation_file) %>%
  mutate(compartment = ifelse(compartment == "cd45-", "stromal", compartment))
color_table <- read_tsv("../data/plotting_colors.tsv")

# update metadata slots
for (so_name in names(so.split.processed)) {
  ca <- cluster_anno %>% filter(compartment == so_name) %>%
    mutate(cluster_label = paste(`Cell-type_label`, Cluster_number, sep = ".")) %>%
    mutate(Cluster_number = factor(Cluster_number))
  
  so.split.processed[[so_name]]$cluster_label <- so.split.processed[[so_name]]@meta.data %>%
    left_join(ca, by = c("seurat_clusters" = "Cluster_number")) %>%
    pull(cluster_label)
  
  so.split.processed[[so_name]]$age_tm <- so.split.processed[[so_name]]@meta.data %>%
    mutate(age_tm = paste(age, treatment, sep = " ")) %>%
    pull(age_tm)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Run clustree analysis
##########################

source("1_clustree.R")

#~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Run proportion tests
#########################

log2FD_threshold = log2(1.5) # log2FD threshold for proportion tests
source("2_proportion_test.R")

#~~~~~~~~~~~~~~~~~~~~~~~
# 3. Downstream analyses
########################

source("3_downstream_analyses.R")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

########################