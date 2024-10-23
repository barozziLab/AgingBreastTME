# This analysis aims to reveal differences in cell populations in the mammary gland of young and old mice

source("TMS_header.R")

#~~~~~~~~~~~~~~~
#setup variables
################

dataFolder <- "../data/"
resultsFolder <- "results/"

#~~~~~~~~~~~~~~~~~~~
#Parameter selection
####################

n_features <- 3000              # number of variable features to use for downstream analyses after SCTransform normalization and variance stabilization
pca_dim_sel <- 50               # number of principal components to use for downstream analyses (e.g. UMAP)
subsample_n <- 0                # number of cells to subsample (for Silhouette score), set to 0 if no subsampling should be performed
run_preAnalysis <- FALSE         # whether to run the pre-analysis pipeline with normalization, PCA, UMAP and FindNeighbors (if there is no corresponding file available, this will be done anyway)
run_silhouette <- FALSE          # whether to run silhouette analysis
run_clustree <- FALSE            # whether to run clustree analysis
run_multiK <- FALSE              # whether to run multiK analysis
n_cores <- 2                    # number of cores available for multi-core processing
log2FD_threshold <- log2(1.5)   # FD threshold for cluster composition analysis

#~~~~~~~~~~~~~~~~~~~
# load seurat object
####################

# The data can be downloaded here: https://figshare.com/ndownloader/files/23872862 and converted to a seurat object
so <- readRDS(file = paste0(dataFolder, "Mammary_Gland_droplet_raw.rds"))

#~~~~~~~~~~~~~~~~~~~~~~
# Pre-Analysis pipeline
#######################

# identify optimal clustering resolutions
source("TMS_clustering_res.R")

# Cell cycle scoring
m.s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", 
                  target_organism = "mmusculus")$ortholog_name

m.g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", 
                    target_organism = "mmusculus")$ortholog_name

so <- CellCycleScoring(so, s.features = m.s.genes, g2m.features = m.g2m.genes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Marker genes by cell type
###########################

Idents(object = so) <- "cell_ontology_class"

markers_cell_ont <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1) %>% 
  filter(p_val_adj <= 0.05) %>%
  mutate(pct.ratio = pct.1 / pct.2, pct.difference = pct.1 - pct.2)

write_tsv(x = markers_cell_ont, file = paste0(resultsFolder, "markers_cell_ontology.tsv"), col_names = T)

#~~~~~~~~~~~~~~~~~~~~~~~~
# Clustering & Signatures
#########################

# clustering resolutions for downstream analyses
clustering_res <- c(0.3, 1.1)

# Clustering, composition analysis (test for significant overrepresantation of certain conditions in each cluster), 
# marker gene identification and stromal gene signature expression

cl <- parallel::makeCluster(n_cores)

doSNOW::registerDoSNOW(cl)
foreach::foreach(i = 1:length(clustering_res), 
                 .export = c("resultsFolder", "pca_dim_sel", "so", "colors_tme", "log2FD_threshold", "pca_dim_sel"), 
                 .packages = c("Seurat", "tidyverse", "patchwork", "RColorBrewer", "readxl", "data.table", "scProportionTest", "ComplexHeatmap", "clusterProfiler"),
                 .verbose = TRUE) %dopar% 
  {
    
    res <- clustering_res[i]
    source("clustering_and_composition.R", local = TRUE)
    NULL # return nothing, all the output is saved
    
  }
parallel::stopCluster(cl)

writeLines(capture.output(sessionInfo()), "sessionInfo_tms.txt")

########################
