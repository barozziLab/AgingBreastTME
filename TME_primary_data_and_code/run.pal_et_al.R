
in_folder <- "./"

in_F <- paste(in_folder, "libs.R", sep = "")
source(in_F)

#~~~~~~~~~~~~~~~~
#Data Preparation
#################

#SeuratObject_TNBCSub.rds
#can be downloaded from figshare
#https://figshare.com/articles/dataset/Data_R_code_and_output_Seurat_Objects_for_single_cell_RNA-seq_analysis_of_human_breast_tissues/17058077


## Load data sets of interest

dois <- c("TNBCSub")

data <- list()

for (doi in dois) {
  in_F <- paste(in_folder, "data_pal_et_al/SeuratObject_", doi, ".rds", sep = "")
  data[[doi]] <- readRDS(in_F)
}

## Load Metadata

in_F <- paste(in_folder, "data_pal_et_al/embj2020107333-sup-0006-tableev4.txt", sep = "")
metadata_v4 <- read_tsv(in_F)
metadata_v4 <- metadata_v4 %>% 
  mutate(group = `Sample Name`) %>% 
  mutate(group = gsub("-", "_", group))

in_F <- paste(in_folder, "data_pal_et_al/embj2020107333-sup-0006-tableev2.txt", sep = "")
metadata_v2 <- read_tsv(in_F)
metadata_v2 <- metadata_v2 %>% 
  mutate(group = `Specimen ID`) %>% 
  mutate(group = gsub("-", "_", group))

## Prepare signatures of interest

# Old vs Young CAFs

ib_sigs <- list()

in_f <- paste(in_folder, "signature.txt", sep = "")
id_tib <- read_tsv(in_f)

#too small, not using it
ib_sigs[["CAFs_old_vs_young_up"]] <- id_tib %>% 
  filter(gs_name == "old_vs_young_up" & !is.na(human_ortholog)) %>% 
  pull(human_ortholog)

ib_sigs[["CAFs_old_vs_young_down"]] <- id_tib %>% 
  filter(gs_name == "old_vs_young_down" & !is.na(human_ortholog)) %>% 
  pull(human_ortholog)

## Add scores

for (doi in dois) {
  for (sig in names(ib_sigs)) {
    data[[doi]] <- AddModuleScore(object = data[[doi]],
                                  features = list(ib_sigs[[sig]]),
                                  ctrl = 10,
                                  name = sig,
                                  seed = 1234,
                                  assay = "RNA")
    w <- colnames(data[[doi]]@meta.data) == paste(sig, "1", sep ="")
    colnames(data[[doi]]@meta.data)[w] <- sig
  }
}

## Rename clusters according to the paper annotation

n <- nrow(data[["TNBCSub"]]@meta.data)
data[["TNBCSub"]]@meta.data$cluster <- rep("None", n)
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 0
data[["TNBCSub"]]@meta.data$cluster[w] <- "T-cells"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 1
data[["TNBCSub"]]@meta.data$cluster[w] <- "TAMs"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 2
data[["TNBCSub"]]@meta.data$cluster[w] <- "Plasma-cells"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 3
data[["TNBCSub"]]@meta.data$cluster[w] <- "CAFs"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 4
data[["TNBCSub"]]@meta.data$cluster[w] <- "T-cells_Cycling"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 5
data[["TNBCSub"]]@meta.data$cluster[w] <- "B-cells"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 6
data[["TNBCSub"]]@meta.data$cluster[w] <- "DCs"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 7
data[["TNBCSub"]]@meta.data$cluster[w] <- "Endothelial"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 8
data[["TNBCSub"]]@meta.data$cluster[w] <- "Pericytes"
w <- data[["TNBCSub"]]@meta.data$seurat_clusters == 9
data[["TNBCSub"]]@meta.data$cluster[w] <- "Myeloid"

## Add metadata

#TNBCSub; add age

meta_add <- as_tibble(data[["TNBCSub"]]@meta.data) %>% 
  select(group) %>% 
  left_join(metadata_v2, by = "group") %>%
  select(group, `Patient Age`) %>%
  dplyr::rename(age = `Patient Age`) %>%
  mutate(age_group = ifelse(age >= 60, "old", "young"))

data[["TNBCSub"]]@meta.data$age <- meta_add$age
data[["TNBCSub"]]@meta.data$age_group <- meta_add$age_group

#################

#~~~~~
#Plots
######

plt_vUp <- VlnPlot(data$TNBCSub, 
                   feature = "CAFs_old_vs_young_up", 
                   group.by = "cluster", 
                   split.by = "age_group", 
                   pt.size = 0)

plt_vDn <- VlnPlot(data$TNBCSub, 
                   feature = "CAFs_old_vs_young_down", 
                   group.by = "cluster", 
                   split.by = "age_group", 
                   pt.size = 0)

#Extract CAFs data alone
caf_tib <- data$TNBCSub@meta.data %>% 
  as_tibble() %>% 
  select(CAFs_old_vs_young_up, CAFs_old_vs_young_down, cluster, age, age_group) %>% 
  filter(cluster == "CAFs") %>%
  mutate(age = as.factor(age))

p <- wilcox.test(caf_tib$CAFs_old_vs_young_up ~ caf_tib$age_group)$p.value
plt_title <- paste0("p = ", formatC(p, 2))
plt_up_1 <- ggplot(caf_tib, aes(x=age_group, y=CAFs_old_vs_young_up)) + 
  geom_boxplot() +
  theme_bw() +
  ggtitle(plt_title)
plt_up_2 <- ggplot(caf_tib, aes(x=age, y=CAFs_old_vs_young_up)) + 
  geom_boxplot() +
  theme_bw() +
  ggtitle("CAFs only")

p <- wilcox.test(caf_tib$CAFs_old_vs_young_down ~ caf_tib$age_group)$p.value
plt_title <- paste0("p = ", formatC(p))
plt_down_1 <- ggplot(caf_tib, aes(x=age_group, y=CAFs_old_vs_young_down)) + 
  geom_boxplot() +
  theme_bw() +
  ggtitle(plt_title)
plt_down_2 <- ggplot(caf_tib, aes(x=age, y=CAFs_old_vs_young_down)) + 
  geom_boxplot() +
  theme_bw() +
  ggtitle("CAFs only")

pdf("signature.pdf", width = 5, height = 6.5)

grid.arrange(plt_vUp, plt_up_1, plt_up_2, ncol = 2, 
             layout_matrix = cbind(c(1,2), c(1,3)), 
             heights=c(2,1), widths=c(1,2))

grid.arrange(plt_vDn, plt_down_1, plt_down_2, ncol = 2, 
             layout_matrix = cbind(c(1,2), c(1,3)), 
             heights=c(2,1), widths=c(1,2))

dev.off()

######


#~~~~~~~
#ML CAFs
########

# ML CAFs young vs old
# using the genes of our signature as input

ml_meta <- data$TNBCSub@meta.data %>% 
  rownames_to_column(var = "cell_id") %>% 
  as_tibble() %>% 
  select(cell_id, cluster, age, age_group) %>% 
  filter(cluster == "CAFs")

goi <- c(ib_sigs$CAFs_old_vs_young_up, ib_sigs$CAFs_old_vs_young_down) %>% unique()
goi_f <- goi[goi %in% rownames(data[[doi]]@assays$RNA)]

coi <- ml_meta$cell_id

ml_data <- data[["TNBCSub"]]@assays$RNA[goi_f,coi] %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>% 
  as_tibble()

#old TRUE/FALSE becomes y
ml_meta_i <- ml_meta %>% 
  mutate(old = ifelse(age_group == "old", TRUE, FALSE)) 

ml_data_i <- ml_data %>% 
  column_to_rownames("gene") %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "cell_id") %>% 
  as_tibble() %>% 
  left_join(ml_meta_i %>% select(cell_id, old), by = "cell_id") %>%
  dplyr::rename(y = old) %>%
  mutate(y = as.factor(y)) %>%
  select(-cell_id)

#delete genes with zero variance
ml_data_i_var <- ml_data_i %>% select(-c(y)) %>% 
  apply(2, var)
gene_bl <- ml_data_i_var[ml_data_i_var == 0] %>% 
  names()
ml_data_i <- ml_data_i %>% 
  select(colnames(ml_data_i)[!colnames(ml_data_i) %in% gene_bl])

#force manual subsampling (400 each)
set.seed(123)
w_pos <- sample(which(ml_data_i$y == TRUE), 400)
set.seed(123)
w_neg <- sample(which(ml_data_i$y == FALSE), 400)
ml_data_i_ss <- ml_data_i[c(w_pos, w_neg), ]

set.seed(123)
ml_split <- initial_split(ml_data_i_ss)
ml_train <- training(ml_split)
ml_test <- testing(ml_split)

ml_recipe <- training(ml_split) %>%
  recipe(y ~.) %>%
  step_center(all_predictors(), -all_outcomes()) %>%
  step_scale(all_predictors(), -all_outcomes()) %>%
  prep()

## RF (simple, no tuning)
#https://rviews.rstudio.com/2019/06/19/a-gentle-intro-to-tidymodels/

ml_testing <- ml_recipe %>%
  bake(testing(ml_split)) 

ml_training <- juice(ml_recipe)

ml_ranger <- rand_forest(trees = 100, mode = "classification") %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(y ~ ., data = ml_training)

my_metrics <- metric_set(accuracy, sens, spec, kap)

rf_confusion_mat <- ml_ranger %>%
  predict(ml_testing) %>%
  bind_cols(ml_testing) %>%
  my_metrics(truth = y, estimate = .pred_class)

#plt <- ml_ranger %>%
#  vip(geom = "col", num_features = 40) + theme_bw()

rf_feature_importance <- vip:::vi(ml_ranger) %>%
  mutate(up = Variable %in% ib_sigs$CAFs_old_vs_young_up) %>%
  mutate(down = Variable %in% ib_sigs$CAFs_old_vs_young_down)

#save conf. matrix and importance
out_F <- paste0("signature.RF.conf_matrix.txt")
write_tsv(x = rf_confusion_mat, file = out_F)
out_F <- paste0("signature.RF.feat_importance.txt")
write_tsv(x = rf_feature_importance, file = out_F)

########
