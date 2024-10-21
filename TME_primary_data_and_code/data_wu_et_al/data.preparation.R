
library("Matrix")
require("tidyverse")
require("Seurat")


#~~~~~~~~~~~~~~~~
#Data Preparation
#################

#This file 
#metadata.csv
#is part of 
#GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078

metadata <- read_csv("metadata.csv")

meta_s1 <- read_tsv("41588_2021_911_MOESM4_ESM.S1.txt")
meta_s1 <- meta_s1 %>%
  rename(orig.ident = Case_ID) %>%
  mutate(orig.ident = paste0("CID", orig.ident))
meta_s1$orig.ident <- gsub("-", "", meta_s1$orig.ident)

#add metadata from table S1 (incl. age)
metadata <- metadata %>% 
  left_join(meta_s1, by = "orig.ident")

##

#These files
#count_matrix_barcodes.tsv; count_matrix_genes.tsv; count_matrix_sparse.mtx 
#are part of 
#GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078

inF <- paste("count_matrix_barcodes.tsv", sep = "")
d_bc <- read_tsv(inF, col_names = FALSE)

inF <- paste("count_matrix_genes.tsv", sep = "")
d_gene <- read_tsv(inF, col_names = FALSE)

inF <- paste("count_matrix_sparse.mtx", sep = "")
d <- readMM(inF)

rownames(d) <- d_gene %>% pull(X1)
colnames(d) <- d_bc %>% pull(X1)

#check
#metadata$cell_id == colnames(d)


#~~~~~~~~~~~~~~~~~~~~~~~
#TNBC - Non-Cancer Cells
########################

#According to the analysis and text
#"UMAP visualization of stromal and immune cells across [..]
#[..] tumors clustered together without batch correction (Extended Data Fig. 1e,f)."
#So I'm extracting and treating TME cells without any correction

#Filter

filt <- metadata$subtype == "TNBC" & 
  metadata$celltype_major != "Cancer Epithelial" & 
  metadata$celltype_major != "Normal Epithelial"
  
metadata_df <- metadata %>%
  column_to_rownames(var = "cell_id") %>%
  as.data.frame()

#Create object

so <- CreateSeuratObject(counts = d[,filt], 
                         min.cells = 10, 
                         min.features = 10, 
                         project = "TNBC", 
                         meta.data = metadata_df[filt,])

so <- NormalizeData(object = so)
so <- ScaleData(object = so)
so <- FindVariableFeatures(object = so)
so <- RunPCA(so, verbose = FALSE)
#plot(so@reductions$pca@stdev)
pca_dim_sel <- 30

so <- FindNeighbors(so, dims = 1:pca_dim_sel)
so <- FindClusters(so, resolution = 0.8)
so <- RunUMAP(object = so, dims = 1:pca_dim_sel)

p1 <- DimPlot(object = so, reduction = 'umap', group.by = 'celltype_major')
#DimPlot(object = so, reduction = 'umap', group.by = 'celltype_minor') 
p2 <- DimPlot(object = so, reduction = 'umap', group.by = 'Age')

pdf("so_TNBC_TME.umap.pdf", width = 10, height = 4)
plot(p1 + p2)
dev.off()

goi <- c("COL1A1", "PECAM1", "SPI1", "JCHAIN", "CD3E", "CD19")
p_vln <- VlnPlot(so, features = c(goi), group.by = 'celltype_major', pt.size = 0)

pdf("so_TNBC_TME.markers.pdf", width = 7, height = 8)
plot(p_vln)
dev.off()

saveRDS(so, "so_TNBC_TME.rds")

