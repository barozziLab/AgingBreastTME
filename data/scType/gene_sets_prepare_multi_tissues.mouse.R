# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
# Modified by Iros Barozzi, April 2022
# Further modified by Stephan Gruener, January 2023 (added line 57 and adapted 18-20, 59+60 to allow for annotating across several celltypes (tissues))
#
# Functions on this page:
# gene_sets_prepare: prepare gene sets and calculate marker sensitivity from input Cell Type tabular file
#
# @params: path_to_db_file - DB file with cell types
# @cell_type - cell type (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
#
# Temporary function that does not use checkGeneSymbols
# Because doesn't seem to work properly with mouse gene symbols (old mappings?)

gene_sets_prepare_multi_bach <- function(path_to_db_file, cell_type){
  
  cell_markers <- read_tsv(path_to_db_file)
  if(!is.null(cell_type)) {
  cell_markers <- cell_markers[cell_markers$tissueType %in% cell_type,] 
  }
  cell_markers$geneSymbolmore1 <- gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 <- gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 <- sapply(1:nrow(cell_markers), function(i){
    
    markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all <- markers_all[markers_all != "NA" & markers_all != ""]
    markers_all <- sort(markers_all)
    
    if(length(markers_all) > 0){
      markers_all <- unique(na.omit(markers_all))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 <- sapply(1:nrow(cell_markers), function(i){
    
    markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all <- markers_all[markers_all != "NA" & markers_all != ""]
    markers_all <- sort(markers_all)
    
    if(length(markers_all) > 0){
      markers_all <- unique(na.omit(markers_all))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 <- gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 <- gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 <- gsub(" ","",cell_markers$geneSymbolmore2)
  
  # need the renaming step here (because no unique labels across tissues)
  cell_markers = cell_markers %>% mutate(comb_name = paste(tissueType, cellName, sep = "_"))
  
  gs <- lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) <- cell_markers$comb_name
  gs2 <- lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) <- cell_markers$comb_name
  
  list(gs_positive = gs, gs_negative = gs2)
}
