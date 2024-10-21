#~~~~~~~~~
# Packages
##########

library(tidyverse)
library(Seurat)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(gprofiler2)
library(clusterCrit)
library(scProportionTest)
library(clustree)
library(data.table)
library(readxl)
library(ComplexHeatmap)
library(clusterProfiler)

#~~~~~~~~~~~~~~~~~
# Custom functions
##################

source("TMS.ds_functions.R")

#source("../../MultiK_parallel_consMem/multiK_conserveMem.fct.R")
#source("../../MultiK_parallel_consMem/clustering.fct.R")

source("multiK_conserveMem.fct.R")
source("clustering.fct.R")

#~~~~~~~~~~~~~~
# Color schemes
###############

colors_tme <- read_tsv("plotting_colors.tsv")