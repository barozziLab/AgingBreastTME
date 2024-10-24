# Generate ortholog resource

library(tidyverse)
library(OmnipathR)
library(liana)
library(magrittr)

# Here, we will convert LIANA's Consensus resource to murine symbols
op_resource <- select_resource("Consensus")[[1]]

# Generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) # mouse

if (!dir.exists("results")) { dir.create("results") }
saveRDS(ortholog_resource, file = "results/ortholog_resource.rds")

########################