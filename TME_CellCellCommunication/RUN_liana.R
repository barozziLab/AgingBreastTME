############################################
## LIANA Cell-cell communication analysis ##
############################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate ortholog resource
############################

source("1_liana_generate_ortholog_resource.R")

#~~~~~~~~~~~~~~~~~~~~
# Run LIANA framework
#####################

source("2_liana_orthologs.R") # This script was executed in a different R environment on an HPC (see sessionInfo_hpc.txt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Results exploration and plotting
##################################

source("liana.fct.lib.R")
source("3_liana_plotting.R")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

########################