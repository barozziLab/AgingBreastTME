# Age-Dependent Differences in Breast Tumor Microenvironment
A repository of code for data analysis/processing generated for the study of age-dependent differences in the breast tumor microenvironment.

The study was performed in collaboration at the European Institute of Oncology (IEO) and the Center for Cancer Research at the Medical University of Vienna (MUW). The dataset contains transcriptional profiles of tumors and the tumor microenvironment in adult (12 months) and young (6-8 weeks) mice, generated by orthotopic injection of 20k 4t1 cells into the mammary fat pad.

Primary findings for this project were published here: (paper link to be added)

The raw data can be accessed at: (GEO links to be added here)

An intermediate object that can be used together with the provided code to reproduce the figures from the publication, but also including further exploratory analyses we conducted is available on zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13960874.svg)](https://doi.org/10.5281/zenodo.13960874)

## Study design
![Approach_Outline](https://github.com/user-attachments/assets/d7645144-6530-48e8-9066-ed9b60b1a3c4)
TT ... Triple Therapy  
C ... Cyclophosphamide  
V ... Vinorelbine  
αPD-1 ... monoclonal antibody targeting PD1  

## Reproducing the analyses

For reproducing the analyses for the publication, we recommend using the intermediate file provided on the zenodo link above (`so.split.processed.rds`), which contains the pre-processed and clustered cells as seurat objects, split by compartment (tumor, immune and stroma). Therefore, download the file to the `data` directory.

### Aging Breast TME

If you wish to additionally redo the analysis from the QC filtered objects to this file, you can do so by additionally downloading the three other files on zenodo. We recommend using our intermediate file, as minimal changes in the clustering results may result in problems with the downstream analyses, since the code contains sections in which manual labels from the publication are added to the object.
In case you still wish to rerun this part, download the three objects (`published_immune.rds`, `stromal.rds`, `tumor_and_immune.rds`) and save them in the following directory: `data/objects_after_QC` and activate the execution of `0_combined_analysis.R` in `RUN_AgingBreastTME_analyses.R`.

All scripts for this part of the analysis can be executed by running `RUN_AgingBreastTME_analyses.R` in `AgingBreastTME_analysis`, which is also where the sessionInfo is provided.

### TME primary data



### TME Cell-cell communication analysis

All analyses can be executed with the `RUN_liana.R` script, however, note that `2_liana_orthologs.R` was executed in a separate R environment on an HPC. The separate sessionInfo file is available in `TME_CellCellCommunication/sessionInfo_hpc.txt`, while `TME_CellCellCommunication/sessionInfo.txt` refers to the environment that all other scripts were executed in.

### TMS analysis

The raw data from the Tabula Muris Senis compendium can be downloaded here: https://figshare.com/ndownloader/files/23872862 and converted to a seurat object, then saved in the `data/` directory.
All analyses can be executed via the run script (`TMS_RUN_young_vs_old.R`), in which parameters and options for the analysis can be selected and adapted.

The `multiK_conserveMem.fct` script contains the code for a memory-efficient and parallelizable reimplementation of multiK (original code see: https://github.com/perou-lab/MultiK), which we used to assess optimal clustering resolutions, together with the silhouette score.

This whole analysis was conducted on an HPC, therefore the used R environment differs from the one used for the other analyses. The respective sessionInfo file can be found in `TMS_analysis/`.

## Directory Structure
