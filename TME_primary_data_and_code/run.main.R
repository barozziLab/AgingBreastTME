
in_folder <- "./"

# Run Pal et al analysis
in_F <- paste(in_folder, "run.pal_et_al.R", sep = "")
source(in_F)

# Prepare Wu et al data
in_F <- paste(in_folder, "data_wu_et_al/data.preparation.R", sep = "")
source(in_F)

# Run Wu et al analysis
in_F <- paste(in_folder, "run.wu_et_al.R", sep = "")
source(in_F)

# Run comparison analyses
in_F <- paste(in_folder, "run.comparison.R", sep = "")
source(in_F)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

########################