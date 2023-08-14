## METADATA ===============================================================
## Description: Exceptional local launch of summary statistics stats
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-13 23:38:56
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

source("./code/2_scripts/source.R")


library(doParallel)
library(foreach)



model <- "M0"
config_n <- 1
expp <- "perms_disp_comp_2000_M0"

# 
# # load parameter table
parms <- read.csv2(file.path(pls$dir_config_gen, expp, paste0(model,"_config_parameters.txt")))
parms[] <- lapply(parms, as.numeric)
# # load landscape file
lc <- readRDS(file.path(pls$dir_env_gen,"landscapes.rds"))



# system variables
output_dir <- file.path(pls$dir_out_zip, expp)
# get dirs
outputs <- list.files(output_dir, recursive = F, full.names = F)
#selection... run only if necessary
selection <- grep(model, outputs, perl=T)
outputs <- outputs[selection]
# order outputs
numbers <- as.numeric(sub(paste0(model,"_(\\d+)\\.zip"), "\\1", outputs))
outputs <- outputs[order(numbers)]

n_outs <- length(outputs)


# determine number of cores of your machines to use
cl <- makeCluster(detectCores()-1)

# register the cluster
registerDoParallel(cl)

# use foreach loop to run simulations in parallel
foreach(i=1:2) %dopar% {
  source("./code/2_scripts/source.R")
  source(file.path(pls$dir_base,"summary_stats/support/load_support_libs.R"))
  call_stats_zip(zf=outputs[o_i], i=output_dir, l=lc, o="3_summary", save=FALSE, t_i=0)
}

stopCluster(cl)
