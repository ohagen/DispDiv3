## METADATA ===============================================================
## Description: Summary statistics
## 
## R version: 4.2.2 for Windows
## Date: 2023-07-27 13:30:26
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

source("./code/2_scripts/source.R")
print(unlist(pls))
source(file.path(pls$dir_base,"summary_stats/support/load_support_libs.R"))

## 
cnfs <- commandArgs(trailingOnly = TRUE)
print(cnfs)
model <- cnfs[1] #e.g. model <- "M0"
config_n <- as.numeric(cnfs[2]) # e.g. config_n <- 1
expp <- cnfs[3] # e.g. expp <- "perms_disp_comp_2000_M0"

# # load parameters table
parms <- load_parameters(dc=pls$dir_config_gen, expp, model)
# # load landscape file
lc <- readRDS(file.path(pls$dir_env_gen,"landscapes.rds"))

# system variables
output_dir <- file.path(pls$dir_out_zip, expp)
# get dirs in order
outputs <- get_ordered_files(output_dir, model, fend=".zip")
# get size
n_outs <- length(outputs)
print(n_outs)
####       CALL STATS!      ###
#### FIX and var TIME #########
for (ti in list(0, 500:0)){ # first for the single time step, then for a sequence of time steps...
  print(paste0("working with: ", outputs[config_n], ", config_n=", config_n, " | at ", length(ti), " time-steps | experiment ", expp))
  # ti <- 0
  call_stats_zip(zf=outputs[config_n], 
                 i=output_dir, 
                 l=lc, 
                 o=file.path(pls$dir_out_zip, "temp_summary", expp), 
                 save=TRUE, # save each as separated .rds
                 t_i=ti) 
}


###################### END ##########################