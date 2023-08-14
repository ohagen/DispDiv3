## METADATA ===============================================================
## Description: Compile single statistics, normaly stored at the folder:
## psl$dir_out_zip/temp_summary and either [t_X] or [t_S-E] being X a 
## single time step, e.g. t_0 and S-E a range of time-steps, e.g. 500-0
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-14 11:02:23
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

source("./code/2_scripts/source.R")
print(pls)
source(file.path(pls$dir_base,"summary_stats/support/load_support_libs.R"))

cnfs <- commandArgs(trailingOnly = TRUE)
print(cnfs)
model <- cnfs[1] #e.g. model <- "M0"
expp <- cnfs[2] # e.g. expp <- "perms_disp_comp_2000_M0"

# # load parameters table
parms <- load_parameters(dc=pls$dir_config_gen, expp, model)
ssfiles <- get_ordered_files(file.path(pls$dir_out_zip, "temp_summary", expp, t_i), model, fend=".rds")
n_ss <- length(ssfiles)
# prepare list
ssl <- list()
for (o_i in 1:n_ss){
  print(paste0("working with: ", ssfiles[o_i], ", o_i=", o_i))
  # dummy for avoiding wrong inheritance by indexes
  vals <- NA
  # load
  vals <- readRDS(file.path(pls$dir_out_zip, "temp_summary", expp, t_i, ssfiles[o_i]))
  #store
  ssl[[o_i]] <- vals
}

# merge single time stats into data frame
strings <- unlist(lapply(strsplit(outputs,"_"), function(x){paste(tail(x,2),collapse="_")}))
names(ssl) <- paste0("config_",gsub(".rds", "", strings))
dbssl <- do.call(rbind, ssl)
dbssl <- as.data.frame(dbssl)
# convert stats to numeric
dbssl[!colnames(dbssl)%in%"finished"]  <- lapply(dbssl[!colnames(dbssl)%in%"finished"], as.numeric)
# get final merged table
mbt <- cbind(parms[rownames(dbssl),], dbssl)
mbt[is.na(mbt)] <- NA
#### SAVE MBT
time_stamp <- format(Sys.time(), "%Y%m%d_%H%M")
saveRDS(mbt, file.path(pls$dir_out_zip, "temp_summary", expp, paste0(time_stamp,"_ti_", expp,".rds")))
# mbt <- readRDS("./3_summary/20230729_2000_ti_permuts_disp_comp_M0.rds")
# mbt <- readRDS("./3_summary/20230730_2316_ti_perms_disp_comp_2000_M0.rds")
# mbt <- readRDS("./3_summary/20230802_2105_ti_perms_disp_comp_2000_M0.rds")
# mbt <- readRDS("./3_summary/20230811_1127_ti_perms_disp_comp_2000_M0.rds")





##### CALL STATS FOR CHANGING TIME


#### CHANGING TIME #########
ssltt <- list()
for (o_i in 1:n_outs){ 
  ssltt[[o_i]] <- call_stats_zip(zf=outputs[config_n], i=output_dir, l=lc, o="3_summary", save=FALSE, t_i=500:0) # save each line as separated .rds
}
# merge single time stats into data frame
names(ssltt) <- paste0("config_",gsub(".zip", "", outputs))
# transform list of matrix to list of data.frames
ssltt <- lapply(ssltt, function(x){
  lapply(x, function(x){
    dbssl <- as.data.frame(x)
    # convert stats to numeric
    dbssl[!colnames(dbssl)%in%"finished"]  <- lapply(dbssl[!colnames(dbssl)%in%"finished"], as.numeric)
    return(dbssl)}
  )})

##### SAVE SSLTT
time_stamp <- format(Sys.time(), "%Y%m%d_%H%M")
saveRDS(ssltt, paste0("./3_summary/", time_stamp,"_tt_", expp,".rds"))




