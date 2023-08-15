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
source(file.path(pls$dir_base,"summary_stats/support/call_summary_sstatstt_time.R"))

cnfs <- commandArgs(trailingOnly = TRUE)
print(cnfs)
model <- cnfs[1] #e.g. model <- "M0"
expp <- cnfs[2] # e.g. expp <- "perms_disp_comp_2000_M0"

# create summary table template
expp_st=list("parms"=NULL, "stats"=list("t"=list(), "tt"=list()))

# # load parameters table
expp_st$parms <- load_parameters(dc=pls$dir_config_gen, expp, model)


#### SINGLE TIME t
# # load summary stats
ssl <- load_temp_summary(expp, model, t_i="t_0", conf="all")
dbssl <- do.call(rbind, ssl)
dbssl <- as.data.frame(dbssl)
# convert stats to numeric
dbssl[!colnames(dbssl)%in%"finished"]  <- lapply(dbssl[!colnames(dbssl)%in%"finished"], as.numeric)
# store
expp_st$stats$t <- dbssl

    # get final merged table
    # mbt <- cbind(parms[rownames(dbssl),], dbssl)
    # mbt[is.na(mbt)] <- NA
    #### SAVE MBT
    # time_stamp <- format(Sys.time(), "%Y%m%d_%H%M")
    # saveRDS(mbt, file.path(pls$dir_out_zip, "temp_summary", expp, paste0(time_stamp,"_ti_", expp,".rds")))
# mbt <- readRDS("./3_summary/20230729_2000_ti_permuts_disp_comp_M0.rds")
# mbt <- readRDS("./3_summary/20230730_2316_ti_perms_disp_comp_2000_M0.rds")
# mbt <- readRDS("./3_summary/20230802_2105_ti_perms_disp_comp_2000_M0.rds")
# mbt <- readRDS("./3_summary/20230811_1127_ti_perms_disp_comp_2000_M0.rds")

##### CALL STATS FOR CHANGING TIME

#### SINGLE TIME t

# load and calc stats individually to avoid large memory requirements
ssl <- list()
ssfiles <- get_ordered_files(file.path(pls$dir_out_zip, "temp_summary", expp, t_i="t_500-0"), model, fend=".rds")
n_ss <- length(ssfiles)
strings <- unlist(lapply(strsplit(ssfiles,"_"), function(x){paste(tail(x,2),collapse="_")}))
strings <- paste0("config_",gsub(".rds", "", strings))
for (conf_i in 1:n_ss){
  # conf_i <- 2
  ssltt <- NULL
  ssltt <- load_temp_summary(expp, model, t_i="t_500-0", conf=conf_i)[[1]]
  ssl[[conf_i]] <- call_sstatstt_time(ssltt)
}
names(ssl) <- strings
dummy <- ssl[[1]]
dummy <- lapply(dummy, function(x)NULL)
n_timephase <- length(dummy)
tpstat <- list()
for (tp_i in 1:n_timephase){
  # tp_i <- 1
  tpstat[[tp_i]] <- lapply(ssl, function(x){
    vals <- x[[tp_i]]
    return(vals)
  })
  
  
}

dbssl <- do.call(rbind, ssl)
dbssl <- as.data.frame(dbssl)


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




