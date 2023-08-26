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

# parameters
if(!interactive()){
  cnfs <- commandArgs(trailingOnly = TRUE)
  print(cnfs)}
model <- if(interactive()){"M0"}else{cnfs[1]}
expp <- if(interactive()){"perms_disp_comp_2000_M0"}else{cnfs[2]}

# create summary table template
expp_st=list("parms"=NULL, "stats"=list("t"=list(), "tt"=list()))

# # load parameters table
expp_st$parms <- load_parameters(dc=pls$dir_config_gen, expp, model)

#### CALL STATS FOR SINGLE TIME t
# # load summary stats
ssl <- load_temp_summary(expp, model, t_i="t_0", conf="all")
dbssl <- do.call(rbind, ssl)
dbssl <- as.data.frame(dbssl)
# convert stats to numeric
dbssl[!colnames(dbssl)%in%"finished"]  <- lapply(dbssl[!colnames(dbssl)%in%"finished"], as.numeric)
# store
expp_st$stats$t <- dbssl

##### CALL STATS FOR CHANGING TIME tt
# load and calc stats individually to avoid large memory requirements
ssl <- list()
ssfiles <- get_ordered_files(file.path(pls$dir_out_zip, "temp_summary", expp, t_i="t_500-0"), model, fend=".rds")
n_ss <- length(ssfiles)
strings <- unlist(lapply(strsplit(ssfiles,"_"), function(x){paste(tail(x,2),collapse="_")}))
strings <- paste0("config_",gsub(".rds", "", strings))
for (conf_i in 1:n_ss){
  # conf_i <- 1
  ssltt <- NULL # dump ssltt in every interaction
  ssltt <- load_temp_summary(expp, model, t_i="t_500-0", conf=conf_i)[[1]]
  ssl[[conf_i]] <- call_sstatstt_time(ssltt)
}

# creat dummy final tt element structure on our expp_st list
names(ssl) <- strings
dummy <- ssl[[1]]
dummy <- lapply(dummy, function(x)NULL)
n_timephase <- length(dummy)


for (conf_i in 1:n_ss){
  for (tp_i in 1:n_timephase){
    dummy[[tp_i]] <- rbind(dummy[[tp_i]], ssl[[conf_i]][[tp_i]])
  }
}

# transform list of matrix to list of data.frames and add correct line label
tt <- lapply(dummy, function(x){
  dbssl <- as.data.frame(x)
  # convert stats to numeric
  dbssl  <- lapply(dbssl, as.numeric)
  dbssl <- do.call(cbind, dbssl)
  rownames(dbssl) <- strings
  return(dbssl)}
)
for (tti in 1:length(tt)){
  tt[[tti]] <- as.data.frame(tt[[tti]])
}
# store
expp_st$stats$tt <- tt
# save file
time_stamp <- format(Sys.time(), "%Y%m%d_%H%M")
saveRDS(expp_st, file.path(pls$dir_out_zip, "temp_summary", expp, paste0(time_stamp,"_summaries_", expp,".rds")))

# remove single stats file from temp_summary/*expp*/t_0/*
# do.call(unlink, list(list.files(file.path(pls$dir_out_zip, "temp_summary", expp, "t_0"), full.names = TRUE)))