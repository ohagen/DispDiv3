## METADATA ===============================================================
## Description: Main call_stats_zip to call set of summary statistics
## from the list of functions: 
## lsstats for a single time-step
## lsstatstt for multiple time-steps
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-02 20:43:29
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

#### CALL STATISTICS FUNCTION
call_stats_zip <- function(zf, i=output_dir, l=landscape, o=dir_summaries, verbouse=FALSE, save=FALSE, t_i=0){
  # zf=zip file name
  # l is the loaded landscape file .rds # we use only the same landscape
  # i = the directory where the zipfiles are
  # o = the directory to save the results
  # t_i is the time-step or the time-steps in the first case, temp_variables_tx is the loaded object , otherwise,
  # its a list of loaded objects, e.g. when t_i=500:0
  if (verbouse){print(paste(zf, ": summaries: "))}
  # reset temp variables to NA in case something is missing... so not to use the last loaded variables
  temp_sgen3sis <<- NA
  temp_phy <<- NA
  temp_human_phy <<- NA
  temp_occ_tx <<- NA
  temp_abundance_tx <<- NA
  temp_traits_tx <<- NA
  temp_l <<- l
  temp_t_i <<- t_i
  # load files from zip
  #sgen3sis # timeless ti=NULL as default!
  temp_sgen3sis <<- read_zip(what="sgen3sis.rds", zf, i=i)
  #phy.nex
  temp_phy <<- read_zip(what="phy.nex", zf, i=i)
  # human phy.txt
  temp_human_phy <<-read_zip(what="phy.txt", zf, i=i)
  
  n_ts <- length(t_i)
  if (n_ts==1){
    ssname <- paste0("t_", t_i)
    runF <- lsstats # specify here the dunctions to be run
    #pa_t_x
    temp_occ_tx <<- read_zip(what="occs/pa_t_x.rds", zf, i=i, ti=t_i)
    #temp_abundance_tx
    temp_abundance_tx <<- read_zip(what="abundance/abundance_t_x.rds", zf, i=i, ti=t_i)
    #temp_traits_tx
    temp_traits_tx <<- read_zip(what="traits/traits_t_x.rds", zf, i=i, ti=t_i)
  } else { # we return a list
    ssname <- paste0("t_", paste0(c(t_i[1],t_i[n_ts]), collapse = "-"))
    runF <- lsstatstt # specify here the dunctions to be run
    temp_occ_tx <<- list()
    temp_abundance_tx <<- list()
    temp_traits_tx <<- list()
    for (txi in 1:n_ts){
      temp_occ_tx[[txi]] <- read_zip(what="occs/pa_t_x.rds", zf, i=i, ti=t_i[txi])
      temp_abundance_tx[[txi]] <- read_zip(what="abundance/abundance_t_x.rds", zf, i=i, ti=t_i[txi])
      temp_traits_tx[[txi]] <- read_zip(what="traits/traits_t_x.rds", zf, i=i, ti=t_i[txi])
    }
    names(temp_occ_tx) <- t_i
    temp_occ_tx <<- temp_occ_tx
    # names(temp_abundance_tx) <- t_i
    temp_abundance_tx <<- temp_abundance_tx
    names(temp_traits_tx) <- t_i
    temp_traits_tx <<- temp_traits_tx
  }
  
  ### RUN STATS ###
  # call runF (either lsstats or lsstatstt) functions
  ss <- lapply(runF, function(x) x())
  ### RUN STATS ###
  #format if ti otherwise keep as is, i.e. list
  if (n_ts==1){ssu <- unlist(ss, use.names = T)}else{ssu <- ss}
  names(ssu) <- gsub("\\.", "_", names(ssu))

  #saving
  if (save){ # if save
    if (!dir.exists(file.path(o, ssname))){
      dir.create(file.path(o, ssname), recursive=T)
    }
    fname <- paste0( "ss_", ssname , "_", gsub(".zip","",tail(strsplit(zf, "/")[[1]],1)),".rds")
    saveRDS(ssu, file.path(o,ssname, fname))
    
    if (verbouse){ # if verbouse print
      print(ssu)
      print(paste("! saving outputs to:", file.path(o,ssname, fname)))
    } # end print
  } # end if save
  return(ssu)
}

#### READ ZIP FILES
read_zip <- function(what="sgen3sis.rds", zf, i=output_dir, ti=NULL){
  # print(paste("reading", what, "at", i))
  # ti, either NULL, a time (numeric)
  what <- gsub("t_x", paste0("t_", ti), what)
  con <- unz(file.path(i,zf), file.path( "output", paste0("config_",gsub(".zip", "", zf)), what))
  fext <- file_ext(what)
  if(fext=="rds"){
    file <- readRDS(gzcon(con)) 
  } else if(fext=="nex"){
    file <- read.nexus(gzcon(con)) 
  } else if(fext=="txt"){
    file <- read.table(unzip(file.path(i,zf), file.path( "output", paste0("config_",gsub(".zip", "", zf)), what)), 
                       header=T, quote="\"", sep="\t") 
  } else {
    close(con)
    stop(paste(fext, "is an unknown file extention! Please try .rds or .nex"))
  }
  close(con)
  return(file)
}


#### CALL ORDERED FILES
get_ordered_files <- function(output_dir, model, fend=".zip"){
  outputs <- list.files(output_dir, recursive = F, full.names = F)
  #selection... doublecheck!
  selection <- grep(model, outputs, perl=T)
  outputs <- outputs[selection]
  # order outputs
  stringsl <- strsplit(outputs,"_")
  strings <- unlist(lapply(stringsl, function(x){x[length(x)]}))
  numbers <- as.numeric(sub(fend, "", strings))
  outputs <- outputs[order(numbers)]
  return(outputs)
}


#### CALL PARAMETERS TABLE
load_parameters <- function(dc=pls$dir_config_gen, expp, model){
  parms <- read.csv2(file.path(dc, expp, paste0(model,"_config_parameters.txt")))
  parms[] <- lapply(parms, as.numeric)
  return(parms)
}
