source(file.path(pls$dir_base,"summary_stats/support/lsstatstt_time.R"))


# LOAD SSI
call_sstatstt_time <- function(ssi){
  # ssi is the summary stat file, e.g. sstatsstt_t_500-0 .rds file
  # return a list of stats for the phases defined at the source.R under patches 
  # clean
  isli <- unlist(lapply(ssi, is.list)) # check if there are sublists
  istrait <- grep( "mean_traits", names(ssi))
  isall <- isli
  isall[istrait] <- TRUE
  ssin <- list()
  ssin <- ssi[!isall] # select only unnested first
  trnames=c("dispersal", "competition", "temp_width", "mean_temp") # make trait a nested list
  temp <- get_traits_colum(ssi[istrait][[1]], traits=trnames)
  names(temp) <- paste0("mean_tr_", trnames)
  ssin <- append(ssin, temp)
  # nested lists
  temp <- ssi[isli][[1]]
  names(temp) <- paste0(names(ssi)[isli], "_", names(temp))
  ssin <- append(ssin, temp) # ssin is now ready for the loop
  timephaz <- patches$`time-phase`
  n_era <- length(timephaz)
  lera <- list()
  for (era_i in 1:n_era){
    # era_i <- 3
    nice_range <- get_timesteps_in_Ma(timephaz[[era_i]])
    print(paste("At time phase:", names(timephaz)[n_era], nice_range))
    # apply time selection
    temp_ssin <- lapply(ssin, function(x){
      x <- x[as.numeric(rownames(x))%in%timephaz[[era_i]],]
      return(x)
    }) #as.numeric(rownames(x))%in%tt
    tempx<- lapply(lssttstt_time, function(x) x())
    tempx <- unlist(tempx)
    names(tempx) <- gsub("\\.", "_", names(tempx))
    lera[[era_i]] <- tempx
  }
  names(lera) <- names(timephaz)
  return(lera)
}

load_temp_summary <- function(expp, model, t_i, conf="all"){
  # load temp summaries as a list
  # if conf="all", then get all files
  # else, conf is the number of the config desired
  ssfiles <- get_ordered_files(file.path(pls$dir_out_zip, "temp_summary", expp, t_i), model, fend=".rds")
  n_ss <- length(ssfiles)
  # prepare list
  ssl <- list()
  strings <- unlist(lapply(strsplit(ssfiles,"_"), function(x){paste(tail(x,2),collapse="_")}))
  strings <- paste0("config_",gsub(".rds", "", strings))
  if (conf=="all"){
    vecfor <- 1:n_ss
    #rename=T
    idx <- vecfor
  }else{
    vecfor <- conf
    #rename=F
    idx <- rep(1,n_ss)
  }
  for (o_i in vecfor){
    print(paste0("working with: ", ssfiles[o_i], ", o_i=", o_i))
    # dummy for avoiding wrong inheritance by indexes
    vals <- NA
    # load
    vals <- readRDS(file.path(pls$dir_out_zip, "temp_summary", expp, t_i, ssfiles[o_i]))
    #store
    ssl[[idx[o_i]]] <- vals
    # ssl[[o_i]] <- vals
  }
  names(ssl) <- strings[1:length(ssl)]
  
  return(ssl)
}
