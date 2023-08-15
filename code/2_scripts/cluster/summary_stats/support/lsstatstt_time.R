## METADATA ===============================================================
## Description: lssttstt_time is a list of functions to be aplied to all
## though times single simulation statistics. This requires call_stats_zip.R 
## to be executed if you do not have the intermediary summary files
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-15 15:38:35
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

get_traits_colum <- function(traitsdf, traits=c("dispersal", "competition", "temp_width", "mean_temp")){
  # transform traits back to list
  n_tr <- length(traits)
  vals <- list()
  for (tri in 1:n_tr){
    # print(tri)
    vals[[tri]] <- traitsdf[, grep(traits[tri],colnames(traitsdf))]
  }
  return(vals)
}
  
trs_change <- function(traitdf){
  # mean trait change
  trcm <- rowMeans(traitdf, na.rm=T)
  # plot(trcm, type="l")
  n_tt <- length(trcm)
  time <- 1:n_tt
  lmtrcm <- lm(trcm~time)
  # lines(predict(lmtrcm, time=1:n_tt))
  slm <- summary(lmtrcm)
  # get saturation of the last 20% time-steps
  l20 <- round(0.2*n_tt,0)
  trcm20 <- tail(trcm, l20)
  t20 <- 1:l20
  lm20 <- lm(trcm20~t20)
  slm20 <- summary(lm20)
  vals <- c("slope"=lmtrcm$coefficients[2],
            "r2"=slm$r.squared,
            "slope_last20"=lm20$coefficients[2],
            "r2_last20"=slm20$r.squared)
  return(vals)
}


lssttstt_time <- list()
{
  lssttstt_time$"speciation"=function(ss=temp_ss){
    sscs <- ss$count_speciation_patch
    tb <- apply(sscs, 2, function(x){
        val1 <- sum(x, na.rm=T)
        val2 <- mean(x, na.rm=T)
        return(c("sum"=val1, "mean"=val2))
      })
    vecvals <- as.vector(tb)
    names(vecvals) <- paste0(rep(colnames(tb), each=2), "_", rownames(tb) )
    return(vecvals)
  }
  
  lssttstt_time$"range"=function(ss=temp_ss){
    sscs <- ss[grep("spatial_sps", names(ss))]
    rm <- lapply(sscs, function(x){
      val <- rowMeans(x, na.rm=T)
      return(val)
    })
    vals <- lapply(rm, function(x){
        val1 <- sum(x, na.rm=T)
        val2 <- mean(x, na.rm=T)
        return(c("sum"=val1, "mean"=val2))
      })
    vals <- unlist(vals)
    names(vals) <- gsub("\\.", "_", names(vals))
    return(vals)
  }
  
  
  lssttstt_time$"trs"=function(ss=temp_ssin, trnames=c("dispersal", "competition", "temp_width", "mean_temp")){
    # traits
    n_trs <- length(trnames)
    ref_index <- grep("mean_tr", names(ss))
    ltrs <- list()
    for (tri in 1:n_trs){
      # tri=1
      temp_listname <- grep(trnames[tri], names(ss)[ref_index], value = T)
      ltrs[[tri]] <- trs_change(ss[[temp_listname]])
    }
    names(ltrs) <- paste0(trnames)
    vals <- unlist(ltrs)
    names(vals) <- gsub("\\.", "_", names(vals1))
    return(vals)
  }
}
