## METADATA ===============================================================
## Description: Sigle statistics functions definition
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-02 20:46:39
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

########### THOUGH TIME STATS tx #########################
{
  # all this functions return objects with rows equal to time steps
  lsstatstt <- list() #MOVE TO TOP TO INCLUDE ABOVE STATISTICS!!!!!!!!!!!!!!!
  
  lsstatstt$"count_speciation_patch" <- function(tt=temp_t_i, hphy=temp_human_phy, occstt=temp_occ_tx, lc=temp_l){
    print("calling count_speciation_patch")
    # tt <- 500:0
    # getting patches names from sourced patches names! Please change this for variations
    ptns <- patches$names
    classes <- c(paste0("between_",ptns), paste0("within_", ptns) )
    df <- matrix(0, ncol=length(classes), nrow=length(tt))
    dimnames(df) <- list(tt, classes)
    
    # #ps is the phylo summary
    # ps <- sg$summary$phylo_summary
    # ps <- ps[as.character(tt),]
    # sps <- ps[,"speciations"]
    # # speciation times index
    # sptindex <- which(sps>0)
    # # e.g. ps[speciation_times,]
    # # speciation times
    # spts <- names(sptindex)
    
    #remove initial species and keep only the ones ithing the range
    #hphy <- hphy[(!hphy$Ancestor==hphy$Descendent)&(hphy$Speciation.Time%in%tt),]
    
    for (ri in which(temp_human_phy$Speciation.Time%in%temp_t_i)){
      # ri <- 1
      # get range of ancestor species prior to speciation
      
      # set previous time step.. this is the same as current in case 
      # of use of the first time-step in the landscape!
      
      #set previous time and check if it's the fisrt time
      t_a_d <- supf$get_previous_time_ancestor_descendant(temp_human_phy, ri, temp_l)
      time <- t_a_d[1]
      ancestor <- t_a_d[2]
      occs_ti <- occstt[[time]]
      patch_ti <- lc$patch[!is.na(lc$patch[,time]), c("x", "y", time)]
      patches_ancestorXtX <- occs_ti[,ancestor]*patch_ti[,3]
      patches_ancestorXtX <- patches_ancestorXtX[!patches_ancestorXtX==0]
      table_patches <- table(patches_ancestorXtX)
      l_table_patches <- length(table_patches)
      table_patches[] <- 1/l_table_patches # set marker
      names(table_patches) <- patches$names[patches$values%in%names(table_patches)]
      
      if (l_table_patches==1){
        df[as.character(temp_human_phy$Speciation.Time[ri]), paste0("within_", names(table_patches))] <- table_patches
      }else{
        df[as.character(temp_human_phy$Speciation.Time[ri]), paste0("between_", names(table_patches))] <- table_patches
      }
    }
    return(df)
  }
  
  lsstatstt$"mean_traits_sps" <- function(tt=temp_t_i, hphy=temp_human_phy, traitstt=temp_traits_tx, lc=temp_l, sg=temp_sgen3sis){
    print("calling mean_traits_sps")
    #remove initial species and keep only the ones withing the time range
    # hphy <- hphy[(!hphy$Ancestor==hphy$Descendent)&(hphy$Speciation.Time%in%tt),]# tt <- 500:0
    traitsn <- sg$parameters$gen3sis$general$trait_names
    species_numb <- unique(paste0(hphy$Ancestor,"_",hphy$Descendent))
    traits_sps <- expand.grid(traitsn, 1:sg$summary$phylo_summary["0","total"])
    classes <- paste0(traits_sps[,1],"_sp_",traits_sps[,2])
    df <- matrix(NA, ncol=length(classes), nrow=length(tt))
    dimnames(df) <- list(tt, classes)
    
    
    # #ps is the phylo summary
    # ps <- sg$summary$phylo_summary
    # ps <- ps[as.character(tt),]
    # sps <- ps[,"speciations"]
    # # speciation times index
    # sptindex <- which(sps>0)
    # # e.g. ps[speciation_times,]
    # # speciation times
    # spts <- names(sptindex)
    
    for (ti in tt){
      # ti <- 0
      # get range of ancestor species prior to speciation
      time <- as.character(ti)
      traits_ti <- traitstt[[time]]
      mean_traits_tisps <- lapply(traits_ti, function(x){(apply(x, 2,mean))})  #quantile, probs=c(0.05,0.5,0.95)
      v_mean_traits <-  unlist(mean_traits_tisps)
      temp_names <- expand.grid(traitsn, names(mean_traits_tisps))
      names(v_mean_traits) <-  paste0(temp_names[,1],"_sp_",temp_names[,2])   
      df[time, names(v_mean_traits)] <- v_mean_traits
    }
    # add trait to the ancestor species to close the phylogeny
    for (ri in which(temp_human_phy$Speciation.Time%in%temp_t_i)){ # add speciation branches
      # ri <- 1
      # get range of ancestor species prior to speciation
      t_a_d <- supf$get_previous_time_ancestor_descendant(temp_human_phy, ri, temp_l)
      time <- t_a_d[1]
      ancestor <- t_a_d[2]
      descendent <- t_a_d[3]
      df[time, colnames(df)%in%paste0(traitsn, "_sp_", descendent) ] <- df[time,paste0(traitsn, "_sp_", ancestor)]
    }
    return(df)
  }
  
  #plot(unlist(lapply(temp_traits_tx, function(x){lapply(x, function(x){if(nrow(x)!=0){return(mean(x[,"temp_width", drop=T]))}else{return(NULL)}})})))
  # unlist(lapply(temp_traits_tx['0'], function(x){if(nrow(x)!=0){return(mean(x[,"temp_width"]))}else{return(NULL)}}))
  
  lsstatstt$"spatial_sps" <- function(tt=temp_t_i, hphy=temp_human_phy, occstt=temp_occ_tx, lc=temp_l, sg=temp_sgen3sis){
    print("calling spatial_sps: i.e. colonization, no change, local extinctions")
    # prepare 3 lists for:
    #   -1 (range decrease, i.e. # of local ext)
    #   0 (range constant, i.e. # of no changes in range)
    #   +1 (range increase, i.e. # of local colonizations)
    
    temp_dummy <- matrix(NA, ncol=sg$summary$phylo_summary["0","total"], nrow=length(tt))
    dimnames(temp_dummy) <- list(tt, 1:sg$summary$phylo_summary["0","total"])
    # declare main table
    vals <- list("-1"=temp_dummy, "0"=temp_dummy, "+1"=temp_dummy)
    
    # move forward in time
    for (ti in tt[-1]){
      print(ti)
      tp <- ti+1
      octp <- occstt[[as.character(tp)]][,-c(1,2)] # get time previous
      range_octp <- colSums(octp) # get previous range for each species, i.e. number of occupied sites
      mask_present <- range_octp>0
      octi <- occstt[[as.character(ti)]][, -c(1,2)] # get time present i
      occdiff <- merge(octp,octi,all=T,by='row.names') # merge
      occdiff[is.na(occdiff)] <- 0
      occdiff <- occdiff[,grep("\\d+\\.(x|y)", colnames(occdiff))] # select only duplicated colummns
      # get before and after matching
      occdiff_tp <- occdiff[,grep("\\d+\\.(x)", colnames(occdiff))]
      occdiff_ti <- occdiff[,grep("\\d+\\.(y)", colnames(occdiff))]
      diffs <- occdiff_ti-occdiff_tp
      tables <- apply(diffs[,mask_present], 2, function(x){c("-1"=sum(x==-1), "0"=sum(x==0), "+1"=sum(x==1))})
      
      # tables <- tables
      for (stat_i in c("-1", "+1")){ # update each time and each sub-list
        temp_stats_tx <- NULL
        temp_stats_tx <- tables[stat_i,]
        vals[[stat_i]][as.character(ti),as.numeric(gsub("^([0-9]+)\\..*$", "\\1", names(temp_stats_tx)))] <- temp_stats_tx
      }
      vals[["0"]][as.character(ti),as.numeric(gsub("^([0-9]+)\\..*$", "\\1", names(temp_stats_tx)))] <- range_octp[mask_present]
    }
    return(vals)
  }
  
  
} # END though time stats