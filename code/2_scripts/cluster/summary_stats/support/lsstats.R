## METADATA ===============================================================
## Description: Sigle statistics functions definition
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-02 20:46:39
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

########## SINGLE TIME STATS tx #########################
{
  lsstats <- list() # single time stats!
  
  lsstats$"finished" <- function(sg=temp_sgen3sis){
    return(sg$flag)  
  }
  
  lsstats$"cpu_time" <- function(sg=temp_sgen3sis){
    return(sg$system$`runtime-hours`)
  }
  
  lsstats$"extinctions_perc" <- function(sg=temp_sgen3sis){
    #ps is the phylo summary
    ps <- sg$summary$phylo_summary
    # lt is the last time, or last row of the ps table
    lt <- as.data.frame(ps[nrow(ps),,drop=F])
    # val is the percentage of extinction events
    val <- round((( lt$"total"- lt$"alive")/ lt$"total"),3)
    return(val)
  }
  
  lsstats$"speciations_perc" <- function(sg=temp_sgen3sis){
    #ps is the phylo summary
    ps <- sg$summary$phylo_summary
    # ft is the first time, meaning first row of the ps table
    ft <- as.data.frame(ps[1,,drop=F])
    # lt is the last time, meaning first row of the ps table
    lt <- as.data.frame(ps[nrow(ps),,drop=F])
    # val is the percentage of extinction events 
    val <- round((( lt$"total"- ft$"total")/ lt$"total"),3)
    return(val)
  }
  
  
  lsstats$"n_sp_alive_t_0" <- function(sg=temp_sgen3sis){
    n_sp_alive <- tail(sg$summary$phylo_summary[,"alive"],1)
    return(unname(n_sp_alive))
  }
  
  lsstats$"turnover" <- function(sg=temp_sgen3sis){
    n_sp_init <- sg$summary$phylo_summary["initial","alive"]
    turn <- sg$summary$phylo_summary[-1,c("speciations", "extinctions")]
    turn <- turn[,1]-turn[,2]
    val <- sum(turn)/n_sp_init #percentage based on the initial abundance
    return(unname(val))
  }
  
  
  lsstats$"speciations_invar" <- function(sg=temp_sgen3sis){
    # Time
    #ps is the phylo summary
    ps <- sg$summary$phylo_summary
    Time <- nrow(ps)-1
    # Temporal speciation invariability: Temporal stability of speciation events
    speciation_invar <- (sum(ps[-1,"speciations"])/Time)/sd(ps[-1,"speciations"])
    val <- round(speciation_invar,3)
    return(val)
  }
  
  
  lsstats$"extinction_invar" <- function(sg=temp_sgen3sis){
    #ps is the phylo summary
    ps <- sg$summary$phylo_summary
    Time <- nrow(ps)-1
    # Temporal speciation invariability: Temporal stability of extinction events
    extinction_invar <- (sum(ps[-1,"extinctions"])/Time)/sd(ps[-1,"extinctions"])
    val <- round(extinction_invar,3)
    return(val)
  }
  
  lsstats$"turnover_invar" <- function(sg=temp_sgen3sis){
    #ps is the phylo summary
    ps <- sg$summary$phylo_summary
    sps <- ps[-1,"speciations"]
    exs <- ps[-1, "extinctions"]
    Time <- nrow(ps)-1
    turnover <- sps-exs
    # Temporal speciation invariability: Temporal stability of extinction events
    turnover_invar <- (sum(turnover)/Time)/sd(turnover)
    val <- round(turnover_invar,3)
    return(val)
  }
  
  lsstats$"gamma" <- function(sg=temp_sgen3sis){
    n_sp_alive <- tail(sg$summary$phylo_summary[,"alive"],1)
    return(unname(n_sp_alive))
  }
  
  lsstats$"mean_abd" <- function(abd=temp_abundance_tx){
    sabds <- unlist(lapply(abd, sum))
    qabd <- quantile(sabds, c(0.05, 0.5, 0.95))
    return(qabd)
  }
  
  lsstats$"mtx" <- function(sg=temp_sgen3sis, phylo=temp_phy, occst=temp_occ_tx, lc=temp_l, tx=temp_t_i){
    # multiple stats here that are , returning a variable large vector of summary statistics
    # remove x y coordinates
    occst <- occst[,-c(1,2)]
    # getting patches names from sourced patches names! Please change this for variations
    ptns <- patches$names
    # getting patches correspondence to numeric
    ptnb <- patches$values
    
    pt_ti <- supf$get_patches_tx(lc=lc, tx=tx) 
    # pt_ti <- lc$patch[!is.na(lc$patch[,as.character(tx)]), as.character(tx)]
    
    # regional scale
    alpha <- rowSums(occst)
    
    all_pts <- c("T", ptns)
    psl <- pt_supf # get patch functions temporarily to use it's names
    psl <- lapply(psl, function(x){
      supf$create_sublists(all_pts) #add empty sublists
      }) 
    # call here each statistics per patch, including T 
    for (pt_x in all_pts){ # for each patch
      patch_mask <- supf$"get_patch_mask"(pt_ti, pt_x)
      psl$mean_alpha_[[pt_x]] <- pt_supf$mean_alpha_(alpha, pt_mask=patch_mask)
      psl$gamma_[[pt_x]] <- pt_supf$gamma_(occst, pt_mask=patch_mask)
      psl$beta_prop_[[pt_x]] <- pt_supf$beta_prop_(mean_alpha=psl$mean_alpha_[[pt_x]], gamma=psl$gamma_[[pt_x]])
      psl$beta_w_[[pt_x]] <- pt_supf$beta_w_(mean_alpha=psl$mean_alpha_[[pt_x]], gamma=psl$gamma_[[pt_x]])
      # psl$comp_alpha_[[pt_x]] <- pt_supf$comp_alpha_(occst, pt_x=pt_x,  pt_ti=pt_ti, gamma=psl$gamma_[["T"]])
    }
    
    vals <- unlist(psl)
    names(vals) <- gsub("\\.", "", names(vals)) # remove points
    
    # call here the stats that use the prepared occst, etc...
    # prepare the occs_pt objects for calculation
    occs_pt <- NULL
    for (cati in 1:length(ptnb)){
      occs_pt <- rbind(occs_pt, as.numeric(as.logical(colSums(occst[pt_ti==ptnb[cati], ]))))
    }
    occs_pt <- rbind(occs_pt, as.numeric(as.logical(colSums(occs_pt)))) # add total col
    dimnames(occs_pt) <- list(c(ptns, "T"), paste0("species",colnames(occst)))
    
    prunedphy <- prune.sample(occs_pt, phylo) # just to make sure! Should allways hold
    occs_pt <- occs_pt[, prunedphy$tip.label] # get same order as phylo
   
    vals2 <- list()
    # call stats for simplified object occs
    vals2[[1]] <- supf$"comp_alpha_p_"(occs_pt=occs_pt, patch_number=ptnb, patch_names=ptns,  gamma=psl$gamma_[["T"]])
    #vals2[[2]] <- supf$"PD_"(occs_pt, prunedphy)
    vals2[[2]] <- supf$"phylo_metrics_"(occs_pt, prunedphy)
    vals2[[3]] <- supf$"community_distance_"(occs_pt, prunedphy)
    
    # get phylogenetic distance matrix
    return(c(vals, unlist(vals2)))
  }
  
  lsstats$"func" <- function(occst=temp_occ_tx, abd=temp_abundance_tx, trs=temp_traits_tx, lc=temp_l, tx=temp_t_i){
    ptns <- patches$names
    # getting patches correspondence to numeric
    ptnb <- patches$values
    pt_ti <- supf$get_patches_tx(lc=lc, tx=tx) 
    occst <- occst[,!colnames(occst)%in%c("x", "y")] # remove, x and y
    # prepare abundance matrix
    for (sp_i in 1:length(abd)){
      occst[occst[,sp_i]==1,sp_i] <- abd[[sp_i]]
    }
    
    abd_pt <- NULL
    for (cati in 1:length(ptnb)){
      abd_pt <- rbind(abd_pt, colSums(occst[pt_ti==ptnb[cati], ]))
    }
    dimnames(abd_pt) <- list(ptns, paste0("",colnames(occst)))
    
    
    # 
    # # getting patches names from sourced patches names! Please change this for variations
    # ptns <- patches$names
    # # getting patches correspondence to numeric
    # ptnb <- patches$values
    # 
    # pt_ti <- supf$get_patches_tx(lc=lc, tx=tx)
    
    
    # prepare traits
    mean_trs <- lapply(trs, colMeans)
    mean_trs <- do.call(rbind, mean_trs)
    
    mask_present <- colSums(abd_pt)>0
    mask_commun <- rowSums(abd_pt)>0
    
    mask_trs_var <- apply(mean_trs,2,function(x){
      rgt <- range(x, na.rm=T)
      diff <- rgt[2]-rgt[1]
      return(diff>0.00001)
    })
    
    stats_FD <- c("FRic", "FEve", "FDiv", "FDis", "RaoQ")
    
    suppressWarnings(vals <- tryCatch(
      FD::dbFD(mean_trs[mask_present,mask_trs_var], abd_pt[mask_commun,mask_present], messages = FALSE)[c("FRic", "FEve", "FDiv", "FDis", "RaoQ")]
                                      , error=function(err) rep(NA,5))) 
    
    names(vals) <- stats_FD # refine results and bypass errors in dbFD
    
    # search for NULLS
    vals <- lapply(vals, function(x){
      if (is.null(x)){
        r <- rep(NA, sum(mask_commun))
      }else{
        r <- x
      }
      return(r)})
    
    # making sure we return allays the same length of output
    dummy <- mask_commun
    dummy[] <- 0
    
    vals <- lapply(vals, function(x){dummy[mask_commun] <- x
    return(dummy)})
    
    
    vals <- unlist(vals)
    names(vals) <- gsub("\\.", "_", names(vals)) # remove points
    
    return(vals)
  }
  
  
  
  lsstats$"trs" <- function(trs=temp_traits_tx){
    # modes
    # get the number of nodes in the data
    trsall <- do.call(rbind, trs)
    modes <- apply(trsall, 2, function(x){Modes(x)$modes}, simplify=FALSE)
    vals_mode <- unlist(lapply(modes, function(x){length(unique(round(x,4) ))}))
    names(vals_mode) <- paste0("modes_", names(vals_mode))
    
    # quantile
    quant <- apply(trsall, 2, function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}, simplify = FALSE)
    vals_quant <- unlist(quant)
    names(vals_quant) <- gsub("\\.", "_", names(vals_quant)) # remove points
    
    # spread
    vals_spread <- apply(trsall, 2,  function(x){
      rg <- range(round(x,4))
      return(rg[2]-rg[1])
    }, simplify = FALSE)
    names(vals_spread) <- paste0("spread_", names(vals_spread))
    
    return(c(vals_mode, vals_quant, vals_spread))
  }
  
  
  
  lsstats$"maxlik_betasplit" <- function(phylo=temp_phy){
    #phylo <- unroot(phylo)
    # get the beta splits for full phylogeny (TF), pruned, (TP), and pruned species 1,2 and 3
    bsplit <- rep(NA,5)
    names(bsplit) <- c("TF", "TP", 1:3)
    tl <- list()
    if (length(phylo$tip.label)>10){
      bsplit["TF"] <- maxlik.betasplit(phylo,confidence.interval="none")$max_lik
      bsplit["TP"] <- maxlik.betasplit(drop.fossil(phylo), confidence.interval="none")$max_lik
      
      list_sis_sps <- list()
      par(mfrow=c(1,3))
      for (sps_i in 1:3){
        # sps_i <- 1
        nms <- getMommy(phylo, N=paste0("species", sps_i))
        n_nms <- length(nms)
        keepit <- rep(FALSE, n_nms)
        sizesis <- rep(NA, n_nms)
        vtips <- list()
        for (nodei in 1:n_nms){
          vtips[[nodei]] <- tips(phylo, nms[nodei])
          sizesis[nodei] <- length(vtips[[nodei]])
          keep <- c("species1", "species2", "species3")%in%vtips[[nodei]]
          # print(keep)
          if ((sum(keep)==1)&(keep[sps_i]==TRUE)&sizesis[nodei]>=10){
            keepit[nodei] <- TRUE
          }
          
        }
        if((sum(keepit)>0)){
          index_f <- (1:length(keepit))[keepit][which.max(sizesis[keepit])]
          #nodef <- nms[index_f]
          vtipsf <- vtips[[index_f]]
          list_sis_sps[[sps_i]] <- keep.tip(phylo, vtipsf)
        }else{
          list_sis_sps[[sps_i]] <- NA
        }
      } # end loop over species
      # run stats
      for (spi in 1:3){
        if (!is.na(list_sis_sps[[spi]][1])){
          bsplit[as.character(spi)] <- tryCatch(
          maxlik.betasplit(drop.fossil(list_sis_sps[[spi]]), confidence.interval="none")$max_lik
          , error=function(err) NA)
        }
      }
    }
    return(bsplit)
  }
  
    lsstats$"maxlik_betasplit" <- function(phylo=temp_phy){
    #phylo <- unroot(phylo)
    # get the beta splits for full phylogeny (TF), pruned, (TP), and pruned species 1,2 and 3
    bsplit <- rep(NA,5)
    names(bsplit) <- c("TF", "TP", 1:3)
    tl <- list()
    if (length(phylo$tip.label)>10){
      bsplit["TF"] <- maxlik.betasplit(phylo,confidence.interval="none", up=30)$max_lik
      bsplit["TP"] <- maxlik.betasplit(drop.fossil(phylo), confidence.interval="none")$max_lik
      
      list_sis_sps <- list()
      par(mfrow=c(1,3))
      for (sps_i in 1:3){
        # sps_i <- 1
        nms <- getMommy(phylo, N=paste0("species", sps_i))
        n_nms <- length(nms)
        keepit <- rep(FALSE, n_nms)
        sizesis <- rep(NA, n_nms)
        vtips <- list()
        for (nodei in 1:n_nms){
          vtips[[nodei]] <- tips(phylo, nms[nodei])
          sizesis[nodei] <- length(vtips[[nodei]])
          keep <- c("species1", "species2", "species3")%in%vtips[[nodei]]
          # print(keep)
          if ((sum(keep)==1)&(keep[sps_i]==TRUE)&sizesis[nodei]>=10){
            keepit[nodei] <- TRUE
          }
          
        }
        if((sum(keepit)>0)){
          index_f <- (1:length(keepit))[keepit][which.max(sizesis[keepit])]
          #nodef <- nms[index_f]
          vtipsf <- vtips[[index_f]]
          list_sis_sps[[sps_i]] <- keep.tip(phylo, vtipsf)
        }else{
          list_sis_sps[[sps_i]] <- NA
        }
      } # end loop over species
      # run stats
      for (spi in 1:3){
        if (!is.na(list_sis_sps[[spi]][1])){
          bsplit[as.character(spi)] <- tryCatch(
          maxlik.betasplit(drop.fossil(list_sis_sps[[spi]]), confidence.interval="none")$max_lik
          , error=function(err) NA)
        }
      }
    }
    return(bsplit)
  }
  
  
  # 
  # # WIP#
  # lsstats$"trs_evol" <- function(phylo=temp_phy, trs=temp_traits_tx, abd=temp_abundance_tx){
  #   #### BROWSER ! ----------
  #   browser()
  #   n_sp <- length(trs)
  #   trsv <- rep(NA, n_sp)
  #   names(trsv) <- paste0("species",names(trs))
  #   for (spi in 1:n_sp){
  #     # spi <- 1
  #     mean_abd <- mean(abd[[spi]])
  #     weight_abd <- /mean_abd
  #     for (ti in trn){
  #       traits[cells_cluster, ti] <- mean(traits[cells_cluster, ti]*weight_abd)
  #     }
  #   }
  #   
  #   
  #   
  #   
  #   
  #   
  #   
  #   trsm <- lapply(trs, function(mean(x$dispersal)))
  #   
  #   #phylo <- unroot(phylo)
  #   # get the beta splits for full phylogeny (TF), pruned, (TP), and pruned species 1,2 and 3
  #   bsplit <- rep(NA,5)
  #   names(bsplit) <- c("TF", "TP", 1:3)
  #   tl <- list()
  #   if (length(phylo$tip.label)>10){
  #     bsplit["TF"] <- maxlik.betasplit(phylo,confidence.interval="none")$max_lik
  #     bsplit["TP"] <- maxlik.betasplit(drop.fossil(phylo), confidence.interval="none")$max_lik
  #     
  #     list_sis_sps <- list()
  #     par(mfrow=c(1,3))
  #     for (sps_i in 1:3){
  #       # sps_i <- 1
  #       nms <- getMommy(phylo, N=paste0("species", sps_i))
  #       n_nms <- length(nms)
  #       keepit <- rep(FALSE, n_nms)
  #       sizesis <- rep(NA, n_nms)
  #       vtips <- list()
  #       for (nodei in 1:n_nms){
  #         vtips[[nodei]] <- tips(phylo, nms[nodei])
  #         sizesis[nodei] <- length(vtips[[nodei]])
  #         keep <- c("species1", "species2", "species3")%in%vtips[[nodei]]
  #         # print(keep)
  #         if ((sum(keep)==1)&(keep[sps_i]==TRUE)&sizesis[nodei]>=10){
  #           keepit[nodei] <- TRUE
  #         }
  #         
  #       }
  #       if((sum(keepit)>0)){
  #         index_f <- (1:length(keepit))[keepit][which.max(sizesis[keepit])]
  #         #nodef <- nms[index_f]
  #         vtipsf <- vtips[[index_f]]
  #         list_sis_sps[[sps_i]] <- keep.tip(phylo, vtipsf)
  #       }else{
  #         list_sis_sps[[sps_i]] <- NA
  #       }
  #     } # end loop over species
  #     # run stats
  #     for (spi in 1:3){
  #       if (!is.na(list_sis_sps[[spi]][1])){
  #         bsplit[as.character(spi)] <- tryCatch(
  #           maxlik.betasplit(drop.fossil(list_sis_sps[[spi]]), confidence.interval="none")$max_lik
  #           , error=function(err) NA)
  #       }
  #     }
  #   }
  #   return(bsplit)
  # }
  # # END WIP
  # 
  
  
  
} # end fix time statistics