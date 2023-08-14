## METADATA ===============================================================
## Description: Support functions for the lsstats and lsstatstt
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-02 20:51:55
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

# metric <- pd.query(tree=phy, matrix=f_sp_matrix, standardize = FALSE, 
#                    null.model="uniform", abundance.weights, reps=1000, seed)


#### Support functions for the summary stats ############
{  
  supf <- list()
  supf$"get_previous_time_ancestor_descendant" <- function(hphy=temp_human_phy, ri=XXX, lc=temp_l){
    # start time is a string
    
    # returns time _ancestor_descendan as names vector
    tti=hphy$Speciation.Time[ri]
    descendent <- as.character(hphy$Descendent[ri])
    #tti is a numetic time,
    
    # get start time from landscape
    start_time <- as.numeric(tail(colnames(lc[[1]]),1))
    if (tti==start_time){
      time <- as.character(tti)
      ancestor <- descendent
    } else{
      time <- as.character(tti+1)
      ancestor <- as.character(hphy$Ancestor[ri])
    }
    return(c("time"=time, "ancestor"=ancestor, "descendent"=descendent))
  }  
  
  supf$"get_patches_tx" <- function(lc=temp_l, tx=temp_t_i){
    #patch_mask <- pt_x==pt_ti
    pt_ti <- lc$patch[!is.na(lc$patch[,as.character(tx)]), as.character(tx)]
    return(pt_ti)
  }
  
  supf$"get_patch_mask"  <- function(pt_ti, pt_x="T"){
    # patch_x is either "T" for all sites or has the value of the patch
    if (pt_x=="T"){
      patch_mask <- !is.na(pt_ti)
    } else {
      # taking data from the source file!
      pt_x <- patches$values[patches$names==pt_x]
      patch_mask <- pt_x==pt_ti
    }
    return(patch_mask)
  }
  
  supf$"create_sublists" <- function(all_pts) {
    sublists <- lapply(all_pts, function(x) NA)
    names(sublists) <- all_pts
    return(sublists)
  }
  
  supf$"comp_alpha_p_" <- function(occs_pt, patch_number=ptnb, patch_names=ptns, gamma){
    # proportional
    #Zeta: species that are common to all assemblages; Zetai corresponds to Zeta diversity sensu Hui and McGeoch (2014)
    #[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8525081/] and DOI: 10.1002/ece3.8096
    # Eta: species that are not unique to an assemblage nor common to all assemblages.
    # Theta: species that are unique to an assemblage.
    
    zeta <- sum(colSums(occs_pt)==length(patch_number))
    names(zeta) <- "zeta_T"
    eta <- sum((colSums(occs_pt)<length(patch_number))&(colSums(occs_pt)>1))
    names(eta) <- "eta_T"
    unique <- colSums(occs_pt)==1
    theta <- rowSums(occs_pt[, unique, drop=F])
    names(theta) <- paste0("theta_", patch_names)
    
    return(c(zeta, eta, theta)/gamma)
  }
  
  supf$"PD_" <- function(occs_pt, prunedphy){
    vals <- pd(occs_pt, prunedphy, include.root=T)
    valsf <- vals[,1]
    names(valsf) <- paste0("PD_", rownames(vals))
    return(vals)
  }
  
  supf$"PD_" <- function(occs_pt, prunedphy){
    vals <- pd(occs_pt, prunedphy, include.root=F)
    valsf <- vals[,1]
    names(valsf) <- paste0("PD_", rownames(vals))
    return(valsf)
  }
  
  supf$"phylo_metrics_" <- function(occs_pt, prunedphy){
    pd_estimate <- pd.query(prunedphy, occs_pt, standardize = T)
    mpd_estimate <- mpd.query(prunedphy, occs_pt, standardize = T)
    mntd_estimate <- mntd.query(prunedphy, occs_pt, standardize = T)
    # get patches names
    ptnames <- rownames(occs_pt)
    names(pd_estimate) <- paste0("PD_S_", ptnames)
    names(mpd_estimate) <- paste0("MPD_S_", ptnames)
    names(mntd_estimate) <- paste0("MNTD_S_", ptnames)
    
    return(c(pd_estimate, mpd_estimate, mntd_estimate))
  }
  
  
  supf$"phylo_community_distance_" <- function(occs_pt, prunedphy){
    # Computes the (standardized) value of the Community Distance measure
    # is the beta diversity version of Mean Pairwise Distance (MPD), giving the average phylogenetic distance between two communities.
    pairs <- cbind(c(1,1,1,2,2,3), c(2,3,4,3,4,4))
    vals <- cd.query(prunedphy, occs_pt, standardize = T,query.matrix =pairs)
    names(vals) <- paste0("CD_S_", paste0(LETTERS[pairs[,1]], LETTERS[pairs[,2]]))
    return(vals)
  }
  
  
}

{  #### Support functions for the summary stats that LOOP OVER TEMPORAL GRAIN ############
  pt_supf <- list()
  
  pt_supf$"mean_alpha_" <- function(alpha, pt_mask){
    # alpha is the full vector or alphas
    # ptvec is the full vector of patch levels
    vals <- mean(alpha[pt_mask])
    return(vals)
  }
  pt_supf$"gamma_" <- function(occs=occst, pt_mask){
    # occs presense/absence matrix siteXspecies
    # alpha is the full vector or alphas
    # ptvec is the full vector of patch levels
    vals <- sum(colSums(occs[pt_mask,])>0)
    return(vals)
  }
  pt_supf$"beta_prop_" <- function(mean_alpha, gamma){
    # Proportional species turnover:  
    #quantifies what proportion of the species diversity in the dataset is not 
    #contained in an average subunit
    vals <- 1-(mean_alpha/gamma)
    return(vals)
  }
  pt_supf$"beta_w_" <- function(mean_alpha, gamma){
    # Whittaker beta div:
    #how many subunits there would be if the total species diversity of the 
    #dataset and the mean species diversity per subunit remained the same, 
    # but the subunits shared no species
    vals <- gamma/mean_alpha
    return(vals)
  }

}

#### RANDOM
# though time statistics
get_gamma_tt <- function(sg=sgen3sis){
  g_ti <- sg$summary$phylo_summary[,"alive"]
  return(g_ti)
}

get_gamma <- function(ss, time=0){
  g_ti <- ss$summary$phylo_summary[as.character(time),"alive"]
  names(g_ti) <- paste0("gamma_t",time)
  return(g_ti)
}