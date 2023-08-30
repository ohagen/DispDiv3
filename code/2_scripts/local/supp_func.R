# lib <- c("raster","fields")
# sapply(lib, require, character.only = TRUE, quietly = TRUE, warn.conflicts = TRUE)
# source(file.path(pls$dir_base,"summary_stats/support/sup_plot_summary.R"))
# source(file.path(pls$dir_base,"summary_stats/support/sup_plot_traits.R"))
# source(file.path(pls$dir_base,"summary_stats/support/call_stats_zip.R"))
# # source(file.path(pls$dir_base,"summary_stats/support/sup_sumstats.R"))
# # source(file.path(pls$dir_base,"summary_stats/support/lsstats.R"))
# # source(file.path(pls$dir_base,"summary_stats/support/lsstatstt.R"))
# lib_ss <- c("tools", "dplyr", "ape", "betapart", "apTreeshape", "picante", "FD", "LaplacesDemon", "RRphylo", "phytools")
# # sapply(lib, require, character.only = TRUE, quietly = TRUE, warn.conflicts = TRUE)
# 

set_par <- function(n_plots, ncols=2){
  par(mfrow=c(ceiling(n_plots/ncols), ncols), mai=c(1,1,1,0.5)*0.3, mgp=c(1.1,0.1,0), tcl = -0.15)
}


plot_lxy_summary <- function(x="dispersal", y="competition", data=lpst, datacol="cpu_time", 
                            absolute=TRUE, colsrap=c("magenta","purple3", "red4"), ramp=FALSE, breaks=10){
  # ramp is either false or a vector of two c(x,y) coordinates
  n_expp <- length(data)
  par(mfrow=c(1, n_expp), mai=c(1,1,1,0.5)*0.3, mgp=c(1,0.2,0), tcl =-0.2)
  cput <- lapply(data, function(x){
    x[datacol]
  })
  
  breaks_cput <- breaks
  
  if (absolute){
    brks <- as.numeric(cut(unlist(cput), breaks = breaks_cput))
    #levels(brks) <- 1:length(levels(brks))
    lbrks <- split(brks,
                   cut(seq_along(brks), length(sss), labels=FALSE))
  }else{
    brks <- lapply(cput, function(x){
      as.numeric(cut(x, breaks = breaks_cput))
    })
    lbrks <- brks
  }

  # lcols <- lapply(lcols, function(x){unlist(x)})
  color_cput=colorRampPalette(colsrap)(breaks_cput)
  for (i in 1:n_expp){
    if (!absolute){
      # add to names the range
      #lapply()
    }
    plot(data[[i]][,x], 
         data[[i]][,y], 
         pch=4, cex=0.4, 
         col=color_cput[lbrks[[i]]], 
         ylab="Competition (l)", 
         xlab="Dispersal (d)", 
         main=names(data)[i])
    title(LETTERS[i], adj=0)
    if (i==3){ # plot colbar
      xpos <- 0.79
      ypos1 <- 0.985
      ypos3 <- 0.975
      colorbar.plot(x=xpos, y=mean(c(ypos1, ypos3)), strip=1:length(color_cput), strip.width = 0.04, strip.length = 0.15, col = color_cput)
      text(x=xpos, y=ypos1, labels = expression('CPU'[time]*' (h)'))
      if (absolute){
        text(x=xpos-0.17, y=ypos3, labels = round(min(unlist(cput)),2))
        text(x=xpos+0.17, y=ypos3, labels = round(max(unlist(cput)),2))
      }
    }
  }
}


add_to_plot_topleft <- function(text=LETTERS[stat_i]){
  mtext(text, 2, adj=2, las=1, padj=-(31-6.7*par()$mfrow[1]))
}

put_brak <- function(x){
  # x is a vector of value to make a range
  str <- get_timesteps_in_Ma(x, convfact = 1, unit = "", space='', rev=F)
  val <- gsub("-", "--", paste0("[",str,"]"))
  return(val)
}


plot_stat_classes_summary <- function(mbt, stats_names, colbar.at=1, limit_val=3, limit_stat=NULL){
  n_stats <- length(stats_names)
  set_par(n_stats)
  for (stat_i in 1:n_stats){
    # stat_i <- 1
    if (any(is.null(limit_val), is.null(limit_stat))){
      mask_mbt <- rep(T,nrow(mbt))
    }else{
      mask_mbt <- mbt[limit_stat]>=limit_val
    }
    
    plot_stat_classes(mbt[mask_mbt,], cats="competition", y=stats_names[stat_i], x="dispersal", plt_type="FALSE")
    title(LETTERS[stat_i], adj=0)
    #add_to_plot_topleft(LETTERS[stat_i])
    if (stat_i==colbar.at){ # plot colbar
      classes <- unique(mbt[,"competition"])
      n_classes <- length(classes)
      #mycol <- colorRampPalette(c("#f72585", "#b5179e", "#3a0ca3", "#4cc9f0" ))
      #cols <- mycol(n_classes)
      cols <- rev(gen3sis::color_richness_non_CVDCBP(n_classes))
      ypos1 <- max(mbt[,stats_names[stat_i]], na.rm=T)
      
      
      width_colbar <- 0.08
      length_colbar <- 0.8
      ypos3 <- ypos1-150-width_colbar
      xpos <- 0.2
      colorbar.plot(x=xpos, y=mean(c(ypos1, ypos3)), strip=1:n_classes, col = cols, 
                    strip.width=width_colbar, strip.length=length_colbar,
                    horizontal = TRUE, adj.y=0.5)
      text(x=xpos, y=ypos1, labels = "Competition (l)")
      text(x=xpos-0.1, y=ypos3, labels = "0.9 \n *")
      text(x=xpos+0.1, y=ypos3, labels = "1 \n neutral", adj=0)
    }
    # mtext(LETTERS[stat_i], side=3, line=-1, text="Here again?", adj=0, outer=T)
    
  }
}







