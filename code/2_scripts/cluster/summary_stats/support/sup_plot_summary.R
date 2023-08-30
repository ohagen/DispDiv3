## METADATA ===============================================================
## Description: plot summaries outputs
## 
## R version: 4.2.2 for Windows
## Date: 2023-07-24 22:27:50
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

plot_stat_classes <- function(mbt, cats="competition", y="gamma", x="dispersal", 
                              plt_type="colbar", ylab=NULL, plotblank=FALSE, ...){
  # plt type is either legend or colbar. If anything else, no leggend/bar is added
  if (is.null(ylab)){
    xlabis=x
    ylabis=y
  } else {
    xlabis=x
    ylabis=ylab
  }
  plot(NULL, xlab=xlabis, ylab=ylabis,
       xlim=range(mbt[,x], na.rm=T), ylim=range(mbt[,y], na.rm=T), bty ="n")#, ...)
  
  classes <- unique(mbt[,cats])
  n_classes <- length(classes)
  # mycol <- colorRampPalette(c("#f72585", "#b5179e", "#3a0ca3", "#4cc9f0" ))
  cols <- rev(gen3sis::color_richness_non_CVDCBP(n_classes)) # mycol(n_classes) #
  for (c_i in 1:n_classes){
    #c_i <- 1
    selection <- mbt[,cats]==classes[c_i]
    if (!plotblank){lines(x=mbt[selection,x],y=mbt[selection,y], col=cols[c_i])}
  }
  if (plt_type=="legend"){
    legend("topleft", legend = classes, col=cols, pch=3, bty="n")
  } else if (plt_type=="colbar"){
    ypos1 <- max(mbt[,y], na.rm=T)
    ypos2 <- ypos1-0.08*ypos1
    ypos3 <- ypos2-0.025*ypos2
    width_colbar <- 0.05
    length_colbar <- 0.4
    xpos <- 0.2
    colorbar.plot(x=xpos, y=ypos2, strip=1:n_classes, col = cols, 
                  strip.width=width_colbar, strip.length=length_colbar,
                  horizontal = TRUE, adj.y=0)
    text(x=xpos, y=ypos1, labels = "competition")
    text(x=xpos-(length_colbar/2.6), y=ypos3, labels = "high")
    text(x=xpos+(length_colbar/2.6), y=ypos3, labels = "neutral")
  }

  
}

plot_time_y <- function(mbtt, time=as.character(500:0), y="gamma", x="divergence_threshold", ...){
  plot(NULL, xlab=x, ylab=y,
       xlim=range(mbtt[,x], na.rm=T), ylim=range(mbtt[,time]), ...)
  
  n_time <- length(time)
  cols <- rev(heat.colors(n_time))
  # cols <- addTrans(cols,200)
  for (t_i in 1:n_time){
    #t_i <- 1
    points(x=jitter(mbtt[,x]),y=jitter(mbtt[,time[t_i]]), col=cols[t_i], pch=4)
  }
  # legend("topleft", legend = classes, col=cols, pch=3, bty="n")
}



plot_trait_phylogeny <- function(traitsdf=df, trait="competition"){
  
  traitsdf <- traitsdf[, grep(trait,colnames(traitsdf))]
  tts <- as.numeric(rownames(traitsdf))
  tts_seq_index <- 1:length(tts)
  nspp <- ncol(traitsdf)
  #ylim=range(traitsdf, na.rm=T)
  plot(tts,  ylim=c(0,1), col=rgb(0,0,0,1), xaxt='n', pch = '', xlab="Ma", ylab=trait, 
       main="Trait evolution")
  axis_lab <- seq(from = 1, to = length(traitsdf[,1]), by = 50)
  axis(1, at = axis_lab, 
       labels = round((tts[axis_lab]/100),2))
  

  cols=rainbow(nspp)
  
  for (line_i in 1:nspp){
    # line_i <- 1
    #no alpha
    lines(traitsdf[,line_i], col="black")#cols[line_i])
    # segments(0:(nrow(traitsdf)-1), head(traitsdf[,line_i], -1), 1:nrow(traitsdf), traitsdf[,line_i][-1])#, color_segmentes[r[,line_i]])
  }  # 
  #   selct <- traitsdf[,line_i]
  #       
  #   segments(0:(nrow(traitsdf)-1), head(traitsdf[,line_i], -1), 1:nrow(traitsdf), traitsdf[,line_i][-1])#, color_segmentes[r[,line_i]])
  # }
  #   #with aplha
  #   alp <- 80# [0,255]
  #   ri <- r[1:time_i,line_i,drop=F]
  #   ri[is.na(ri)] <- 6
  #   
  #   segments(0:(nrow(ti)-1), head(ti[,line_i], -1), 1:nrow(ti), ti[,line_i][-1]
  #            , col=rgb(
  #              veccol[1,ri], 
  #              veccol[2,ri], 
  #              veccol[3,ri], 
  #              alpha=rep(alp,3), 
  #              maxColorValue = 255)
  #   )
  # }
  # legend("bottomleft", new_legend[c(5,1,2,3,4)], col=color_segmentes[c(5,1,2,3,4)], pch=15, bty='n')
  # abline(v=time_i)
  
}

plot_space_changes <- function(space=ssltt[[1]]$spatial_sps$`+1`,  na.val=0, ...){
  spdx <- rowMeans(space, na.rm=T)
  spdx[is.na(spdx)] <- 0
  plot(spdx, type='l', ...)
}

plot_stat_classes_p <- function(mbt=lt[[3]], 
                                y = "gamma", 
                                x="range_spatial_sps_0_mean",
                                cats="trs.competition_50%", 
                                plt_type=NULL,
                                xposbar=NULL,
                                yposbar=NULL,
                                ...){
  classes <- unique(mbt[,cats])
  n_classes <- length(classes)
  cols <- rev(gen3sis::color_richness_non_CVDCBP(n_classes))
  #mycol <- colorRampPalette(c("#f72585", "#b5179e", "#3a0ca3", "#4cc9f0" ))
  #cols <- mycol(n_classes)
  plot(x=mbt[,x],y=mbt[,y], col=cols[as.factor(mbt[,cats])], bty ="n", pch=16, cex=0.5, ...)
  if (plt_type=="legend"){
    legend("topleft", legend = classes, col=cols, pch=3, bty="n")
  } else if (plt_type=="colbar"){
    ypos1 <- if (is.null(yposbar)){0.8*max(mbt[,y])}else{yposbar}
    ypos2 <- ypos1-0.08*ypos1
    ypos3 <- ypos2-0.025*ypos2
    width_colbar <- 0.05
    length_colbar <- 0.3
    xpos <- if (is.null(xposbar)){0.8*max(mbt[,x])}else{yposbar}
    
    colorbar.plot(x=xpos, y=ypos2, strip=1:n_classes, col = cols, 
                  strip.width=width_colbar, strip.length=length_colbar,
                  horizontal = TRUE, adj.y=0)
    text(x=xpos, y=ypos1, labels = "competition")
    text(x=xpos-(length_colbar*xpos), y=ypos3, labels = "high")
    text(x=xpos+(length_colbar*xpos), y=ypos3, labels = "neutral")
  }
}

