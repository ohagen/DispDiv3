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
                              plt_type="colbar", ylab=NULL, xlab=NULL, plotblank=FALSE,
                              ylim=NULL, ...){
  # plt type is either legend or colbar. If anything else, no leggend/bar is added
  if (is.null(ylab)){
    xlabis=x
    ylabis=y
  } else {
    xlabis=xlab
    ylabis=ylab
  }
  if (is.null(ylim)){
    ylimbis <- range(mbt[,y], na.rm=T)
  } else {
    ylimbis <- ylim
  }
  
  plot(NULL, xlab=xlabis, ylab=ylabis,
       xlim=range(mbt[,x], na.rm=T), ylim=ylimbis, bty ="n")#, ...)
  
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



plot_trait_phylogeny <- function(traitsdf=df, trait="competition", type="absolute", maint="", ylab=NULL){

  # type can be absolute,  relative or a given min\max for the y axis
  traitsdf <- traitsdf[, grep(trait,colnames(traitsdf))]
  tts <- as.numeric(rownames(traitsdf))
  tts_seq_index <- 1:length(tts)
  nspp <- ncol(traitsdf)
  #ylim=range(traitsdf, na.rm=T)
  
  if (length(type)==2&is.numeric(type)){
    yrangecust <- type
  } else{
    yrangecust <- if(type=="absolute"){c(0,1)}else{range(traitsdf, na.rm = T)}
  }

  plot(tts,  ylim=yrangecust, col=rgb(0,0,0,1), xaxt='n', pch = '', xlab="Ma", ylab=if(is.null(ylab)){trait}else{ylab}, 
       main=maint)
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
                                cats="trs_competition_50%", 
                                plt_type=NULL,
                                xposbar=NULL,
                                yposbar=NULL,
                                cex_p=0.5,
                                pch_p=16,
                                mycols=NULL,
                                labcust=list("Competition", "high", "zero"),
                                ...){
  classes <- unique(mbt[,cats])
  n_classes <- length(classes)
  if (is.null(mycols)){
    cols <- rev(gen3sis::color_richness_non_CVDCBP(n_classes))
  } else{
    #mycols <- colorRampPalette(c("#f72585", "#b5179e", "#3a0ca3", "#4cc9f0" ))
    cols <- mycols(n_classes)
  }
  
  #mycol <- colorRampPalette(c("#f72585", "#b5179e", "#3a0ca3", "#4cc9f0" ))
  #cols <- mycol(n_classes)
  plot(x=mbt[,x],y=mbt[,y], col=cols[as.factor(mbt[,cats])], bty ="n", pch=pch_p, cex=cex_p, ...)
  if (plt_type=="legend"){
    legend("topleft", legend = classes, col=cols, pch=pch_p, bty="n")
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
    text(x=xpos, y=ypos1, labels = labcust[[1]])
    text(x=xpos-(length_colbar*xpos), y=ypos3, labels = labcust[[2]])
    text(x=xpos+(length_colbar*xpos), y=ypos3, labels = labcust[[3]])
  }
}

# stats_symbol_lib <- list(
#   "log10_mean_aplha"= bquote(log[10](bar(alpha))),
#   "trs_dispersal_50%"=expression("Mean Dispersal  ("~bar("d"['t'['0']])~")"),
#   "mtx_beta_prop_T"= expression(beta~'%'),
#   "gamma"=  expression(gamma),
#   "mtx_zeta_T"=  expression("Zeta "~zeta),
#   "mtx_eta_T"=  expression("Eta "~eta),
#   "mtx_MNTD_S_T"=  expression('MNTD'['S']),
#   "speciations_perc"=  "Speciation %",
#   "extinctions_perc"=   "Extinction %",
#   "mean_abd_50%"=   bquote(bar(B)),
#   "mean_range"= bquote(bar(range)~Km^2),
#   "change_prop"=  "Range change %",
#   "maxlik_betasplit_TF"=expression(beta['max split']~"Total"),
#   "maxlik_betasplit_TF"=expression(beta['max split']~"Total"),
#   "func_FDiv_D"=expression('F'['unc']~'D'['iv']~" D"),
#   "trs_dispersal_slope_time"=expression('Slope'~"dispersal"~'T'['start']/'T'['end']),
#   "trs_temp_width_slope_time"=expression('Slope'~omega~'T'['start']/'T'['end']),
#   )

stats_symbol_lib <- list(
  "log10_mean_aplha"= c(bquote(log[10](bar(alpha))), description="Log 10 of mean alpha diversity (all sites) at final timestep."),
  "mean_range"= c(bquote(bar(occup)~Km^2), description="Mean regional occupancy calculated in Km2 according to the mean number of sites occupied by all species through time."),
  "change_prop"=  "Occupancy change %",
  "trs_dispersal_slope_time"=expression('Slope'~'d t'['start']~"~"~'t'['end']),
  "trs_competition_slope_time"=expression('Slope'~"l"~'t'['start']~"~"~'t'['end']),
  "trs_mean_temp_slope_time"=expression('Slope T t'['start']~"~"~'t'['end']),
  "trs_temp_width_slope_time"=expression('Slope'~omega~'t'['start']~'t'['end']),
  "dispersal"=c(pt=expression("Dispersal (d)"~" t"['start']), description="Dispersal trait value at the initial timestep." ),             
  "divergence_threshold"=c(pt=expression(phi) , description="Speciation threshold. Time of isolation until reproductive isolation." ),  
  "competition"=c(pt=expression("Competition (l)"~"t"['start']) , description="Competition trait value at the initial timestep."),           
  "finished"=c(pt="Finished" , description="Status of final simulation, if OK, simulation was finished." ),              
  "cpu_time"=c(pt=expression("CPU time (h)") , description="CPU time in hours for simulation." ),              
  "extinctions_perc"=c(pt="Extinction prop." , description="Proportion of extinction events in a simulation over final gamma diversity. " ),      
  "speciations_perc"=c(pt="Speciation prop.", description="Proportion of speciation events in a simulation over final gamma diversity. " ),      
   "n_sp_alive_t_0"=c(pt=expression(gamma) , description="Final gamma diversity, or final regional taxonomic richness (number of species alive at the last timestep)." ),        
   "turnover"=c(pt="Turnover" , description=" Turnover between speciation and extinction. I.e. the sum of speciation - extinction for all timesteps divided by the number of initial species." ),              
   "speciations_invar"=c(pt="Speciation invar." , description="Temporal speciation invariability, also refered to as temporal stability of speciation events." ),
  "extinction_invar"=c(pt="Extinction invar." , description="Temporal extinction invariability, also refered to as temporal stability of extinction events." ),
   "turnover_invar"=c(pt="Turnover invar." , description="Temporal turnover invariability, also refered to as temporal turnover stability." ),        
   "gamma"=c(pt=expression(gamma) , description="Final gamma diversity, or final regional taxonomic richness (number of species alive at the last timestep)." ),                  
   "mean_abd_5%"=c(pt=expression(bar("N"['5% t'['end']])) , description="Mean population size of its 5% quantile at final timestep." ),             
   "mean_abd_50%"=c(pt=expression(bar("N"['50%'])[' t'['end']]) , description="Mean population size at final timestep." ),     
   "mean_abd_95%"=c(pt=expression(bar("N"['95% t'['end']])), description="Mean population size of its 95% quantile at final timestep." ),          
   "mean_PD_alpha"=c(pt=expression(bar("PD"[alpha])) , description="The Mean Faith's Phylogenetic Diversity (PD) at the site level, quantifies the total branch lengths and was calculates with the function PD from the picante pacakge." ),         
   "mtx_mean_alpha_T"=c(pt=expression(bar(alpha)) , description="Mean alpha diversity of the entire arquipelago." ),      
   "mtx_mean_alpha_A"=c(pt=expression(bar(alpha['A'])) , description="Mean alpha diversity of sites on island A." ),      
   "mtx_mean_alpha_B"=c(pt=expression(bar(alpha['B'])) , description="Mean alpha diversity of sites on island B." ),     
   "mtx_mean_alpha_C"=c(pt=expression(bar(alpha['C'])) , description="Mean alpha diversity of sites on island C." ),      
   "mtx_mean_alpha_D"=c(pt=expression(bar(alpha['D'])) , description="Mean alpha diversity of sites on island D." ),     
   "mtx_gamma_T"=c(pt=expression(gamma) , description="Final gamma diversity, or final regional taxonomic richness (number of species alive at the last timestep)." ),           
   "mtx_gamma_A"=c(pt=expression(gamma['A'][" t"['end']]) , description="Final gamma diversity of island A, or final regional taxonomic richness (number of species alive at the last timestep)." ),           
   "mtx_gamma_B"=c(pt=expression(gamma['B'][" t"['end']]) , description="Final gamma diversity of island B." ),           
   "mtx_gamma_C"=c(pt=expression(gamma['C'][" t"['end']]) , description="Final gamma diversity of island C." ),          
   "mtx_gamma_D"=c(pt=expression(gamma['D'][" t"['end']]) , description="Final gamma diversity of island D." ),           
   "mtx_beta_prop_T"=c(pt=expression(beta~'%') , description="Proportional species turnover, i.e. 1-mean(alpha/gamma) quantifies what proportion of the species diversity in the dataset that is not contained in an average site" ),       
   "mtx_beta_prop_A"=c(pt=expression(beta['A']~'%') , description="Proportional species turnover in island A." ),       
   "mtx_beta_prop_B"=c(pt=expression(beta['B']~'%') , description="Proportional species turnover in island B." ),      
   "mtx_beta_prop_C"=c(pt=expression(beta['C']~'%') , description="Proportional species turnover in island C." ),       
   "mtx_beta_prop_D"=c(pt=expression(beta['D']~'%') , description="Proportional species turnover in island D." ),        
   "mtx_beta_w_T"=c(pt=expression(beta~'W') , description="Whittaker beta diversity, i.e. gamma/mean(alpha), how many subunits there would be if the total species diversity of the mean species diversity per subunit remained the same, but the subunits shared no species." ),          
   "mtx_beta_w_A"=c(pt=expression(beta['A']~'W') , description="Whittaker beta diversity of island A." ),          
   "mtx_beta_w_B"=c(pt=expression(beta['B']~'W') , description="Whittaker beta diversity of island B." ),          
   "mtx_beta_w_C"=c(pt=expression(beta['C']~'W') , description="Whittaker beta diversity of island C." ),           
   "mtx_beta_w_D"=c(pt=expression(beta['D']~'W') , description="Whittaker beta diversity of island D." ),           
   "mtx_zeta_T"=c(pt=expression("Zeta "~zeta) , description="Proportion of species that are common to all assemblages. Zeta diversity sensu Hui and McGeoch (2014)." ),            
   "mtx_eta_T"=c(pt=expression("Eta "~eta) , description="Proportion of species that are not unique to an assemblage nor common to all assemblages." ),             
   "mtx_theta_A"=c(pt=expression("Theta "~theta['A']) , description="Theta: species that are unique to an assemblage." ),           
   "mtx_theta_B"=c(pt=expression("Theta "~theta['B']) , description="" ),           
   "mtx_theta_C"=c(pt=expression("Theta "~theta['C']) , description="" ),           
   "mtx_theta_D"=c(pt=expression("Theta "~theta['D']) , description="" ),           
   "mtx_theta_T"=c(pt=expression("Theta "~theta) , description="" ),           
   "mtx_PD_S_A"=c(pt=expression("PD"['A'['S']]) , description="Standardized value of the unrooted Phylogenetic Diversity measure for species in island A." ),            
   "mtx_PD_S_B"=c(pt=expression("PD"['B'['S']]) , description="Standardized value of the unrooted Phylogenetic Diversity measure for species in island B." ),          
   "mtx_PD_S_C"=c(pt=expression("PD"['C'['S']]) , description="Standardized value of the unrooted Phylogenetic Diversity measure for species in island C." ),             
   "mtx_PD_S_D"=c(pt=expression("PD"['D'['S']]) , description="Standardized value of the unrooted Phylogenetic Diversity measure for species in island D." ),             
   "mtx_PD_S_T"=c(pt=expression("PD"[''['S']]) , description="Standardized value of the unrooted Phylogenetic Diversity measure for species on the arquipelago, calculated using the pd function from the picante package." ),             
   "mtx_MPD_S_A"=c(pt=expression("MPD"['A'['S']]) , description="Standardized value of the Mean Pairwise Distance measure in island A." ),           
   "mtx_MPD_S_B"=c(pt=expression("MPD"['B'['S']]) , description="Standardized value of the Mean Pairwise Distance measure in island B." ),          
   "mtx_MPD_S_C"=c(pt=expression("MPD"['C'['S']]) , description="Standardized value of the Mean Pairwise Distance measure in island C." ),         
   "mtx_MPD_S_D"=c(pt=expression("MPD"['D'['S']]) , description="Standardized value of the Mean Pairwise Distance measure in island D." ),          
   "mtx_MPD_S_T"=c(pt=expression("MPD"[''['S']]) , description="Standardized value of the Mean Pairwise Distance measure on the arquipelago, calculated using the mpd function from the picante package." ),            
   "mtx_MNTD_S_A"=c(pt= expression('MNTD'['A']) , description="Standardized value of the mean nearest taxon distance measure for the island A." ),         
   "mtx_MNTD_S_B"=c(pt= expression('MNTD'['B']) , description="Standardized value of the mean nearest taxon distance measure for the island B." ),          
   "mtx_MNTD_S_C"=c(pt= expression('MNTD'['C']) , description="Standardized value of the mean nearest taxon distance measure for the island C." ),         
   "mtx_MNTD_S_D"=c(pt= expression('MNTD'['D']) , description="Standardized value of the mean nearest taxon distance measure for the island D." ),          
   "mtx_MNTD_S_T"=c(pt= expression('MNTD'['T']) , description="Standardized value of the mean nearest taxon distance measure for the entire archipelago."),          
   "mtx_CD_S_AB"=c(pt=expression('CD'['A-B']) , description="Standardized Community Distance between island A and B. It is the beta diversity version of Mean Pairwise Distance (MPD), giving the average phylogenetic distance between two communities." ),           
   "mtx_CD_S_AC"=c(pt=expression('CD'['A-C']) , description="Standardized Community Distance between island A and C." ),           
   "mtx_CD_S_AD"=c(pt=expression('CD'['A-D']) , description="Standardized Community Distance between island A and D." ),         
   "mtx_CD_S_BC"=c(pt=expression('CD'['B-C']) , description="Standardized Community Distance between island B and C." ),          
   "mtx_CD_S_BD"=c(pt=expression('CD'['B-D']) , description="Standardized Community Distance between island B and D." ),          
   "mtx_CD_S_CD"=c(pt=expression('CD'['C-D']) , description="Standardized Community Distance between island C and D." ),        
   "trs_modes_dispersal"=c(pt="Mode d" , description="Number of modes at the dispersal trait distribution at final timestep, i.e. no modes (such as in a uniform distribution), or more. Calculated with the function 'Mode' from package LaplacesDemon." ),   
   "trs_modes_mean_temp"=c(pt=expression("Mode"~" T") , description="Number of modes at the thermal optimum distribution at final timestep." ),   
   "trs_modes_temp_width"=c(pt=expression("Mode"~" "~omega) , description="Number of modes at the thermal range trait distribution at final timestep." ),  
   "trs_modes_competition"=c(pt=expression("Mode"~" l") , description="Number of modes at the competition trait (i.e. tolerance to other species) distribution at final timestep." ), 
   "trs_dispersal_5%"=c(pt=expression(bar("d"['5% t'['end']])) , description="Mean dispersal trait of its 5% quantile at final timestep." ),       
   "trs_dispersal_50%"=c(pt=expression(bar("d"['50%'])~'t'['end']) , description="Mean dispersal trait at final timestep." ),     
   "trs_dispersal_95%"=c(pt=expression(bar("d"['95% t'['end']])) , description="Mean dispersal trait of its 95% quantile at final timestep." ),     
   "trs_mean_temp_5%"=c(pt=expression(bar("T"['5% t'['end']])) , description="Mean thermal optimum of its 5% quantile at final timestep." ),  
   "trs_mean_temp_50%"=c(pt=expression(bar("T"['50% t'['end']])) , description="Mean thermal optimum at final timestep." ),  
   "trs_mean_temp_95%"=c(pt=expression(bar("T"['95% t'['end']])) , description="Mean thermal optimum of its 95% quantile at final timestep." ),  
   "trs_temp_width_5%"=c(pt=expression(bar(omega['5% t'['end']])) , description="Mean thermal range trait of its 5% quantile at final timestep." ),    
   "trs_temp_width_50%"=c(pt=expression(bar(omega['50% t'['end']])) , description="Mean thermal range trait at final timestep." ),   
   "trs_temp_width_95%"=c(pt=expression(bar(omega['95% t'['end']])) , description="Mean thermal range trait of its 95% quantile at final timestep." ),     
   "trs_competition_5%"=c(pt=expression(bar("l"['5% t'['end']])) , description="Mean tolerance to other species trait of its 5% quantile at final timestep." ),      
   "trs_competition_50%"=c(pt=expression(bar("l"['50% t'['end']])) , description="Mean tolerance to other species trait at final timestep." ),   
   "trs_competition_95%"=c(pt=expression(bar("l"['95% t'['end']])) , description="Mean tolerance to other species trait of its 95% quantile at final timestep." ),   
   "trs_spread_dispersal"=c(pt=expression("Spread d t"['end']) , description="The range between the lowest and highest dispersal trait values at the final timestep." ),  
   "trs_spread_mean_temp"=c(pt=expression("Spread T t"['end']) , description="The range between the lowest and highest thermal optimum values at the final timestep." ), 
   "trs_spread_temp_width"=c(pt=expression("Spread "~omega~" t"['end']) , description="The range between the lowest and highest thermal range trait values at the final timestep." ), 
   "trs_spread_competition"=c(pt=expression("Spread l t"['end']) , description="The range between the lowest and highest tolerance to other species trait values at the final timestep." ), 
   "maxlik_betasplit_TF"=c(pt=expression(beta['max split']) , description="The complete phylogenetic beta value derived from the Maximum Likelihood estimation within the Beta-splitting model. This is computed using the maxlik.betasplit function from the apTreeshape package."),   
   "maxlik_betasplit_TP"=c(pt=expression(beta['max split']~"Pruned") , description="The pruned phylogeny beta value, excluding extinct branches, derived from the Maximum Likelihood estimation within the Beta-splitting model." ),
   "maxlik_betasplit_1"=c(pt=expression(beta['max split']~"Sp"['1']) , description="The pruned phylogeny beta value of species 1, excluding extinct branches." ),
   "maxlik_betasplit_2"=c(pt=expression(beta['max split']~"Sp"['2']) , description="The pruned phylogeny beta value of species 2, excluding extinct branches." ),
   "maxlik_betasplit_3"=c(pt=expression(beta['max split']~"Sp"['3']) , description="The pruned phylogeny beta value of species 3, excluding extinct branches." )
)




