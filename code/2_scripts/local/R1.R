library(gen3sis)
library(shape)
library(ape)


source("./code/2_scripts/source.R")
temp_folder <- "C:/temp/dispdiv3/book/slice"

times_plotted <- c(500:489, seq(490,250, by=-20),seq(230,190, by=-10), seq(180,0, by=-1))


plot_ranges_start <- function(species_list, landscape, disturb=0, max_sps=9) {
  #browser()
  disturb=abs(disturb)
  max_sps <- abs(max_sps)
  #plot landscape
  #oldpar <- par(no.readonly = TRUE)
  #on.exit(par(oldpar))
  #layout( matrix(c(1,1,2),nrow=1, byrow =TRUE)  )
  #layout.show(2)
  #par(mar=c(4,3,3,7), oma=c(0.1,0.8,0.3,0.8))
  # par(xpd = FALSE)
  raster::image(raster::rasterFromXYZ(cbind(landscape$coordinates,1)), main="species ranges", col="navajowhite3", asp = 1)
  n_species <- length(species_list)
  alive <- unlist(lapply(species_list, function(x){length(x$abundance)}))
  alive <- alive>0
  #visual combinations
  sp_cols <- c("#FF0000", "#001AFF", "#FFE500", 
               "#CCFF00", "green3", "#CC00FF", 
               "#FF0099", "#7F00FF", "#00FFB2") #FIX
  sp_pchs <- 1:9 #FIX
  
  #limit plot
  omitted <- 0
  n_sps_max <- sum(alive)
  if (n_sps_max>max_sps){
    omitted <- n_sps_max-max_sps
    n_sps_max <- max_sps
  }
  # set to case
  cols <- rep(sp_cols, length.out=n_species)
  pchs <- rep(sp_pchs, length.out=n_species)
  #par(xpd = TRUE)
  for (i in 1:n_sps_max){
    sp_i <- (1:n_species)[alive][i]
    img <- cbind(landscape[["coordinates"]][names(species_list[[sp_i]]$abundance),,drop=FALSE], species_list[[sp_i]]$id)
    df <- as.data.frame(img)
    plot_diturbance <- sample(seq(-disturb, disturb, by=0.01), 1)
    points(x=as.numeric(df$x)+plot_diturbance, y=as.numeric(df$y)+plot_diturbance, pch=pchs[sp_i], col=cols[sp_i])
  }
  # legend plotted species
  # empty plot
  #plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE, ann=FALSE)
  # legend
  #par(xpd=TRUE)
  legend("top", inset=c(-0.15,0), title=paste(n_sps_max, "species", paste0("\n[", omitted, ' omitted]')), legend=(1:n_species)[alive][1:n_sps_max], pch=pchs[alive][1:n_sps_max], col=cols[alive][1:n_sps_max], bty="n")
}


add_cut_mark <- function(l=1){
  mtext(" x                                                                                         x                 x                     x                                                                        x                                                                                                            x", 
        line=if(l==1){-35.5}else{-27})
  
}





# sims H
parms


m0 <- sss$M0$parms
for (config_n in 1:nrow(parms)){
  # config_n <- 1
  print(parms[config_n,,drop=F])
  fname1 <- paste0(c(round(parms[config_n,c( "dispersal", "competition")],2),"H",rownames(parms[config_n,]), ".png"), collapse="_")
  png(file.path(pls$dir_out, "plots", fname1), width = 880, height = 680, pointsize = 22)
  plot_summary(readRDS(file.path(pls$dir_out, rownames(parms[config_n,,drop=F]), "sgen3sis.rds")))
  dev.off()
  match <- m0[abs(m0$dispersal-parms[config_n,]$dispersal)==min(abs(m0$dispersal-parms[config_n,]$dispersal))&abs(m0$competition-parms[config_n,]$competition)<=min(abs(m0$competition-parms[config_n,]$competition)),]
  # print(paste("Matched:")
  print((match))
  name_match <- gsub("config_", "", rownames(match))
  # read matched sim
  con <- unz(file.path(pls$dir_out_zip,"perms_disp_comp_2000_M0", paste0(name_match,".zip")), file.path("output", rownames(match),"sgen3sis.rds"))
  msim <- readRDS(gzcon(con))
  
  fname2 <- paste0(c(round(match[c( "dispersal", "competition")],2),"R",rownames(match), ".png"), collapse="_")
  png(file.path(pls$dir_out, "plots", fname2), width = 880, height = 680, pointsize = 22)
  plot_summary(msim)
  dev.off()
}


sss$M0$parms$
  
  
  
  
  mbt <- get_parm_stats(parms = sss[[1]]$parms, stat = sss[[1]]$stats$t)
main_stats_names <- c("gamma", "mtx_beta_prop_T","mtx_mean_alpha_T", "mtx_zeta_T","mtx_eta_T","maxlik_betasplit_TF","speciations_perc", "extinctions_perc", "mtx_MPD_S_T", "mtx_MNTD_S_T")#
y_is <- "gamma"  
plot_stat_classes(mbt[,], cats="competition", y=y_is, 
                    x="dispersal",
                    xlab="Dispersal (d)",
                    ylab=stats_symbol_lib[[y_is]],
                    plt_type="FALSE")
slices <- c(0.63,0.84 ) # 0.73,
match <- rep(NA, length(slices))

for (s_i in 1:length(slices)){
  # s_i <- 1
  abline(v=slices[s_i])
  text(slices[s_i]-0.028,max(mbt[y_is])-30, paste0("d=", slices[s_i]), srt=90)
}

configs_Poster <- NA
for (s_i in 1:length(slices)){
  mask <- abs(mbt$dispersal-slices[s_i])==min(abs(mbt$dispersal-slices[s_i]))
  mbt[mask,]$dispersal
  match[s_i] <- mbt[mask,]$dispersal[1]
  print(paste("dispersal=",slices[s_i]))
  sims_t <- mbt[mask & (mbt$competition==0.9 | mbt$competition==1), ]
  #sims_t <- sims_t[,c("seed", "dispersal", "competition",  "cpu_time", "n_sp_alive_t_0", "extinctions_perc" )]
  configs_Poster <- rbind(configs_Poster, sims_t)
  # sims_t[!is.na(sims_t)]
}
# very dirty way of cleaning.. I know...
configs_Poster <- configs_Poster[-1,]
# add poster numbering
configs_Poster$poster_n <- paste("config_M0_P", 1:4,sep="_")
print(configs_Poster)

# create dir folders

dir.create(paste0(temp_folder, 1))
dir.create(paste0(temp_folder, 3))
{
  ############## BIG LOOP ##################
}

# loop over each book/slice for printing... this is done per dispersal value here

for (s_i in c(1,3)){
  # s_i <- 3
  page <- 0 # set page to zero for the new book!
  {
  jpeg(file = paste0("C:/temp/dispdiv3/book/slice", s_i, "/", formatC(page, width = 6, format = "d", flag = "0") ,"_", s_i, ".jpeg"), width = 2480, height = 3508, units = "px", 
       pointsize = 40)
  par(mfrow=c(2,2), oma=c(0,0,0,0), mai=c(5,0,5,0.5))
  
  for (pos_i in c(s_i, s_i+1)){
    # pos_i <- s_i 
    configs_slice <- configs_Poster$poster_n[pos_i]
    print(paste0("using: ",configs_slice))
    
    plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    add_cut_mark()     
    plot_stat_classes(mbt[,], cats="competition", y=y_is, 
                      x="dispersal",
                      xlab="Dispersal (d)",
                      ylab=stats_symbol_lib[[y_is]],
                      plt_type="FALSE")
    title(if(selected_table$competition[pos_i]<0.999){"high competition"}else{"no competition"}, font.main = 1)
    dispersal_slice <-round(configs_Poster$dispersal[configs_Poster$poster_n==configs_slice],2)
    abline(v=dispersal_slice)
    text(dispersal_slice-0.028,max(mbt[y_is])-30, paste0("d=", dispersal_slice), srt=90)
    points(dispersal_slice, y=selected_table$n_sp_alive_t_0[pos_i], cex=10)
    points(dispersal_slice, y=selected_table$n_sp_alive_t_0[pos_i], cex=3)
    points(dispersal_slice, y=selected_table$n_sp_alive_t_0[pos_i], cex=1)
    Arrows(code=1, x0 =dispersal_slice , y0 = selected_table$n_sp_alive_t_0[pos_i], 
           x1 = 0, y1 = max(mbt[y_is]), 
           cex=3, lwd = 2,
            arr.adj = 1)  
    legend("topleft", title=paste0(gsub("P_", "P", configs_slice), ".R"), 
           legend="", selected_table$n_sp_alive_t_0, cex=1.1, text.col="white", bg ="grey55")
  }
  dev.off()
  }
  ##### PLOT ANIMATION timesteps
  for (t_i in times_plotted){
    page <- page+1
    {
    jpeg(file = paste0("C:/temp/dispdiv3/book/slice", s_i, "/", formatC(page, width = 6, format = "d", flag = "0") ,"_", s_i, ".jpeg"), width = 2480, height = 3508, units = "px", 
         pointsize = 40)
    # par(mfrow=c(2,2), oma=c(0,0,0,0), mai=c(5,0,5,0.5))
      par(mfrow=c(4,2), oma=c(0,0,0,0), mai=c(0.5,0,0.5,0.5))
    # s_i <- 1
    for (pos_i in c(s_i, s_i+1)){
      # pos_i <- s_i 
      configs_slice <- configs_Poster$poster_n[pos_i]
      print(paste0("using: ",configs_slice))
      sps <- readRDS(file.path(pls$dir_out, configs_slice, "species", paste0("species_t_",t_i,".rds")))
      phy <- read.nexus(file.path(pls$dir_out, configs_slice, "phylogeny", paste0("phylogeny_t_",t_i,".nex")))
      lc <- readRDS(file.path(pls$dir_out, configs_slice, "landscapes", paste0("landscape_t_",t_i,".rds")))
      
      plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      title(formatC(page, width = 6, format = "d", flag = "0"), line=-1)
      formated_year <- paste0(formatC(round(t_i/100, 2),  width=2, flag = "0", format="f", digits=2), " Ma")
      if (length(phy$tip.label)>9){
        plot_richness(sps, lc)
        legend("bottomright", legend =formated_year, cex=1.5, text.col="white", bg ="black")
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
        add_cut_mark() 
        plot(phy, show.tip.label=F)
      }else{
        plot_ranges_start(sps, lc)
        legend("bottomright", legend = formated_year, cex=1.5, text.col="white", bg ="black")
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
        add_cut_mark(l=2) 
        plot(phy)
        
      }
      
      #
      
      # # title(if(selected_table$competition[pos_i]<0.999){"high competition"}else{"no competition"}, font.main = 1)
      # dispersal_slice <-round(configs_Poster$dispersal[configs_Poster$poster_n==configs_slice],2)
      # abline(v=dispersal_slice)
      # text(dispersal_slice-0.028,max(mbt[y_is])-30, paste0("d=", dispersal_slice), srt=90)
      # points(dispersal_slice, y=selected_table$n_sp_alive_t_0[pos_i], cex=10)
      # points(dispersal_slice, y=selected_table$n_sp_alive_t_0[pos_i], cex=3)
      # points(dispersal_slice, y=selected_table$n_sp_alive_t_0[pos_i], cex=1)
      # Arrows(code=1, x0 =dispersal_slice , y0 = selected_table$n_sp_alive_t_0[pos_i], 
      #        x1 = 0, y1 = max(mbt[y_is]), 
      #        cex=3, lwd = 2,
      #        arr.adj = 1)  
      # legend("topleft", title=paste0(gsub("P_", "P", configs_slice), ".R"), 
      #        legend="", selected_table$n_sp_alive_t_0, cex=1.1, text.col="white", bg ="grey55")
    }
    dev.off()
    }
  }
  
  
  
  
  
  
  
  
  
  
  
} # end loop per book




