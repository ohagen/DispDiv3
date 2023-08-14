## METADATA ===============================================================
## Description: Summary statistics
## 
## R version: 4.2.2 for Windows
## Date: 2023-07-27 13:30:26
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

source("./code/2_scripts/source.R")
print(unlist(pls))
source(file.path(pls$dir_base,"summary_stats/support/load_support_libs.R"))

## 
cnfs <- commandArgs(trailingOnly = TRUE)
print(cnfs)
model <- cnfs[1] #e.g. model <- "M0"
config_n <- cnfs[2] # e.g. config_n <- 1
expp <- cnfs[3] # e.g. expp <- "perms_disp_comp_2000_M0"

# # load parameters table
parms <- load_parameters(dc=pls$dir_config_gen, expp, model)
# # load landscape file
lc <- readRDS(file.path(pls$dir_env_gen,"landscapes.rds"))

# system variables
output_dir <- file.path(pls$dir_out_zip, expp)
# get dirs in order
outputs <- get_ordered_files(output_dir, model, fend=".zip")
# get size
n_outs <- length(outputs)

####       CALL STATS!      ###
#### FIX and var TIME #########
for (ti in list(0, 500:0)){ # first for the single time step, then for a sequence of time steps...
  print(paste0("working with: ", outputs[config_n], ", config_n=", config_n, " | at ", length(ti), " time-steps | experiment ", expp))
  call_stats_zip(zf=outputs[config_n], i=output_dir, l=lc, o=file.path(pls$dir_out_zip, "temp_summary", expp), save=TRUE, t_i=ti) # save each line as separated .rds
}


###################### END ##########################


# plot!
for (o_i in 1:n_outs){
par(mfrow=c(2,2))
  for (trt_i in c("dispersal", "competition", "mean_temp", "temp_width")){
    plot_trait_phylogeny(ssltt[[o_i]][[2]],trait =trt_i)
  }
}

gen3sis::plot_summary(read_zip("sgen3sis.rds",zf=outputs[o_i], i=output_dir))



########## WIP ######### 
source("./2_scripts/source.R" )
source("./2_scripts/5_support_summary/load_support_libs.R")
o_i <- 6

ssl[[o_i]] <- call_stats_zip(zf=outputs[o_i], i=output_dir, l=lc, o="3_summary", save=FALSE, t_i=0)
########## END ######### 


plot_stat_classes(mbt, cats="competition", y="extinctions_perc", x="dispersal")
plot_stat_classes(mbt, cats="competition", y="n_sp_alive_t_0", x="dispersal")
plot_stat_classes(mbt, cats="competition", y="turnover", x="dispersal")
plot_stat_classes(mbt, cats="competition", y="speciations_invar", x="dispersal")
plot_stat_classes(mbt, cats="competition", y="extinction_invar", x="dispersal")
plot_stat_classes(mbt, cats="competition", y="cpu_time", x="dispersal")







plot_space_changes(space=ssltt[[1]]$spatial_sps$`-1`,  na.val=0, ylab="# -1", main="absolute local extinctions")
plot_space_changes(space=ssltt[[1]]$spatial_sps$`-1`/ssltt[[1]]$spatial_sps$`0`,  na.val=0, ylab="# -1", main="proportional local extinctions")


plot_space_changes(space=ssltt[[1]]$spatial_sps$`+1`,  na.val=0, ylab="# +1")
plot_space_changes(space=ssltt[[1]]$spatial_sps$`+1`/ssltt[[1]]$spatial_sps$`0`,  na.val=0, ylab="# +1", main="proportional colonization extinctions")

plot_space_changes(space=ssltt[[1]]$spatial_sps$`+1`-ssltt[[1]]$spatial_sps$`-1`,  na.val=0, ylab="net range change", main="Range turnover")
plot_space_changes(space=(ssltt[[1]]$spatial_sps$`+1`-ssltt[[1]]$spatial_sps$`-1`)/ssltt[[1]]$spatial_sps$`0`,  na.val=0, ylab="proportional net range change", main="Range turnover")


hist(ssltt[[1]]$spatial_sps$`0`["0",], main="ranges of species at t0")


plot_space_changes(space=ssltt[[1]]$spatial_sps$`0`,  na.val=0)


gen3sis::plot_summary(read_zip("sgen3sis.rds",zf=outputs[1], i=output_dir))


# stop
st

# low good disp, high comp
mbt[(mbt$dispersal>0.36)&(mbt$dispersal<0.42)&(mbt$competition<0.93),]
ldhc <- which((mbt$dispersal>0.36)&(mbt$dispersal<0.42)&(mbt$competition<0.93))

# low good disp, low comp
mbt[(mbt$dispersal>0.36)&(mbt$dispersal<0.42)&(mbt$competition>0.98),]
ldlc <- which((mbt$dispersal>0.36)&(mbt$dispersal<0.42)&(mbt$competition>0.98))


# high  disp, high comp
mbt[(mbt$dispersal>0.75)&(mbt$dispersal<0.85)&(mbt$competition<0.93),]
hdhc <- which((mbt$dispersal>0.36)&(mbt$dispersal<0.42)&(mbt$competition>0.98))

# high  disp, low comp
mbt[(mbt$dispersal>0.75)&(mbt$dispersal<0.85)&(mbt$competition>0.98),]
hdlc <- which((mbt$dispersal>0.36)&(mbt$dispersal<0.42)&(mbt$competition>0.98))



hhdLc <- 400
hhdHC <- 20

for (o_i in hdlc){
  par(mfrow=c(2,2))
  for (trt_i in c("dispersal", "competition", "mean_temp", "temp_width")){
    plot_trait_phylogeny(ssltt[[o_i]][[2]],trait =trt_i)
  }
}


gen3sis::plot_summary(read_zip("sgen3sis.rds",zf=outputs[20], i=output_dir))


plot_stat_classes(mbt, cats="competition", y="n_sp_alive_t_0", x="dispersal")
abline(v=c(seq(0,1, length.out=10)))
seq(0.9,1, length.out=5)




dbssl <- do.call(rbind, ssl)
dbssl <- as.data.frame(dbssl)
# convert stats to numeric
dbssl[!colnames(dbssl)%in%"finished"]  <- lapply(dbssl[!colnames(dbssl)%in%"finished"], as.numeric)
# get final merged table
mbt <- cbind(parms, dbssl)


#plot summary
gen3sis::plot_summary(read_zip("sgen3sis.rds",zf=outputs[21], i=output_dir))

#plot mbt
plot(mbt$cpu_time, mbt$n_sp_alive_t_0, col=rev(viridis(50))[mbt$divergence_threshold-49], pch=19)
plot(mbt$cpu_time, mbt$speciations_perc, col=rev(viridis(50))[mbt$divergence_threshold-49], pch=19)
# save ssl

saveRDS(mbt, file = "3_summary/permut_M0_all.rds")

plot(mbt$dispersal, mbt$competition, col=heat.colors(max(mbt$'n_sp_alive_t_0'))[mbt$n_sp_alive_t_0], pch=3)

text(mbt$dispersal, mbt$competition, labels=mbt$n_sp_alive_t0)



### WIP ########
# plot richness

plot(rasterFromXYZ(ssl[[1]]$summary$`richness-final`))
plot(rasterFromXYZ(ssl[[80]]$summary$`richness-final`))

# plot spi
ti <- 0
spi <- 3
occi <- readRDS(file.path(oid,"occs",paste0("pa_t_", ti,".rds")))
plot(rasterFromXYZ(occi[,c("x", "y", as.character(spi))]), main=paste0(spi, " t=",ti))
plot(rasterFromXYZ(cbind(occi[,c("x", "y")], rowSums(occi[,-c(1,2)]))), main="richness", col=rainbow(15))
plot(phyi, show.node.label = TRUE)
#### END WIP


gammas <- lapply(ssl, get_gamma)

gammas <- unlist(gammas)

names(gammas) <- outputs



# merge summary with parameters

# merged big table
mbt <- cbind(parms, "gamma"=gammas)
plot_stat_classes(mbt, cats="competition", y="n_sp_alive_t_0", x="dispersal")

plot_stat_classes(mbt, cats="competition", y="n_sp_alive_t_0", x="dispersal")
plot_stat_classes(mbt, cats="competition", y="extinctions_perc", x="dispersal")

plot(mbt$dispersal, mbt$competition, col=gen3sis::color_richness_non_CVDCBP(max(mbt$n_sp_alive_t0))[mbt$n_sp_alive_t0], pch=15)


plot(mbt[,"divergence_threshold"],mbt[,"gamma"])
head(mbt)


gammas_tt <- lapply(ssl, get_gamma_tt )
gammas_tt <- lapply(gammas_tt, function(x){
  x <- x[-1]
})
gammas_tt_table <- do.call(rbind, gammas_tt)
# merge table over time
mbtt <- cbind(parms, gammas_tt_table)
plot_time_y(mbtt, bty="n")


#### TAXONOMIC DIVERSITY
# gamma


#### Phylogenetic DIVERSITY
# phylogenetic diversity
occnxy <- occi[,-c(1,2)]
# alternative occnxy just by patch
mdf <- merge(occi, lc$patch[,c("x", "y", "0")])
mdff <- mdf[1:4,-c(1,2,length(mdf))]
for (p_i in 1:4){
  # p_i <- 1
  mdff[p_i, ] <- as.numeric(as.logical(colSums(mdf[mdf[,"0"]==p_i,-c(1,2,length(mdf))])))
}
colnames(mdff) <- paste0("species", colnames(mdff))

pd(mdff, phyi, include.root = FALSE)

colnames(occnxy) <- paste0("species", colnames(occnxy))
pd_estimate <- pd.query(phyi, occnxy, standardize = T)
plot(rasterFromXYZ(cbind(occi[,c(1,2)], pd_estimate)))
mpd_estimate <- mpd.query(phyi, occnxy, standardize = T)
plot(rasterFromXYZ(cbind(occi[,c(1,2)], mpd_estimate)))
mntd_estimate <- mntd.query(phyi, occnxy,standardize = T)
plot(rasterFromXYZ(cbind(occi[,c(1,2)], mntd_estimate)))

# Computes the (standardized) value of the Community Distance measure
# is the beta diversity version of Mean Pairwise Distance (MPD), giving the average phylogenetic distance between two communities.
cd_estimate <- cd.query(phyi, occnxy, standardize = T)
plot(rasterFromXYZ(cbind(occi[,c(1,2)], cd_estimate[,120])))

cd_estimate <- cd.query(phyi, mdff, standardize = F)


#### Functional Diversity



# plotting

plot(prunedphy, show.tip.label = FALSE, main = "x")
tiplabels(tip = which(samp[i, ] > 0), pch = 19, cex = 2)