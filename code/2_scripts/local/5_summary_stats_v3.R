## METADATA ===============================================================
## Description: summary statistics
## 
## R version: 4.2.2 for Windows
## Date: 2023-07-27 13:30:26
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

source("./2_scripts/source.R")
source("./2_scripts/5_support_summary/load_support_libs.R")
source("./2_scripts/6_support_plot/5_support_plot_summary.R")


model <- "M0"
expp <- "perms_disp_comp_2000_M0"


# load parameter table
parms <- read.csv2(paste0("./1_gen3sis_formalization/mx_config/", expp,"/",model,"_config_parameters.txt"))
# parms <- read.csv2("./1_gen3sis_formalization/mx_config/st50-100_M1/M1_config_parameters.txt")
lc <- readRDS("./1_gen3sis_formalization/mx_space/final_elev_2/landscapes.rds")
output_dir <- paste0("outputs_eve/", expp,"/")

parms[] <- lapply(parms, as.numeric)
# get dirs

outputs <- list.files(output_dir, recursive = F, full.names = F)

#selection... run only if necessary
selection <- grep(model, outputs, perl=T)
outputs <- outputs[selection]

# order outputs
numbers <- as.numeric(sub(paste0(model,"_(\\d+)\\.zip"), "\\1", outputs))
outputs <- outputs[order(numbers)]


# system variables
n_outs <- length(outputs)
# finished <- rep(1,n_outs)
# ssl <- list()
# phyl <- list()
# occl <- list()
# traitsl <- list()
# cputl <- list()


# # unzip all
# for (o_i in 1:n_outs){
#   unzip(file.path(output_dir, outputs[o_i]), exdir=output_dir)
#   print(paste("done with", outputs[o_i]))
# }

# prepare summary stats output table


##### CALL STATS FOR FIX TIME


#### FIX TIME #########
ssl <- list()
for (o_i in 1:n_outs){
  # o_i <- 1
  print(paste0("working with: ", outputs[o_i], ", o_i=", o_i))
  ssl[[o_i]] <- call_stats_zip(zf=outputs[o_i], i=output_dir, l=lc, o="3_summary", save=FALSE, t_i=0) # save each line as separated .rds
}
# merge single time stats into data frame
names(ssl) <- paste0("config_",gsub(".zip", "", outputs))
dbssl <- do.call(rbind, ssl)
dbssl <- as.data.frame(dbssl)
# convert stats to numeric
dbssl[!colnames(dbssl)%in%"finished"]  <- lapply(dbssl[!colnames(dbssl)%in%"finished"], as.numeric)
# get final merged table
mbt <- cbind(parms[rownames(dbssl),], dbssl)
mbt[is.na(mbt)] <- NA
#### SAVE MBT
time_stamp <- format(Sys.time(), "%Y%m%d_%H%M")
saveRDS(mbt, paste0("./3_summary/", time_stamp,"_ti_", expp,".rds"))
# mbt <- readRDS("./3_summary/20230729_2000_ti_permuts_disp_comp_M0.rds")
# mbt <- readRDS("./3_summary/20230730_2316_ti_perms_disp_comp_2000_M0.rds")
mbt <- readRDS("./3_summary/20230802_2105_ti_perms_disp_comp_2000_M0.rds")


#### TODO, SAVE ALSO SEPARATE ONLY TO DIFERENTIATE THE STATS




##### CALL STATS FOR CHANGING TIME


#### CHANGING TIME #########
ssltt <- list()
for (o_i in 1:n_outs){ 
  ssltt[[o_i]] <- call_stats_zip(zf=outputs[o_i], i=output_dir, l=lc, o="3_summary", save=FALSE, t_i=500:0) # save each line as separated .rds
}
# merge single time stats into data frame
names(ssltt) <- paste0("config_",gsub(".zip", "", outputs))
# transform list of matrix to list of data.frames
ssltt <- lapply(ssltt, function(x){
  lapply(x, function(x){
  dbssl <- as.data.frame(x)
  # convert stats to numeric
  dbssl[!colnames(dbssl)%in%"finished"]  <- lapply(dbssl[!colnames(dbssl)%in%"finished"], as.numeric)
  return(dbssl)}
  )})

##### SAVE SSLTT
time_stamp <- format(Sys.time(), "%Y%m%d_%H%M")
saveRDS(ssltt, paste0("./3_summary/", time_stamp,"_tt_", expp,".rds"))

## WIP ###
ssltt <- list(readRDS("C:/temp/dispdiv3/outputs_eve/selection_500t/ss_t_500-0_M2_1119.rds")) # resonable
ssltt <- list(readRDS("C:/temp/dispdiv3/outputs_eve/selection_500t/ss_t_500-0_M1_1119.rds")) # resonable
outputs <- c("M0_1119.zip", "M1_1119.zip", "M2_1119.zip")

gen3sis::plot_summary(read_zip("sgen3sis.rds",zf=outputs[1], i="C:/temp/dispdiv3/outputs_eve/selection_500t"))

# ssltt <- list(readRDS("C:/temp/dispdiv3/outputs_eve/selection_500t/ss_t_500-0_M2_1714.rds")) # this one starts with a very low
#### END WIP



# plot!
for (o_i in 1:n_outs){
par(mfrow=c(4,1))
  for (trt_i in c("dispersal", "competition", "mean_temp", "temp_width")){
    plot_trait_phylogeny(ssltt[[o_i]][[2]],trait =trt_i, type="ass")
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


plot_stat <- colnames(dbssl)[-1] # remove finished!
n_stats <- length(plot_stat)

for (stat_i in 1:n_stats){
  # stat_i <- 1
  print(paste("Statistics: ", plot_stat[stat_i]))
  plot_stat_classes(mbt, cats="competition", y=plot_stat[stat_i], x="dispersal", main=plot_stat[stat_i])
                      # expression("t"[0]))
                      #substitute(paste("CO"[2],'=',CO2_level,sep='')))
                      #expression(paste(parse(plot_stat[stat_i]),"t"[0])))
                      #paste(plot_stat[stat_i], "t0"))
}




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
