## METADATA ===============================================================
## Description: Plot summaries
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-11 16:22:06
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##


source(file.path(pls$dir_base, "source.R"))
source(file.path(pls$dir_base,"summary_stats/support/plots/load_plot_libs.R"))






# mbt <- readRDS("./3_summary/20230811_1127_ti_perms_disp_comp_2000_M0.rds")


mbt <- get_parm_stats(parms = sss[[1]]$parms, stat = sss[[1]]$stats$t)
plot_main <- c("gamma", "speciations_perc", "extinctions_perc", "mtx.mean_alpha_T", "mtx.beta_prop_T",
               "mtx.zeta_T","mtx.eta_T")#
# plot_stat <-colnames(sss[[1]]$stats$t)[!colnames(sss[[1]]$stats$t)%in%
#                              c("finished",
#                                "func.FRic_A","func.FRic_B","func.FRic_C","func.FRic_D",
#                                "func.FEve_A","func.FEve_B", "func.FEve_C", "func.FEve_D",
#                                "func.FDiv_A", "func.FDiv_B", "func.FDiv_C", "func.FDiv_D",
#                                "func.FDis_A", "func.FDis_B", "func.FDis_C", "func.FDis_D",
#                                "func.RaoQ_A","func.RaoQ_B","func.RaoQ_B","func.RaoQ_C",
#                                "trs.NA")] # remove finished!
plot_stat_patch <- c("mtx.theta_","mtx.MNTD_S_")
plot_stat_between_patches <- c("mtx.CD_S")
plot_stat_traits <- c("trs.dispersal_50%", "tr.mean_temp_50%", "trs.temp_width_50%","trs.competition_50")
plot_phylo_trees <- c("maxlik_betasplit_.1",  "maxlik_betasplit_.2",  "maxlik_betasplit_.3")

plot_stat <- plot_main
n_stats <- length(plot_stat)

mask_mbt <- mbt$n_sp_alive_t_0>10
for (stat_i in 1:n_stats){
  # stat_i <- 1
  print(paste("Statistics: ", plot_stat[stat_i]))
  plot_stat_classes(mbt[mask_mbt,], cats="competition", y=plot_stat[stat_i], x="dispersal", 
                    main='paste(plot_stat[stat_i], "final time-step")')
}


plot_stat_classes()





####### plot on side
##

LOAD sss


plot_stat <- #plot_phylo_trees
n_stats <- length(plot_stat)
n_expp <- length(sss)
par(mfrow=c(n_stats,n_expp))
for (stat_i in 1:n_stats){
    # stat_i <- 1
  print(paste("Statistics: ", plot_stat[stat_i]))
  for (m_i in 1:n_expp){
    mbt <- get_parm_stats(parms = sss[[m_i]]$parms, stat = sss[[m_i]]$stats$t)
    mask_mbt <- mbt$n_sp_alive_t_0>10
    plot_stat_classes(mbt[mask_mbt,], cats="competition", y=plot_stat[stat_i], x="dispersal", 
                      main='paste(plot_stat[stat_i], "final time-step")')
  }
}








################# WIP ##############
par(mfrow=c(1, n_expp), mai=c(1,1,1,1)*0.3, mgp=c(1,0.2,0), tcl = -0.2)
# cols <- 
for (i in 1:n_expp){
  plot(sss[[i]]$parms$dispersal, sss[[i]]$parms$competition, pch=4, cex=0.4, col=rgb(0,0,0,0.5,1), 
       ylab="Dispersal (d)", xlab="Competition (l)",
       
       main=names(sss)[i])
}
