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






mbt <- readRDS("./3_summary/20230811_1127_ti_perms_disp_comp_2000_M0.rds")



plot_stat <- colnames(dbssl)[-1] # remove finished!
n_stats <- length(plot_stat)

for (stat_i in 1:n_stats){
  # stat_i <- 1
  print(paste("Statistics: ", plot_stat[stat_i]))
  plot_stat_classes(mbt, cats="competition", y=plot_stat[stat_i], x="dispersal", main=paste(plot_stat[stat_i], "final time-step"))
}


plot_stat_classes()