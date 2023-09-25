## METADATA ===============================================================
## Description: Load relevant libraries for the summary stats
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-02 20:57:54
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

lib <- c("tools", "dplyr", "ape", "betapart", "apTreeshape", "picante", "FD", 
         "LaplacesDemon", "RRphylo", "phytools", "diversitree", "lessR")
sapply(lib, require, character.only = TRUE, quietly = TRUE, warn.conflicts = TRUE)
source(file.path(pls$dir_base,"summary_stats/support/call_stats_zip.R"))
source(file.path(pls$dir_base,"summary_stats/support/sup_sumstats.R"))
source(file.path(pls$dir_base,"summary_stats/support/lsstats.R"))
source(file.path(pls$dir_base,"summary_stats/support/lsstatstt.R"))

