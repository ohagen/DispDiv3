## METADATA ===============================================================
## Description: Load relevant libraries for the simulation experiment
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-02 20:57:54
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

lib <- c("gen3sis","utils")
sapply(lib, require, character.only = TRUE, quietly = TRUE, warn.conflicts = TRUE)