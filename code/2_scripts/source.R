## METADATA ===============================================================
## Description: Source file of project Dispersal Diversity
## 
## R version: 4.2.2 for Windows
## Date: 2023-07-06 23:59:31
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================#


#### DEFINE PATCHES
# define color code for patches
patches_center_coords <- matrix(NA, ncol=2, nrow=4)
colnames(patches_center_coords) <- c("x","y")
patches_center_coords[1,] <- c(9,51)
patches_center_coords[2,] <- c(51,51)
patches_center_coords[3,] <- c(9,10)
patches_center_coords[4,] <- c(51,10)
patches <- list(
  "names"=LETTERS[1:4],
  "values"=1:4,
  "colors"=c("#1b9e77","#d95f02","#7570b3","#e7298a"),
  "center_coords"=patches_center_coords,
  "center_ids"=c("490","532", "2950", "2992")
)

# decleare paths list 
pls <- list(
  "env_vars"="CLUSTER_IRRELEVANT",
  "sea_level"="CLUSTER_IRRELEVANT",
  "dir_base"=      if(interactive()){"./code/2_scripts/cluster"}else{"."},
  "dir_out"=       if(interactive()){"c:/temp/dispdiv3/output"}else{"output"},
  "dir_env_gen"=   if(interactive()){"c:/temp/dispdiv3/mx_space/ddl"}else{"landscapes/ddl"},
  "dir_config_gen"=if(interactive()){"./code/1_gen3sis_formalization/config"}else{"configs/"},
  "dir_out_zip"=   if(interactive()){"c:/temp/dispdiv3/outputs_eve"}else{"/data/idiv_onstein/gen3sis/output"}# either NULL, no zipping and moving zipped file
  ) 


###### Functions #########

# get ranges used for scaling!
get_range <- function(range=list("min"=env_vars$min_temp, "max"=env_vars$max_temp)){
  min <- lapply(range$min, function(x){
    cellStats(x, stat='range', na.rm=TRUE, asSample=TRUE)
  })
  max <- lapply(range$max, function(x){
    cellStats(x, stat='range', na.rm=TRUE, asSample=TRUE)
  })
  range <- range(c(unlist(min), unlist(max)))
  rt <- c("min"=range[1], "max"=range[2])
  return(rt)
}

