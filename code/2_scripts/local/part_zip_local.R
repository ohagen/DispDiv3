## METADATA ===============================================================
## Description: partial loop to zip simulations run locally
## 
## R version: 4.2.2 for Windows
## Date: 2023-12-11 19:03:07
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##


for (config_n in 2:nrow(parms)){
  config_name <- paste0(model,"_", config_n)
  print("zipping...")
  z_file_name <- file.path(pls$dir_out,paste0( config_name, ".zip"))
  zip(zipfile=z_file_name, 
      files=file.path(pls$dir_out,paste0("config_",config_name)))
  # move it
  final_zip_path <- file.path(pls$dir_out_zip, expp)
}