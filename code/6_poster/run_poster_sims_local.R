## METADATA ===============================================================
## Description: Run single LOCAL poster session IBS, Prague 2024
## 
## R version: 4.2.2 for Windows
## Date: 2024-01-07 10:00:00
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

source("./code/2_scripts/source.R")
source(file.path(pls$dir_base,"simulations/load_sim_libs.R"))

## 
model <- "M0"
config_n <- config_n <- 1
expp <- "poster_IBS_2024"

# run simulation and store summary obejct to sim
config_name <- paste0(model,"_P_", config_n)
print(getwd())
sim <- run_simulation(config = file.path(pls$dir_config_gen, expp, paste0("configs_",model), paste0("config_",config_name,".R")), 
                      landscape = pls$dir_env_gen,
                      output_directory = pls$dir_out, verbose = 1)























# remove the output landscape folder to save space
print("unlinking...")
unlink(file.path(pls$dir_out,paste0("config_", config_name), "landscapes"),recursive=TRUE)

if (!is.null(pls$dir_out_zip)){
  # zip it
  print("zipping...")
  z_file_name <- file.path(pls$dir_out,paste0( config_name, ".zip"))
  zip(zipfile=z_file_name, 
      files=file.path(pls$dir_out,paste0("config_",config_name)))
  # move it
  final_zip_path <- file.path(pls$dir_out_zip, expp)
  print(paste("moving...", z_file_name, "to", final_zip_path))
  if(!dir.exists(final_zip_path)){
    dir.create(final_zip_path)
  }
  file.copy(z_file_name, final_zip_path)
  # remove from sources
  print("removing zip from work dir")
  unlink(z_file_name)
  # print(paste("removing", config_name, "from",pls$dir_out ))
  
}

# cat results
print("done run this single gen3sis simulation ;)")
