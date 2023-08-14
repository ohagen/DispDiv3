## METADATA ===============================================================
## Description: Create config files and cluster batch files as well as
## parameters from experiments
## 
## R version: 4.0.0 for Windows
## Date: 2021-07-07 17:32:50
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

#install.packages("randtoolbox", repos="https://stat.ethz.ch/CRAN/")
library(randtoolbox)



#### FUNCTIONS ----------
# parameter values convertion
linMap <- function(x, from, to) {(x - min(x)) / max(x - min(x)) * (to - from) + from}

# config creation
create_configs <- function(comb_matrix, model){
  print("creating configs")
  for(i in 1:nrow(comb_matrix)){ #nrow(comb_matrix)
    # i <- 1
    
    comb_vector <- comb_matrix[i,]
    if( !dir.exists(paste0('configs_',model))){
      dir.create(paste0('configs_',model))
    }
    output_dir <- paste0('../../output/sd', model, '/config_',model,'_', i)
    
    config_i <- readLines(paste0('config_template_', model,'.R'))
    config_i <- attribute_variables(comb_vector, config_i, model)
    # config_i <- gsub('output_directory', as.character(output_dir), config_i)
    
    writeLines(config_i, paste0('configs_',model,'/config_',model,'_', i, '.R'))
    
    cat(i, '\n')
    
  }
}

# atribute variables
attribute_variables <- function(comb_vector, config_i, model){
  variables <- names(comb_vector)
  # variables <- variables[!variables%in%c("starting_time")]
  for (variable_i in variables){
    # variable_i <- variables[1] 
    config_i <- gsub(paste0('*.comb_vector\\$',variable_i), comb_vector[[variable_i]], config_i)
  }
  return(config_i)
}

# bat creation (need modification to adapt for the WSL cluster)
create_bat <- function(comb_matrix, model){
  print("creating bat")
  #create bat file
  run_head <- '@runAsMultiple, @Node_NODE21' #only relevant for a specific HPC
  script_name <- 'run_gen3sis.R'
  r_version <- "/_shared/R3.6.1/r-with-tools.bat"
  config_dir <- paste0('-c config/sd/configs_', model,'/config_', model,'_', 1:nrow(comb_matrix), '.R')
  input_dir <- paste0('-i ', 'E:/PhD/simulations/input/WorldMap200-0Ma_scotese_multiple_4d')
  output_dir <- paste0('-o output/', model)
  other_par <- ''#'-s all'
  run_body <- paste(r_version, script_name, config_dir, input_dir, output_dir, other_par)
  write(c(run_head, run_body), file=paste0('../run_', model, '_',Sys.Date(),'_4d.bat'))
}

#make experiment, call all previous functions
make_experiment <- function(model, parameters, n=NULL){
  # n is either NULL or the number of sobol samples between the range
  # if n is NULL, we dont use sobol and create the parameters simply based on the vector of parameters
  # if 
  if(is.null(n)){
     comb_matrix <- data.frame(do.call(cbind, parameters))
     print(paste0("experiment model: ", model, " (n=",length(parameters[[1]]),") taken from range of parameters"))
  }else{
    comb_matrix <- data.frame(sobol(n, length(parameters), init = T))
    print(paste0("experiment model: ", model, " (n=",n,")"))
  }
  
 
  parnames <- names(parameters)
  colnames(comb_matrix) <- parnames
  
  print(paste("n configs", nrow(comb_matrix)))
  
  # add seed
  comb_matrix <- cbind(seed=round(runif(nrow(comb_matrix))*100000,0),comb_matrix)
  
  #update parameters
  comb_matrix <- as.data.frame(comb_matrix)
  
  if (!is.null(n)){ # in case we use ranging for all parameters
    #update parameters
    for (pi in parnames[parnames!="starting_time"]){
      # pi="dispersal"
      ddd <- lapply(strsplit(as.character(parameters[[pi]][1]),"\\."), nchar)[[1]][2]
      if (is.na(ddd)){
        ddd <- 0
      }
      comb_matrix[,pi] <- round(linMap(comb_matrix[,pi], from=parameters[[pi]][1], to=parameters[[pi]][2]),ddd)
    }
  }

  
  #save comb_matrix
  rownames(comb_matrix) <- paste('config', model, 1:nrow(comb_matrix), sep='_')
  write.table(comb_matrix, paste0(model, '_config_parameters.txt'), row.names=T, col.names=T, sep=';')
  
  # write config and bat
  create_configs(comb_matrix, model)
  create_bat(comb_matrix, model)
}

#### SCRIPT ----------

# n=20
# 
# stime <- sample(c(1200, 960, 840), n, replace=T) #equivalent to 200, 160 and 140 Ma

### ancestor all world ---------------

# grid exploration method 26.07.2023, full 
# model <- "M0"
# # extremes
# 
# # add extremes
# parameters <- list("dispersal"=seq(0,1, length.out=20),
#                    "divergence_threshold"=65, #10 and 50 as 0.5 Myrs or 100-500 kyrs
#                    "competition"=seq(0.9,1,length.out=20))
# permuts <- expand.grid(parameters)
# plot(permuts$dispersal, permuts$competition)
# 
# permuts <- as.list(permuts)
# parameters <- permuts[1:length(parameters)]
# 
# current_wd <- getwd()
# setwd("./1_gen3sis_formalization/mx_config")
# make_experiment(model, parameters, n=NULL)
# setwd(current_wd)

### PERMUTS FULL RANGE DISPERSAL 1-50 weibull scale
model <- "M0"
# extremes

# add extremes
parameters <- list("dispersal"=seq(0,1, length.out=50),
                   "divergence_threshold"=65, #10 and 50 as 0.5 Myrs or 100-500 kyrs
                   "competition"=seq(0.9,1,length.out=40))
permuts <- expand.grid(parameters)
dim(permuts)
plot(permuts$dispersal, permuts$competition)

permuts <- as.list(permuts)
parameters <- permuts[1:length(parameters)]

current_wd <- getwd()
setwd("./1_gen3sis_formalization/mx_config")
make_experiment(model, parameters, n=NULL)
setwd(current_wd)




### SOBOL
# 
# model <- "M0"
# # extremes
# 
# # add extremes
# parameters <- list("dispersal"=rep(seq(0.001,0.999, length.out=90),3),
#                    "divergence_threshold"=rep(65,270), #10 and 50 as 0.5 Myrs or 100-500 kyrs
#                    "competition"=c(rep(0.9,90), rep(0.95, 90), rep(1,90))
#                    )
# current_wd <- getwd()
# setwd("./1_gen3sis_formalization/mx_config")
# make_experiment(model, parameters, n=NULL)
# setwd(current_wd)

# MOVE FILES MANUALY, i.e. configs and parameters_table

model <- "M0"

parameters <- list("dispersal"=c(0.01, 0.99), 
                   "divergence_threshold"=c(20,50), #20 and 50 as 0.5 Myrs or 100-500 kyrs
                   "competition"=c(0.91,1.00))

current_wd <- getwd()
setwd("./1_gen3sis_formalization/mx_config/")
make_experiment(model, parameters, n=100)
setwd(current_wd)


model <- "M1"

parameters <- list("dispersal"=rep(0.4, 30), ### CHECK AND CHANGE THIS VALUE!!!!!
                   "divergence_threshold"=round(seq(50,100,length.out=30),0), #20 and 50 as 0.5 Myrs or 100-500 kyrs
                   "competition"=rep(0.91,30))

current_wd <- getwd()
setwd("./1_gen3sis_formalization/mx_config/")
make_experiment(model, parameters, n=NULL)
setwd(current_wd)
# move things to experiment folder

