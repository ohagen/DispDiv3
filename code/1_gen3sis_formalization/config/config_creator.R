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

model <- "M0"

parameters <- list("dispersal"=seq(0,1, length.out=50),
                   "divergence_threshold"=65,
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


model <- "M1"

parameters <- list("dispersal"=seq(0,1, length.out=50),
                   "divergence_threshold"=65,
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

model <- "M2"

parameters <- list("dispersal"=seq(0,1, length.out=50),
                   "divergence_threshold"=65,
                   "competition"=seq(0.9,1,length.out=40))

permuts <- expand.grid(parameters)
dim(permuts)
### APPLY TRADE-OFF FOR all initial conditions!
# I go though all of this because I really want to have this trade off function
# all declared in one place! If the template is false, so should be the definition of
# the input variables and the simulations... so it;s easier to see if a mistake 
# is present.
templatef  <- readLines("./1_gen3sis_formalization/mx_config/config_template_M2.R") # read template function to get the apply_trade_off function
function_line <- grep("^apply_trs_tradeoff <- function", templatef)
# Find the end of the function declaration
brace_count <- 0
for (i in function_line:length(templatef)) {
  brace_count <- brace_count + sum(grepl("\\{", templatef[i])) - sum(grepl("\\}", templatef[i]))
  if (brace_count == 0) {
    end_line <- i
    break
  }
}
function_string <- paste(templatef[function_line:end_line], collapse = "\n")
function_expr <- parse(text = function_string)
eval(function_expr) # load the apply_trs_tradeoff function!
# end 
permuts <- apply_trs_tradeoff(permuts) # apply it

plot(permuts$dispersal, permuts$competition)

permuts <- as.list(permuts)
parameters <- permuts[1:length(parameters)]

current_wd <- getwd()
setwd("./1_gen3sis_formalization/mx_config/")
make_experiment(model, parameters, n=NULL)
setwd(current_wd)


