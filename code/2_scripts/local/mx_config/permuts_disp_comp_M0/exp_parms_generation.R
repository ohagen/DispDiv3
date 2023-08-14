# grid exploration method, full 

model <- "M0"
# extremes

# add extremes
parameters <- list("dispersal"=seq(0,1, length.out=20),
                   "divergence_threshold"=65, #10 and 50 as 0.5 Myrs or 100-500 kyrs
                   "competition"=seq(0.9,1,length.out=20))
permuts <- expand.grid(parameters)
plot(permuts$dispersal, permuts$competition)

permuts <- as.list(permuts)
parameters <- permuts[1:length(parameters)]

current_wd <- getwd()
setwd("./1_gen3sis_formalization/mx_config")
make_experiment(model, parameters, n=NULL)
setwd(current_wd)