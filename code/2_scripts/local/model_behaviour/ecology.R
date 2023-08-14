## METADATA ===============================================================
## Description: Ecology function
## 
## R version: 4.2.2 for Windows
## Date: 2023-07-11 12:48:21
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##
ss_eff_emp <- list("dispersal"=NA, "dispersal_success"=NA, "environmental_filter"=NA, "competition_c"=NA, "competition_l"=NA)
ss_eff <- list()

# prepare input from traits list (t_l) and landscape list (l_l)
# return traits and landscape
prep_imp <- function(t_l,l_l){
  traits <- matrix(unlist(t_l), ncol=length(t_l))
  colnames(traits) <- names(t_l)
  landscape <- matrix(unlist(l_l), ncol=length(l_l))
  colnames(landscape) <- names(l_l)
  return(list("landscape"=landscape, "traits"=traits))
}

set_exp <- function(abundance, t_l, l_l, t_s){
  prep <- prep_imp(t_l = t_l, l_l = l_l)
  cols <- rainbow(length(abundance))
  t_s_l <- length(t_s)
  abt <- matrix(NA,ncol=length(abundance), nrow=t_s_l)
  colnames(abt) <- names(abundance)
  abt[1,] <- abundance
  ss_eff <- list()
  assign("ss_eff", ss_eff, envir = .GlobalEnv)
  return(list("abt"=abt, "t_s_l"=t_s_l, "prep"=prep, "cols"=cols, "t_s"=t_s))
}

plot_exp <- function(abt, epx, ...){
  plot(NULL, xlab="", ylab="",
       xlim=c(0, epx$t_s_l), ylim=c(0, max(abt, na.rm=T)), ...)
  for (spi in names(abundance)){
    lines(abt[,spi], col=epx$cols[as.numeric(spi)])
  }
  legend("topleft", legend = colnames(epx$abt), col=epx$cols, pch=3, bty="n")
}
# for arbitrary real constants a, b and non zero c. 
# It is named after the mathematician Carl Friedrich Gauss. 
# The graph of a Gaussian is a characteristic symmetric "bell curve" shape. 
# The parameter a is the height of the curve's peak, b is the position of the center of the peak 
# and c (the standard deviation, sometimes called the Gaussian RMS width) controls the width of the "bell". 
# this functions scales the abundance and with ns, meaning, a higher ns means a stricter niche strength,
# with ns=0 as no niche selection at all

fg <- function(x,a,b,c,ns=1){
  #### BROWSER ! ----------
  # browser()
  a <- a/c
  c <- c/ns # scale width for larg ns... ns<<<0.0001 ns lim to zero aproxes a straight line
  v <- a*exp(-((x-b)^2/(2*c^2)))
  return(v)
}
plot(fg(x=0.5, a=1, b=seq(0,1,0.05), c=.1), type='l', xaxt="n", ylab = "Environmental fitness") 
axis(1, at=1:length(seq(0,1,0.05)), labels=seq(9,26,length.out=length(seq(0,1,0.05)) ))


s_e_f_tt <- function(abundance, traits, landscape, ss_eff_i=ss_eff){
  ns <- length(abundance)
  # set temp_id
  temp_id_size <- ns+2
  temp_id <- matrix(NA, ncol=temp_id_size, nrow=length(ss_eff_emp))
  colnames(temp_id) <- c("id", "patch", names(abundance))
  rownames(temp_id) <- names(ss_eff_emp)
  temp_id[,"id"] <- 1
  temp_id[,"patch"] <- 1
  mask_patch <- rep(T,temp_id_size)
  mask_patch[1:2] <- FALSE
  # set env niche
  env_min_fg <- fg(x=landscape[,"min_temp"], a=1, b=traits[, "mean_temp"], c=traits[, "temp_width"])
  env_max_fg <- fg(x=landscape[,"max_temp"], a=1, b=traits[, "mean_temp"], c=traits[, "temp_width"])
  # set growth rate
  g <- .1
  # abundance_tii first is only what the env. determines to be the new abundances
  r_f <- g*sqrt(env_min_fg*env_max_fg) # geometric mean
  temp_id["environmental_filter",mask_patch] <- r_f
  # Competition
  c_c <- traits[,"competition_c"] # intra competition
  c_l <- traits[,"competition_l"]
  # vector c_c
  v_c_c <- (1-c_c)*abundance # the larger the c_c the better for the species
  temp_id["competition_c",mask_patch] <- v_c_c
  # abundance_tii accounts now to the reduction of intra competition c_c 
  # abundance_tc <- abundance-v_c_c
  # remove conspecifics
  # prepare matrx
  abd_l <- matrix(abundance, nrow=ns, ncol=ns)
  diag(abd_l) <- 0
  # apply c_l
  v_c_l <- (1-c_l)*colSums(abd_l)
  temp_id["competition_l",mask_patch] <- v_c_l
  # abundance_tii accounts now to the reduction of inter competition c_l
  # abundance_tl <- abundance-v_c_l
  temp_id["dispersal",mask_patch] <- as.numeric(abundance==1)
  abundance_tii <- abundance*(r_f-v_c_c-v_c_l)
  #abundance thhreashold
  # abundance_tii[abundance_tii<0.001] <- 0
  temp_id["dispersal_success",mask_patch] <- as.numeric(abundance>0)*temp_id["dispersal",mask_patch]
  ss_eff_i[[length(ss_eff_i)+1]] <- temp_id
  return(list("abundance"=abundance_tii, "ss_eff"=ss_eff_i))
}

# EQUILIBRIUM SOLUTION
s_e_f_t <- function(abundance, traits, landscape){
  # traits <- epx$prep$traits
  # landscape <- epx$prep$landscape
  ns <- length(abundance)
  #### get rf, here r_f is the per capita growth rate of biomass that depends on the local site conditions 
  # set env niche
  env_min_fg <- fg(x=landscape[,"min_temp"], a=1, b=traits[, "mean_temp"], c=traits[, "temp_width"])
  env_max_fg <- fg(x=landscape[,"max_temp"], a=1, b=traits[, "mean_temp"], c=traits[, "temp_width"])
  # set growth rate
  g <- .1
  # abundance_tii first is only what the env. determines to be the new abundances
  r_f <- g*sqrt(env_min_fg*env_max_fg) # geometric mean
  
  ###### get (a_ff) = same species interaction coefficient and (afh)= heterospecific interaction coefficient 
  # get traits Competition
  c_c <- traits[,"competition_c"] # intra competition
  c_l <- traits[,"competition_l"]
  
  # set a_ff and a_fh
  a_ff <- 1-c_c
  a_fh <- 1-c_l
  
  # check if conditions are met in order to continue
  if (any(a_ff<=a_fh)){
    stop("a_ff has to be bigger than a_fh for this equilibrium function to be used! Check your intial and evolutionary conditions of the competition traits.")
  }
  
  ####### get K_f = the carrying capacity of species f (that is the equilibrium for the case without heterospecific biomass
  K_f <- r_f/a_ff
  ###### get J = the total biomass J* of the community in equilibrium
  J <- get_J(a_ff, a_fh, K_f)
  wistop <- FALSE
  keep_on_while <- rep(TRUE, ns)
  while(wistop==FALSE){
    shall_live <- (a_ff*K_f)>(a_fh*J)
    if (all(shall_live[keep_on_while])){
      B_f <- ((a_ff*K_f)-(a_fh*J))/(a_ff-a_fh)
      B_f[!shall_live] <- 0
      wistop <- TRUE
    } else{
      print(paste("here we go: "))
      #### BROWSER ! ----------
      a_ff[!shall_live] <- 0
      a_fh[!shall_live] <- 0
      K_f[!shall_live] <- 0
      keep_on_while <- shall_live
      J <- get_J(a_ff, a_fh, K_f)
    }
  }
  names(B_f) <- names(abundance)
  B_f[B_f<0] <- 0
  return(B_f)
}

get_J <- function(a_ff, a_fh, K_f){
  J <- sum((a_ff*K_f)/(a_ff-a_fh), na.rm=T)/(1+sum((a_fh/(a_ff-a_fh)), na.rm=T)) # new
  return(J)
}



run_exp <- function(epx, ss_eff){
  for (ti in epx$t_s[-1]){
    prev_abd <- epx$abt[ti-1,]
    dt <- s_e_f_tt(abundance=prev_abd, 
                   traits=epx$prep$traits, 
                   landscape=epx$prep$landscape, 
                   ss_eff_i=ss_eff)
    epx$abt[ti,] <- prev_abd+dt$abundance
    epx$abt[ti,epx$abt[ti,]<0] <- 0
    ss_eff <<- dt$ss_eff
  }
  return(list("epx"=epx, "ss_eff"=ss_eff))
}



#### RUN EXPERIMENT 
es <- 0
t_l <- list("mean_temp"=c(.5,.5,.5,rep(.6,es)), 
            "temp_width"=c(0.3,0.2,0.1,rep(0.23,es)), 
            "competition_c"=c(0.8,0.8,0.8, rep(0.8,es)), 
            "competition_l"=c(1,1,1,rep(1,es)))
l_l <- list("min_temp"=0.43, "max_temp"=0.54)
abundance <- c(1, 1, 1)
names(abundance) <- 1:(es+3)
epx <- set_exp(abundance, t_l, l_l, t_s=1:30)
out <- run_exp(epx, ss_eff)
plot_exp(out$epx$abt, out$epx)






# set for random species
es <- 300 # community size
grid_size <- 1000
range_c_c <- seq(0.8,0.8,length.out=grid_size)
range_c_l <- seq(0.85,1,length.out=grid_size)
range_mean_temp <- seq(0.4,0.5,length.out=grid_size)
range_temp_width <- seq(0.2,0.8,length.out=grid_size)
t_l <- list("mean_temp"=sample(range_mean_temp, size=es,replace = T), 
            "temp_width"=sample(range_temp_width, size=es,replace = T), 
            "competition_c"=sample(range_c_c, size=es,replace = T), 
            "competition_l"=sample(range_c_l, size=es,replace = T))
l_l <- list("min_temp"=0.44, "max_temp"=0.46)
abundance <- rep(1,es)
names(abundance) <- 1:es
epx <- set_exp(abundance, t_l, l_l, t_s=1:50)
out <- run_exp(epx, ss_eff=list())
par(mfrow=c(2,1))
plot_exp(out$epx$abt, out$epx, 
         main=paste("alive%=",(sum(out$epx$abt[epx$t_s_l,]>0)/es)*100, "total_biomass=", round(sum(out$epx$abt[epx$t_s_l,]),0)))
hist(out$epx$abt[epx$t_s_l,out$epx$abt[epx$t_s_l,]>0])



#### RUN EXPERIMENT TO COMPARE AND TEST EQUILIBRIUM EQUATION

es <- 0
t_l <- list("mean_temp"=c(.5,.5,.5,rep(.6,es)), 
            "temp_width"=c(0.3,0.2,0.1,rep(0.23,es)), 
            "competition_c"=c(0.8,0.8,0.8, rep(0.8,es)), 
            "competition_l"=c(0.90, 0.92,0.91,rep(1,es)))
l_l <- list("min_temp"=0.43, "max_temp"=0.54)
abundance <- c(1, 1, 1, rep(1,es))
names(abundance) <- 1:(es+3)
epx <- set_exp(abundance, t_l, l_l, t_s=1:50)
out <- run_exp(epx, ss_eff)
plot_exp(out$epx$abt, out$epx)
abd_eq <- s_e_f_t(abundance, epx$prep$traits, epx$prep$landscape)
points(y=abd_eq, x=rep(epx$t_s_l, length(abd_eq)), col=epx$cols)




es <- 400 # community size
grid_size <- 1000
range_c_c <- seq(0.8,0.8,length.out=grid_size)
range_c_l <- seq(0.9,1,length.out=grid_size)
range_mean_temp <- seq(0.4,0.6,length.out=grid_size)
range_temp_width <- seq(0.3,0.6,length.out=grid_size)
t_l <- list("mean_temp"=sample(range_mean_temp, size=es,replace = T), 
            "temp_width"=sample(range_temp_width, size=es,replace = T), 
            "competition_c"=sample(range_c_c, size=es,replace = T), 
            "competition_l"=sample(range_c_l, size=es,replace = T))
l_l <- list("min_temp"=0.4, "max_temp"=0.6)
abundance <- rep(0.1, es)# c(rep(0.1,20),rep(1,length.out=(es-20)))
names(abundance) <- 1:es
epx <- set_exp(abundance, t_l, l_l, t_s=1:520)
out <- run_exp(epx, ss_eff)
par(mfrow=c(2,1))
plot_exp(out$epx$abt, out$epx, main=paste("alive%=",(sum(out$epx$abt[epx$t_s_l,]>0)/es)*100, "total_biomass=", round(sum(out$epx$abt[epx$t_s_l,]),0), "init_sp=",es))
# get equilibria
abd_eq <- s_e_f_t(abundance, traits=epx$prep$traits, landscape=epx$prep$landscape)
points(y=abd_eq, x=rep(epx$t_s_l, length(abd_eq)), col=epx$cols)
colmeans <- colMeans(out$epx$abt[(out$epx$t_s_l-10):out$epx$t_s_l,])
plot(jitter(abd_eq), jitter(colmeans), col=epx$cols, ylab="mean last 10 time-steps ODE", xlab="Estimated equilibria")
abline(h=0)
abline(v=0)
title(paste(sum(colmeans<=0&abd_eq>0), "Equilibria dubious alive,", sum(colmeans>0&abd_eq<=0), "dubious Equilibria dead,"))
lines(0:2, 0:2)
