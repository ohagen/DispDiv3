## METADATA ===============================================================
## Description: Create landscape
## 
## R version: 4.2.2 for Windows
## Date: 2023-06-30 00:21:44
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##
library(gen3sis)

source("./code/2_scripts/source.R") # includes pls list

# load
env_vars <- readRDS(pls$env_vars)
# to get time formated
sc <- read.csv(pls$sea_level)

names(env_vars)

elev_range <- get_range(range=list("min"=env_vars$elevation, "max"=env_vars$elevation))
# run once again to get finally the correct temp range after setting the maximum!
temp_range <- get_range(range=list("min"=env_vars$min_temp, "max"=env_vars$max_temp))

# Ranges
print(c("temperature"=temp_range))
# 9 and 26
print(c("elevation"=elev_range))
#-115 and 803

################## PLOT ENV ######################
# plot mean temp
par(mfrow=c(2,1))
plot(rev(unlist(lapply(env_vars$mean_temp, function(x){mean(x[], na.rm=T)}))), type="l", 
     xaxt="n", xlab="Ma", ylab="Mean Temp", main="")
axis(1, seq(0,501,100), labels=rev(sc[,1])[seq(1,501,100)], lwd=0.5)

# plot patches size
# get patches are change over time
df_patches <- lapply(env_vars$patch, function(x){
  df <- as.data.frame(x)
  ft <- c("1"=0,"2"=0,"3"=0,"4"=0)
  tdf <- table(df)
  ft[names(tdf)] <- tdf
  return(ft)
  })
df_patches <- do.call(rbind, df_patches)
df_patches <- df_patches[nrow(df_patches):1, ]
plot(df_patches[,1], ylim=range(df_patches),
     type='l', col=patches_color[1], ylab="Patch size", xlab="Ma", xaxt='n')
axis(1, seq(0,501,100), labels=rev(sc[,1])[seq(1,501,100)]   )
for(li in 4:2){
  lines(df_patches[,li], col=patches_color[li])
}
text(x = -5, y=c(df_patches[1,]+20), col=patches_color, labels=LETTERS[1:4])

# plot patches size
# get patches temp+change over time
sm <- list("A"=list(), "B"=list(), "C"=list(), "D"=list())
for (ti in 1:length(env_vars[[1]])){
  patches_mask <- env_vars$patch[[ti]]
  for (p_i in 1:4){
    sm[[p_i]][[ti]] <- c(
      "mean"=mean(env_vars$mean_temp[[ti]][patches_mask==p_i]),
      "5%"=quantile(env_vars$min_temp[[ti]][patches_mask==p_i],0.05),
      "95"=quantile(env_vars$max_temp[[ti]][patches_mask==p_i],0.95)
    )
  }
}
# do summary
smm <- lapply(sm, function(x){
  do.call(rbind,x)
})
smm <- lapply(smm, function(x){
  x[nrow(x):1,]
})
# get all
melted <-do.call(c, sm)
melted <- unlist(melted)
mm <- range(melted, na.rm=T)
plot(1, type="n", xlab="Ma", ylab="Temperature", xlim=c(1, 501), ylim=mm, xaxt='n')
for(li in 1:4){

  lines(smm[[li]][,1], col=patches_color[li])
  lines(smm[[li]][,2], col=patches_color[li])
  lines(smm[[li]][,3], col=patches_color[li])
}

plot(df_patches[,1], ylim=range(df_patches),
     type='l', col=patches_color[1], ylab="Patch size", xlab="Ma", xaxt='n')

df_patches <- df_patches[nrow(df_patches):1, ]
plot(df_patches[,1], ylim=range(df_patches),
     type='l', col=patches_color[1], ylab="Patch size", xlab="Ma", xaxt='n')
axis(1, seq(0,501,100), labels=rev(sc[,1])[seq(1,501,100)]   )
for(li in 4:2){
  lines(df_patches[,li], col=patches_color[li])
}
text(x = -5, y=c(df_patches[1,]+20), col=patches_color, labels=LETTERS[1:4])


################## END PLOT ENV ######################


# table(env_vars$patch[[1]][]) # example to gent patches area

#### create gen3sis input ------------

cost_function <- function(source, habitable_src, dest, habitable_dest){
  if (!all(habitable_src, habitable_dest)) {
    return(4)
  } else {
    # cy is the hypothenusa from Pitagoras at 1 degree scale, i.e. 100Km, thus :
    #-115 and 802 for min and max temp
    ab <- c(1,abs((source["elevation"]+115)/10000 -(dest["elevation"]+115)/10000) )
    cy <- sqrt(x = sum(ab^2))
    return(cy)
  }
}
cost_function_elev <- function(source, habitable_src, dest, habitable_dest){
  if (!all(habitable_src, habitable_dest)) {
    return(4)
  } else {
    # get elevation source and destination
    elevs <- c(source["elevation"], dest["elevation"])
    # avoid negatives
    elevs <- elevs+115
    # add 0.1 increase per 100m slope increase
    ab <- c(1, abs(elevs[1]-elevs[2])/1000)
    return(sum(ab))
  }
}
create_input_landscape(landscapes = env_vars, cost_function = cost_function_elev, directions=8,
                       output_directory = "./1_gen3sis_formalization/mx_space/final_elev_2", #timesteps = string_time_step,
                       calculate_full_distance_matrices = T, verbose=T)


plot(rasterFromXYZ(tit$elevation[,c(1,2,3)]))