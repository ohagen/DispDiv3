## METADATA ===============================================================
## Description: apply sea_levels
## 
## R version: 4.2.2 for Windows
## Date: 2023-06-22 15:27:45
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##
library(raster)
library(plot3D)
library(viridis)
library(utils)

x <- seq(1, 60)
y <- x

lapserate <- 0.01 #Celsius decrease per m in altitude increase

#functions
plot_persp <- function(x, y, z, sea_level=sl_i, breaks_persp=seq(-700, 800, 5), main="") {
  ## getting the value of the midpoint
  #zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
  
  ## calculating the breaks
  breaks <- breaks_persp# hist(zz, plot=FALSE)$breaks
  ## cutting up zz
  # cols <- lcol(length(breaks)-1)
  land_color <- rev(heat.colors(sum(breaks>=sea_level)))
  water_palette <- colorRampPalette(c("dodgerblue4", "dodgerblue", "lightblue"))   
  water_color <- water_palette(sum(breaks<sea_level))  
  color <- c(water_color, land_color)
  # zzz <- cut(zz, breaks=breaks, labels=cols)
  
  ## plotting
  persp3D(x, y, z, col=color,
          phi=45, theta=98, ltheta = 105,
          expand = 0.5, border = NA, box = FALSE, shade = 0.75, main=main, 
          clab="sea                                            land      ", 
          colkey=list("line.clab"=1, adj.clab=0.5,side.clab=2 ))
  ## return breaks 
  list(breaks=breaks)#, colors=cols) #...and colors for the legend
  
}

# create patches names
attribute_patch_name <- function(raster_i, pn=LETTERS[1:4], doplot=TRUE){
  # raster_i <- env_vars$mean_temp[[13]]
  # pn vector of names
  # returns a 
  df <- as.data.frame(raster_i, xy=TRUE)
  df$layer[!is.na(df$layer)] <- 1
  df$layer[is.na(df$layer)|df$layer=="<NA>"] <- NA
  df[(df$x>=36)&(df$y>=30)&!is.na(df$layer),3] <- 2
  df[(df$x<36)&(df$y<30)&!is.na(df$layer),3] <- 3
  df[(df$x>=36)&(df$y<30)&!is.na(df$layer),3] <- 4
  r <- rasterFromXYZ(df)
  if (doplot){
    plot(r, col=rainbow(4), main=paste0("Patches: ",paste0(pn, collapse=", ")))
  }
  return(r)
}

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



lc <- readRDS("./2_scripts/20230515_landscapes_l.rds") # artificial islands
sc <- read.csv("../data/raw/sea_level_change_art_501.csv") # artifical sea level change made with 1_create_landscapes_v2.R
bgtemp <- readRDS("./2_scripts/background_temperature.rds") # background temperature
tc <- read.csv("../data/raw/temp_global.csv") # temperature change


e <- extent(c(0, 60, 0, 60))

# create empty vector to store rasters for each variable
env_vars <- list("elevation"=list(), 
                 "mean_temp"=list(), 
                 "min_temp"=list(), 
                 "max_temp"=list(),
                 "patch"=list())

breaks_m_temp <- seq(10, 26, by = 2)
m_temp_cols <- rev(heat.colors(length(breaks_m_temp) - 1))
breaks_r_temp <- seq(0, 5, by = 1)
r_temp_cols <- rev(heat.colors(length(breaks_r_temp) - 1))

for (i in 501:1){ #501
  # i <- 1
  issl <- lc[[i]]
  sl_i <- sc[i,"sea_level"]
  tc_i <- tc[i, "temp_reduction_c"]
  
  issl[issl<sl_i] <- NA
  #temp_i <- issl
  r_i <- raster(issl)
  extent(r_i) <- e
  # plot(r_i, main=paste(sc[i,"Ma"], "Ma"))
  temp_i <- bgtemp-tc_i
  temp_i <- temp_i-(lapserate*r_i) # apply lapse rate
  #min and max
  #temp_min_i <- temp_i
  # set min temp based on a higher sd for higher elevation
  temp_min_i <- temp_i-abs(rnorm(length(temp_i[]), 0, sd=(r_i[]+115)/(803+115)))
  temp_max_i <- temp_i+abs(rnorm(length(temp_i[]), 0, sd=(r_i[]+115)/(803+115)))
  
  env_vars$elevation[[i]] <- r_i
  env_vars$mean_temp[[i]] <- temp_i
  env_vars$min_temp[[i]] <- temp_min_i
  env_vars$max_temp[[i]] <- temp_max_i
  env_vars$patch[[i]] <- attribute_patch_name(temp_i, doplot=FALSE)
  {
    png(file = paste0("./4_figures/dynamic_landscape/",501-i, "_", sc[i,"Ma"], "Ma.png"),
        width = 980, height = 450, units = "px", pointsize = 23,
        bg = "white", res = NA)
    par(mfrow=c(1,3),mar=c(0.1, 0.1, 2, 0.1), mai=c(0.1,0.2,0.25,0.1)*2)
    plot_persp(x,y,lc[[i]], 
               main=paste0("[",formatC(sc[i,"Ma"], digits=2, format = "f"), "Ma]   ", "Elevation (m)"), 
               sea_level=sc[i,"sea_level"])
    #par(mai=c(0.1,0.1,0.1,0.1)*3)
    
    # pmat <- persp(x, y, lc[[i]], phi=30, 
    #               theta=-30,nticks=8,ticktype="detailed", add=F)
    text(trans3d(22,50,800,pmat), "A")
    text(trans3d(72,50,600,pmat), "B")
    text(trans3d(0,20,500,pmat), "C")
    if(i<340){
      text(trans3d(50,21,100,pmat), "D", cex=1.1, col=rgb(1,1,1,0.7))
      text(trans3d(50,20,100,pmat), "D", cex=1.1, col=rgb(1,1,1,0.7))
      text(trans3d(49.5,21,99,pmat), "D", cex=1.1, col=rgb(1,1,1,0.7))
      text(trans3d(59.5,20,98,pmat), "D", cex=1.1, col=rgb(1,1,1,0.7))
      text(trans3d(50,20,100,pmat), "D")
    }
    text(30,sc[i,"sea_level"], "-----")
    
    plot(temp_i, main=paste0("Mean Temperature (", intToUtf8(176), "C)"), col=m_temp_cols, breaks=breaks_m_temp, axes=FALSE, box=FALSE, colNA="black", legend=FALSE)
    #par(mar=c(0.1, 0.1, 2, 0.1), mai=c(0.1,0.1,1,0.1))
    plot(temp_i, legend.only=TRUE, legend.shrink=1, legend.width=1.7, zlim=c(10, 26),
         axis.args=list(at=breaks_m_temp, labels=breaks_m_temp),
         legend.args=list(text="", side=4, font=2, line=2),
         col=m_temp_cols, breaks=breaks_m_temp)
    plot(temp_max_i-temp_min_i, main=paste0("Temperature Range (", intToUtf8(176), "C)"), col=r_temp_cols, breaks=breaks_r_temp, axes=FALSE, box=FALSE, colNA="black", legend=FALSE)
    plot(temp_max_i-temp_min_i, legend.only=TRUE, legend.shrink=1, legend.width=1.7, zlim=c(0, 5),
         axis.args=list(at=breaks_r_temp, labels=breaks_r_temp),
         legend.args=list(text="", side=4, font=2, line=2),
         col=r_temp_cols, breaks=breaks_r_temp)
    dev.off()
  }
  # issl[issl<=sl_i] <- sl_i
  
  
  
  
  #Sys.sleep(0.2)
}

# save
saveRDS(env_vars, "./1_gen3sis_formalization/mx_space/env_vars.rds")

# load
env_vars <- readRDS("./1_gen3sis_formalization/mx_space/env_vars.rds")

names(env_vars)

elev_range <- get_range(range=list("min"=env_vars$elevation, "max"=env_vars$elevation))
# run once again to get finally the correct temp range after setting the maximum!
temp_range <- get_range(range=list("min"=env_vars$min_temp, "max"=env_vars$max_temp))
print(c("temperature"=temp_range))
# 10 and 26
print(c("elevation"=elev_range)) 
#-115 and 802

# plot mean temp
plot(rev(unlist(lapply(env_vars$mean_temp, function(x){mean(x[], na.rm=T)}))), type="l")
table(env_vars$patch[[1]][]) # example to gent patches area

#### create gen3sis input ------------

cost_function <- function(source, habitable_src, dest, habitable_dest){
  if (!all(habitable_src, habitable_dest)) {
    return(4)
  } else {
    return(1)
  }
}
create_input_landscape(landscapes = env_vars, cost_function = cost_function, directions=8,
                       output_directory = "./1_gen3sis_formalization/mx_space/temp_tit", #timesteps = string_time_step,
                       calculate_full_distance_matrices = T, verbose=T)


plot(rasterFromXYZ(tit$elevation[,c(1,2,3)]))
