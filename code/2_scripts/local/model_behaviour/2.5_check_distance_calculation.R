
wd <- "C:/Documents and Settings/am92guke/Documents/iDiv/Main_papers/DiversityDispersal/code/1_gen3sis_formalization/mx_space/"

setwd(wd)


# set i TODO{
d_m <- readRDS("./final_elev_2/distances_full/distances_full_13.rds")
# C:\VITAL LOCAL\Meus Documentos\ETH PhD\Code\R\package\Gen3sis\inst\extdata\SouthAmerica\landscape\distances_full\cost_function_null\distances_full
lc <- readRDS("final_elev_2/landscapes.rds")
   
lc_ti <-lc$mean_temp[dimnames(d_m)[[1]], c('x', 'y', 15)]    
#}
# set site
{
  pdf("../../4_figures/dynamic_landscape/resistance/resitance_t13.pdf", width=14, height=16)
  par(mfrow=c(2,2), mai=c(2,2,3,2))
  for (s_i in 1:4){
    # s_i <- 1
    dm_ti <- cbind(lc_ti[,c('x', 'y')], cost=as.numeric(d_m[,patches$center_ids[s_i]]))
    r_ti_si <- rasterFromXYZ(dm_ti)
    plot(r_ti_si, main=patches$names[s_i], axes=FALSE, box=FALSE, colNA="black")
    points(x=patches$center_coords[,"x"]+2, y=patches_center_coords[,"y"]+3, pch=c("A", "B", "C", "D"), cex=0.9, col="grey6")
    points(patches$center_coords[s_i,"x"],patches$center_coords[s_i,"y"] , pch=17, cex=0.9)
    # plot costs
    destin_patches <- c(1:4)[!(1:4%in%s_i)]
    destin_patches_center_dist <- d_m[patches$center_ids[destin_patches],patches$center_ids[s_i]]
    text(x=patches$center_coords[destin_patches,"x"],y=patches$center_coords[destin_patches,"y"]-2, 
         round(destin_patches_center_dist,0), cex=0.65)
  }
  dev.off()
}


abline(h=51)
abline(h=9)
abline(v=10)
abline(v=51)

dm_ti[(dm_ti$x==51.5|dm_ti$x==9.5)&(dm_ti$y==10.5|dm_ti$y==51.5),]

plot(r_ti_si, col=viridis(100))
# points(patches_center_coords[1,1],patches_center_coords[1,2])
points(x=patches$center_coords[,"y"], y=patches_center_coords[,"x"], pch=c("A", "B", "C", "D"))
abline(h=51)
abline(h=9)
abline(v=10)
abline(v=51)

#### END WIP



dist_water_t0 <- rasterFromXYZ(dist_water_t0_mat)

maxcost <- 6200
mincost <- 0
cost_breaks <- seq(mincost, maxcost, by=20)
cost_colors <- rev(gray(seq(0.03, 0.9, length.out=length(cost_breaks)-1)))

par(mar=c(1,2,1,2))
layout(matrix(c(1,1,1,1,2), ncol=1))

#par(mar=c(1,1,1,2))
image(dist_water_t0, col=cost_colors, breaks=cost_breaks)
#plot(dist_water_t0, col=cost_colors, breaks=cost_breaks, axes=F, box=F, legend.args = list(text = 'connection cost', side = 2,
#         font = 2, line = 0.5, cex = 0.8))
arrows(-61.5, -30.5, -61.5, -64.5, lwd=3, col='red')
text(-55, -40, labels=paste('connection cost =', round(dist_water_t0_mat$cost[dist_water_t0_mat$x==(-61.5) & dist_water_t0_mat$y==(-64.5)], 0)), cex=1.5, font=2)
