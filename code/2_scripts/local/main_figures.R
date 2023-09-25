## METADATA ===============================================================
## Description: Main figures plot
## 
## R version: 4.2.2 for Windows
## Date: 2023-08-28 22:25:27
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##

library(here)
library(fields)
library(RRphylo)
library(corrplot)
# library(rcartocolor)
dir_base <- "code/2_scripts"
here::i_am("./code/2_scripts/source.R")
source(here(dir_base,"source.R"))

# bypass is.interactive() trigered by quarto
pls <- list(
  "dir_base"=file.path(dir_base, "cluster"),
  "dir_out"="c:/temp/dispdiv3/output",
  "dir_env_gen"="c:/temp/dispdiv3/mx_space/ddl",
  "dir_config_gen"="../../1_gen3sis_formalization/config",
  "dir_out_zip"="c:/temp/dispdiv3/outputs_eve"
)
source(here(pls$dir_base,"../local/supp_func.R"))
source(here(pls$dir_base,"summary_stats/support", "sup_plot_summary.R"))
source(here(pls$dir_base,"summary_stats/support", "sup_plot_traits.R"))
source(here(pls$dir_base,"summary_stats/support", "call_stats_zip.R"))
# get summary files for all model!
ss_f <-get_ordered_files(file.path(pls$dir_out_zip, "temp_summary"), model=".rds", fend=".rds")
ss_f <- ss_f[1:3] # select only main stats
sss <- list()
for (i in 1:length(ss_f)){
  sss[[i]] <- readRDS(file.path(pls$dir_out_zip, "temp_summary", ss_f[i]))
}
names(sss) <- unlist(lapply(strsplit(ss_f, "_"), function(x){
  r <- x[length(x)] # TODO expp names hare
  r <- gsub(".rds", "", r)
  return(r)
}))

for (li in 1:length(sss)){
  sss[[li]]$stats$tt <- lapply(sss[[li]]$stats$tt, as.data.frame)
}

# Plot failed
failed <- lapply(sss, function(x){
  apply(x$stats$t, 2, function(y){
    sum(is.na(y))
  })
})
set_par(length(failed),1)
lapply(failed, lollipoPlot)

#### FIGURE 1 PATTERNS ############
mbt <- get_parm_stats(parms = sss[[1]]$parms, stat = cbind(sss[[1]]$stats$t,sss[[1]]$stats$tt$Total ))
change <- mbt$`range_spatial_sps_+1_mean`-mbt$`range_spatial_sps_-1_mean`
mbt$mean_range <- change+mbt$`range_spatial_sps_0_mean`
mbt$change_prop <- change/mbt$`range_spatial_sps_0_mean`
mbt$log10_mean_aplha <- log(mbt$mtx_mean_alpha_T, base=10)

f1_stats_names <- c("log10_mean_aplha",  
                    "mtx_beta_w_T",
                    "gamma"
                    #"mtx_eta_T",
                    #"maxlik_betasplit_.TF", 
                    #"mtx_MPD_S_T"
                    #"mtx_MNTD_S_T"
)#
f1_pross_names <- c("speciations_perc",
                    "extinctions_perc",
                    #"mean_abd_50%",  
                    "mean_range"
                    #"change_prop"
)#

stats_symbol_lib <- list(
  "log10_mean_aplha"= bquote(log[10](bar(alpha))), 
  "mtx_beta_prop_T"= expression(beta~'%'),
  "gamma"=  expression(gamma),
  "mtx_zeta_T"=  expression(zeta),
  "mtx_eta_T"=  expression(eta), 
  "mtx_MNTD_S_T"=  expression('MNTD'['S']),
  "speciations_perc"=  "Speciation %",
  "extinctions_perc"=   "Extinction %",
  "mean_abd_50%"=   bquote(bar(B)),
  "mean_range"= bquote(bar(range)~Km^2),
  "change_prop"=  "Range change %",
  "maxlik_betasplit_TF"=expression(beta['max split']),
  "func_FDiv_D"=expression('F'['unc']~'D'['iv']~" D"))

pnames <- c(rbind(f1_stats_names, f1_pross_names))
n_stats <- length(pnames)


mask_mbt <- mbt$n_sp_alive_t_0>=3
lp <- mbt[mask_mbt,c("dispersal", "competition", pnames)]


# plot_stat_classes_summary(lp, stats_names = pnames)

{
  pdf(file.path(pls$dir_out, "f1_1.pdf"), width = 6, height = 5.5)
  set_par(n_stats, 2)
  for (stat_i in 1:n_stats){
    # stat_i <- 1
    
    plot_stat_classes(lp, cats="competition", y=pnames[stat_i], x="dispersal", plt_type="FALSE", ylab=stats_symbol_lib[[pnames[stat_i]]])
    title(LETTERS[stat_i], adj=0)
    if (stat_i==1){ # plot colbar
      classes <- unique(mbt[,"competition"])
      n_classes <- length(classes)
      cols <- rev(gen3sis::color_richness_non_CVDCBP(n_classes))
      ypos1 <- max(mbt[,f1_stats_names[stat_i]], na.rm=T)
      
      
      width_colbar <- 0.08
      length_colbar <- 0.8
      ypos3 <- ypos1-width_colbar
      xpos <- 0.2
      colorbar.plot(x=xpos, y=mean(c(ypos1, ypos3)), strip=1:n_classes, col = cols, 
                    strip.width=width_colbar, strip.length=length_colbar,
                    horizontal = TRUE, adj.y=0.5)
      text(x=xpos, y=ypos1, labels = "Competition (l)")
      text(x=xpos-0.1, y=ypos3, labels = "0.9 \n *")
      text(x=xpos+0.1, y=ypos3, labels = "1 \n neutral", adj=0)
      
    }
    #if (stat_i==2){
    abline(v=unlist(patches$disprange), lty=3, col="grey", lwd=2)
    #}
    
  }
  dev.off()
  
}

### F1_2
stats_names <- c("speciation_between_A_mean", "speciation_between_B_mean", "speciation_between_C_mean", "speciation_between_D_mean", "speciation_within_A_mean", "speciation_within_B_mean","speciation_within_C_mean", "speciation_within_D_mean")
disprange <- patches$disprange
names(disprange) <- gsub("--", "-", unlist(lapply(disprange, put_brak)))
n_dispr <- length(disprange)
## START PLOT  
# set_par(n_dispr,ncols = 3)
lvals <- list()
# set_par(3,1)
for (disprange_i in 1:n_dispr) {
  # disprange_i <- 1
  # print(disprange_i)
  temp_mbt <-mbt[mbt["dispersal"]>=disprange[[disprange_i]][1]&mbt["dispersal"]<disprange[[disprange_i]][2],c("gamma", stats_names)]
  gamma<-temp_mbt$gamma
  temp_mbt<-temp_mbt[,-1]
  lvals[[disprange_i]]<-apply(temp_mbt,2, function(x){
    val<-mean(x/gamma)
  })
  pb<-proportions(lvals[[disprange_i]][grep("between", names(lvals[[disprange_i]]))])
  pw<-proportions(lvals[[disprange_i]][grep("within", names(lvals[[disprange_i]]))])
  pdf(file.path(pls$dir_out, paste0(disprange_i,"dougnut.pdf")), height = 2, width = 2)
  if (!is.nan(sum(pb))){
    PieChart(pb, fill=patches$colors, col="black", lwd=0.5, main=names(disprange)[disprange_i], labels="")  
  }
  if (!is.nan(sum(pw))){
    # PieChart(pw, col=patches$colors, add=T)
    PieChart(pw, fill=patches$colors, hole=0.02, radius=0.65, col="black", lwd=0.5, main=names(disprange)[disprange_i], labels="")
  }
  dev.off()
}






barplot(lvals[[1]])
tips=lvals[[1]]
ggplot(tips, aes(x= day,  group=sex)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="day") +
  facet_grid(~sex) +
  scale_y_continuous(labels = scales::percent)



# dataplot <- unlist(lvals[[disprange_i]])

lollipoPlot(lvals[[3]])
lollipoPlot(dataplot, 
            col=patches$colors, 
            pt.col=patches$colors,
            pch=rep(c(17,16),each=length(patches$names)),
            xaxt="n",
            lwd=1,
            ylim=c(0,7),
            ylab="Mean speciation per Myr",
            lty=rep(c(1,2,3), each=c(2*length(patches$names))))
text(x=c(1,9,17), y=-0.75, labels=names(sss[[1]]$stats$tt), adj=0)
title(paste0("d",dispnames[disprange_i]))
title(LETTERS[disprange_i], adj=0)
if (disprange_i==2){  
  legend("top", bty="n", c("between patches", "withing patches"), pch=c(17,16),title="Speciation")
}
if (disprange_i==1){
  legend("top", bty="n", LETTERS[1:4], text.col = patches$colors, title="Patch", title.col="black" )
}
if (disprange_i==3){
  legend("top", bty="n", names(patches$`time-phase`), lty=c(1:3), title="Time phase")
}
}





###### FIGURE 2 MODEL COMPARISON

mbtl <- lapply(sss, function(x){
  mbt<-get_parm_stats(parms = x$parms, stat = cbind(x$stats$t,x$stats$tt$Total ))
  change <- mbt$`range_spatial_sps_+1_mean`-mbt$`range_spatial_sps_-1_mean`
  mbt$mean_range <- change+mbt$`range_spatial_sps_0_mean`
  mbt$change_prop <- change/mbt$`range_spatial_sps_0_mean`
  mbt$log10_mean_aplha <- log(mbt$mtx_mean_alpha_T, base=10)
  return(mbt)
})

f2_stats_names <- c("log10_mean_aplha",  
                    "mtx_beta_prop_T",
                    "gamma",
                    "mtx_eta_T",
                    "maxlik_betasplit_TF", 
                    #"mtx_MPD_S_T"
                    "mtx_MNTD_S_T"
)#
f2_pross_names <- c("speciations_perc",
                    "extinctions_perc",
                    "mean_abd_50%",  
                    "mean_range",
                    "change_prop")
pnames<-c(f2_stats_names)#, f2_pross_names)
n_stats<-length(pnames)

pdf(file.path(pls$dir_out, "f2_1.pdf"), width = 6, height = 5.5)
set_par(round(n_stats*3,0), length(mbtl))
ita<<-0
for (stat_i in 1:n_stats){
  # stat_i<-1
  lapply(mbtl, function(x){
    ita<<-ita+1
    plot_stat_classes_p(mbt=x, y=pnames[stat_i], x="trs_dispersal_50%", cat="trs_competition_50%", 
                        ylab=stats_symbol_lib[[pnames[stat_i]]], xlab="Dispersal (d)", plt_type="goodie",
                        main=if (stat_i==1){names(mbtl)[ita]}else{""}, xlim=c(0,1))
  })
}
dev.off()



### REP
pnames<-c(f2_pross_names)#, f2_pross_names)
n_stats<-length(pnames)

pdf(file.path(pls$dir_out, "f2_2.pdf"), width = 6, height = 5.5)
set_par(round(n_stats*3,0), length(mbtl))
ita<<-0
for (stat_i in 1:n_stats){
  # stat_i<-1
  lapply(mbtl, function(x){
    ita<<-ita+1
    plot_stat_classes_p(mbt=x, y=pnames[stat_i], x="trs_dispersal_50%", cat="trs_competition_50%", 
                        ylab=stats_symbol_lib[[pnames[stat_i]]], xlab="Dispersal (d)", plt_type="goodie",
                        main=if (stat_i==1){names(mbtl)[ita]}else{""}, xlim=c(0,1))
  })
}
dev.off()
## END REP


