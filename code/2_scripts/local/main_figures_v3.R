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


# load fg
fg <- function(x,a,b,c,ns=1){
  a <- a/c
  v <- a*exp(-((x-b)^2/(2*c^2)))
  return(v)
}

###### fig pre
# load config eco:
plot(fg(seq(0,1,0.01), a=1, b=0.23,c=0.04), type="l", lwd=2, lty=1)
lines(fg(seq(0,1,0.01), a=1, b=0.6,c=0.1), lwd=2, lty=2)
hist(rweibull(10000, shape = 2, scale = 30), lty=1, col=rev(gen3sis:::color_richness(10))[1], freq=T, breaks = 80)
hist(rweibull(10000, shape = 2, scale = 10), lty=2, col=rev(gen3sis:::color_richness(10))[9], freq=T)  #from 1 to 50)




#### FIGURE 1 PATTERNS ############
lmbt <- lapply(sss, function(x){
  mbt <- get_parm_stats(parms = x$parms, stat = cbind(x$stats$t,x$stats$tt$Total))
  change <- mbt$`range_spatial_sps_+1_mean`-mbt$`range_spatial_sps_-1_mean`
  mbt$mean_range <- change+mbt$`range_spatial_sps_0_mean`
  mbt$change_prop <- change/mbt$`range_spatial_sps_0_mean`
  mbt$log10_mean_aplha <- log(mbt$mtx_mean_alpha_T, base=10)
  mbt$log10_mean_PD_aplha <- log(mbt$mean_PD_alpha, base=10)
  return(mbt)
})
lopt <- lapply(sss, function(x){
  mbt <- get_parm_stats(parms = x$parms, stat = cbind(x$stats$t,x$stats$tt$Total))
  mask_opt_space <- get_optima_space_mask(mbt)
  ordm <- order(mbt[mask_opt_space,'trs_dispersal_50%'])
  return(list("mask_opt_space"=mask_opt_space, "ordm"=ordm))
})


stats_names <- c("gamma", 
                    "mtx_beta_prop_T",
                 "log10_mean_aplha",
                    #"mtx_eta_T",
                    #"maxlik_betasplit_.TF", 
                    "log10_mean_PD_aplha",
                    "speciations_perc",
                    "extinctions_perc"
                    #"mean_abd_50%",  
                    
                    #"change_prop"
)

pnames <- stats_names
n_stats <- length(pnames)


# mask_mbt <- mbt$n_sp_alive_t_0>=3
# lp <- mbt[mask_mbt,c("dispersal", "competition", pnames)]
# 
# # M2 tradeoff
# m1T <- get_parm_stats(parms = sss$M1$parms, stat = cbind(sss$M1$stats$t,sss$M1$stats$tt$Total ))
# mask_opt_spacem1T <- get_optima_space_mask(m1T)
# ordm1T <- order(m1T[mask_opt_spacem1T,'trs_dispersal_50%'])
# 
# # M2 tradeoff
# m2 <- get_parm_stats(parms = sss$M2$parms, stat = cbind(sss$M2$stats$t,sss$M2$stats$tt$Total ))
# mask_opt_space <- get_optima_space_mask(m2)
# ord <- order(m2[mask_opt_space,'trs_dispersal_50%'])



# plot_stat_classes_summary(lp, stats_names = pnames)

{
  pdf(file.path(pls$dir_out, "f1_1_v5.pdf"), width = 6.5, height = 6.5)
  set_par(8, 2)
  for (stat_i in 1:n_stats){
    # stat_i <- 1
    # stat_i <- 3

      plot_stat_classes(lmbt$M0, cats="competition", y=pnames[stat_i], x="dispersal", 
                      plt_type="FALSE", 
                      xlab=if (stat_i==n_stats){"Dispersal (d)"}else{" "},
                      ylab=stats_symbol_lib[[pnames[stat_i]]])
                      #ylim=if(pnames[stat_i]=="extinctions_perc"){c(0,0.18)}else{NULL})
    
      
    # add ablines and title
    abline(v=unlist(patches$disprange), lwd=3, lty=1, col="grey")
    title(LETTERS[stat_i], adj=0) 
    ltys <- c(3,2,1)
      # PLOT M0T
    for (i in 1:length(lopt)){
      smoothingSpline <- smooth.spline(y=lmbt[[i]][lopt[[i]]$mask_opt_space, pnames[stat_i]][lopt[[i]]$ordm],
                                         x=lmbt[[i]][lopt[[i]]$mask_opt_space, 'trs_dispersal_50%'][lopt[[i]]$ordm])
      lines(smoothingSpline, lwd=1.8, lty=ltys[i])
    }
    
   
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
      text(x=xpos-0.1, y=ypos3, labels = "0.9; aff=0.1\n *")
      text(x=xpos+0.1, y=ypos3, labels = "1; aff=0 \n **", adj=0)

      legend("topleft", title="Trade-off optima", legend=c("M0", "ME", "MET"),lty=ltys, lwd=c(2,2), bty = "n")

    }
    #if (stat_i==2){
    
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
  # disprange_i <- 3
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
    PieChart(pb, fill=patches$colors, col="white", lwd=0.5, main=names(disprange)[disprange_i], labels="")  
  }
  if (!is.nan(sum(pw))){
    # PieChart(pw, col=patches$colors, add=T)
    PieChart(pw, fill=patches$colors, hole=0.02, radius=0.65, col="white", lwd=0.5, main=names(disprange)[disprange_i], labels="")
  }
  dev.off()
  pdf(file.path(pls$dir_out, paste0(disprange_i,"bar.pdf")), height = 5, width = 5)
  mtbar <- barplot(c(pw,pb), col=patches$colors, ylim=c(0,0.5), ylab="Speciation events %", border=NA, names="")
  lines(y=rep(-0.00,2), x=c(0.5,4.5), col="grey50", lwd=2)
  lines(y=rep(-0.00,2), x=c(5,9.5), col="grey50", lwd=2)
  text(y=rep(-.025,2),x=c(2.5,7), c("Within", "Between"), outter=T, adj=0.5,xpd=NA, cex=3)
  text(mtbar, c(pw,pb)+0.025, rep(patches$names, 2), adj=0.5,xpd=NA, cex=3)
  dev.off()
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
names(mbtl)[3] <- "M1T"

f2_stats_names <- c("gamma",
                    "mtx_beta_prop_T"
                    # "log10_mean_aplha",
                    #"mtx_eta_T",
                    #"maxlik_betasplit_TF", 
                    #"mtx_MPD_S_T"
                    # "mtx_MNTD_S_T"
)#
f2_pross_names <- c(#"speciations_perc",
                    #"extinctions_perc",
                    "mean_abd_50%",  
                    "mean_range",
                    "change_prop")
pnames<-c(f2_stats_names, f2_pross_names)
pnames<-c(f2_pross_names)
n_stats<-length(pnames)

pdf(file.path(pls$dir_out, "f2_1.pdf"), width = 6, height = 5.5)
set_par(round(n_stats*2,0), length(mbtl))
ita<<-0
for (stat_i in 1:n_stats){
  # stat_i<-1
  ylim_temp <- range(unlist(lapply(mbtl, function(x){x[[pnames[stat_i] ]]})))
  lapply(mbtl, function(x){
    ita<<-ita+1
    plot_stat_classes_p(mbt=x, y=pnames[stat_i], x="trs_dispersal_50%", cat="trs_competition_50%", 
                        ylab=if(ita%in%seq(1,n_stats*3,by=3)){stats_symbol_lib[[pnames[stat_i]]]}else{" "}, 
                        xlab=if(ita>((n_stats*3)-3)){stats_symbol_lib[["trs_dispersal_50%"]]}else{" "}, 
                        plt_type="goodie",
                        main=if(stat_i==1){names(mbtl)[ita]}else{" "}, xlim=c(0,1), ylim=ylim_temp)
    # TODO ADD FUNCTION TO APPLY TRADEOFF
    
    # lo <- loess(x[tradeoff_mask,pnames[stat_i]][ord]~x[tradeoff_mask,'trs_dispersal_50%'][ord], degree=0)
    #permuts <- apply_trs_tradeoff(x[,c("dispersal","competition")])
    # tradeoff_mask <- (permuts$dispersal==m0$dispersal)&(permuts$competition==m0$competition)
    
    mask_opt_space <- get_optima_space_mask(x)
    #print(table(mask_opt_space))
    
    ord <- order(x[mask_opt_space,'trs_dispersal_50%'])
    smoothingSpline <- smooth.spline(y=x[mask_opt_space, pnames[stat_i]][ord],
                                     x=x[mask_opt_space,'trs_dispersal_50%'][ord])
    
    
    
    print("OK")
    lines(smoothingSpline, lwd=2, lty=1)
    
    # lines(predict(lo), x=x[tradeoff_mask,'trs_dispersal_50%'][ord])
    
    
    #toffline <- aggregate(x = x[pnames[stat_i]], by = list(disp=x$"trs_dispersal_50%"), FUN = "mean")
    #lines(toffline)
  })
}
dev.off()

m0 <- mbtl[[1]]



m0t <- m0[tradeoff_mask,]
plot_stat_classes_p(mbt=m0, y=pnames[stat_i], x="trs_dispersal_50%", cat="trs_competition_50%", 
                    ylab=if(ita%in%seq(1,n_stats*3,by=3)){stats_symbol_lib[[pnames[stat_i]]]}else{" "}, 
                    xlab=if(ita>((n_stats*3)-3)){stats_symbol_lib[["trs_dispersal_50%"]]}else{" "}, 
                    plt_type="goodie",
                    main=if(stat_i==1){names(mbtl)[ita]}else{" "}, xlim=c(0,1), ylim=ylim_temp)




lines(toffline)
######## FIGURE sp/ext

lt <- lapply(sss, function(x){
  get_parm_stats(parms = x$parms,
                 stat = cbind(x$stats$t, x$stats$tt$Total))# cbind(x$stats$t, x$stats$tt$Total)
})

lt <- lapply(lt, function(x){
  x <- x[x$gamma>=20,]
})


set_par(length(lt),3)
#stats <- grep("_slope_time", names(sss$M0$stats$tt$Total), value = T)
y_n <- "extinctions_perc"
x_n <- "speciations_perc"
c_n <- "trs_competition_50%"
yliml <- range(unlist(lapply(lt, function(x){x[y_n]})), na.rm=T)
#yliml[1] <- 0.85
xliml <- range(unlist(lapply(lt, function(x){x[x_n]})), na.rm=T)
#yliml <- range(c(yliml, xliml))
#xliml <- yliml
#yliml[1] <- 0.8
ita <<-0
lapply(lt, function(x){
  ita <<-ita+1
  plot_stat_classes_p(x, y=y_n, x=x_n, cats=c_n, ylim=yliml, xlim=xliml,  
                      ylab=stats_symbol_lib[[y_n]],
                      cex_p=0.5+x$`trs_dispersal_50%`*2,
                      pch_p=1,
                      xlab=stats_symbol_lib[[x_n]],
                      plt_type=if(ita==1){"colbar"}else{"NULL"}, 
                      yposbar = mean(yliml)+0.03,
                      xposbar = 0.1)
  title(names(lt)[ita])
  title(LETTERS[ita], adj=0)
  lines(y=c(1,1), x=c(0,0.6), col="grey")
  if (ita==1){
    #legend("bottomright", legend=c(seq(0,1,0.25)), pch=1, pt.cex=c(0.5+seq(0,1,0.25)*2), bty="n", main=)
  }
})
# dev.off()