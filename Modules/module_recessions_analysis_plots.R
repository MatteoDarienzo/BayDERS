


###############################################################################################################
plot.extracted.recessions <- function(curves.h, Data_h, hmin_grid, hmax_grid, dir.extraction, 
                                      station.name, data.period, chi) {
###############################################################################################################

  # Plot extracted curves in time series:
  #****************************************
  rec.time.plot <- ggplot(Data_h, aes(x=curves.h$t.real, y=curves.h$Qcurve, color=curves.h$RCi)) + 
    geom_point(size = 0.1) + 
    #geom_point(aes(x=curves$tpeak, y=curves$Qpeak, fill= curves$RCi), size = 4) +
    scale_color_gradientn(colours = rainbow(5)) +
    scale_x_continuous(name = "Time [day]", expand = c(0,0)) +
    scale_y_continuous(name = "Stage h [cm]", limits = c(hmin_grid,hmax_grid), expand = c(0,0)) +
    #ggtitle(paste("All recession curves",station.name, data.period, sep=" ")) +
    xlab("Time [day]") + 
    ylab("Stage [cm]") + 
    labs(color = "Recession index") +
    theme_bw(base_size = 15)+
    theme( axis.text=element_text(size=15)
           ,axis.title=element_text(size=20,face="plain")
           ,panel.grid.major=element_blank()
           ,panel.grid.minor=element_blank()
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="none"
           ,plot.margin=unit(c(0.5, 0.5, 0.5, 2),"cm"))
           #,axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l =0 )))
  pdf(paste0(dir.extraction,"/Figure1.pdf"), 6, 3 ,useDingbats=F)
  print(rec.time.plot)
  dev.off()
  ggsave(rec.time.plot, filename =paste0(dir.extraction,"/time_recessions_",chi,"cm.png"), 
         bg = "transparent", width = 6, height =3, dpi = 200)
  
  

  
  
  
  # Plot extracted curves against recession time t*:
  #*************************************************
  rec.curves.plot <- ggplot(Data_h, aes(x=curves.h$tcurve, y=curves.h$Qcurve, color=curves.h$RCi)) + 
    geom_point(size = 0.1) + labs(color = "# of the recession")+
    scale_color_gradientn(colours = rainbow(5)) +
    scale_x_continuous(name = "Recession time [day]", limits =c(0,80), expand = c(0,0)) +
    scale_y_continuous(name = "Stage h [cm]", limits = c(hmin_grid, hmax_grid), expand = c(0,0)) +
    #ggtitle(paste("All recession curves",station.name, data.period, sep=" ")) +
    theme_bw(base_size = 15)+
    theme( axis.text=element_text(size=15)
           ,axis.title=element_text(size=20,face="plain")
           ,panel.grid.major=element_blank()
           ,panel.grid.minor=element_blank()
           ,legend.text=element_text(size=15)
           ,legend.title=element_text(size=20)
           ,legend.key.size=unit(0.5, "cm")
           ,legend.position = c(.99, .99)
           ,legend.justification = c("right", "top")
           ,legend.box.just = "right"
           #,legend.margin = margin(6, 6, 6, 6)
           ,plot.margin=unit(c(0.5,0.5, 0.5, 2),"cm")
           )
  pdf(paste0(dir.extraction,"/Figure2.pdf"), 8, 4 ,useDingbats=F)
  print(rec.curves.plot)
  dev.off()
  ggsave(rec.curves.plot, filename =paste0(dir.extraction,"/all_curves_",chi,"cm.png"), 
         bg = "transparent", width = 8, height =4, dpi = 300)
  
  #************************************************************************************
  plot4paper.extraction =  plot_grid( 
                            rec.time.plot, rec.curves.plot,
                            labels = c('a)', 'b)'), 
                            label_size = 30,
                            label_fontface = "plain",
                            ncol = 2, nrow=1)  
  
  pdf(paste0(dir.extraction,"/Figure3.pdf"), 10, 4 ,useDingbats=F)
  print(plot4paper.extraction)
  dev.off()
  ggsave(filename = paste0(dir.extraction,"/Figure4paper_extraction.png"), 
         plot=plot4paper.extraction,
         width = 16, height =6, dpi = 300) 
  
  
  return(list(rec.time.plot = rec.time.plot,
              rec.curves.plot=rec.curves.plot,
              plot4paper.extraction=plot4paper.extraction))
}



































#################################################################################################################
recession.sensitivity.chi.for.paper  <- function(t_limni, h_limni, 
                                                 dir.case_study, dir_code, dir_exe,
                                                 station.name, data.period, 
                                                 chi.to.test, uh.rec, Nburn.rec, Nmin.rec, tmin.rec,  
                                                 tgood, delta.t.max,  delta.t.min, rec.model, 
                                                 tmax.rec, hmax.rec,
                                                 stage.scale.shift) {
#################################################################################################################
  #EXTRACT CURVES:
  dir.create(paste0(dir.case_study,"/Results/segmentation_recessions/curves_extraction"))
  dir.extraction <- paste0(dir.case_study,"/Results/segmentation_recessions/curves_extraction")
  dir.create(paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression"))
  dir.regression <- paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression")
  rec.plot.test =NULL
  
  
  
  
  
  #iterate for the tests of chi:
  for (jj in 1:length(chi.to.test)){
  #*********************************************************************************
  min_loc.h <- localmin(t_limni, 
                        h_limni*100)  #stage in centimeters !!!
  recess.h <- rec(min_loc.h$tmin, 
                  min_loc.h$Qmin, 
                  chi= chi.to.test[jj],
                  delta.t.max = delta.t.max)
  #extract all the recession curves:
  curves.h <- extract_curve(trec = recess.h$trec, 
                            Qrec = recess.h$Qrec)
  hmin = min(curves.h$Qcurve); 
  hmax = max(curves.h$Qcurve); 
  hmin_grid = 0; 
  hmax_grid = 0
  hgrid = seq(-10000,100000,100); 
  tgridd = seq(0,10000,10);
  tmax = max(curves.h$tcurve); 
  tmax_grid = 0
  for (i in 1:length(hgrid)) {
    if((hmin >= hgrid[i]) & (hmin <= hgrid[i+1])) {
      hmin_grid = hgrid[i]
    }
  }
  for (i in 1:length(hgrid)) {
    if((hmax >= hgrid[i]) & (hmax <= hgrid[i+1])) {
      hmax_grid = hgrid[i+1]
    }
  }
  for (i in 1:length(tgridd)) {
    if((tmax >= tgridd[i]) & (tmax <= tgridd[i+1])) {
      tmax_grid = tgridd[i+1]
    }
  }
  Data_h = data.frame(RCi    = curves.h$RCi, 
                      hcurve = curves.h$Qcurve, 
                      tcurve = curves.h$tcurve, 
                      treal  = curves.h$t.real, 
                      tpeak  = curves.h$tpeak,
                      hpeak  = curves.h$Qpeak, 
                      tfinal = curves.h$tfinal,
                      hfinal = curves.h$Qfinal)
  #detect all time final of the recession curves:
  tf.h = 0 ; Qf.h = 0 ; index.h = 0 ; icurve.h = 0; Qpeak.h = 0; tpeak.h =0;
  for (ii in 1:length(curves.h$Qcurve)) {
    if (curves.h$Qfinal[ii] != 0) {
      icurve.h = icurve.h + 1 
      tf.h[icurve.h] = curves.h$tfinal[ii]      # Final time of a recession
      Qf.h[icurve.h] = curves.h$Qfinal[ii]      # Final Q of a recession
      Qpeak.h[icurve.h] = curves.h$Qpeak[ii]    # Qpeak of a recession
      tpeak.h[icurve.h] = curves.h$tpeak[ii]    # tpeak of a recession
      index.h[icurve.h] = ii                  # curve index
    } else {
      icurve.h = icurve.h
    }
  }
  # Saving results of recession curves extraction (in stage h) :
  write.table(Data_h, file=paste0(dir.extraction,"/Extract_rec_curves.csv"), sep=";")
  #initialisation:
  asymptote.h = 0; curve_good.h = 0;hpeakgood.h = 0; t.real.good.h =0; index.good.h = 0; 
  asym.h.maxpost = 0; asym.h.stdev= 0; asym.h.Q10= 0; asym.h.Q90 = 0; asym.h.mean = 0;
  asym.h.Q2 = 0; asym.h.Q95 =0;  Resultss =NULL; quantiles.rec = NULL;
  theta1.maxpost = 0; theta1.stdev= 0; theta1.Q10= 0; theta1.Q90 = 0; theta1.mean = 0;
  theta2.maxpost = 0; theta2.stdev= 0; theta2.Q10= 0; theta2.Q90 = 0; theta2.mean = 0;
  theta3.maxpost = 0; theta3.stdev= 0; theta3.Q10= 0; theta3.Q90 = 0; theta3.mean = 0;
  theta4.maxpost = 0; theta4.stdev= 0; theta4.Q10= 0; theta4.Q90 = 0; theta4.mean = 0;
  theta5.maxpost = 0; theta5.stdev= 0; theta5.Q10= 0; theta5.Q90 = 0; theta5.mean = 0;
  theta6.maxpost = 0; theta6.stdev= 0; theta6.Q10= 0; theta6.Q90 = 0; theta6.mean = 0;
  results.regression = 0; asymptote.df <- data.frame(NULL); asym.df.temp =NULL
  d.h = NULL
  d.h.selected=NULL; d.h.selected.with.true.time=NULL;
  ncurves.h = icurve.h
  Nrec_longer_100=0
  
  for (i in 2:ncurves.h) {         
  #*************************
    d.h[[i]] = data.frame(t=curves.h$tcurve[(index.h[i-1]+1):(index.h[i]-1)],
                          h=curves.h$Qcurve[(index.h[i-1]+1):(index.h[i]-1)],
                          uh = abs(0/100*curves.h$Qcurve[(index.h[i-1]+1):(index.h[i]-1)]))
    if ((tail(d.h[[i]]$t,1) >= tgood) & (length(d.h[[i]]$t) >= Nmin.rec) ) {  
      # only curves with period > tgood!!!
      setwd(dir_exe)
      curve_good.h = curve_good.h + 1
      print(curve_good.h)
      # dir.create(paste(dir.regression,"/C",curve_good.h, sep=""))
      # dir.create(paste(dir.regression,"/C",curve_good.h,"/BaM", sep=""))
      # dir.rec <- paste(dir.regression,"/C",curve_good.h,"/BaM", sep="")
      dir.BaM.rec <- paste(dir_code,"/BaM_exe/Recession_h",sep="")  
      hpeakgood.h[curve_good.h] = Qpeak.h[i]; 
      t.real.good.h[curve_good.h] = tpeak.h[i];
      #tpeakgood.h[curve_good.h] = tpeak.h[i];
      index.good.h[curve_good.h] = index.h[i];
      #tcc = tcurve[(index[i-1]+1):index[i]]
      
      #start from the data 
      d.h.selected[[curve_good.h]]  = data.frame(t =curves.h$tcurve[(index.h[i-1]+ Nburn.rec):(index.h[i]-1)],
                                                   h = curves.h$Qcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)],
                                                   uh = uh.rec) 
      d.h.selected.with.true.time[[curve_good.h]]  = data.frame(t =  curves.h$tcurve[(index.h[i-1] + Nburn.rec):(index.h[i]-1)],
                                                                  h =  curves.h$Qcurve[(index.h[i-1] + Nburn.rec):(index.h[i]-1)],
                                                                  uh = uh.rec,
                                                                  treal = curves.h$tcurve[(index.h[i-1]+ Nburn.rec):(index.h[i]-1)] +
                                                                          t.real.good.h[curve_good.h],
                                                                  ind.rec = rep(curve_good.h, length(curves.h$tcurve[(index.h[i-1]+ Nburn.rec):(index.h[i]-1)]))
                                                                  )
        #abs(uh.rec/100*curves.h$Qcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)]))
  
      # if ((rec.model == "1expWithAsymptNorm") |(rec.model == "2expWithAsymptNorm")) {
      #   d.h.selected[[curve_good.h]]  = data.frame(curves.h$tcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)],
      #                                              curves.h$Qcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)]/curves.h$Qcurve[(index.h[i-1]+Nburn.rec)],
      #                                              uh.rec/100*curves.h$Qcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)])
      # }
        
      eliminate.index =0
      for (ccc in 2:length( d.h.selected[[curve_good.h]]$t)) {
        if ((d.h.selected[[curve_good.h]]$t[ccc] - d.h.selected[[curve_good.h]]$t[ccc-1]) < delta.t.min) {
          eliminate.index = c(eliminate.index,ccc)
        }
      }
      if (eliminate.index > 0) {
      d.h.selected[[curve_good.h]]  = d.h.selected[[curve_good.h]][-eliminate.index,]
      d.h.selected.with.true.time[[curve_good.h]]  = d.h.selected.with.true.time[[curve_good.h]][-eliminate.index,]
      }
      #nobs.h <- index.h[i] - 1 - index.h[i-1]
      nobs.h =length(d.h.selected[[curve_good.h]]$t)
      curve_data.h = "Recession_h/Curves_Data.txt"
      write.table(d.h.selected[[curve_good.h]] , 
                  file=curve_data.h, append = FALSE, sep = "\t", eol = "\n", 
                  na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh"))
      if (tail(d.h.selected.with.true.time[[curve_good.h]]$t, 1) > 80) {
        Nrec_longer_100 = Nrec_longer_100 +1
      } 
    }
  }
  df.curves= bind_rows( d.h.selected, .id = "column_label")
  df.curves.t = bind_rows( d.h.selected.with.true.time, .id = "column_label")
  Nrec = length(df.curves.t$t)
  Nk = tail(df.curves.t$ind.rec,1)

  #plot:
  rec.plot.test[[jj]] = plot.extracted.recessions.paper(Data_rec = df.curves.t,
                                                        hmin_grid, 
                                                        hmax_grid, 
                                                        dir.extraction, 
                                                        station.name, 
                                                        data.period, 
                                                        chi = chi.to.test[jj],
                                                        tmax.rec,
                                                        hmax.rec, 
                                                        nobs = Nrec,
                                                        Nperiods = Nk,
                                                        Nrec_longer_100 = Nrec_longer_100,
                                                        stage.scale.shift = stage.scale.shift)
  }
  
  bind.plots =  plot_grid( rec.plot.test[[1]]$plot4paper.extraction, 
                           rec.plot.test[[2]]$plot4paper.extraction, 
                           rec.plot.test[[3]]$plot4paper.extraction,   
                           labels = c('a)', 'b)', 'c)'), 
                           label_size = 30,
                           label_fontface = "plain",
                           ncol = 1, nrow=3)  
  
  pdf(paste0(dir.extraction,"/Figure.pdf"), 12, 14 ,useDingbats=F)
  print(bind.plots)
  dev.off()
}































###############################################################################################################
plot.extracted.recessions.paper <- function(Data_rec, 
                                            stage.record,
                                            hmin_grid,
                                            hmax_grid, 
                                            dir.extraction, 
                                            chi, 
                                            tmax.rec, 
                                            hmax.rec,
                                            nobs, 
                                            Nperiods, 
                                            Nrec_longer_100,
                                            stage.scale.shift) {
###############################################################################################################
  
  # Plot extracted curves in time series:
  #**************************************
  Data_rec_not_shifted   = Data_rec
  Data_rec_not_shifted$h = Data_rec$h - stage.scale.shift
  hmin_grid_not_shifted  = hmin_grid  - stage.scale.shift
  hmax_grid_not_shifted  = hmax_grid  - stage.scale.shift
  
  rec.time.plot <- ggplot(Data_rec_not_shifted, aes(x = treal, y= h, color=ind.rec)) + 
    geom_line(data = stage.record,  aes(x = t_limni, y = h_limni*100), color = "gray90", size =0.2) +
    geom_point(size = 1.2) + 
    coord_cartesian(clip = "off")+
    #geom_point(aes(x=curves$tpeak, y=curves$Qpeak, fill= curves$RCi), size = 4) +
    scale_color_gradientn(colours = rainbow(5)) +
    scale_x_continuous(name = "Time [day]",   expand = c(0,0)) +
    scale_y_continuous(name = "Stage h [cm]", expand = c(0,0))  + #limits = c(hmin_grid_not_shifted, hmax.rec)) +
    annotate("text", 
             x=0,
             y=Inf,
             hjust=0,vjust=1,
             label = TeX((paste0("$\\chi$ =", chi, " cm"))), 
             color = "black", size=7, hjust=0) +
    xlab("Time [day]") +   ylab("Stage [cm]") +  labs(color = "Recession index") +
    theme_bw(base_size = 15) +
    theme(  axis.text        = element_text(size=15)
           ,axis.title       = element_text(size=15, face="plain")
           ,panel.grid.major = element_blank()
           ,panel.grid.minor = element_blank()
           ,legend.text      = element_text(size=20)
           ,legend.title     = element_text(size=30)
           ,legend.key.size  = unit(1.5, "cm")
           ,legend.position  = "none"
           ,plot.margin      = unit(c(0.5, 0.5, 0.5, 2),"cm"))
    # ggsave(rec.time.plot, filename =paste0(dir.extraction,"/time_recessions_",chi,"cm.png"), 
    #       bg = "transparent", width = 6, height =3, dpi = 200)

  
  

  
  # Plot extracted curves against recession time t*:
  #*************************************************
  rec.curves.plot <-   ggplot(Data_rec_not_shifted, aes(x = t, y = h, color = ind.rec)) + 
                       geom_point(size = 1.2) +
                       geom_text( aes(x=0,
                                      y=Inf,
                                      hjust=0,vjust=1,
                                      label=paste0("Nrec = ", nobs, "\nNk = ", Nperiods, "\nNk (t > 80 days) =", Nrec_longer_100)),
                                  color ="black", size  = 4) +
                       # annotate("text", 
                       #          x = 0,
                       #          y = c(hmax.rec - 50, hmax.rec - 100, hmax.rec - 150), 
                       #          label =  c(TeX((paste0("$N_{rec} =$", nobs))), 
                       #                     TeX((paste0("$N_{k} =$", Nperiods))), 
                       #                     TeX((paste0("$N_{k}(t \\; > 80 \\; days) =$", Nrec_longer_100 )))), 
                       #          color = "blue", 
                       #          size  = 5,
                       #          hjust = 0) +
                       labs(color = "# recession,  k")+
                       coord_cartesian(clip = "off")+
                       scale_color_gradientn(colours = rainbow(5)) +
                       scale_x_continuous(name = "Recession time [day]",  expand = c(0,0)) + #, limits =c(0, tmax.rec)) +
                       scale_y_continuous(name = "Stage h [cm]", expand = c(0,0)) + #, limits = c(hmin_grid_not_shifted, hmax.rec) +
                       theme_bw(base_size = 15) +
                       theme( axis.text             = element_text(size=15)
                              ,axis.title           = element_text(size=15)
                              ,panel.grid.major     = element_blank()
                              ,panel.grid.minor     = element_blank()
                              ,legend.text          = element_text(size=10)
                              ,legend.title         = element_text(size=15)
                              ,legend.key.size      = unit(0.3, "cm")
                              ,legend.position      = c(.99, .99)
                              ,legend.justification = c("right", "top")
                              ,legend.box.just      = "right"
                              #,legend.margin       = margin(6, 6, 6, 6)
                              ,plot.margin          = unit(c(0.5,0.5, 0.5, 2),"cm"))
  
  # ggsave(rec.curves.plot, filename =paste0(dir.extraction,"/all_curves_",chi,"cm.png"), 
  #        bg = "transparent", width = 8, height =4, dpi = 300)
  #***********************************************************
  plot4paper.extraction =  plot_grid(rec.time.plot, 
                                     rec.curves.plot,
                                     label_fontface = "plain",
                                     ncol = 2, nrow=1)  
  
  
  #***********************************************************
  return(list(rec.time.plot         = rec.time.plot,
              rec.curves.plot       = rec.curves.plot,
              plot4paper.extraction = plot4paper.extraction))
}











































################################################################################################
plot.recessions =function(dir.case_study, 
                          rec.model, 
                          curves,
                          which.recession,
                          df.h.infinity, 
                          BayesianOption , 
                          limits.y, 
                          limits.x){
  ################################################################################################  
  dir.create(paste(dir.case_study,"/Results/segmentation_recessions/curves_regression", sep=""))
  dir.regression <- paste(dir.case_study,"/Results/segmentation_recessions/curves_regression", sep="")
  colfunc = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))   
  
  
  if (BayesianOption ==1) {
    ####################################### 
    Ncurves= length(curves[[1]])
    # plot of recessions:
    #********************
    if (rec.model == "2expWithAsympt") {
      cat(c("a1","b1","a2","b2","a3"), file = output_file_h, sep=";", append = FALSE)
      cat("", file = output_file_h, sep="\n", append = TRUE)
      reg.plot <- ggplot() + 
        scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]), 
                           expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
        scale_y_continuous(name = "Stage h (cm)", limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        theme_bw(base_size = 15)
      ggsave(reg.plot, filename =paste(dir.regression,"/regression.png", sep=""), 
             bg = "transparent", width = 8, height =4, dpi = 400)
    } else if(rec.model == "3expWithAsympt") {
      cat(c("a1","b1","a2","b2","a3","b3","a4"), file = output_file_h, sep=";", append = FALSE)
      cat("", file = output_file_h, sep="\n", append = TRUE)
      reg.plot <- ggplot() + 
        theme_bw(base_size = 15)+
        scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
                           expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
        scale_y_continuous(name = "Stage h (cm)", limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) 
      #axis.line = element_line(colour = "black"))
      
      ggsave(reg.plot, filename =paste(dir.regression,"/regression.png", sep=""), 
             bg = "transparent", width = 12, height =7, dpi = 200)
    } else if (rec.model == "1expWithAsymptNorm") {
      cat(c("a1","b1"), file = output_file_h, sep=";", append = FALSE)
      cat("", file = output_file_h, sep="\n", append = TRUE)
      reg.plot <- ggplot() + 
        scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]), 
                           expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
        scale_y_continuous(name = "Stage h (cm)", limits = c(0,1), expand = c(0,0)) +
        theme_bw(base_size = 15)+
        theme(text=element_text(size=16,  family="Serif"))
    } else if (rec.model == "2expWithAsymptNorm") { 
      cat(c("a1","b1","a2","b2"), file = output_file_h, sep=";", append = FALSE)
      cat("", file = output_file_h, sep="\n", append = TRUE)
      reg.plot <- ggplot() + 
        scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]), 
                           expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.y[3]))+
        scale_y_continuous(name = "Stage h (cm)", limits = c(0,1), expand = c(0,0)) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        theme_bw(base_size = 20)
    }
    
    temp=NULL
    temp2 =NULL
    for (r in which.recession) {
      temp[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C", r,"/BaM/temp.txt"),header=TRUE)
      temp2[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C", r,"/BaM/temp2.txt"),header=TRUE)
      reg.plot =  reg.plot +
        # geom_point(data = temp , aes(x = x, y= y), color=colfunc(80)[curve_good.h] , size= 0.5)+
        geom_line(data = temp2[[r]], aes(x=xx, y=yy), color=colfunc(Ncurves)[r], size = 0.1)+
        # geom_ribbon(data = temp2, aes(x = xx, ymin = zz ,ymax = kk), 
        # fill = colfunc(80)[curve_good.h], alpha = 0.3) 
        geom_point(data = temp[[r]] , aes(x = x, y= y), color= colfunc(Ncurves)[r], size= 2)+
        geom_ribbon(data = temp2[[r]], aes(x = xx, ymin = zz , ymax = kk), 
                    fill = colfunc(Ncurves)[r], alpha = 0.2)
    }
    ggsave(reg.plot, filename =paste0(dir.regression,"/regression.png"),
           bg = "transparent", width = 12, height =7, dpi = 200)
    
    #***********************************
    # plot of the asymptote time series:
    #***********************************
    df.h.infinity =df.h.infinity[which.recession,]
    asympt.plot =  ggplot()+
      theme_bw(base_size = 15) +
      scale_x_continuous(name="time [day]",expand = c(0,0), limits = c(0,tail(t_limni,1))) +
      scale_y_continuous(name="Asymptotic stage [cm]",
                         limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
      xlab("Time [day]")+ ylab("Asymptotic stage [cm]")+
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      #axis.line = element_line(colour = "black"))
      geom_point(data = df.h.infinity, aes(x = t, y = X50.), colour = "black",size = 2) +
      geom_errorbar(data = df.h.infinity, aes(x= t, ymin = X2.5., ymax =X97.5.), width=50, size = 0.3)
    ggsave(asympt.plot, filename =paste0(dir.regression,"/asymptote.png"),
           bg = "transparent", width = 12, height =6, dpi = 200)#"
    
  } else if (BayesianOption ==2) {
    ####################################  
    ncurves.pooling= length(which.recession)
    if (rec.model == "2expWithAsympt") {
      # for (param in 1:7) {
      #   quantiles.rec[[param]] <- quantile(x = Results_mcmc.pool[,param], p = c(0.025, 0.1, 0.5, 0.9, 0.975)) 
      # }
      col.asympt = seq((2*ncurves.pooling+2  + 1), (2*ncurves.pooling+2 +  ncurves.pooling), 1)
      asymptote.h = Results_summary.pool[16,col.asympt]
      u.asymptote.h = Results_summary.pool[11,col.asympt]
      asymptote.df = data.frame(
        time = t.real.good.h[1], 
        h2 = quantile(Results_mcmc.pool[,col.asympt[1]], 0.025),
        h97 = quantile(Results_mcmc.pool[,col.asympt[1]], 0.975),
        hMAP = asymptote.h[[1]]
      )
      for (f in 2:ncurves.pooling) {
        asymptote.df = rbind(asymptote.df, data.frame(
          time = t.real.good.h[f], 
          h2 = quantile(Results_mcmc.pool[,col.asympt[f]], 0.025),
          h97 = quantile(Results_mcmc.pool[,col.asympt[f]], 0.975),
          hMAP = asymptote.h[[f]]
        ))
      }
      Data.segm.rec <- data.frame("t" = t.real.good.h[1], 
                                  "Y" = asymptote.h[[1]], 
                                  "uY"= u.asymptote.h[[1]]
      )
      for (f in 2:ncurves.pooling) {
        Data.segm.rec = rbind(Data.segm.rec, data.frame("t" = t.real.good.h[f], 
                                                        "Y" = asymptote.h[[f]], 
                                                        "uY"= u.asymptote.h[[f]]
        ))
      } 
      #******************************************
    } else if (rec.model == "3expWithAsympt") {
      #******************************************
      
      
      
      #*****************************************************************************
    }
    ggplot(data=df.h.infinity)+ 
      geom_point(aes(x=t, y=maxpost)) + 
      geom_errorbar(aes(x=t, ymin =`2.5%`, ymax =`97.5%`)) +
      theme_bw() + 
      scale_y_continuous(limits = c(-150, 50))
    
  }
}



































####################################################################################################
#plot the BIC
####################################################################################################
plot.BIC.rec <- function(BIC.rec.df, dir_code) {
  bic.rec.plot <- ggplot(BIC.rec.df) +   
    theme_light(base_size = 10)+
    geom_line(aes(x = y, y = BIC.rec), size = 1, color ="gray") +
    geom_point(aes(x = y, y = BIC.rec), size = 2) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill ="transparent", color = NA),
          panel.background = element_rect(fill ="transparent"), 
          axis.line = element_line(colour = "black"))+
    xlab("Number of segments nS") + ylab("BIC recessions")
  
  ggsave(bic.rec.plot, filename =paste0(dir_code,"/Results/segmentation_recessions/segmentation/BIC.png"), 
         bg = "transparent", device = "png", width = 16, height =8, dpi = 400, units = "in")
}























#######################################################################
plot.segmentation.results.rec = function(dir.results, 
                                         df, 
                                         df.shifts, 
                                         ts.res.before, 
                                         ts.res,
                                         ts.real,
                                         ts.res.plus,  
                                         df.mu, 
                                         mu.res, 
                                         mcmc.segment, 
                                         nS,
                                         x.name, 
                                         limits.X, 
                                         y.name, 
                                         limits.Y,
                                         error.type, 
                                         plot.axis.x, 
                                         plot.axis.y,
                                         plot.title,
                                         plot.mu.uncertainty,
                                         plot.tau.uncertainty,
                                         plot.shift.dates,
                                         plot.pdf,
                                         points.size,
                                         plot.file.name) {
  #######################################################################  
  df$title <- plot.title 
  if (limits.Y[1] != "automatic"){
    error.max = 0
    error.min = 0
    for (iii in 1:length(df$t)) {
      error.min[iii] = max((df$h[iii] - 2*df$uh[iii])) #,  limits.Y[1])
      error.max[iii] = min((df$h[iii] + 2*df$uh[iii])) #,  limits.Y[2])
    }
    df$error.max = error.max
    df$error.min = error.min
  } else {
    df$error.max = df$h + 2*df$uh
    df$error.min = df$h - 2*df$uh
  }  
  
  if (length(mu.res) >1) {
  ########################
    df.shift = data.frame( ts.res.before = ts.res.before, 
                           ts.res.plus   = ts.res.plus,  
                           mu.res        = mu.res )
    df.shift2 = data.frame(ts.res  = ts.res,
                           ts.real = as.numeric(ts.real))
    col.ts.distrib <- rainbow((nS-1)) 
    X1 = mcmc.segment
    X2 = X1[,(nS+1):(2*nS-1)]
    
    if (nS ==2) {
      X = X2
    } else {
      X = data.frame(time    = X2[,1],
                     ord     = rep(1, length(X2[,1])) , 
                     colorrr = rep(col.ts.distrib[1], length(X2[,1])))
      for (orderr in 2:(nS-1)) {
        X =rbind(X,    data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1])), 
                                  colorrr = rep(col.ts.distrib[orderr], length(X2[,1]))) )
      }
    }

    
    # plot 1:
    ##########
    p= ggplot(data = df)
    if (plot.mu.uncertainty ==TRUE){
      p=p+ annotate("rect", 
                    xmin = ts.res.before, 
                    xmax = ts.res.plus, 
                    ymin = df.mu$mu.q2,
                    ymax = df.mu$mu.q97, 
                    fill = "red", alpha=0.2)
    }
    if (plot.tau.uncertainty ==TRUE){
      p=p+  annotate("rect",
                     xmin= df.shifts$tau.q2, 
                     xmax= df.shifts$tau.q97, 
                     ymin=-Inf, 
                     ymax=Inf, 
                     fill="blue", alpha=0.2)
    }
    p = p+
      geom_point( aes(x = t, y = h), size = points.size, color = "black") +
      # geom_line(data = df, aes(x = t, y = Y), size = 0.3, color = "blue")+
      scale_x_continuous(name   = x.name, 
                         expand = c(0,0),  
                         limits = c(limits.X[1], limits.X[2]), 
                         breaks = floor(df.shift2$ts.real))
    if (limits.Y[1] =="automatic"){
      p = p+ 
        scale_y_continuous(name = y.name, 
                           expand = c(0,0))
    } else {
      p = p+ 
        scale_y_continuous(name = y.name, 
                           expand = c(0,0))
                           #limits = c(limits.Y[1], limits.Y[2]),
                           #breaks = seq(limits.Y[1], limits.Y[2], limits.Y[3] ))
    }
    p = p+
      geom_segment(data= df.shift, mapping = aes(x    = ts.res.before , 
                                                 y    = mu.res, 
                                                 xend = ts.res.plus,
                                                 yend = mu.res), col="red") +
      coord_cartesian(clip = 'off')+
      geom_vline( xintercept = df.shift2$ts.res,  col="blue",  lwd =0.2, linetype= "dashed")+
      geom_vline( xintercept = df.shift2$ts.real, col="blue", lwd =0.5, linetype= "solid") +
      # annotate("rect",xmin= ts.res.before, xmax=ts.res.plus, ymin=(mu.res),
      #         ymax=(mu.res+10), fill="red", alpha=1) +
      theme_light(base_size = 15) +
      theme(axis.title.x      = element_blank()
           ,axis.text.x       = element_text(size = 8, color ="blue")
           #,panel.grid.major = element_line(size=1.2)
           #,panel.grid.minor = element_line(size=0.8)
           ,legend.text       = element_text(size=10)
           ,legend.title      = element_text(size=10)
           ,legend.key.size   = unit(1.5, "cm")
           ,legend.position   = "none"
           ,plot.margin       = unit(c(0.5,0.5,0, 0.5),"cm")
           ,axis.title.y      = element_text(size =15, margin = margin(t = 0, r = 25, b = 0, l = 0))
           ,panel.grid        = element_blank())
    
    ############
  } else {
    ############
    df.shift = data.frame( mu.res = mu.res )
    p = ggplot(data =df)+
      geom_point(data =df, aes(x = t, y = h), size = points.size, color = "black") +
      # geom_line(data = df, aes(x = t, y = Y), size = 0.3, color = "blue")+
      geom_segment(data =df.shift, mapping=aes(x =df$t[1] , y = mu.res, xend = tail(df$t,1),
                                               yend = mu.res[1]), col="red") +
      scale_x_continuous(name = x.name, 
                         expand = c(0,0),  
                         limits =c(limits.X[1], limits.X[2]))
    if (limits.Y[1] =="automatic"){
      p = p+ 
        scale_y_continuous(name = y.name,   expand = c(0,0))
    } else {
      p = p+ 
        scale_y_continuous(name = y.name, 
                           expand = c(0,0)) 
                           #limits = c(limits.Y[1], limits.Y[2]),
                           #breaks = seq(limits.Y[1], limits.Y[2], limits.Y[3] ))
    }
    p = p+
      coord_cartesian(clip = 'off')+
      theme_light(base_size = 15) +
      theme(axis.title.x    = element_blank()
           #,panel.grid.major=element_line(size=1.2),panel.grid.minor=element_line(size=0.8)
           ,legend.text     = element_text(size=10)
           ,legend.title    = element_text(size=10)
           ,legend.key.size = unit(1.5, "cm")
           ,legend.position = "none"
           ,plot.margin     = unit(c(0.5, 0.5, 0, 0.5),"cm")
           ,axis.title.y    = element_text(size =15, margin = margin(t = 0, r = 25, b = 0, l = 0))
           ,panel.grid      = element_blank())
    if (plot.mu.uncertainty ==TRUE){
      p = p+ annotate("rect", 
                      xmin = limits.X[1] , 
                      xmax = limits.X[2], 
                      ymin = df.mu$mu.q2,
                      ymax = df.mu$mu.q97, 
                      fill ="red", alpha=0.2)
    }
  }
  if ((error.type == 2) | (error.type == 3)) {
    p=p + geom_errorbar(data = df, 
                        aes(x = t, 
                            ymin = error.min,
                            ymax = error.max), 
                        size = 0.2, width = 1, col= "black")
  }
  if (plot.axis.y == FALSE) {
    p=p + theme(  axis.text.y  = element_blank()             
                 ,axis.ticks.y = element_blank()
                 ,axis.title.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 0)))
  }
  if (plot.axis.x == FALSE) {
    p=p + theme( axis.text.x  = element_blank()             
                ,axis.ticks.x = element_blank())
  }
  p=p  +
    facet_grid(. ~ title)+
    theme(strip.text = element_text(size=15))
  
  
  
  
  
  
  #--------------------------------------------------------------------------------
  if ((plot.pdf == TRUE) & (!is.null(ts.res))) {
  #-------------------------------------------------------------------------------- plot 2
    init.ts.dens = ggplot()+
      scale_x_continuous(name   = x.name, 
                         expand = c(0,0),  
                         limits = c(limits.X[1], limits.X[2])) +
      scale_y_continuous(name =  "Scaled pdf", 
                         expand = c(0,0),
                         breaks = c(0,1)) +
      theme_light(base_size = 15) +
      theme( axis.title       = element_text(size=15)
            ,axis.title.y     = element_blank()
            ,axis.text.y      = element_blank()
            ,axis.ticks.y     = element_blank()
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,legend.text      = element_text(size=20)
            ,legend.title     = element_text(size=30)
            ,legend.key.size  = unit(1.5, "cm")
            ,legend.position  = "none"
            ,plot.margin      = unit(c(0.3, 0.5, 0.2, 2.9),"cm"))
    
    if (length(ts.res)==1) {
      init.ts.dens = init.ts.dens + geom_density(aes(x= X, ..scaled..),
                                                 fill= "blue", 
                                                 colour="blue", alpha=0.3)
    } else {
      
      init.ts.dens = init.ts.dens + 
        geom_density(aes(x= X$time, ..scaled.., group =X$ord),  
                     fill= "blue", colour="blue", alpha=0.3)
    }
    # merge the two plots:     #needs cowplot "package"
    initial.ts.plot2 = plot_grid( p, 
                                  init.ts.dens, 
                                  ncol = 1, 
                                  nrow = 2, 
                                  rel_heights = c(1, 0.4))  
    ggsave(initial.ts.plot2, filename =paste0(dir.results,"/", plot.file.name,".png"), bg = "transparent", width = 8, height =5, dpi = 300)
    return(initial.ts.plot2)
    ######
  } else {
    ######
    ggsave(p, filename =paste0(dir.results, "/", plot.file.name, ".png"),bg = "transparent", width = 8, height =5, dpi = 300)
    return(p)
  }
}







































##############################################################################################################
segment.asymptote.plot <- function(dir.segm.recessions, df.asymptote, data.annotate.off, seg.rec, df.limni) {
##############################################################################################################
  asympt.plot = ggplot()+ 
    geom_line(data = df.asymptote ,aes(x = df.asymptote$t, y = df.asymptote$hinfinity), size = 0.1, color = "darkgray") +
    geom_errorbar(data = df.asymptote, aes(x = df.asymptote$t, ymin= df.asymptote$h2, ymax = df.asymptote$h97), 
                  size = 0.2, width=60, color = "black") +
    geom_point(data= df.asymptote, aes(x = df.asymptote$t, y = df.asymptote$hinfinity), size = 1)+
    theme_light(base_size = 15)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,7000)) +    
    scale_y_continuous(name = TeX("$ Asymptotic \\; stage , \\; h_{\\infty} \\; \\left[ cm \\right] $"), 
                       expand = c(0,0),limits =c(-100,20))+ 
    coord_cartesian(clip = 'off') +
    annotate("rect",xmin= seg.rec[[1]]$ts.q2, xmax=seg.rec[[1]]$ts.q97,ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2) +
    annotate("rect",xmin= c(0,seg.rec[[1]]$t.MAP.before), xmax=c(seg.rec[[1]]$t.MAP.plus, tail(df.limni$t_limni,1)), 
             ymin=seg.rec[[2]]$q10, ymax=seg.rec[[2]]$q90, fill="gray", alpha=0.5) +
    geom_segment(mapping= aes(x =c(0,seg.rec[[1]]$t.MAP.before) , y = seg.rec[[2]]$MAP, 
                              xend = c(seg.rec[[1]]$t.MAP.plus,tail(df.limni$t_limni,1)), 
                              yend = seg.rec[[2]]$MAP), 
                 color = "black", size = 1) +
    geom_vline(xintercept = seg.rec[[1]]$t.real, col="red", lwd =0.5, linetype = "dashed")+
    geom_vline(xintercept = seg.rec[[1]]$t.MAP, col="blue", lwd =0.5, linetype = "solid")+           
    theme(plot.background = element_rect(fill ="transparent", color = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill ="transparent"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"))+
    theme(plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"))+
    geom_point(data=data.annotate.off, aes(x = x, y=-100),shape=4, size=3, color = "black", stroke = 1.5)
  ggsave(asympt.plot, filename= paste0(dir.segm.recessions,"/asympt.png"), bg = "transparent",width = 7, height =3.5, dpi = 400)
  
}
























##################################################################################################################
plot.time.shifts.recess <- function(dir.case_study, data.annotate.off, data.annotate.gaug, data.annotate.gaug.adjust, data.annotate.recess, 
                                    data.annotate.recess.adjust, data.annotate.step1, colo, t_Gaug, h_Gaug,
                                    df.limni, limni.labels
) {
  ##################################################################################################################
  
  # grouping the gaugings per period (considering the final subdivision):
  color.step1 = 0; Period.step1 = 0; t.shift.with.gaug.step1 =0
  if (!is.null(data.annotate.step1$t.step1[1])) {
    for (i in 1:length(t_Gaug)) { 
      if(t_Gaug[i] <= data.annotate.step1$t.step1[1]) {
        color.step1[i] = colo[1]
        Period.step1[i] = 1
        t.shift.with.gaug.step1[i] = data.annotate.step1$t.step1[1]
      }
    }
    
    for (j in 2:(length(data.annotate.step1$t.step1) +1)) {
      for (i in 1:length(t_Gaug)) { 
        if ((t_Gaug[i] <= tail(df.limni$t_limni,1)) & (t_Gaug[i] > data.annotate.step1$t.step1[j-1])) {
          #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
          color.step1[i] = colo[j]
          Period.step1[i] = j
          t.shift.with.gaug.step1[i] = data.annotate.step1$t.step1[j]
        }
      }
    }
  } else {
    for (i in 1:length(t_Gaug)) {
      #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
      color.step1[i] = colo[1]
      Period.step1[i] = 1
    }
  }
  
  
  #plot:
  #########################
  t.plot <- ggplot()
  if (is.null(df.limni)==FALSE) {
    t.plot = t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "blue",size = 0.2)
  }
  
  if (exists("data.annotate.off")==TRUE)  {
    if (is.null("data.annotate.off")==FALSE) {
      t.plot <- t.plot +
        geom_point(data = data.annotate.off, aes(x = x, y = start), color= "black", size = 4, shape =4, stroke=2 )
    }
  }
  t.plot <- t.plot + geom_point( aes(x = t_Gaug , y= h_Gaug), color =color.step1, size = 2.5)
  
  t.plot <- t.plot + 
    #Step1:
    annotate("rect", xmin= data.annotate.gaug$q2, xmax=data.annotate.gaug$q97, ymin=data.annotate.gaug$start, 
             ymax=data.annotate.gaug$finish, fill="blue", alpha=0.1) +
    geom_segment(data = data.annotate.gaug, aes(x = MAP, y = start, yend = finish, xend= MAP), size = 0.7, color ="blue")+
    geom_segment(data = data.annotate.gaug.adjust, aes(x = t.adj, xend = t.adj, y = start, yend = finish), size = 0.7, color ="blue")+
    annotate("rect", xmin= data.annotate.recess$q2, xmax=data.annotate.recess$q97, ymin=data.annotate.recess$start, 
             ymax=data.annotate.recess$finish, fill="green", alpha=0.1) +
    geom_segment(data = data.annotate.recess, aes(x = MAP, y = start, yend = finish, xend = MAP), size = 0.7, color ="green") +
    geom_segment(data = data.annotate.recess.adjust, aes(x = t.adj, y = start, yend = finish, xend = t.adj), size = 0.7, color ="green") +
    geom_segment(data = data.annotate.step1, aes(x = t.step1, y = start, yend = finish, xend = t.step1), size = 0.7, color ="red")+
    scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1))) +
    scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(-4,5), breaks=seq(-1,5,1))+
    xlab(limni.labels[1])+ 
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    #theme_bw(base_size = 10)+theme(panel.grid = element_blank())+
    geom_hline(yintercept = c(-1, -1.5, -2 , -2.5, -3, -3.5, -4), color="darkgray", linetype="dashed", size = 0.5)+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none"
          , plot.margin = margin(t=1, b=1, l=0.1, r=0, "cm"))
  #******************
  annot <- ggplot() +
    scale_x_continuous(name=element_blank(),  limits = c(0,2))+
    scale_y_continuous(name=element_blank(), expand = c(0,0), limits = c(-4,5))+
    annotate(geom="text", x=0, y=-3.75, size= 7, label="Official dates of RC update",  color="black",  hjust = 0)+
    annotate(geom="text", x=0, y=-1.25, size= 7, label="1. Segmentation of gaugings",   color="blue", hjust = 0)+
    annotate(geom="text", x=0, y=-1.75, size= 7, label="2. Adjusted to real events",   color="blue",  hjust = 0) +
    annotate(geom="text", x=0, y=-2.25, size= 7, label="3. Segmentation of recessions",   color="green",  hjust = 0)+
    annotate(geom="text", x=0, y=-2.75, size= 7, label="4. Adjusted to real events",   color="green",  hjust = 0) +
    annotate(geom="text", x=0, y=-3.25, size= 7, label="5. Results of Step 1",   color="red",  hjust = 0)+
    #annotate(geom="text", x=1, y=-0.75, size= 7, label="Different used methods:",   color="black")+
    #theme_bw(base_size = 10)+theme(panel.grid = element_blank())+
    geom_hline(yintercept = c(-1, -1.5, -2 , -2.5, -3, -3.5,-4), color="gray", size = 0.5)+
    geom_text() +
    theme_void()
  
  library(cowplot)
  tt.plot = plot_grid( t.plot, annot, ncol = 2, rel_widths = c(1, 0.32), align = 'h')
  #library(egg)      
  #tt.plot = ggarrange(t.plot, annot, ncol = 2)
  ggsave(tt.plot, filename =paste0(dir.case_study,"/Results/time_series_step1.png"), 
         device = "png", width = 16, height =8, dpi = 600, units = "in")
  
  
}





































##################################################################################################### plot asymptote time series
plot.seg.asympt.time <- function(dir.seg.rec, ts.res.rec, asym.df, asym.summer.df, asym.winter.df, 
                                 Q10.ts.rec, Q90.ts.rec, ts.res.before.rec, ts.res.plus.rec, Q10.mu.res.rec,
                                 Q90.mu.res.rec, mu.res.rec, b2.from.SPD, q10, q90, t.real.good.h, asym.h.maxpost) { 
  ################################################################################################################
  #par(mar = c(1,1,1,1))
  if (!is.null(ts.res.rec[1])) {
    ppp <- ggplot()+
      geom_line(data = asym.df ,aes(x = t, y =maxpost), size = 0.4, color = "darkgray") +
      geom_errorbar(data = asym.df ,aes(x = t, ymin= q2, ymax = q97), 
                    size = 0.2, width=60, color = "black") +
      geom_point(data = asym.summer.df ,aes(x = t, y = h), size = 1.5, color = "black") +
      geom_point(data = asym.winter.df ,aes(x = t, y = h), size = 1.5,shape = 21, colour = "black", fill = "white") +
      theme_light(base_size = 15)+
      #geom_vline(xintercept = officialShiftsTime, size=0.5, color = "red")+    
      theme(plot.title = element_text(hjust = 0.5))+
      # scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,tail(t.real.good.h,1)+1000)) +
      scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,5000)) +    
      #scale_x_date(name = "time [days]", limits = c(as.Date("2001-01-01"), as.Date("2015-01-01"))) + 
      # scale_y_continuous(name = TeX("Stage asymptote ($h_\\infty $ [cm])"), limits =c(min(asym.h.maxpost),max(asym.h.maxpost)))+ 
      scale_y_continuous(name = TeX("$ Asymptotic \\; stage , \\; h_{\\infty} \\; \\left[ cm \\right] $"), 
                         expand = c(0,0),limits =c(-100,50))+ 
      coord_cartesian(clip = 'off') +
      # annotate("rect",xmin= Q10.ts.rec, xmax=Q90.ts.rec,
      #           ymin=-Inf, ymax=Inf, fill="forestgreen", alpha=0.2) + 
      # geom_vline(xintercept = ts.res.rec, col="forestgreen", lwd =0.7)+
      
      # annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=(mu.res.rec),
      #          ymax=(mu.res.rec +1), fill="red", alpha=1) +
      annotate("rect",xmin= Q10.ts.rec, xmax=Q90.ts.rec,
               ymin=-Inf, ymax=Inf, fill="black", alpha=0.3) +
      annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=Q10.mu.res.rec,
               ymax=Q90.mu.res.rec, fill="gray", alpha=0.5) +
      #geom_segment(data = offic.times, aes(x = t, y = -Inf, yend = b2.from.SPD.before , xend = t),
      #             linetype = "dotdash", size = 0.2, color ="black")+
      #geom_segment(data= seg.results, aes(x = t, y = -Inf, yend = h, xend = t),linetype = "longdash", size = 0.2)+
      # annotate("rect",xmin= winter, xmax=summer,
      #          ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2) +
      # annotate("rect",xmin= summer.plus, xmax=winter.plus,
      #          ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
      geom_vline(xintercept = ts.res.rec, col="black", lwd =0.5, linetype = "longdash")+
      #geom_vline(xintercept = offic.times$t, col="black", lwd =0.3, linetype = "dotdash")+
      annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=(mu.res.rec),
               ymax=(mu.res.rec +0.5), fill="gray45", alpha=1) +
      annotate("rect",xmin= offic.times.before , xmax=offic.times.plus, ymin= b2.from.SPD,
               ymax=b2.from.SPD +1, fill="black", alpha=1) +
      theme(plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"))+
      theme(plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"))+
      geom_point(data=offic.times, aes(x = t, y=h),shape=4, size=3, color = "black", stroke = 1.5)
    ggsave(ppp, filename =paste(dir.seg.rec,"/asymptote_segmentation2.png", sep=""), bg = "transparent", width = 7, height =5, dpi = 800)
  } else {  
    ppp <- ggplot()+
      geom_point(data = asym.df ,aes(x = t, y = h), size = 1.2) +
      geom_errorbar(data = asym.df, aes(x=t, ymin= (q10), ymax = (q90)), size = 0.4, width=30, col= "black") +
      #annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=(mu.res.rec),
      #         ymax=(mu.res.rec +1), fill="red", alpha=1) +
      theme_light(base_size = 10)+
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,tail(t.real.good.h,1))) +
      scale_y_continuous(name = TeX("$h_\\infty [cm]$"), expand = c(0,0), limits =c(min(asym.h.maxpost),max(asym.h.maxpost)))+
      theme(plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"), 
            axis.line = element_line(colour = "black"))
    ggsave(ppp, filename =paste(dir.seg.rec,"/asymptote_segmentation.png", sep=""), bg = "transparent", width = 6, height =3.5, dpi = 800)
  }  
}























########################################################################################
# Plot stage h time series  with shift times (first set + second set +official times) :
#****************************************************************************************
plot.ts.recession <- function(officialShiftsTime, times.of.shift, times.of.shift.rec, ht.df, gaug.df,
                              t.q10, t.q90, Q10.ts.rec, Q90.ts.rec, dir_code ) {
  data.annotate.off <- data.frame(start = -2, finish = -1, x = officialShiftsTime)
  data.annotate.gaug <- data.frame(start = -2, finish = -1.5, x = times.of.shift[1:(length(times.of.shift)-1)] )
  data.annotate.rec <- data.frame(start = -1.5, finish = -1, x = times.of.shift.rec[1:(length(times.of.shift.rec))])
  
  tplot.rec <- ggplot()+
    #geom_vline(xintercept = times.of.shift[1:(length(times.of.shift)-1)], col= "blue", lwd = 0.3)+
    #geom_vline(xintercept = times.of.shift.rec[1:(length(times.of.shift.rec))], col= "green3", lwd = 0.3)+
    #geom_vline(xintercept = officialShiftsTime, col= "red", lwd = 0.2,lty =1)   +
    geom_line(data = ht.df, aes(x = t_limni, y = h_limni), color = "darkgrey", size = 0.1)+
    geom_point(data = gaug.df, aes(x = t_Gaug, y = h_Gaug), color = "black", size = 0.1) +
    geom_segment(data = data.annotate.off, aes(x = x, y = start, yend = finish, xend = x, 
                                               color = "Official segment."),size = 0.2)+
    geom_segment(data = data.annotate.gaug, aes(x = x, y = start, yend = finish, xend = x, 
                                                color = "Gaugings segment."),size = 0.2)+
    geom_segment(data = data.annotate.rec, aes(x = x, y = start, yend = finish, xend = x, 
                                               color = "Recession segment."),size = 0.2)+
    theme_light(base_size = 10) +  
    scale_x_continuous(name="time (day)",expand = c(0,0)) +
    scale_y_continuous(name="Stage h (m)",limits = c(-2,max(h_Gaug)))+
    xlab("Time [day]")+ ylab("Stage h [cm]")+
    # annotate("rect", xmin= (ts.res.Q10.P2), xmax=(ts.res.Q90.P2),
    #          ymin=-Inf, ymax=Inf, fill="blue", alpha=0.1)+
    #geom_vline(xintercept = ts.res.P2, col= "blue", lwd = 0.6,lty =1)+
    # annotate("rect", xmin= (Q10.ts), xmax=(Q90.ts),
    #          ymin=-Inf, ymax=Inf, fill="blue", alpha=0.1)+
    #geom_vline(xintercept = ts.res, col= "blue", lwd = 0.6, lty =1)+  
    
    # annotate("rect", xmin= (officialShiftsTime), 
    #       xmax=(officialShiftsTime+10),ymin=0, ymax=1, 
    #       fill="red", alpha=1)+
    annotate("rect", xmin= t.q10, xmax=t.q90,
             ymin=-2, ymax=-1.5, fill="blue", alpha=0.2) +
    # annotate("rect", xmin= times.of.shift[1:(length(times.of.shift)-1)], 
    #          xmax=times.of.shift[1:(length(times.of.shift)-1)]+10,
    #       ymin=0, ymax=0.5, fill="blue", alpha=1) +
    annotate("rect", xmin= (Q10.ts.rec), xmax=(Q90.ts.rec),
             ymin=-1.5, ymax=-1, fill="green", alpha=0.2)+
    # annotate("rect", xmin= (times.of.shift.rec[1:(length(times.of.shift.rec))]), 
    #          xmax=(times.of.shift.rec[1:(length(times.of.shift.rec))]+10),ymin=0.5, ymax=1, 
    #          fill="green", alpha=1)+
    labs(colour="Legend")+
    theme(plot.title = element_text(hjust = 0.5),
          #plot.background = element_rect(fill ="transparent", color = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #panel.background = element_rect(fill ="transparent"), 
          axis.line = element_line(colour = "black"),
          legend.position = "bottom",#c(1, 1), legend.justification = c(1, 1), 
          legend.key.size = unit(5, "mm"), 
          legend.text = element_text(size = 4), 
          legend.title=element_text(size=4, face ="bold"))+
    #legend.background = element_rect(#fill="lightblue",
    #  size=0.1, linetype="solid", 
    #  colour ="black")) +
    scale_colour_manual(values=c("blue","red","green"))
  #guides(fill = guide_legend(order = 2, override.aes = list(shape = 15, size = 9)), shape = guide_legend(order = 1))
  #annotation_logticks(sides="l")
  ggsave(tplot.rec, filename =paste(dir_code,"/Results/segmentation_recessions/time_series.png", sep=""), bg = "transparent",
         width = 4, height =2, dpi = 300)
}












# ########################################################################################
# # Plot stage h time series  with shift times (firste set + second set +official times) :
# #****************************************************************************************
# t2.plot <- ggplot()+
#    geom_vline(xintercept = times.of.shift[1:(length(times.of.shift)-1)], col= "blue", lwd = 0.8)+
#    geom_vline(xintercept = times.of.shift.rec[1:(length(times.of.shift.rec))], col= "green3", lwd = 0.8)+
#    #geom_point(data = CdT, aes(x = t_Gaug, y = Q_Gaug), colour = "black",size = 2.5) +
#    geom_line(data = ht.df, aes(x = t_limni, y = h_limni), color = "darkgrey",size = 0.1)+
#    theme_light(base_size = 20) 
#    for (i in 1:length(RCmaxpost)) {
#       if (!is.null(RCmaxpost[[i]])) {
#          t2.plot <- t2.plot + geom_point(data = gaug[[i]], aes(x = X.tP. , y= X.h.), size = 0.4, color =colo[i])
#       }
#    }
#    t2.plot <- t2.plot + scale_x_continuous(name="time (day)",expand = c(0,0)) +
#    scale_y_continuous(name="Stage h (m)",limits = c(min(h_Gaug),max(h_Gaug)),expand = c(0,2)) +
#    #ggtitle("Discharge time series and shift detection, Ardeche River at Meyras (2000-2014)") +
#    xlab("Time [day]")+ ylab("Discharge Q [m3/s]")+
#    theme(plot.title = element_text(hjust = 0.5))+
#    #geom_vline(xintercept = ts.final, col= "blue", lwd = 1.4)+
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#    annotate("rect", xmin= t.q10, xmax=t.q90,
#            ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2)+
#    annotate("rect", xmin= t.q10.rec, xmax=t.q90.rec,
#               ymin=-Inf, ymax=Inf, fill="green3", alpha=0.2)
#    ggsave(t2.plot, filename =paste(dir_code,"/Results/segmentation_recessions/time_series.png", sep=""), 
#        bg = "transparent", width = 6, height =3.5, dpi = 800)
# 
# p <- ggplot()+
#      geom_vline(xintercept = times.of.shift[1:(length(times.of.shift)-1)], col= "blue", lwd = 0.8)+
#      geom_vline(xintercept = times.of.shift.rec[1:(length(times.of.shift.rec))], col= "green3", lwd = 0.8)+
#      geom_line(data = ht.df, aes(x = t_limni, y = h_limni), color = "darkgrey",size = 0.1)+
#      geom_point(data = CdT, aes(x = t_Gaug, y = h_Gaug), colour = "black",size = 0.3) +
#      theme_light(base_size = 10) +  
#      scale_x_continuous(name="time (day)",expand = c(0,0)) +
#      scale_y_continuous(name="Stage h (m)",limits = c(min(h_Gaug),max(h_Gaug)),expand = c(0,2))+
#      xlab("Time [day]")+ ylab("Stage h [cm]")+
#      theme(plot.title = element_text(hjust = 0.5))+
#      theme(plot.background = element_rect(fill ="transparent", color = NA),
#            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#            panel.background = element_rect(fill ="transparent"), 
#            axis.line = element_line(colour = "black"))+
#      # annotate("rect", xmin= (ts.res.Q10.P2), xmax=(ts.res.Q90.P2),
#      #          ymin=-Inf, ymax=Inf, fill="blue", alpha=0.1)+
#      geom_vline(xintercept = ts.res.P2, col= "blue", lwd = 0.6,lty =1)+
#      # annotate("rect", xmin= (Q10.ts), xmax=(Q90.ts),
#      #          ymin=-Inf, ymax=Inf, fill="blue", alpha=0.1)+
#      geom_vline(xintercept = ts.res, col= "blue", lwd = 0.6, lty =1)+
#      # annotate("rect", xmin= (Q10.ts.rec), xmax=(Q90.ts.rec),
#      #          ymin=-Inf, ymax=Inf, fill="red", alpha=0.1)+
#      geom_vline(xintercept = ts.res.rec, col= "red", lwd = 0.6,lty =1)+
#      geom_vline(xintercept = officialShiftsTime, col= "green", lwd = 0.6,lty =1)
#      #annotation_logticks(sides="l")
#    ggsave(p, filename ="final_ts.png", bg = "transparent",
#           width = 4, height =2, dpi = 800)
#    
#   











#############################################################################################
# plot asymptote time series with segmentation
plot.asympt.segment <- function(officialShiftsTime, asym.df, t.real.good.h, b2.from.SPD, b2.from.SPD.before) {
  #############################################################################################  
  offic.times <- data.frame(t = officialShiftsTime, h=-100)
  offic.times.before <- c(0, officialShiftsTime[1:(length(officialShiftsTime))])
  offic.times.plus <- c(officialShiftsTime[1:(length(officialShiftsTime))],tail(t.real.good.h,1))
  b2.from.SPD <- c(-4,-5,-47,-50,-58,-70)
  b2.from.SPD.before <- c(-4,-5,-47,-50,-58)
  
  # summer <- c(seq(144,5619,365))     # with summer starting from april:
  # summer.plus <- c(seq(144,5254,365))
  # winter <- c(0,seq(327,5437,365))
  # winter.plus <- c(seq(327,5437,365))
  summer <- c(seq(206,5681,365))       # with summer starting from june:
  summer.plus <- c(seq(206,5316,365))
  winter <- c(0,seq(327,5437,365))
  winter.plus <- c(seq(327,5437,365))
  
  kkkk = 0; tttt = 0; winter.index = 0; summer.index = 0;
  for (iiii in 1:length(asym.df$t)) {
    for (jjjj in 1:length(winter)) {
      if ((asym.df$t[iiii] > winter[jjjj]) & (asym.df$t[iiii] <= summer[jjjj])) {
        kkkk = kkkk+1
        winter.index[kkkk] = iiii
      }
    }
  }
  for (iiii in 1:length(asym.df$t)) {
    for (jjjj in 1:length(summer)) {
      if ((asym.df$t[iiii] > summer[jjjj]) & (asym.df$t[iiii] <= winter.plus[jjjj])) {
        tttt = tttt+1
        summer.index[tttt] = iiii
      }
    }
  }
  
  seg.results <- data.frame(t = ts.res.rec, h= mu.res.rec[1:(length(mu.res.rec)-1)])
  # asym.winter.df <- data.frame(t = asymptote.pool$t[winter.index], h = asymptote.pool$MAP[winter.index], 
  #                              sd = asymptote.pool$stdev[winter.index], q10 = asymptote.pool$q10[winter.index], 
  #                              q90 = asymptote.pool$q90[winter.index])
  # asym.summer.df <- data.frame(t = asymptote.pool$t[summer.index], h = asymptote.pool$MAP[summer.index], 
  #                              sd = asymptote.pool$stdev[summer.index], q10 = asymptote.pool$q10[summer.index], 
  #                              q90 = asymptote.pool$q90[summer.index])
  asym.winter.df <- data.frame(t = asym.df$t[winter.index], h = asym.df$maxpost[winter.index], 
                               uh = asym.df$stdev[winter.index], q2 = asym.df$q2[winter.index], 
                               q97 = asym.df$q97[winter.index])
  asym.summer.df <- data.frame(t = asym.df$t[summer.index], h = asym.df$maxpost[summer.index], 
                               uh = asym.df$stdev[summer.index], q2 = asym.df$q2[summer.index], 
                               q97 = asym.df$q97[summer.index])
  
  ###################################################################### plot asymptote time series
  #par(mar = c(1,1,1,1))
  if (!is.null(ts.res.rec[1])) {
    ppp <- ggplot()+
      geom_line(data = asym.df ,aes(x = t, y =maxpost), size = 0.4, color = "darkgray") +
      geom_errorbar(data = asym.df ,aes(x = t, ymin= q2, ymax = q97), 
                    size = 0.2, width=60, color = "black") +
      geom_point(data = asym.summer.df ,aes(x = t, y = h), size = 1.5, color = "black") +
      geom_point(data = asym.winter.df ,aes(x = t, y = h), size = 1.5,shape = 21, colour = "black", fill = "white") +
      theme_light(base_size = 15)+
      #geom_vline(xintercept = officialShiftsTime, size=0.5, color = "red")+    
      theme(plot.title = element_text(hjust = 0.5))+
      # scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,tail(t.real.good.h,1)+1000)) +
      scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,5000)) +    
      #scale_x_date(name = "time [days]", limits = c(as.Date("2001-01-01"), as.Date("2015-01-01"))) + 
      # scale_y_continuous(name = TeX("Stage asymptote ($h_\\infty $ [cm])"), limits =c(min(asym.h.maxpost),max(asym.h.maxpost)))+ 
      scale_y_continuous(name = TeX("$ Asymptotic \\; stage , \\; h_{\\infty} \\; \\left[ cm \\right] $"), 
                         expand = c(0,0),limits =c(-100,50))+ 
      coord_cartesian(clip = 'off') +
      # annotate("rect",xmin= Q10.ts.rec, xmax=Q90.ts.rec,
      #           ymin=-Inf, ymax=Inf, fill="forestgreen", alpha=0.2) + 
      # geom_vline(xintercept = ts.res.rec, col="forestgreen", lwd =0.7)+
      
      # annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=(mu.res.rec),
      #          ymax=(mu.res.rec +1), fill="red", alpha=1) +
      annotate("rect",xmin= Q10.ts.rec, xmax=Q90.ts.rec,
               ymin=-Inf, ymax=Inf, fill="black", alpha=0.3) +
      annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=Q10.mu.res.rec,
               ymax=Q90.mu.res.rec, fill="gray", alpha=0.5) +
      #geom_segment(data = offic.times, aes(x = t, y = -Inf, yend = b2.from.SPD.before , xend = t),
      #             linetype = "dotdash", size = 0.2, color ="black")+
      #geom_segment(data= seg.results, aes(x = t, y = -Inf, yend = h, xend = t),linetype = "longdash", size = 0.2)+
      # annotate("rect",xmin= winter, xmax=summer,
      #          ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2) +
      # annotate("rect",xmin= summer.plus, xmax=winter.plus,
      #          ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
      geom_vline(xintercept = ts.res.rec, col="black", lwd =0.5, linetype = "longdash")+
      #geom_vline(xintercept = offic.times$t, col="black", lwd =0.3, linetype = "dotdash")+
      annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=(mu.res.rec),
               ymax=(mu.res.rec +0.5), fill="gray45", alpha=1) +
      annotate("rect",xmin= offic.times.before , xmax=offic.times.plus, ymin= b2.from.SPD,
               ymax=b2.from.SPD +1, fill="black", alpha=1) +
      theme(plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"))+
      theme(plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"))+
      geom_point(data=offic.times, aes(x = t, y=h),shape=4, size=3, color = "black", stroke = 1.5)
    ggsave(ppp, filename =paste(dir.seg.rec,"/asymptote_segmentation2.png", sep=""), bg = "transparent", width = 7, height =5, dpi = 800)
  } else {  
    ppp <- ggplot()+
      geom_point(data = asym.df ,aes(x = t, y = h), size = 1.2) +
      geom_errorbar(data = asym.df, aes(x=t, ymin= (q10), ymax = (q90)), size = 0.4, width=30, col= "black") +
      #annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=(mu.res.rec),
      #         ymax=(mu.res.rec +1), fill="red", alpha=1) +
      theme_light(base_size = 10)+
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,tail(t.real.good.h,1))) +
      scale_y_continuous(name = TeX("$h_\\infty [cm]$"), expand = c(0,0), limits =c(min(asym.h.maxpost),max(asym.h.maxpost)))+
      theme(plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"), 
            axis.line = element_line(colour = "black"))
    ggsave(ppp, filename =paste(dir.seg.rec,"/asymptote_segmentation.png", sep=""), bg = "transparent", width = 6, height =3.5, dpi = 800)
  }  
}










  
  
  
   




































############################################################################################################
plot.segmentation.recession = function(dir.segm.recessions,
                                       read.reg.rec,
                                       read.res.rec, 
                                       df.limni,
                                       limits.X, limits.Y, 
                                       x.name, y.name,
                                       limni.labels,
                                       grid_limni.ylim,
                                       obs.uncertainty.y,
                                       plot.gamma.uncertainty,
                                       plot.b.from.gaugings,
                                       #**********************
                                       dir.case_study,
                                       rec.model, 
                                       curves,  
                                       which.recession,
                                       which.recession.long,
                                       df.h.infinity,
                                       BayesianOption,
                                       limits.y.recess,  # stage is in "cm"
                                       limits.x.recess,
                                       stage.scale.shift ,
                                       which.recession.realtime,
                                       starting.time
                                       ) {
############################################################################################################
  # Quantiles for tot uncertainty:
  Ysim=list(); 
  data.tmp= list();
  for (jj in 1:read.res.rec$nS.ok) {
    Ysim[[jj]] =0
    for (ii in 1: length(read.res.rec$mcmc.segment[,jj])) {
      Ysim[[jj]][ii] = read.res.rec$mcmc.segment[ii,jj] + rnorm(1, 
                                                                 mean = 0, 
                                                                 sd = read.res.rec$mcmc.segment[,2*read.res.rec$nS.ok])
    }
    data.tmp[[jj]] = quantile(Ysim[[jj]], probs=c(0.025,0.975) ) #, na.rm=TRUE)
  } 
  # preparing vector with pdf of shift times:
  if (!is.null(read.res.rec$pdf.ts.rec)) {
    X1=read.res.rec$pdf.ts.rec
    if (ncol(read.res.rec$pdf.ts.rec)==1) {
      X = X1
    } else {
      X =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(read.res.rec$pdf.ts.rec))) {
        X =rbind(X, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
    }
  }
  # preparing  the dataframe for riverbed estimated through gaugings:
  if (plot.b.from.gaugings == TRUE){
     dir.sed.transp = paste0(dir.case_study,"/Results/segmentation_sed_transp")
     dir.sed.transp.SPD = paste0(dir.sed.transp,"/SPD")
     bt1.df = read.table(paste0(dir.sed.transp.SPD,"/bt1_df.txt"), header =TRUE)
     bt2.df = read.table(paste0(dir.sed.transp.SPD,"/bt2_df.txt"), header =TRUE)
     ts.gaug =  read.table(paste0(dir.sed.transp.SPD,"/STEP1_shift_times.txt"), header =TRUE)
     ts.before.gaug =  sort(c(0, ts$t.adj))
     ts.plus.gaug = sort(c(ts$t.adj, tail(df.limni$t_limni,1)))
  }
  
  
  #********************************************************************************
  # First plot: time series of asymptote of recessions with the shift times and b(t)
  if ( read.res.rec$nS.ok >1) {
    seg.rec.plot <- ggplot()
    
    if (obs.uncertainty.y == TRUE) {
      seg.rec.plot = seg.rec.plot +  
        geom_errorbar(data = read.reg.rec$df.h.infinity,
                      aes(x=t, 
                          ymin= (mean - 2*stdev), 
                          ymax = (mean+2*stdev)),
                      size = 0.3, width=40, col= "black")
    }
    time.adjust.before = c(0, read.res.rec$data.annotate.recess.adjust$t.adj)
    time.adjust.plus   = c(read.res.rec$data.annotate.recess.adjust$t.adj, limits.X[2])
    
    if (plot.gamma.uncertainty ==TRUE){
      for (jj in 1:read.res.rec$nS.ok){
        seg.rec.plot = seg.rec.plot + 
                        annotate("rect", 
                                 xmin = time.adjust.before[jj], 
                                 xmax = time.adjust.plus[jj], 
                                 ymin = data.tmp[[jj]][1],
                                 ymax = data.tmp[[jj]][2],
                                 fill = "pink", alpha = 0.3)
      }
    }
    seg.rec.plot = seg.rec.plot + 
                   annotate("rect", 
                            xmin = time.adjust.before, 
                            xmax = time.adjust.plus, 
                            ymin = read.res.rec$Q2.mu,
                            ymax = read.res.rec$Q97.mu, 
                            fill ="red", alpha=0.3) +
      #geom_line(data = Data.segm.rec, aes(x = t, y = Y), size = 0.2, color = "blue")+
      geom_point(data = read.reg.rec$df.h.infinity, aes(x = t, y = mean), size = 3, pch=21, fill="white") +
      theme_light(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(name = x.name, expand = c(0,0),  limits =limits.X) +
      scale_y_continuous(name = y.name, expand = c(0,0),  limits =limits.Y) +
      # annotate("rect",
      #          xmin= read.res.rec$Q2.ts, 
      #          xmax=read.res.rec$Q97.ts, 
      #          ymin=-Inf, 
      #          ymax=Inf, fill="blue", alpha=0.2) +
      geom_segment(mapping=aes(x = time.adjust.before , 
                               y = read.res.rec$mu.results.df$mu.mean,
                               # y = read.res.rec$mu.res,
                               xend = time.adjust.plus, 
                               yend = read.res.rec$mu.results.df$mu.mean
                               #yend = read.res.rec$mu.res
                               ), col="red") +
      #geom_vline(xintercept = read.res.rec$ts.res, col="black", lwd =0.5, linetype= "dashed")+
      geom_vline(aes(xintercept = read.res.rec$data.annotate.recess.adjust$t.adj), 
                 linetype= "solid", lwd = 1, color="red") +
      #annotate("rect",xmin= ts.res.before, xmax=ts.res.plus, ymin=(mu.res),
      #         ymax=(mu.res+10), fill="red", alpha=1) +
      # annotate(geom="text", x=ts.res, y=limits.Y[1] + (limits.Y[2]-limits.Y[1])/100, 
      #          label=signif(ts.res, digits = 4),color="black", size=1.5)+
      theme(
        axis.text=element_text(size=15)
        ,axis.title=element_text(size=20,face="plain")
        # ,panel.grid.major=element_line(size=1.2),panel.grid.minor=element_line(size=0.8)
        ,legend.text=element_text(size=10)
        ,legend.title=element_text(size=10)
        ,legend.key.size=unit(1.5, "cm")
        ,legend.position="none"
        ,plot.margin = margin(1, 2, 1, 1, "cm")
        #,axis.text=element_blank()
        #,axis.ticks=element_blank()
        #,panel.grid =element_blank()
      )
  } else {
    #===========================================================================================
    seg.rec.plot <- ggplot()+
      geom_point(data = read.reg.rec$df.h.infinity, aes(x = t, y = mean), size = 2, pch=1) +
      geom_line(data = read.reg.rec$df.h.infinity, aes(x = t, y = mean), size = 0.2, color = "blue")+
      theme_light(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5))
    
    if (plot.gamma.uncertainty ==TRUE){
      p = p+ annotate("rect",  
                      xmin = read.reg.rec$df.h.infinity$t[1], 
                      xmax = tail(read.reg.rec$df.h.infinity$t,1),
                      ymin = data.tmp[[1]][1],
                      ymax = data.tmp[[1]][2],
                      fill = "pink", alpha=0.3)
    }
      annotate("rect", 
               xmin= read.reg.rec$df.h.infinity$t[1], 
               xmax= tail(read.reg.rec$df.h.infinity$t,1),
               ymin= read.res.rec$Q2.mu,
               ymax= read.res.rec$Q97.mu, fill="red", alpha=0.1) +
        
      geom_segment(mapping=aes(x = read.reg.rec$df.h.infinity$t[1] , 
                                 y = read.res.rec$mu.res, 
                                 xend = tail(read.reg.rec$df.h.infinity$t,1), 
                                 yend = read.res.rec$mu.res[1]), col="red") +
      scale_x_continuous(name = x.name, expand = c(0,0), limits =limits.X,
                         breaks = seq(limits.X[1], limits.X[2]+1, 2)) +
      scale_y_continuous(name = y.name, expand = c(0,0), limits =limits.Y) +
      theme(axis.text=element_text(size=15)
            ,axis.title=element_text(size=20,face="plain")
            # ,panel.grid.major=element_line(size=1.2),panel.grid.minor=element_line(size=0.8)
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin = margin(1, 2, 1, 1, "cm")
      )
    #axis.line = element_line(colour = "black"))
    if (obs.uncertainty.y == TRUE) {
      seg.rec.plot = seg.rec.plot + 
        geom_errorbar(data = read.reg.rec$df.h.infinity, aes(x=t, 
                                                ymin= (mean - 2*stdev), 
                                                ymax = (mean+2*stdev)),
                      size = 0.5, width=0.3, col= "black")
    }
  }
  
  
  if (plot.b.from.gaugings == TRUE) {
        seg.rec.plot = seg.rec.plot + 
        geom_segment(mapping= aes(x =ts.before.gaug , 
                                  y = bt2.df$mean, 
                                  xend = ts.plus.gaug, 
                                  yend = bt2.df$mean), 
                                  color = "red", size = 1,
                                  linetype = "dashed") +
        geom_rect(mapping = aes(xmin = ts.before.gaug, 
                                xmax = ts.plus.gaug, 
                                ymin = bt2.df$X2.5.,
                                ymax = bt2.df$X97.5.), 
                                fill="red", alpha=0.3)+
        geom_vline(aes(xintercept = ts.gaug$t.adj), color = "black",
                   lwd =0.7, linetype = "dotdash")
  }
  
  ggsave(seg.rec.plot, filename =paste0(dir.segm.recessions,"/segmentation_recession.png"),
         bg = "transparent", width = 8, height =5, dpi = 300)
  
  #****************************************************************************************
  # deltab.rec = 0
  # for (i in 1:(nS.rec-1)) {
  # deltab.rec[i] = (mu.res.rec[i]-mu.res.rec[i+1])/100
  # }
  
  #plots:
  #asym.df <- theta7.df
  # ppp <- ggplot()+
  #   geom_line(data = asym.df ,aes(x = t, y =maxpost), size = 0.4, color = "darkgray") +
  #   geom_errorbar(data = asym.df ,aes(x = t, ymin= q2, ymax = q97),
  #                 size = 0.2, width=60, color = "black") +
  #   geom_point(data = asym.summer.df ,aes(x = t, y = h), size = 1.5, color = "black") +
  #   geom_point(data = asym.winter.df ,aes(x = t, y = h), size = 1.5,shape = 21, colour = "black", fill = "white") +
  #   theme_light(base_size = 15)+
  #   theme(plot.title = element_text(hjust = 0.5))+
  #   scale_x_continuous(name = "time [days]", expand = c(0,0), limits =c(0,5000)) +   
  #   scale_y_continuous(name = TeX("$ Asymptotic \\; stage , \\; h_{\\infty} \\; \\left[ cm \\right] $"),
  #                      expand = c(0,0),limits =c(-100,50))+
  #   coord_cartesian(clip = 'off') +
  #   annotate("rect",xmin= Q10.ts.rec, xmax=Q90.ts.rec,
  #            ymin=-Inf, ymax=Inf, fill="black", alpha=0.3) +
  #   annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=Q10.mu.res.rec,
  #            ymax=Q90.mu.res.rec, fill="gray", alpha=0.5) +
  #   geom_vline(xintercept = ts.res.rec, col="black", lwd =0.5, linetype = "longdash")+
  #   annotate("rect",xmin= ts.res.before.rec, xmax=ts.res.plus.rec, ymin=(mu.res.rec),
  #            ymax=(mu.res.rec +0.5), fill="gray45", alpha=1) +
  #   annotate("rect",xmin= offic.times.before , xmax=offic.times.plus, ymin= b2.from.SPD,
  #            ymax=b2.from.SPD +1, fill="black", alpha=1) +
  #   theme(plot.background = element_rect(fill ="transparent", color = NA),
  #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.background = element_rect(fill ="transparent"),
  #         axis.line = element_line(colour = "black"),
  #         axis.ticks = element_line(colour = "black"))+
  #   theme(plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"))+
  #   geom_point(data=offic.times, aes(x = t, y=h),shape=4, size=3, color = "black", stroke = 1.5)
  # ggsave(ppp, filename =paste(dir.seg.rec,"/asymptote_segmentation2.png", sep=""), bg = "transparent", width = 7, height =5, dpi = 800)
  #
  
  # plot.seg.asympt.time(dir.seg.rec=dir.segm.recessions, ts.res.rec, asym.df, asym.summer.df, asym.winter.df,
  #                      Q10.ts.rec, Q90.ts.rec, ts.res.before.rec, ts.res.plus.rec, Q10.mu.res.rec,
  #                      Q90.mu.res.rec, mu.res.rec, b2.from.SPD, q10, q90, t.real.good.h, asym.h.maxpost)
  
  
  
  #***************************************
  #Identification of the real shift times:
  #***************************************
  #looking inside the u95% interval of the shift times.
  #Options:  1) if there is a flood, assign the shift time to the max peak !
  #          2) if the shift is due to other causes (vegetation, works, ice ...)
  #             then assign the shift time to the maxpost
  #          3) if the shift time is known then fix it.
  # interval.rec = NULL; ts.real.rec = NULL; flood.rec =NULL;
  # ts.morpho.real.rec = NULL; ts.morpho.MAP.rec = NULL; ts.morpho.q2.rec = NULL; ts.morpho.q97.rec = NULL;
  # times.of.shifts.rec =NULL
  # if ( nS.rec > 1) {
  #   for (i in 1:length(ts.res.rec)) {
  #     #look for the shift i:
  #     interval.rec[[i]] = which((t_limni >= min(Q2.ts.rec[i],ts.res.rec[i])) &
  #                             (t_limni <= max(ts.res.rec[i],Q97.ts.rec[i])))
  #     flood.rec[i] = max(h_limni[interval.rec[[i]]])
  #     if (flood.rec[i] < flood.threshold) {  #  ===> To improve with a sediment transport analysis !!!!!!!!!!!!!!!!
  #       ts.real.rec[i] = ts.res.rec[i] 
  #      
  #       #(...)
  #      
  #     } else {
  #       ts.real.rec[i] = t_limni[which(h_limni==flood.rec[i],1)]
  #       ts.morpho.real.rec= c(ts.morpho.real.rec, ts.real.rec[i])
  #       ts.morpho.MAP.rec = c(ts.morpho.MAP.rec, ts.res.rec[i])
  #       ts.morpho.q2.rec = c(ts.morpho.q2.rec, Q2.ts.rec[i])
  #       ts.morpho.q97.rec = c(ts.morpho.q97.rec, Q97.ts.rec[i])
  #     }
  #   }
  #   plot(t_limni, h_limni, type = "l")
  #   abline(v = ts.real.rec, col="red")
  #   abline(v = ts.res.rec, col="green")
  # }
  # times.of.shifts.rec <- c(times.of.shifts.rec, ts.real.rec)
  # times.of.shifts.rec <- sort(times.of.shifts.rec)
  # times.of.shifts.rec <- times.of.shifts.rec[c(TRUE, !times.of.shifts.rec[-length(times.of.shifts.rec)]
  #                                      == times.of.shifts.rec[-1])]
  # df.seg.rec = data.frame(t.MAP = ts.res.rec,
  #                         t.MAP.before = ts.res.before.rec,
  #                         t.MAP.plus = ts.res.plus.rec,
  #                         t.real =ts.real.rec,
  #                         ts.MAP = ts.morpho.MAP.rec,
  #                         ts.q2 = ts.morpho.q2.rec,
  #                         ts.q97 = ts.morpho.q97.rec)
  # df.mean = data.frame(MAP= mu.res.rec,
  #                      q10 = Q10.mu.res.rec,
  #                      q90 =Q90.mu.res.rec
  #                      )
  
  
  # plot(df.limni$t_limni,df.limni$h_limni, type ="l", col="gray")
  # abline(v =df.seg.rec$t.MAP, col = "red" , lwd = 3)
  # abline(v= df.seg.rec$t.real, col="green", lwd=3)
  
  #return(list(df.seg.rec, df.mean))
  
  #***************************************************************
  t.plot <- ggplot()+
            scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.X) +
            scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), 
                               breaks=seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3])) +
            ylab(limni.labels[2]) +
            coord_cartesian(clip = 'off')
  if ( read.res.rec$nS.ok >1) {
  if (plot.gamma.uncertainty == TRUE){
    for (jj in 1:read.res.rec$nS.ok){
      t.plot = t.plot + 
        annotate("rect", 
                 xmin = time.adjust.before[jj], 
                 xmax = time.adjust.plus[jj], 
                 ymin = data.tmp[[jj]][1]/100,
                 ymax = data.tmp[[jj]][2]/100,
                 fill = "pink", alpha = 0.3)
    }
  }
  #
  t.plot = t.plot + 
    annotate("rect", 
             xmin = time.adjust.before, 
             xmax = time.adjust.plus, 
             ymin = read.res.rec$Q2.mu/100,
             ymax = read.res.rec$Q97.mu/100, 
             fill ="red", alpha=0.2)
  } else {
    if (plot.gamma.uncertainty ==TRUE){
      t.plot = t.plot+ annotate("rect",  
                      xmin = read.reg.rec$df.h.infinity$t[1], 
                      xmax = tail(read.reg.rec$df.h.infinity$t,1),
                      ymin = data.tmp[[1]][1]/100,
                      ymax = data.tmp[[1]][2]/100,
                      fill = "pink", alpha=0.3)
    }
    t.plot = t.plot + 
    annotate("rect", 
             xmin= read.reg.rec$df.h.infinity$t[1], 
             xmax= tail(read.reg.rec$df.h.infinity$t,1),
             ymin= read.res.rec$Q2.mu/100,
             ymax= read.res.rec$Q97.mu/100, fill="red", alpha=0.1)
  }
  t.plot = t.plot +   
           geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray40", size = 0.3) +
           geom_point(data=read.res.rec$gaugings.df, 
                      aes(x = t , y= h), size = 3, pch =21, fill= read.res.rec$gaugings.df$color)+
           #geom_vline(aes(xintercept = data.annotate.recess.adjust$t.adj), color = "blue", lwd =0.3, linetype = "dashed")+
           geom_segment(mapping=aes(x = time.adjust.before , 
                                    y = read.res.rec$mu.res/100, 
                                    xend = time.adjust.plus, 
                                    yend = read.res.rec$mu.res/100), 
                                    col="red")+
          theme_bw(base_size=20) +
          theme( axis.text=element_text(size=15)
                ,axis.title=element_text(size=20,face="plain")
                ,panel.grid.major=element_blank()
                ,panel.grid.minor=element_blank()
                ,legend.text=element_text(size=20)
                ,legend.title=element_text(size=30)
                ,legend.key.size=unit(1.5, "cm")
                ,legend.position="none"
                ,plot.margin=unit(c(1, 0.5, 0, 2),"cm")
                ,axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0))
                ,axis.text.x=element_blank()
                ,axis.line.x = element_blank())
  # #***************************************************************
  if (!is.null(read.res.rec$dates)){
    t.plot <- t.plot +
      annotate("text", 
               x = read.res.rec$data.annotate.recess.adjust$t.adj,
               y = grid_limni.ylim[1] + 0.1, 
               label = read.res.rec$dates, color = "red", size=5)
  }
  # #*******************************************************************
  t.plot2 <- ggplot() 
  if (ncol(read.res.rec$pdf.ts.rec)==1) {
    t.plot2 = t.plot2 + 
      geom_density(aes(x= X, ..scaled.. /2),
                   fill="blue", 
                   colour=NA, alpha=0.3)
  } else {
    
    t.plot2 = t.plot2 +
      geom_density(aes(x= X$time, ..scaled.. /2, group =X$ord),  
                   fill= "blue", 
                   colour=NA, alpha=0.2)
  }
  t.plot2 = t.plot2 + 
    # annotate("rect", xmin= data.annotate.gaug$q2, xmax=data.annotate.gaug$q97, ymin=0, 
    #         ymax=0.5, fill="blue", alpha=0.1) +
    #geom_segment(data = data.annotate.gaug, aes(x = MAP, y = 0, yend =0.5, xend= MAP), size = 0.7, color ="blue")+
    geom_segment(data = read.res.rec$data.annotate.recess.adjust, 
                 aes(x = t.adj, y = -0.5, yend =0, xend=t.adj ),
                 size = 1, color ="red")+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,0.5))+
    xlab(limni.labels[1])+ 
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c( -0.5,0, 0.5), color="darkgray", linetype="dashed", size = 0.5)+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=15)
          ,axis.title=element_text(size=20,face="plain")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none"
          ,plot.margin=unit(c(0,0.5,0.2, 2),"cm")
          ,axis.title.y = element_text(margin = margin(t = 0, r = 40, b = 0, l =0 ))
          ,axis.text.y=element_blank()
          ,axis.ticks.y = element_blank()) +
    scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.X)
  
  if (is.null(read.res.rec$data.annotate.off)==FALSE) {
      t.plot2 <- t.plot2 +
        geom_point(data = read.res.rec$data.annotate.off, 
                   aes(x = xeffect, y = -1), color= "black", size = 4, shape =4, stroke=2 )
      #geom_point(data = data.annotate.off, aes(x = xpotent, y = -1), color= "green", size = 4, shape =4, stroke=2 )
      
  }
  seg.rec.plot = seg.rec.plot+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=15)
          ,axis.title=element_text(size=20,face="plain")
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin=unit(c(2, 0.5, 2, 2),"cm")
          ,axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l =0 ))
          ,axis.title.x = element_blank())
  #**********************************************************************
  t.plot3 = plot_grid(seg.rec.plot, t.plot, t.plot2, 
                      ncol = 1, nrow = 3, rel_heights = c(1, 1, 0.5)) #needs cowplot "package"
  ggsave(t.plot3, 
         filename =paste0(dir.segm.recessions,"/time_series_segm_recession.png"), 
         device = "png", width = 16, height =13, dpi = 400, units = "in")
  
  

  
  #add the other plots (for the figure of article Darienzo et al., 2020):
  #######################################################################
  colfunc = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))
  dir.recess = paste0(dir.case_study,"/Results/segmentation_recessions")
  
  if (BayesianOption ==1) {
    Ncurves= length(curves[[1]])
    # plot of recessions:
    reg.plot <- ggplot() + 
                      theme_bw(base_size = 20)+
                      scale_x_continuous(name = "Time [day]", 
                                         limits =c(limits.x.recess[1], limits.x.recess[2]), 
                                         expand = c(0,0), 
                                         breaks = seq(limits.x.recess[1], limits.x.recess[2], limits.x.recess[3])) +
                      scale_y_continuous(name = "Stage h [cm]", 
                                         limits = c(limits.y.recess[1], limits.y.recess[2]), 
                                         expand = c(0,0)) +
                      theme( axis.text=element_text(size=15)
                            ,axis.title=element_text(size=20,face="plain")
                            ,panel.grid.major=element_blank()
                            ,panel.grid.minor=element_blank()
                            ,legend.text=element_text(size=20)
                            ,legend.title=element_text(size=30)
                            ,legend.key.size=unit(1.5, "cm")
                            ,legend.position="none"
                            ,plot.margin=unit(c(1.5,0.5, 2, 2),"cm")
                            ,axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l =0 )))

    temp=NULL;  temp2 =NULL
    for (r in which.recession) {
      temp[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C", r,"/BaM/temp.txt"),header=TRUE)
      temp2[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C", r,"/BaM/temp2.txt"),header=TRUE)
      reg.plot =  reg.plot +
        # geom_point(data = temp , aes(x = x, y= y), color=colfunc(80)[curve_good.h] , size= 0.5)+
        geom_line(data = temp2[[r]], aes(x=xx, y=yy), color=colfunc(Ncurves)[r], size = 0.1)+
        # geom_ribbon(data = temp2, aes(x = xx, ymin = zz ,ymax = kk), 
        # fill = colfunc(80)[curve_good.h], alpha = 0.3) 
        geom_point(data = temp[[r]] , aes(x = x, y= y), color= colfunc(Ncurves)[r], size= 1)+
        geom_ribbon(data = temp2[[r]], aes(x = xx, ymin = zz , ymax = kk), 
                    fill = colfunc(Ncurves)[r], alpha = 0.1)
    }
    #********************************************************************************************************
    Ncurves.long = length(which.recession.long)
    # plot of recessions:
    reg.long.plot <- ggplot() + 
                     theme_bw(base_size = 20)+
                     scale_x_continuous(name = "Time [day]", 
                         limits =c (limits.x.recess[1], limits.x.recess[2]), 
                         expand = c(0,0), 
                         breaks = seq(limits.x.recess[1], limits.x.recess[2], limits.x.recess[3])) +
                     scale_y_continuous(name = "Stage h [cm]", 
                         limits = c(limits.y.recess[1], limits.y.recess[2]), 
                         expand = c(0,0)) +
                     theme(axis.text=element_text(size=15)
                           ,axis.title=element_text(size=20,face="plain")
                           ,panel.grid.major=element_blank()
                           ,panel.grid.minor=element_blank()
                           ,legend.text=element_text(size=20)
                           ,legend.title=element_text(size=30)
                           ,legend.key.size=unit(1.5, "cm")
                           ,legend.position="none"
                           ,plot.margin=unit(c(1.5,0.5,2, 2),"cm")
                           ,axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l =0 )))
    temp=NULL;  temp2 =NULL
    for (rlong in which.recession.long) {
        temp[[rlong]]  = read.table(file=paste0(dir.recess,"/curves_regression/C", rlong,"_long/BaM/temp.txt"),header=TRUE)
        temp2[[rlong]] = read.table(file=paste0(dir.recess,"/curves_regression/C", rlong,"_long/BaM/temp2.txt"),header=TRUE)
        reg.long.plot  = reg.long.plot +
                         # geom_point(data = temp , aes(x = x, y= y), color=colfunc(80)[curve_good.h] , size= 0.5)+
                         geom_line(data = temp2[[rlong]], aes(x=xx, y=yy), color=colfunc(Ncurves.long)[rlong], size = 0.1)+
                         geom_point(data = temp[[rlong]] , aes(x = x, y= y), color= colfunc(Ncurves.long)[rlong], size= 2)+
                         geom_ribbon(data = temp2[[rlong]], aes(x = xx, ymin = zz , ymax = kk), 
                                     fill = colfunc(Ncurves.long)[rlong], alpha = 0.1)
    }
    
    
    ggsave(reg.long.plot, filename =paste0(dir.segm.recessions,"/Long_recessions_estimation.png"),
           bg = "transparent", width = 12, height =7, dpi = 200)
  #########################
  } else {
  #########################
    # Quantiles for tot uncertainty:
    dir.rec.pool.test = paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression/Pooling/test3")
    data.MCMC.cooked=as.matrix(read.table(paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt")  #### Cooked MCMC
                                          ,header=TRUE,dec=".", sep=""))
    data.MCMC.MaxPost = as.numeric(read.table(paste(dir.rec.pool.test,"/Results_Summary.txt",sep="")
                                              ,row.names=1,dec=".",sep="", skip = 16))
    nsample = length(data.MCMC.cooked[,1])
    tgrid=seq(limits.x.recess[1],  limits.x.recess[2],  0.5) 
    Ncurves.pool = length(which.recession)
    #
    # recess model with tot uncertainty:
    Recession.model = function(theta, t){ # Rating curve with discharge error propagation:
      h=0*t
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
      res.h = sapply(h, function(h, theta){h + 
          rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
          theta=c(theta[8],theta[9]))
      return(res.h)
    }
    # recess model for maxpost:
    Recess.Maxpost.model = function(theta,t){  # Exponential recession (maximum posterior):
      h=0*t
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
      return(h)
    }
    #Initialisation:
    b1.var =FALSE
    MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period 
    for(i in 1:Ncurves.pool){
      if (b1.var == TRUE){
        col.num = c( i, Ncurves.pool+i,                      #a1, b1,
                     Ncurves.pool*2 +1, Ncurves.pool*2 +2,   #a2, b2
                     Ncurves.pool*2 +3, Ncurves.pool*2 +4,   #a3, b3,
                     Ncurves.pool*2 +4 + i,                  #a4 (asymptotic level parameter)
                     Ncurves.pool*3 + 5, Ncurves.pool*3 +6)  #gamma1, gamma2
      } else {
        col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                     Ncurves.pool +2, Ncurves.pool +3,       #a2, b2
                     Ncurves.pool +4, Ncurves.pool +5,       #a3, b3,
                     Ncurves.pool +5 + i,                    #a4 (asymptotic level parameter)
                     Ncurves.pool*2 + 5 +1, Ncurves.pool*2 +5+2)  #gamma1, gamma2
      }
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    }
    Rec.Post <-   apply(MCMC.save,    MARGIN=1, Recession.model,       t=tgrid) #apply Rec model for tgrid with structural error:
    Rec.MaxPost = apply(MaxPost.save, MARGIN=1, Recess.Maxpost.model,  t=tgrid) # Maximum posterior computing
    ############
    # Quantiles for Figure RC: 
    message("- Plotting all Recession curves !!!  Wait ... "); flush.console()
    List.Rec.quants = list(NULL)
    for(i in 1:Ncurves.pool){
      data.tmp = apply(Rec.Post[,(nsample*(i-1)+1):(nsample*i)], 
                       MARGIN=1, quantile, probs=c(0.025,0.975) ) #, na.rm=TRUE)
      List.Rec.quants[[i]] = data.frame(cbind(tgrid, t(data.tmp), Rec.MaxPost[,i]-100))
      colnames(List.Rec.quants[[i]]) = c("t", "inf", "sup", "maxpost")
    }
    inter.per=seq(1,Ncurves.pool,1)   # inter.per=c(25:35)
    nRec = length(inter.per)
    palette.per = colfunc(Ncurves.pool)
    palette.per = palette.per[1:Ncurves.pool]
    data.plot.Rec = data.frame(do.call("rbind", 
                                       List.Rec.quants[inter.per]), 
                               Period=rep(inter.per,each=length(tgrid)))
    inter.null=which(data.plot.Rec$maxpost==0)
    # Initialise the plot:
    #err.per=0.1; color.inter="#F4A582"; color.interBorder="#F4A582"; alpha.inter=0.5; color.gaugings="#0571B0"
    # na.rm=TRUE; shape.gaugings=18; size.gaugings=4; color.RC="#CA0020"; size.RC=2; size.wind=20; size.axis=20
    # face.wind="bold"; hydr.op=NULL; color.measurement="#92C5DE"; na.rm.meas=TRUE; shape.meas=18; size.meas=4
    # error.bar.op=NULL; width.errorbar=0.005; time.unit="s"; inter.hydr=NULL;blank=FALSE; width.wind=20
    # height.wind=20; param.plot=FALSE; color.param.inter="#FEE0B6"; color.param.interBorder="#FEE0B6"
    # alpha.param.inter=0.5; hydr.op=1;time.unit="h"; blank=TRUE; alpha.inter=0.5; width.wind=40; height.wind=20
    # param.plot = TRUE; color.param.inter = "#B2ABD2"; color.param.interBorder = "#B2ABD2";
    # vect.break=seq(1,nperiod,1) #seq(1,9,1)
    # vect.emp=rep("",22)
    # breaks.wind=c(vect.break*0.01,vect.break*0.1,vect.break,vect.break*10,vect.break*100,vect.break*1000)
    # labels.wind=c(0.01,vect.emp,0.1,vect.emp,1,vect.emp,10,vect.emp,100,vect.emp,1000,vect.emp)
    data.rec.obs = read.table(paste0(dir.rec.pool.test,"/Curves_Data.txt"),header=TRUE,dec=".",sep="") # Gaugings loading
    ylim.wind = c(stage.limits[1], stage.limits[2])
    xlim.wind =c(limits.x.recess[1], limits.x.recess[2])
    pos.num=function(x.int){inter=which(x.int==data.rec.obs$period);return(inter)}
    inter.rec.obs =unlist(sapply(inter.per, pos.num), recursive = TRUE, use.names = TRUE)
    # data for the plot:
    data.Rec = data.plot.Rec #[-inter.null,]
    data.Rec$Period=factor(data.Rec$Period)
    data.rec.obs$Period = as.factor(data.rec.obs$Period)
    write.table( data.Rec, paste0(dir.rec.pool.test,"/Rec_SPD_env.txt"), sep ="\t", row.names=FALSE)
    
    
    # RC ggplot::
    reg.pool.plot = ggplot(data.Rec)+
      geom_smooth(aes(x=t, y=maxpost, ymax=sup-100, ymin=inf-100, group=Period, fill=Period), 
                  size=1, stat='identity', alpha=0.1)  + #alpha=0.1
      geom_path(aes(x=t,   y=maxpost, group=Period, colour=Period), size=1) +
      ### recession curve obs:
      geom_linerange(aes(x=time, ymax= (h + 2*uh)-100, ymin=(h-2*uh)-100, colour=Period), data=data.rec.obs, size=0.3)+
      geom_point(aes(x=time, y=h-100, colour=Period), data=data.rec.obs, shape=16, size=1.5)+
      ### Labels
      xlab(expression(paste("Recession time [day]",sep="")))+
      ylab(expression(paste("Stage h [cm]")))+
      labs(colour = "Period")+
      scale_colour_manual(values = palette.per)+
      scale_fill_manual(values = palette.per)+
      #scale_y_continuous(breaks=breaks.wind,labels=labels.wind)+
      coord_cartesian(ylim=ylim.wind, xlim=xlim.wind)+
      scale_x_continuous(limits = xlim.wind , expand=c(0,0))+
      scale_y_continuous(expand=c(0,0)) +
      ### Theme
      theme_bw(base_size=15)+
      theme( axis.text=element_text(size=15)
             ,axis.title=element_text(size=20,face="bold")
             ,panel.grid.major=element_blank() #element_line(size=1.2)
             ,panel.grid.minor=element_blank()
             ,legend.text=element_text(size=20)
             ,legend.title=element_text(size=30)
             #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
             ,legend.key.size=unit(1.5, "cm")
             ,legend.position="none"
             ,plot.margin=unit(c(1.5,0.5,2, 2),"cm")
             ,axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l =0 )))
    ggsave(reg.pool.plot, filename =paste(dir.segm.recessions,"/regression_pool.png", sep=""),
           bg = "transparent", width = 12, height =7, dpi = 200)     
    
  }
  
  ############################################################################################
  starting.time  = starting.time #days
  tot.length.rec = length(curves[[1]][[which.recession.realtime]][,1])
  
  Ncurves.realtime =  tot.length.rec - starting.time
  # plot of recessions:
  reg.realtime.plot <- ggplot() + 
              theme_bw(base_size = 20)+
              scale_x_continuous(name = "Time [day]", 
                       limits =c(limits.x.recess[1], limits.x.recess[2]), 
                       expand = c(0,0), 
                       breaks = seq(limits.x.recess[1], limits.x.recess[2], limits.x.recess[3])) +
              scale_y_continuous(name = "Stage h [cm]", 
                       limits = c(limits.y.recess[1], limits.y.recess[2]), 
                       expand = c(0,0)) +
            theme( axis.text=element_text(size=15)
           ,axis.title=element_text(size=20,face="plain")
           ,panel.grid.major=element_blank()
           ,panel.grid.minor=element_blank()
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="none"
           ,plot.margin=unit(c(1.5,0.5, 2, 2),"cm")
           ,axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l =0 )))
           
           dir.regression.rt <- paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression/Real_Time")
           temp.rt =NULL;  temp2.rt =NULL
           for (r in 1:Ncurves.realtime) {
               dir.rec.rt= paste0(dir.regression.rt,"/C",r,"/BaM")
               temp.rt[[r]] = read.table(file=paste0(dir.rec.rt,"/temp.txt"),header=TRUE)
               temp2.rt[[r]] = read.table(file=paste0(dir.rec.rt,"/temp2.txt"),header=TRUE)
               #
               reg.realtime.plot =  reg.realtime.plot +
                                    geom_point(data = temp.rt[[r]] , aes(x = x, y = y), color= "black", size= 2) +
                                    geom_ribbon(data = temp2.rt[[r]], aes(x = xx,   ymin = zz ,    ymax = kk), 
                                                fill = "blue", alpha=0.1) 
           }
  reg.realtime.plot = reg.realtime.plot+
     theme( axis.text=element_text(size=15)
           ,axis.title=element_text(size=20,face="plain")
           ,panel.grid.major=element_blank()
           ,panel.grid.minor=element_blank()
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="none"
           ,plot.margin=unit(c(1.5,0.5, 2, 2),"cm")
           ,axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l =0 )))

  
  # plot4paper1 =  plot_grid( 
  #                          reg.long.plot,
  #                          reg.plot,
  #                          labels = c('a)', 'b)'), 
  #                          label_size = 30,
  #                          label_fontface = "plain",
  #                          ncol = 2, nrow=1)  
  plot4paper1 =  plot_grid( 
                           reg.pool.plot,
                           reg.realtime.plot,
                           labels = c('a)', 'b)'), 
                           label_size = 30,
                           label_fontface = "plain",
                           ncol = 2, nrow=1)  
  plot4paper = plot_grid(plot4paper1, 
                         seg.rec.plot,
                         t.plot, 
                         t.plot2, 
                         labels = c('',  'c)', 'd)', '' ), 
                         label_size = 30,
                         label_fontface = "plain",
                         ncol = 1, nrow =4,
                         rel_heights = c(1, 1, 1, 0.5)) 
  ####################################################
  ggsave(filename=paste0(dir.segm.recessions,"/Figure4paper_recessions.png"), 
         plot=plot4paper,
         width = 16, height =20, dpi = 400) 
  ####################################################
}













































#######################################################################################
# Stage h time series  with  shift times:
initial.ts.plot.rec <- function(CdT.P ,
                                df.limni , 
                                tshift ,
                                limni.labels,
                                grid_limni.ylim,
                                dir.seg.gaug,
                                seg.iter, 
                                t_Gaug, h_Gaug, 
                                mcmc.segment, 
                                nS) {
#######################################################################################
  # reading data:
  Utot.times = "90% total uncertainty of change point times"
  MAPtimes ="Change point times (MAP)"
  col.tflood ="Time of flood"
  color.segm = c("90% total uncertainty of change point times"="blue")
  col.segm = c("Change point times (MAP)"="blue", "Time of flood"="red")

  ts.res = tshift$tau.MAP;
  tflood = tshift$tflood;
  Q2.ts = tshift$tau.q2;
  Q97.ts = tshift$tau.q97;
  
  col.ts.distrib <- rainbow((nS-1)) 
  
  X1=mcmc.segment
  X2=X1[,(nS+1):(2*nS-1)]
  if (nS ==2) {
    X = X2
  } else {
    X =data.frame(time= X2[,1], ord = rep(1, length(X2[,1])) , colorrr = rep(col.ts.distrib[1], length(X2[,1])))
    for (orderr in 2:(nS-1)) {
      X =rbind(X,    data.frame(time= X2[,orderr], ord = rep(orderr, length(X2[,1])),  colorrr = rep(col.ts.distrib[orderr], length(X2[,1]))) )
    }
  }
  
  
  # plot 1:
  ###############################
  initial.ts.plot <- ggplot() 
  if (is.null(df.limni)==FALSE) {
    initial.ts.plot= initial.ts.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "darkgray",size = 0.2)+
      scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0, tail(df.limni$t_limni,1))) +
      scale_y_continuous(name=limni.labels[2], limits = grid_limni.ylim[1:2], expand = c(0,0))+
      coord_cartesian(clip = 'off')
  } else {
    initial.ts.plot= initial.ts.plot + 
      scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0, tail(t_Gaug,1))) +
      scale_y_continuous(name=limni.labels[2], limits = c(min(h_Gaug), max(h_Gaug)), expand = c(0,0))+
      coord_cartesian(clip = 'off')
  }
  initial.ts.plot= initial.ts.plot +
    geom_point(aes(x = t_Gaug, y =h_Gaug), fill = "gray60", size = 1.5, pch =21) +
    geom_point(data = CdT.P, aes(x = tP, y = hP), fill = "black", size = 1.5, pch=21) +
    theme_classic(base_size = 15) +
    ylab(limni.labels[2]) +
    geom_rect(mapping= aes(xmin= Q2.ts, xmax=Q97.ts ,ymin=-Inf, ymax=Inf), fill=col.ts.distrib, alpha=0.2 ) +
    geom_vline(aes(xintercept = ts.res), col=col.ts.distrib, lwd =0.5, linetype = "solid") +
    scale_fill_manual(name=element_blank(), values=color.segm) +
    scale_colour_manual(name = element_blank(), values=col.segm, breaks=c(MAPtimes),labels = c(MAPtimes),
                        guide = guide_legend(override.aes = list(linetype = c("solid","longdash"),shape = c(NA, NA)))) +
    theme(text = element_text(size=10),
          plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"),
          plot.margin=unit(c(0.3,0.3,0,0.5),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0.3)),
          axis.text.x=element_blank(),
          axis.line = element_line(colour = "black", 
                                   size = 0.4, linetype = "solid"),
          axis.ticks.x = element_line(colour = "black", 
                                      size = 0.4, linetype = "solid"),
          axis.ticks.y = element_line(colour = "black", 
                                      size = 0.4, linetype = "solid"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position ="none",
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.background = element_rect(colour = "transparent", fill = "transparent"))
  if (is.null(tflood)==FALSE) {
    initial.ts.plot= initial.ts.plot +
      geom_vline(aes(xintercept = tflood, col=col.tflood), linetype="longdash", lwd =0.4)
  }
  
  # plot 2:
  #####################################
  init.ts.dens = ggplot()
  if (is.null(df.limni)==FALSE) {
    init.ts.dens= init.ts.dens + 
      scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0, tail(df.limni$t_limni,1)))
  } else {
    init.ts.dens= init.ts.dens + 
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits = c(0, tail(t_Gaug,1)))
  }
  init.ts.dens= init.ts.dens +  
    ylab("Scaled pdf")+
    theme_classic(base_size = 15)+ theme(text = element_text(size=10),
                                         plot.title = element_text(hjust = 0.5),
                                         #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
                                         plot.background = element_rect(fill = "transparent", color = NA),
                                         panel.background = element_rect(fill = "transparent"),
                                         plot.margin=unit(c(0.3,0.3,0,0.05),"cm"),
                                         axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
                                         axis.line = element_line(colour = "black", 
                                                                  size = 0.4, linetype = "solid"),
                                         axis.ticks.x = element_line(colour = "black", 
                                                                     size = 0.4, linetype = "solid"),
                                         axis.ticks.y = element_line(colour = "black", 
                                                                     size = 0.4, linetype = "solid"),
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank(),
                                         legend.position ="none",
                                         legend.key = element_rect(colour = "transparent", fill = "transparent"),
                                         legend.background = element_rect(colour = "transparent", fill = "transparent"))
  
  if (nS == 2) {
    init.ts.dens= init.ts.dens + geom_density(aes(x= X, ..scaled..),
                                              fill=col.ts.distrib[1], 
                                              colour=NA, alpha=0.3)
  } else {
    
    init.ts.dens= init.ts.dens + geom_density(aes(x= X$time, ..scaled.., group =X$ord,  
                                                  fill= X$colorrr, colour=X$colorrr), 
                                              alpha=0.3)
  }
  # combined plot:
  ######################
  initial.ts.plot2 = plot_grid( initial.ts.plot, 
                                init.ts.dens, ncol = 1, nrow = 2, rel_heights = c(1, 0.5))   #needs cowplot "package"
  return(initial.ts.plot2)
}




















































############################################################################################################
plot.segmentation.recession.all.methods = function(dir.segm.recessions,
                                                   df.limni,
                                                   limits.X, limits.Y, limits.Y.alpha, limits.Y.lambda,
                                                   x.name, y.name,
                                                   limni.labels,
                                                   grid_limni.ylim,
                                                   obs.uncertainty.y,
                                                   plot.gamma.uncertainty,
                                                   plot.b.from.gaugings,
                                                   #***********************
                                                   dir.case_study,
                                                   rec.model, 
                                                   which.recession,
                                                   BayesianOption,
                                                   limits.y.recess,   # stage is in "cm"
                                                   limits.x.recess,
                                                   stage.scale.shift,
                                                   starting.time,
                                                   chi.test,
                                                   plot.dates,
                                                   model.names,
                                                   model.titles,
                                                   data.annotate.off,
                                                   bt.from.gaugingsss,
                                                   pdf.ts.gaugings,
                                                   plot_ts_gaugings
) {
  ##############################################################################################################
  # Recession models studied :
  #***************************
  # 1) 1exp =>                 h(t) = alpha1 *exp(-lambda1* t) + beta
  # 2) 2exp =>                 h(t) = alpha1 *exp(-lambda1* t) + alpha1 *exp(-lambda1* t) + beta
  # 3) 3exp =>                 h(t) = alpha1 *exp(-lambda1* t) + alpha2 *exp(-lambda2* t) + alpha3*exp(-lambda3* t) + beta
  # 4) double exp (horton) =>  h(t) = alpha1 *exp(-lambda1* t^n1) + beta
  # 5) hyperb   =>             h(t) = alpha1 /(1+lambda1*t)^2
  # 6) Boussinesq (1903) =>    h(t) = alpha1 *(1 + lambda1*t)^(-2)
  # model.names = c("expexp", "hyperb", "expexp", "hyperb", "expexp", "hyperb")
  
  #directories:
  dir.recess = paste0(dir.case_study,"/Results/segmentation_recessions")
  dir.regression = paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression")
  dir.rec.pool = paste0(dir.regression,"/Pooling")
  dir.BaM.rec.pool = paste0(dir_code,"/BaM_exe/Recession_h_pooling")
  # preparing  the dataframe for riverbed estimated through gaugings :
  # if (plot.b.from.gaugings == TRUE){
  #   dir.sed.transp = paste0(dir.case_study,"/Results/segmentation_sed_transp")
  #   dir.sed.transp.SPD = paste0(dir.sed.transp,"/SPD")
  #   bt1.df = read.table(paste0(dir.sed.transp.SPD,"/bt1_df.txt"), header =TRUE)
  #   bt2.df = read.table(paste0(dir.sed.transp.SPD,"/bt2_df.txt"), header =TRUE)
  #   ts.gaug =  read.table(paste0(dir.sed.transp.SPD,"/STEP1_shift_times.txt"), header =TRUE)
  #   ts.before.gaug =  sort(c(0, ts$t.adj))
  #   ts.plus.gaug = sort(c(ts$t.adj, tail(df.limni$t_limni,1)))
  # }
  
  # recess model with tot uncertainty:
  #*******************************************
  Recession.model = function(theta, model, t){ 
    #*******************************************
    h=0*t
    if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[8],theta[9]))
    } else if (model =="2expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="2expWithAsympt_bis"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="2expWithAsympt_rel"){
      h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
          theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="1expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[4],theta[5]))
    } else if ((model =="expexp")|(model =="expexp_bis")){
      h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } else if ((model =="hyperb")|(model =="hyperb_bis")){
      h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } 
    return(res.h)
  }
  
  # recess model for maxpost:
  #************************************************
  Recess.Maxpost.model = function(theta, model, t){
    #************************************************
    h=0*t
    if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
    }  else if (model =="2expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
    }  else if (model =="2expWithAsympt_bis"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
    }  else if (model =="2expWithAsympt_rel"){
      h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
        theta[5] 
    } else if ((model =="1expWithAsympt")|(model =="1expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) + theta[3] 
    } else if ((model =="expexp")|(model =="expexp_bis")){
      h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
    } else if ((model =="hyperb")|(model =="hyperb_bis")){
      h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
    } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
    }
    return(h)
  }
  

  #initialisation of the lists of plot objects:
  reg.pool.plot=NULL
  reg.pool.plot2=NULL  
  seg.rec.plot=NULL
  t.plot=NULL
  t.plot2 =NULL
  t.plot3 =NULL
  t.plot4 = NULL
  title.model=NULL
  dir.rec.segm.test.param=NULL;
  tau.results.df.rec =NULL; mu.results.df.rec=NULL; gamma.results.df.rec=NULL; df.shift.times.rec=NULL
  df.shift.times.plus.rec=NULL; ts.res.before.rec=NULL; ts.res.rec=NULL; ts.res.plus.rec=NULL; Q2.ts.rec=NULL
  Q97.ts.rec=NULL; Q2.mu.rec=NULL; mu.res.rec=NULL; Q97.mu.rec=NULL; Data.segm.rec=NULL; nS.ok.rec=NULL;
  mcmc.seg.rec=NULL; pdf.ts.rec=NULL; gamma_segm_recess=NULL; gaugings.df.recess=NULL;
  shift.times.recessions=NULL; data.annotate.recess=NULL; data.annotate.recess.adjust=NULL; 
  X=NULL; X1=NULL; data.tmp= NULL;data.tmp.2= NULL; Ysim=NULL; time.adjust.before =NULL; time.adjust.plus=NULL; parameters =NULL; parameters.names =NULL 
  model.title =0
  
  
  
  
  
  
  
  
  # do a specific subfigure column for the results of each model:  
  ###########################################################################################
  for (mod in 1:length(model.names)){
    ###########################################################################################
    print(paste0("model ",mod,"/",length(model.names)," ********************"))
    #looking for the directories of results of regression estimation for the model:
    dir.create(paste0(dir.rec.pool,"/test_" ,model.names[mod]))
    dir.rec.pool.model = paste0(dir.rec.pool,"/test_", model.names[mod])   
    dir.create(paste0(dir.rec.pool.model,"/chi_", chi.test))
    dir.rec.pool.test = paste0(dir.rec.pool.model,"/chi_",chi.test)
    # directories where to save the results of var parameters for the segmentation:
    dir.create(paste0(dir.segm.recessions,"/",model.names[mod]))
    dir.segm.recessions.model = paste0(dir.segm.recessions,"/", model.names[mod])
    dir.create(paste0(dir.segm.recessions.model,"/chi_",chi.test))
    dir.segm.recessions.model.param = paste0(dir.segm.recessions.model,"/chi_",chi.test)
    # read data mcmc:  
    data.MCMC.cooked=as.matrix(read.table(paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"), header=TRUE,dec=".", sep=""))
    data.MCMC.MaxPost = as.numeric(read.table(paste0(dir.rec.pool.test,"/Results_Summary.txt"), row.names=1,dec=".",sep="", skip = 16))
    colfunc = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))
    summary.rec   = read.table(file=paste0(dir.rec.pool.test,"/Results_Summary.txt"), header=TRUE)
    mcmc.rec      = read.table(file=paste0(dir.rec.pool.test,"/Results_MCMC_cooked.txt"), header=TRUE)
    residuals.rec = read.table(file=paste0(dir.rec.pool.test,"/Results_Residuals.txt"), header=TRUE)
    curves.data.rec = read.table(file=paste0(dir.rec.pool.test,"/Curves_Data.txt"), header=TRUE)
    
    
    
    
    
    ###########################################################################################
    # Recessions regression plot:
    ###########################################################################################
    print("plot 1")
    nsample = length(data.MCMC.cooked[,1])
    tgrid=seq(limits.x.recess[1],  limits.x.recess[2],  0.5) 
    Ncurves.pool = tail(curves.data.rec$Period,1)
    
    
    
    #Initialisation:
    if (model.names[mod] =="3expWithAsympt"){
    ########################################
      model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2 e^{-\\lambda_2 t} +  \\alpha_3 e^{-\\lambda_3 t} + \\beta^{(k)}$")
      b1.var =FALSE
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                                  #a1, b1,
                     Ncurves.pool +2, Ncurves.pool + 3,                    #a2, b2
                     Ncurves.pool +4, Ncurves.pool + 5,                      #a3, b3,
                     Ncurves.pool +5 + i,                            #a4 (asymptotic level parameter)
                     Ncurves.pool*2 + 5+1, Ncurves.pool*2 + 5+2)  #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
   } else if (model.names[mod] =="3expWithAsympt_bis"){
   #####################################################
      model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2^{(k)} e^{-\\lambda_2 t} + \\alpha_3 e^{-\\lambda_3 t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                               #a1(var), b1,
                     Ncurves.pool +1 + i, 2*Ncurves.pool + 2,         #a2(var), b2
                     2*Ncurves.pool + 3, 2*Ncurves.pool + 4,          #a3, b3,
                     2*Ncurves.pool + 4 + i,                          #a4(var) (asymptotic level parameter)
                     3*Ncurves.pool + 4 + 1, 3*Ncurves.pool + 6)      #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }  
      
      
    }  else if (model.names[mod] =="2expWithAsympt_bis"){
    #####################################################
      model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2^{(k)} e^{-\\lambda_2 t}  + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=8) # 7 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=8)      # 7 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                     Ncurves.pool + 1 + i, Ncurves.pool*2 +1 + 1,   #a2, b2
                     Ncurves.pool*2 +1+1+i,                       #a3,
                     Ncurves.pool*3 + 2 + 1, Ncurves.pool*3 +2+2)  #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
      
    }  else if ((model.names[mod] =="2expWithAsympt")|(model.names[mod] =="2expWithAsympt_rel")){
    ##############################################################################################
      model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2 e^{-\\lambda_2  t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=8) # 7 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=8)      # 7 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                             #a1, b1,
                     Ncurves.pool +2, Ncurves.pool +3,              #a2, b2
                     Ncurves.pool + 3 + i,                           #a3,
                     Ncurves.pool*2 + 3 + 1, Ncurves.pool*2 +3 +2)  #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
       
    } else if (model.names[mod] =="1expWithAsympt"){
    ################################################
      model.title[mod] = TeX("$h(t,k) = \\alpha^{(k)} \\; e^{-\\lambda t} + \\beta^{(k)}$")
      MCMC.save    = matrix(NA, nrow=Ncurves.pool*nsample, ncol=6) # 5 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=6)      # 5 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                     Ncurves.pool +1+ i,                      #a2
                     Ncurves.pool*2 +1+1, Ncurves.pool*2 +1 + 2)  #gamma1, gamma2
         MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
         MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
      
    } else if ((model.names[mod] =="expexp")|(model.names[mod] =="hyperb")|(model.names[mod] =="Coutagne")) {
    ##########################################################################################################
      model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t}^{\\eta t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool + 1,                          # a1(k), b1,
                     Ncurves.pool + 2,                             # n1
                     Ncurves.pool + 2 + i,                         # a2(k) (asymptotic level parameter)
                     Ncurves.pool*2 + 2 +1, Ncurves.pool*2 +2 +2)  # gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
      
    } else if ((model.names[mod] =="expexp_bis")|(model.names[mod] =="hyperb_bis")|(model.names[mod] =="Coutagne_bis")){
    ####################################################################################################################
      model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t}^{\\eta^{(k)} t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool + i,                     # a1(k), b1(k),
                     Ncurves.pool*2 + 1,                      # n1
                     Ncurves.pool*2 + 1 + i,                  # a2 (asymptotic level parameter)
                     Ncurves.pool*3 + 2, Ncurves.pool*3 +3)   # gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      } 
    }
    
    #######################################################################################################
    # Apply recession model (total uncertainty and maxpost):
    Rec.Post    = apply(MCMC.save,    MARGIN=1, Recession.model,      model=model.names[mod],  t=tgrid) #add structural error:
    Rec.MaxPost = apply(MaxPost.save, MARGIN=1, Recess.Maxpost.model, model=model.names[mod],  t=tgrid) # Maximum posterior 
    
    # Quantiles: 
    message("- Plotting all Recession curves !!!  Wait ... "); flush.console()
    List.Rec.quants = list(NULL)
    for(i in 1:Ncurves.pool){
      data.tmp = apply(Rec.Post[,(nsample*(i-1)+1):(nsample*i)], 
                       MARGIN=1, quantile, probs=c(0.025, 0.975),  na.rm=TRUE)
      List.Rec.quants[[i]] = data.frame(cbind(tgrid, 
                                              t(data.tmp - stage.scale.shift),
                                              Rec.MaxPost[,i] - stage.scale.shift))  #!!!!!!!!!!!!!! 
      colnames(List.Rec.quants[[i]]) = c("t", "inf", "sup", "maxpost")
    }
    
    #prepare plot:
    inter.per=seq(1,Ncurves.pool,1) 
    nRec = length(inter.per)
    palette.per = colfunc(Ncurves.pool)
    palette.per = palette.per[1:Ncurves.pool]
    data.plot.Rec = data.frame(do.call("rbind", 
                                       List.Rec.quants[inter.per]), 
                               Period=rep(inter.per,each=length(tgrid)))
    inter.null=which(data.plot.Rec$maxpost==0)
    data.rec.obs = read.table(paste0(dir.rec.pool.test,"/Curves_Data.txt"),
                              header=TRUE,dec=".",sep="") # Gaugings loading
    ylim.wind = c(limits.y.recess[1], limits.y.recess[2])
    xlim.wind =c(limits.x.recess[1], limits.x.recess[2])
    pos.num =function(x.int){
             inter=which(x.int==data.rec.obs$period);
    return(inter)}
    inter.rec.obs =unlist(sapply(inter.per, pos.num), recursive = TRUE, use.names = TRUE)
    data.Rec = data.plot.Rec #[-inter.null,]
    data.Rec$Period = factor(data.Rec$Period)
    data.rec.obs$Period = as.factor(data.rec.obs$Period)
    write.table( data.Rec, paste0(dir.rec.pool.test,"/Rec_SPD_env.txt"), sep ="\t", row.names=FALSE)
    
    #recession regression plot:
    PltData.rec <- data.Rec[data.Rec$inf > -200,]
    ##############################################
    reg.pool.plot[[mod]] = ggplot(PltData.rec)+
      # geom_smooth(aes(x=t,
      #                 y=maxpost,
      #                 ymax=sup, 
      #                 ymin=inf, 
      #                 group=Period, 
      #                 fill=Period), 
      #             size=0.1, stat='identity', alpha=0.1)  + #alpha=0.1
      geom_path(aes(x=t,   
                    y=maxpost, 
                    group=Period, 
                    colour=Period), size=0.1) +
      ### recession curve obs:
      geom_linerange(aes(x=time, 
                         ymax= (h + 2*uh) - stage.scale.shift, 
                         ymin=(h-2*uh) - stage.scale.shift, 
                         colour=Period),
                     data=data.rec.obs, size=0.1)+
      geom_point(aes(x=time, 
                     y=h-stage.scale.shift,
                     colour=Period), 
                 data=data.rec.obs, shape=16, size=1)+
      ### Labels
      xlab(expression(paste("Recession time [days]",sep="")))+
      ylab(expression(paste("Stage h [cm]")))+
      labs(colour = "Period")+
      ggtitle( model.titles[mod]) +
      scale_colour_manual(values = palette.per)+
      scale_fill_manual(values = palette.per)+
      #scale_y_continuous(breaks=breaks.wind,labels=labels.wind)+
      coord_cartesian(ylim=ylim.wind, xlim=xlim.wind)+
      scale_x_continuous(limits = xlim.wind , expand=c(0,0))+
      scale_y_continuous(expand=c(0,0)) +
      ### Theme
      theme_light(base_size=15)+
      theme(axis.text=element_text(size=10)
            ,axis.title=element_text(size=15,face="plain")
            ,plot.title =element_text(size=25, face="bold", hjust = 0.5)
            ,panel.grid.major=element_blank() #element_line(size=1.2)
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(2, 0.5, 0, 1),"cm")
            ,axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
    # add the title with the name of the model:
    title.model[[mod]] <- ggdraw() + 
      draw_label( model.title[mod],
                 x = 0, hjust = 0,
                 size = 15) +
      theme(plot.margin = margin(0, 0, 0, 7))
    reg.pool.plot2[[mod]] =  plot_grid(title.model[[mod]],
                                       reg.pool.plot[[mod]] ,  
                                       ncol = 1,
                                       # rel_heights values control vertical title margins
                                       rel_heights = c(0.2, 1))
    
    pdf(paste0(dir.segm.recessions,"/Figure_recession_model_",model.names[mod],"_chi",chi.test,".pdf"), 15,17 ,useDingbats=F)
    print(reg.pool.plot[[mod]] )
    dev.off()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Segmentation for the specific model mod and all parameters:
    #####################################################################################
    print("plot 2")
    dir.rec.segm.test.param[[mod]] = list(); tau.results.df.rec[[mod]] = list(); mu.results.df.rec[[mod]] =list(); 
    gamma.results.df.rec[[mod]] =list();  df.shift.times.rec[[mod]] = list(); df.shift.times.plus.rec[[mod]] =list(); ts.res.before.rec[[mod]] =list(); 
    ts.res.rec[[mod]] =list(); ts.res.plus.rec[[mod]] =list();  Q2.ts.rec[[mod]] =list(); Q97.ts.rec[[mod]] =list();   Q2.mu.rec[[mod]] =list();
    mu.res.rec[[mod]] =list(); Q97.mu.rec[[mod]] =list(); Data.segm.rec[[mod]] =list(); nS.ok.rec[[mod]] =list(); mcmc.seg.rec[[mod]] =list();
    pdf.ts.rec[[mod]] =list(); gamma_segm_recess[[mod]] =list(); gaugings.df.recess[[mod]] =list(); shift.times.recessions[[mod]] = list(); 
    data.annotate.recess[[mod]] =list(); data.annotate.recess.adjust[[mod]] =list(); X[[mod]] =list(); X1[[mod]] =list(); data.tmp.2[[mod]] = list(); Ysim[[mod]] =list();
    time.adjust.before[[mod]]  =list(); time.adjust.plus[[mod]] =list(); 
    seg.rec.plot[[mod]] =list();
    
    if (model.names[[mod]] == "3expWithAsympt"){
      # parameters = c("a1", "a2", "a3", "a4")
      # parameters.names = c(TeX("$\\alpha_{1}$"), TeX("$\\alpha_{2}$"), TeX("$\\alpha_{3}$"), TeX("$\\beta$"))
      parameters[[mod]] = c("a1", "a4")
      parameters.names[[mod]] = c(TeX("$\\alpha_{1}$"), TeX("$\\beta$"))
      
    } else if (model.names[[mod]] =="hyperb"){  #Drogue model
      parameters[[mod]] = c("a1", "a2")
      parameters.names[[mod]] = c(TeX("$\\alpha$"),  TeX("$\\beta$"))
      
    } else if (model.names[[mod]] =="expexp"){  #Horton model
      parameters[[mod]] = c("a1", "a2")
      parameters.names[[mod]] = c(TeX("$\\alpha$"), TeX("$\\beta$"))
      
    } else if (model.names[[mod]] =="2expWithAsympt_bis"){
      parameters[[mod]] = c("a1", "a2", "a3")
      parameters.names[[mod]] = c(TeX("$\\alpha_{1}$"), TeX("$\\alpha_{2}$"), TeX("$\\beta$"))
      
    } else if (model.names[[mod]] =="2expWithAsympt"){
      parameters[[mod]] = c("a1", "a3")
      parameters.names[[mod]] = c(TeX("$\\alpha_{1}$"), TeX("$\\beta$"))
      
    } else if ((model.names[[mod]] =="1expWithAsympt")){
      parameters[[mod]] = c("a1", "a2")
      parameters.names[[mod]] = c(TeX("$\\alpha$"), TeX("$\\beta$"))
      
    } else if ((model.names[[mod]] =="2expWithAsympt_rel")){
      parameters[[mod]] = c("a1", "a3")
      parameters.names[[mod]] = c(TeX("$\\alpha_{1}$"), TeX("$\\beta$"))
      
    } else if (model.names[[mod]] =="Coutagne"){
      parameters[[mod]] = c("a1", "a2")
      parameters.names[[mod]] = c(TeX("$\\alpha$"),  TeX("$\\beta$"))
      
    } else if (model.names[[mod]] == "3expWithAsympt_bis"){
      parameters[[mod]] = c("a1", "a2", "a4")
      parameters.names[[mod]] = c(TeX("$\\alpha_{1}$"),TeX("$\\alpha_{2}$"), TeX("$\\beta$"))
      
    } else if ((model.names[[mod]] == "Coutagne_bis")|(model.names[[mod]] == "expexp_bis")|(model.names[[mod]] == "hyperb_bis")){
      parameters[[mod]] = c("a1", "b1", "a2")
      parameters.names[[mod]] = c(TeX("$\\alpha$"), TeX("$\\lambda$"),  TeX("$\\beta$"))
    }
      
  
    ########################################################
    for (param in 1:length(parameters[[mod]])) {
    ########################################################
      dir.rec.segm.test.param[[mod]][[param]] = paste0(dir.segm.recessions.model.param,"/",parameters[[mod]][param])
      mu.results.df.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/mu.results.df.txt"), header=TRUE)
      gamma.results.df.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/gamma.results.df.txt"), header=TRUE)
      Q2.mu.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/Q2.mu.txt"), header=TRUE)
      mu.res.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/mu.res.txt"), header=TRUE)
      Q97.mu.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/Q97.mu.txt"), header=TRUE)
      Data.segm.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/Data.segm.rec.txt"), header=TRUE)
      mcmc.seg.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]], "/mcmc_segmentation.txt"), header=TRUE)
      gamma_segm_recess[[mod]][[param]]= c(mean(mcmc.seg.rec[[mod]][[param]]$Y1_gamma1), std(mcmc.seg.rec[[mod]][[param]]$Y1_gamma1))
      #
      #
      if (length(mu.results.df.rec[[mod]][[param]]$mu.MAP) >1) {
        tau.results.df.rec[[mod]][[param]]      = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/tau.results.df.txt"), header=TRUE)
        df.shift.times.rec[[mod]][[param]]      = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/df.shift.times.txt"), header=TRUE)
        df.shift.times.plus.rec[[mod]][[param]] = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/df.shift.times.plus.txt"), header=TRUE)
        ts.res.before.rec[[mod]][[param]]       = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/ts.res.before.txt"), header=TRUE)
        ts.res.rec[[mod]][[param]]              = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/ts.res.txt"), header=TRUE)
        ts.res.plus.rec[[mod]][[param]]         = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/ts.res.plus.txt"), header=TRUE)
        Q2.ts.rec[[mod]][[param]]               = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/Q2.ts.txt"), header=TRUE)
        Q97.ts.rec[[mod]][[param]]              = read.table(file=paste0(dir.rec.segm.test.param[[mod]][[param]] , "/Q97.ts.txt"), header=TRUE)
        nS.ok.rec[[mod]][[param]]               = nrow(df.shift.times.rec[[mod]][[param]])+1
        pdf.ts.rec[[mod]][[param]]              = mcmc.seg.rec[[mod]][[param]][, (nS.ok.rec[[mod]][[param]]+1):(2*nS.ok.rec[[mod]][[param]]-1)]
        shift.times.recessions[[mod]][[param]]  = read.table(file =paste0(dir.rec.segm.test.param[[mod]][[param]] ,
                                                                  "/df.shift.times.txt"), header = TRUE)
        data.annotate.recess[[mod]][[param]] = data.frame(t     = shift.times.recessions[[mod]][[param]]$ts.res,
                                                          tstart= shift.times.recessions[[mod]][[param]]$Q2.ts,
                                                          tend  = shift.times.recessions[[mod]][[param]]$Q97.ts)
        #data.annotate.recess.adjust[[mod]][[param]] = data.frame(t.adj = shift.times.recessions[[mod]][[param]]$ts.real)
        data.annotate.recess.adjust[[mod]][[param]] = data.frame(t.adj = shift.times.recessions[[mod]][[param]]$ts.res)
      } else {
        nS.ok.rec[[mod]][[param]]                   = 1
        shift.times.recessions[[mod]][[param]]      = list()
        pdf.ts.rec[[mod]][[param]]                  = list()
        data.annotate.recess[[mod]][[param]]        = list()
        data.annotate.recess.adjust[[mod]][[param]] = list()
      } 
      

      
      # plot dates of shifts:
      if (plot.dates == TRUE){
        dates.recess = 0
        for (dat in 1:length(data.annotate.recess.adjust[[mod]][[param]]$t.adj)) {
          dates.recess[dat] <-  data.annotate.recess.adjust$t.adj[[mod]][[param]][dat] +  Gaugings$X.3[1]
        }
        RealDates.recess <- as.Date(dates.recess,  origin = "1899-12-30")
        RealDates.recess.new <- strptime(as.character(RealDates.recess), "%Y-%m-%d" )
        RealDates.recess.newnew <- format( RealDates.recess.new, "%d/%m/%Y")
      } else {
        RealDates.recess.newnew = NULL
      }
      # gaugings to plot in the stage record with the periods derived from segmentation:
      gaugings.df.recess[[mod]][[param]] = read.table(file =paste0(dir.rec.segm.test.param[[mod]][[param]], 
                                                            "/data_with_periods.txt"), header = TRUE)
      # Quantiles for tot uncertainty for each parameter:
      Ysim[[mod]][[param]] = list();
      data.tmp.2[[mod]][[param]] = list()
      for (jj in 1:nS.ok.rec[[mod]][[param]]) {
        Ysim[[mod]][[param]][[jj]] =0
        for (ii in 1: length(mcmc.seg.rec[[mod]][[param]][,jj])) {
          Ysim[[mod]][[param]][[jj]][ii] = mcmc.seg.rec[[mod]][[param]][ii,jj] + 
                                           rnorm(1, 
                                                 mean = 0, 
                                                 sd = mcmc.seg.rec[[mod]][[param]][, 2*nS.ok.rec[[mod]][[param]]])
        }
        if (param == length(parameters[[mod]])) {
            data.tmp.2[[mod]][[param]][[jj]] = quantile(Ysim[[mod]][[param]][[jj]], probs=c(0.025,0.975) ) #, na.rm=TRUE)
        } else {
            delete.mcmc = which(Ysim[[mod]][[param]][[jj]] <0)
            data.tmp.2[[mod]][[param]][[jj]] = quantile(Ysim[[mod]][[param]][[jj]][-delete.mcmc], probs=c(0.025,0.975) ) #, na.rm=TRUE)
        }
      } 
      # Preparing vector with pdf of shift times :
      if (length(pdf.ts.rec[[mod]][[param]]) > 0) {
        X1[[mod]][[param]] = pdf.ts.rec[[mod]][[param]]
        if (nS.ok.rec[[mod]][[param]]==2) {
          X[[mod]][[param]] = X1[[mod]][[param]]
        } else {
          X[[mod]][[param]] = data.frame(time= X1[[mod]][[param]][,1], 
                                         ord = rep(1, length(X1[[mod]][[param]][,1])))
          for (orderr in 1:(ncol(pdf.ts.rec[[mod]][[param]]))) {
            X[[mod]][[param]] = rbind(X[[mod]][[param]], 
                                data.frame(time= X1[[mod]][[param]][,orderr],
                                           ord = rep(orderr, length(X1[[mod]][[param]][,1]))))
          }
        }
      }  
      

      # Second  plots: time series of parameters of the regression model with the shift times and b(t)
      #print(param)
      if (param == length(parameters[[mod]])) {
          inddd = which((Data.segm.rec[[mod]][[param]]$X97.5. < limits.Y[2]) & 
                       (Data.segm.rec[[mod]][[param]]$X2.5. > limits.Y[1]))
      } else if (parameters[[mod]][param] == "b1") {
          inddd = which((Data.segm.rec[[mod]][[param]]$X97.5. < limits.Y.lambda[2]) & 
                       (Data.segm.rec[[mod]][[param]]$X2.5. > limits.Y.lambda[1]))
      } else {
          inddd = which((Data.segm.rec[[mod]][[param]]$X97.5. < limits.Y.alpha[2]) & 
                        (Data.segm.rec[[mod]][[param]]$X2.5. > limits.Y.alpha[1]))
      } 
      dddf = Data.segm.rec[[mod]][[param]][inddd,]
      
      #############################
      if ( nS.ok.rec[[mod]][[param]] >1) {
      ##############################
        seg.rec.plot[[mod]][[param]]  <- ggplot()
        if (obs.uncertainty.y == TRUE) {
          seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +  
            # geom_errorbar(data = Data.segm.rec[[mod]][[param]],
            #               aes(x=t, 
            #                   ymin= (mean - 2*stdev), 
            #                   ymax = (mean+2*stdev)),
            #               size = 0.1, width=20, col= "black")
            geom_errorbar(data = dddf,
                          aes(x    = t, 
                              ymin = X2.5., 
                              ymax = X97.5.),
                          size = 0.2, width=80, col= "gray30")
        }
        time.adjust.before[[mod]][[param]] = c(0, 
                                               data.annotate.recess.adjust[[mod]][[param]]$t.adj)
        time.adjust.plus[[mod]][[param]]   = c(data.annotate.recess.adjust[[mod]][[param]]$t.adj, 
                                               limits.X[2])
        if (plot.gamma.uncertainty ==TRUE){
          for (jj in 1:nS.ok.rec[[mod]][[param]]){
            seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] + 
              annotate("rect", 
                       xmin = time.adjust.before[[mod]][[param]][jj], 
                       xmax = time.adjust.plus[[mod]][[param]][jj], 
                       ymin = data.tmp.2[[mod]][[param]][[jj]][1],
                       ymax = data.tmp.2[[mod]][[param]][[jj]][2],
                       fill = "red", alpha = 0.2)
          } 
        }
        
        seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
              # annotate("rect",
              #      xmin = time.adjust.before[[mod]][[param]],
              #      xmax = time.adjust.plus[[mod]][[param]],
              #      ymin = Q2.mu.rec[[mod]][[param]]$x,
              #      ymax = Q97.mu.rec[[mod]][[param]]$x, 
              #      fill ="red", alpha=0.2) +

          # geom_point(data = Data.segm.rec[[param]], aes(x = t, y = mean), size = 2, pch=1) +
          # geom_line(data = Data.segm.rec[[param]], aes(x = t, y = mean), size = 0.1, color = "blue")+
          geom_point(data = dddf, aes(x = t, y = mean), size = 1, pch=21, fill="gray30", color= "gray30")
          if (param == length(parameters[[mod]])) {
                seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
                                                annotate("rect", 
                                                         xmin = time.adjust.before[[mod]][[param]],
                                                         xmax = time.adjust.plus[[mod]][[param]],
                                                         ymin = mu.res.rec[[mod]][[param]]$x, 
                                                         ymax = mu.res.rec[[mod]][[param]]$x + 3, 
                                                         fill="red", alpha=1) 

          } else if (parameters[[mod]][param] == "b1") {
               seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
                       annotate("rect", 
                       xmin = time.adjust.before[[mod]][[param]],
                       xmax = time.adjust.plus[[mod]][[param]],
                       ymin = mu.res.rec[[mod]][[param]]$x, 
                       ymax = mu.res.rec[[mod]][[param]]$x + 0.5, 
                       fill="red", alpha=1) 
          } else {
               seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
                       annotate("rect", 
                       xmin = time.adjust.before[[mod]][[param]],
                       xmax = time.adjust.plus[[mod]][[param]],
                       ymin = mu.res.rec[[mod]][[param]]$x, 
                       ymax = mu.res.rec[[mod]][[param]]$x + 12, 
                       fill="red", alpha=1) 
          } 
          seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
          theme_light(base_size = 10) +
          theme(plot.title = element_text(hjust = 0.5))+
          scale_x_continuous(name = x.name, expand = c(0,0), limits =limits.X)+
          # geom_segment(x = time.adjust.before[[mod]][[param]],
          #              xend = time.adjust.plus[[mod]][[param]],
          #              y = mu.res.rec[[mod]][[param]]$x,
          #              yend = mu.res.rec[[mod]][[param]]$x, 
          #              col="red") +
          #geom_vline(xintercept = read.res.rec$ts.res, col="black", lwd =0.5, linetype= "dashed")+
          geom_vline(aes(xintercept = data.annotate.recess.adjust[[mod]][[param]]$t.adj), 
                     linetype= "solid", lwd = 0.9, color="blue") +
          coord_cartesian(clip = "off") +
          theme_light(base_size=15)+
          theme(axis.text=element_text(size=10)
                ,axis.title=element_text(size=15, face="plain")
                ,panel.grid.major=element_blank() #element_line(size=1.2)
                ,panel.grid.minor=element_blank()
                ,legend.text=element_text(size=20)
                ,legend.title=element_text(size=30)
                #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
                ,legend.key.size=unit(1.5, "cm")
                ,legend.position="none"
                ,plot.margin=unit(c(2, 0.5, 0, 1),"cm")
                ,axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20, face="bold"))
        
        ##############
      } else {
        ##############
          seg.rec.plot[[mod]][[param]] <- ggplot()
          if (obs.uncertainty.y == TRUE) {
              seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +  
              # geom_errorbar(data = Data.segm.rec[[mod]][[param]],
              #               aes(x=t, 
              #                   ymin= (mean - 2*stdev), 
              #                   ymax = (mean+2*stdev)),
              #               size = 0.3, width=40, col= "black")
              geom_errorbar(data = dddf,
                            aes(x    = t, 
                                ymin = X2.5., 
                                ymax = X97.5.),
                            size = 0.2, width=80, col= "gray30")
          }
          seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +  
          #geom_line(data = Data.segm.rec[[param]], aes(x = t, y = mean), size = 0.1, color = "blue")+
          geom_point(data = Data.segm.rec[[mod]][[param]], aes(x = t, y = mean),  size = 1, pch=21, fill="gray30", color= "gray30") +
          theme_light(base_size = 10) +  
          coord_cartesian(clip = "off") +
          theme(plot.title = element_text(hjust = 0.5))
        
          if (plot.gamma.uncertainty ==TRUE){
          seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]]+
                                          annotate("rect",  
                                                   xmin = Data.segm.rec[[mod]][[param]]$t[1], 
                                                   xmax = tail(Data.segm.rec[[mod]][[param]]$t,1),
                                                   ymin = data.tmp.2[[mod]][[param]][[1]][1],
                                                   ymax = data.tmp.2[[mod]][[param]][[1]][2],
                                                   fill = "red", alpha=0.2)
         }
         
          # annotate("rect", 
          #          xmin= Data.segm.rec[[mod]][[param]]$t[1], 
          #          xmax= tail(Data.segm.rec[[mod]][[param]]$t,1),
          #          ymin= Q2.mu.rec[[mod]][[param]]$x,
          #          ymax= Q97.mu.rec[[mod]][[param]]$x, 
          #          fill="red", alpha=0.2) 
          
          if (param == length(parameters[[mod]])) {
            seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
              annotate("rect",
                       xmin= Data.segm.rec[[mod]][[param]]$t[1],
                       xmax= tail(Data.segm.rec[[mod]][[param]]$t,1),
                       ymin= mu.res.rec[[mod]][[param]]$x,
                       ymax= mu.res.rec[[mod]][[param]]$x + 3,
                       fill="red", alpha=1)
            
          } else if (parameters[[mod]][param] == "b1") {
            seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
              annotate("rect",
                       xmin= Data.segm.rec[[mod]][[param]]$t[1],
                       xmax= tail(Data.segm.rec[[mod]][[param]]$t,1),
                       ymin= mu.res.rec[[mod]][[param]]$x,
                       ymax= mu.res.rec[[mod]][[param]]$x + 0.5,
                       fill="red", alpha=1)
          } else {
            seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
              annotate("rect",
                       xmin= Data.segm.rec[[mod]][[param]]$t[1],
                       xmax= tail(Data.segm.rec[[mod]][[param]]$t,1),
                       ymin= mu.res.rec[[mod]][[param]]$x,
                       ymax= mu.res.rec[[mod]][[param]]$x + 12,
                       fill="red", alpha=1)
          } 

          
          # geom_segment(x     = Data.segm.rec[[mod]][[param]]$t[1], 
          #              xend  = tail(Data.segm.rec[[mod]][[param]]$t,1),
          #              y     = mu.res.rec[[mod]][[param]]$x, 
          #              yend  = mu.res.rec[[mod]][[param]]$x, 
          #              color = "red") +
          seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
          scale_x_continuous(name = x.name, expand = c(0,0), limits =limits.X)+
          theme_light(base_size=15)+
          theme(axis.text=element_text(size=10)
                ,axis.title=element_text(size=15, face="plain")
                ,panel.grid.major=element_blank() #element_line(size=1.2)
                ,panel.grid.minor=element_blank()
                ,legend.text=element_text(size=20)
                ,legend.title=element_text(size=30)
                #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
                ,legend.key.size=unit(1.5, "cm")
                ,legend.position="none"
                ,plot.margin=unit(c(2, 0.5, 0, 1),"cm")
                ,axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20, face="bold"))
        # # add the title with the name of the model:
        #   theme(axis.text=element_text(size=10)
        #         ,axis.title=element_text(size=10,face="plain")
        #         ,axis.title.y=element_text(size=20,face="plain")
        #         # ,panel.grid.major=element_line(size=1.2),
        #         # panel.grid.minor=element_line(size=0.8)
        #         ,legend.text=element_text(size=20)
        #         ,legend.title=element_text(size=30)
        #         ,legend.key.size=unit(1.5, "cm")
        #         ,legend.position="none"
        #         ,plot.margin = margin(1, 2, 1, 1, "cm")
        #         ,axis.line = element_line(colour = "black"))
      }
      
      ####################################################################
      if (param == length(parameters[[mod]])) {
        if (plot.b.from.gaugings == TRUE) {
          seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] + 
            # geom_segment(mapping= aes(x = c(limits.X[1], bt.from.gaugingsss$t.shift.for.b$treal), 
            #                           y = bt.from.gaugingsss$bt1.df$maxpost*100, 
            #                           xend = c(bt.from.gaugingsss$t.shift.for.b$treal, limits.X[2]),  
            #                           yend = bt.from.gaugingsss$bt1.df$maxpost*100), 
            #              color = "black", 
            #              size = 0.5,
            #              linetype = "dashed") +
            # geom_segment(mapping= aes(x = c(limits.X[1], bt.from.gaugingsss$t.shift.for.b$treal), 
            #                           y = bt.from.gaugingsss$bt2.df$maxpost*100, 
            #                           xend = c(bt.from.gaugingsss$t.shift.for.b$treal, limits.X[2]),  
            #                           yend = bt.from.gaugingsss$bt2.df$maxpost*100), 
            #              color = "blue", 
            #              size = 0.5,
            #              linetype = "dashed")+
            theme(axis.text=element_text(size=10)
                  ,axis.title=element_text(size=15, face="plain")
                  ,panel.grid.major=element_blank() #element_line(size=1.2)
                  ,panel.grid.minor=element_blank()
                  ,legend.text=element_text(size=20)
                  ,legend.title=element_text(size=30)
                  #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
                  ,legend.key.size=unit(1.5, "cm")
                  ,legend.position="none"
                  ,plot.margin=unit(c(2, 0.5, 0, 1),"cm")
                  ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=20, face="bold"))
            
            # geom_rect(mapping = aes(xmin = ts.before.gaug, 
            #                         xmax = ts.plus.gaug, 
            #                         ymin = bt2.df$X2.5.,
            #                         ymax = bt2.df$X97.5.), 
            #           fill="red", alpha=0.3)+
            # geom_vline(aes(xintercept = ts.gaug$t.adj), color = "black",
            #            lwd =0.7, linetype = "dotdash")
        }
        # fix the axis scale y for the asymptotic level:
        seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] + 
                                        scale_y_continuous(name =  parameters.names[[mod]][param], 
                                                           limits =c(limits.Y[1], limits.Y[2]), 
                                                           expand = c(0.01, 0.01))
      #########
      }  else if (parameters[[mod]][param] == "b1") { 
        seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
          scale_y_continuous(name =  parameters.names[[mod]][param],
                             limits =c(limits.Y.lambda[1], limits.Y.lambda[2]),  ##############    change this !!!!!!!!!!!!!!!!!!!!!!!!
                             expand = c(0.01,0.01))
        
        
      }  else {
      #########
        # automatic y scale axis:
        seg.rec.plot[[mod]][[param]] =  seg.rec.plot[[mod]][[param]] +
                                        scale_y_continuous(name =  parameters.names[[mod]][param],
                                                           limits =c(limits.Y.alpha[1], limits.Y.alpha[2]),  ##############    change this !!!!!!!!!!!!!!!!!!!!!!!!
                                                           expand = c(0.01,0.01))
      }
      # TO VERIFY::
      
      # seg.rec.plot = seg.rec.plot+
      #   coord_cartesian(clip = 'off')+
      #   theme_bw(base_size=20)+
      #   theme(axis.text=element_text(size=15)
      #         ,axis.title=element_text(size=20,face="plain")
      #         ,panel.grid.major=element_blank()
      #         ,panel.grid.minor=element_blank()
      #         ,legend.text=element_text(size=20)
      #         ,legend.title=element_text(size=30)
      #         ,legend.key.size=unit(1.5, "cm")
      #         ,legend.position="none"
      #         ,plot.margin=unit(c(2, 0.5, 2, 2),"cm")
      #         ,axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l =0 ))
      #         ,axis.title.x = element_blank())
      pdf(paste0(dir.segm.recessions,"/Figure_segment_model_",model.names[mod],"_chi",chi.test,"par_",parameters[[mod]][param],".pdf"), 15,10 ,useDingbats=F)
      print(seg.rec.plot[[mod]][[param]])
      dev.off()
    }
    
    # combine the plots with the parameters time series:
    if ((any(model.names =="2expWithAsympt_bis"))|(any(model.names =="3expWithAsympt_bis"))|(any(model.names =="expexp_bis"))|(any(model.names =="hyperb_bis"))|(any(model.names =="Coutagne_bis"))){
    if (length(parameters[[mod]])==2){
        # t.plot3[[mod]] = plot_grid(seg.rec.plot[[mod]][[1]], NULL, NULL, seg.rec.plot[[mod]][[2]],
        #                         ncol = 1, nrow = 4, rel_heights = c(1, 1, 1, 1))
        t.plot3[[mod]] = plot_grid(seg.rec.plot[[mod]][[1]],
                                   NULL,
                                   seg.rec.plot[[mod]][[2]],
                                   ncol = 1, nrow = 3, rel_heights = c(1, 1, 1))
    } else if (length(parameters[[mod]])==3){
        t.plot3[[mod]] = plot_grid(seg.rec.plot[[mod]][[1]],
                                   seg.rec.plot[[mod]][[2]],
                                   seg.rec.plot[[mod]][[3]],
                                   ncol = 1, nrow = 3, rel_heights = c(1, 1, 1))
    }
    } else { 
       t.plot3[[mod]] = plot_grid(seg.rec.plot[[mod]][[1]],
                                  seg.rec.plot[[mod]][[2]],
                                  ncol = 1, nrow = 2, 
                                  rel_heights = c(1, 1))
       pdf(paste0(dir.segm.recessions,"/Figure_segment_model_",model.names[mod],"_chi",chi.test,"all_param.pdf"), 15,15 ,useDingbats=F)
       print(t.plot3[[mod]]) 
       dev.off()
    } 
    
    
    
    
    
    
    
    
    
    ###########################################################################################
    #PLot of stage record with the shifts:
    print("plot 3")
    t.plot[[mod]] <- ggplot()+
                     scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.X) +
                     scale_y_continuous(name=limni.labels[2], expand = c(0,0), 
                                        limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), 
                                        breaks=seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3])) +
                     ylab(limni.labels[2]) +
                     coord_cartesian(clip = 'off') +
                     geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray60", size = 0.2) + # limni
                     geom_point(data=gaugings.df.recess[[mod]][[length(parameters[[mod]])]],                 # gaugings
                                  aes(x = t , y= h), size = 1, color = gaugings.df.recess[[mod]][[length(parameters[[mod]])]]$color) +
                       theme_light(base_size=15)+
                       theme(axis.text=element_text(size=10)
                             ,axis.title=element_text(size=15,face="plain")
                             ,panel.grid.major=element_blank()
                             ,panel.grid.minor=element_blank()
                             ,legend.text=element_text(size=20)
                             ,legend.title=element_text(size=30)
                             ,legend.key.size=unit(1.5, "cm")
                             ,legend.position="none"
                             ,plot.margin=unit(c(2, 0.5, 0, 1),"cm")
                             ,axis.title.y = element_text(margin = margin(t = 0, r = 23, b = 0, l = 0))
                             ,axis.text.x=element_blank()
                             ,axis.line.x = element_blank())

    ###################
    if ( nS.ok.rec[[mod]][[length(nS.ok.rec[[mod]])]] >1) {
      if (plot.gamma.uncertainty == TRUE){
        for (jj in 1:nS.ok.rec[[mod]][[length(nS.ok.rec[[mod]])]]){
          t.plot[[mod]]  = t.plot[[mod]]  + 
            annotate("rect", 
                     xmin = time.adjust.before[[mod]][[length(nS.ok.rec[[mod]])]][jj], 
                     xmax = time.adjust.plus[[mod]][[length(nS.ok.rec[[mod]])]][jj], 
                     ymin = data.tmp.2[[mod]][[length(nS.ok.rec[[mod]])]][[jj]][1]/100,
                     ymax = data.tmp.2[[mod]][[length(nS.ok.rec[[mod]])]][[jj]][2]/100,
                     fill = "red", alpha = 0.2)
        }
      }
      t.plot[[mod]]  = t.plot[[mod]]  + 
        # annotate("rect", 
        #          xmin = time.adjust.before[[mod]][[length(nS.ok.rec[[mod]])]], 
        #          xmax = time.adjust.plus[[mod]][[length(nS.ok.rec[[mod]])]], 
        #          ymin = unlist(Q2.mu.rec[[mod]][[length(nS.ok.rec[[mod]])]])/100,
        #          ymax = unlist(Q97.mu.rec[[mod]][[length(nS.ok.rec[[mod]])]])/100, 
        #          fill ="red", alpha=0.2)+
        #geom_vline(aes(xintercept = data.annotate.recess.adjust$t.adj), color = "blue", lwd =0.3, linetype = "dashed")+
        geom_segment(mapping=aes(x = time.adjust.before[[mod]][[length(parameters[[mod]])]] , 
                                 y = unlist(mu.res.rec[[mod]][[length(parameters[[mod]])]])/100, 
                                 xend = time.adjust.plus[[mod]][[length(parameters[[mod]])]], 
                                 yend = unlist(mu.res.rec[[mod]][[length(parameters[[mod]])]])/100), col="red")
      #################
    } else {
      #################  
      if (plot.gamma.uncertainty ==TRUE){
        t.plot[[mod]]  = t.plot[[mod]] + annotate("rect",  
                                                  xmin = Data.segm.rec[[mod]][[length(nS.ok.rec[[mod]])]]$t[1], 
                                                  xmax = tail(Data.segm.rec[[mod]][[length(nS.ok.rec[[mod]])]]$t,1),
                                                  ymin = data.tmp.2[[length(nS.ok.rec[[mod]])]][[1]][1]/100,
                                                  ymax = data.tmp.2[[length(nS.ok.rec[[mod]])]][[1]][2]/100,
                                                  fill = "red", alpha=0.2)
      }
      # t.plot[[mod]]  = t.plot[[mod]]  + 
      #   annotate("rect", 
      #            xmin= Data.segm.rec[[mod]][[length(nS.ok.rec[[mod]])]]$t[1], 
      #            xmax= tail(Data.segm.rec[[mod]][[length(nS.ok.rec[[mod]])]]$t,1),
      #            ymin= unlist(Q2.mu.rec[[mod]][[length(nS.ok.rec[[mod]])]])/100,
      #            ymax= unlist(Q97.mu.rec[[mod]][[length(nS.ok.rec[[mod]])]])/100, fill="red", alpha=0.1)
    }

    
    if (plot.dates == TRUE){
      if (!is.null(dates[[length(parameters[[mod]])]])){
        t.plot[[mod]]  <- t.plot[[mod]]  +
          annotate("text", 
                   x = read.res.rec$data.annotate.recess.adjust$t.adj,
                   y = grid_limni.ylim[1] + 0.1, 
                   label = read.res.rec$dates, color = "red", size=5)
      }
    }
    # ##########################################################################################
    # # PLot of the pdf of the shifts:
    # t.plot2[[mod]]  <- ggplot() +
    #   scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,0.5))+
    #   xlab(limni.labels[1])+ 
    #   coord_cartesian(clip = 'off') + 
    #   geom_hline(yintercept = c( -0.5,0, 0.5), color="darkgray", linetype="dashed", size = 0.5)+
    #   theme_light(base_size=15)+
    #   theme(axis.text=element_text(size=10)
    #         ,axis.title=element_text(size=15,face="plain")
    #         ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
    #         ,legend.text=element_text(size=20),legend.title=element_text(size=30)
    #         ,legend.key.size=unit(1.5, "cm"),legend.position="none"
    #         ,plot.margin=unit(c(0, 0.5, 0, 1),"cm")
    #         ,axis.title.y = element_text(margin = margin(t = 0, r = 43, b = 0, l = 0))
    #         ,axis.text.y=element_blank()
    #         ,axis.ticks.y = element_blank()) +
    #   scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.X)
    # 
    # if (is.null(data.annotate.off)==FALSE) {
    #   t.plot2[[mod]]  <- t.plot2[[mod]]  +
    #     geom_point(data = data.annotate.off, 
    #                aes(x = xeffect, y = -1), color= "black", size = 3, shape =4, stroke=1.5 )
    # }
    # if (!is.null(pdf.ts.rec[[mod]][[length(parameters[[mod]])]])) {
    #   if (nS.ok.rec[[mod]][[length(nS.ok.rec[[mod]])]] ==2) {
    #     t.plot2[[mod]]  = t.plot2[[mod]]  + 
    #       geom_density(aes(x= X[[mod]][[length(parameters[[mod]])]], ..scaled.. /2),
    #                    fill="blue", 
    #                    colour=NA, alpha=0.4)
    #   } else {
    #     t.plot2[[mod]]  = t.plot2[[mod]]  +
    #       geom_density(aes(x= X[[mod]][[length(parameters[[mod]])]]$time, 
    #                        ..scaled.. /2, 
    #                        group =X[[mod]][[length(parameters[[mod]])]]$ord),  
    #                    fill= "blue", 
    #                    colour=NA, alpha=0.4)
    #   }
    #   t.plot2[[mod]]  = t.plot2[[mod]]  + 
    #     geom_segment(data = data.annotate.recess.adjust[[mod]][[length(parameters[[mod]])]], 
    #                  aes(x = t.adj, y = -0.5, yend =0, xend=t.adj ),
    #                  size = 0.9, color ="blue", linetype = "solid")
    # }
    # 
    
    ##########################################################################################
    # PLot of the pdf of the shifts:
    t.plot2[[mod]]  <- ggplot() +
      scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 0.5))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off') + 
      geom_hline(yintercept = c( 0.5), color="darkgray", linetype="dashed", size = 0.5)+
      theme_light(base_size=15)+
      theme(axis.text=element_text(size=10)
            ,axis.title=element_text(size=15,face="plain")
            ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20),legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm"),legend.position="none"
            ,plot.margin=unit(c(0, 0.5, 0, 1),"cm")
            ,axis.title.y = element_text(margin = margin(t = 0, r = 45, b = 0, l = 0))
            ,axis.text.y=element_blank()
            ,axis.ticks.y = element_blank()) +
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.X)
    if (!is.null(pdf.ts.rec[[mod]][[length(parameters[[mod]])]])) {
      if (nS.ok.rec[[mod]][[length(nS.ok.rec[[mod]])]] ==2) {
        t.plot2[[mod]]  = t.plot2[[mod]]  + 
          geom_density(aes(x= X[[mod]][[length(parameters[[mod]])]], ..scaled.. /2),
                       fill="blue", 
                       colour=NA, alpha=0.4)
      } else {
        t.plot2[[mod]]  = t.plot2[[mod]]  +
          geom_density(aes(x= X[[mod]][[length(parameters[[mod]])]]$time, 
                           ..scaled.. /2, 
                           group =X[[mod]][[length(parameters[[mod]])]]$ord),  
                       fill= "blue", 
                       colour=NA, alpha=0.4)
      }
    }
    
    
    ###########################################################################################
    #cowplot "package", combine plot stage record h(t) with the plot of pdf of the shifts:
    if (plot_ts_gaugings==TRUE){
      if (mod==1){
        t.plot4[[mod]] = plot_grid(reg.pool.plot[[mod]], t.plot3[[mod]], t.plot[[mod]] , t.plot2[[mod]] , 
                                   ncol = 1, nrow = 4, rel_heights = c(1.1, 2, 0.8, 0.15),  labels = c('a) ', 
                                                                                                      'b)' ,
                                                                                                      'c)', 
                                                                                                      ''), 
                                   label_size = 25, label_fontface = "bold")
      } else {
        t.plot4[[mod]] = plot_grid(reg.pool.plot[[mod]], t.plot3[[mod]], t.plot[[mod]] , t.plot2[[mod]] , 
                                   ncol = 1, nrow = 4, rel_heights = c(1.1, 2, 0.8, 0.15),  labels = c('', '' ,'', ''), 
                                   label_size = 25, label_fontface = "plain")
      }

    } else {
      t.plot[[mod]]  = t.plot[[mod]] +
                       scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.X)+
                       theme(axis.text.x = element_text(size=10))
      
        
      if (mod==1){
        t.plot4[[mod]] = plot_grid(reg.pool.plot[[mod]], t.plot3[[mod]], t.plot[[mod]], 
                                   ncol = 1, nrow = 3, rel_heights = c(1.1, 2, 0.8, 0.3),  labels = c('a) ', 
                                                                                                      'b)' ,
                                                                                                      'c)'), 
                                   label_size = 25, label_fontface = "bold")
      } else {
        t.plot4[[mod]] = plot_grid(reg.pool.plot[[mod]], t.plot3[[mod]], t.plot[[mod]], 
                                   ncol = 1, nrow = 3, rel_heights = c(1.1, 2, 0.8),  labels = c('', '' ,''), 
                                   label_size = 25, label_fontface = "plain")
      }
      
    }
    pdf(paste0(dir.segm.recessions,"/Figure_article_recession_model_", model.names[mod], ".pdf"), 10, 16 ,useDingbats=F)
    print( t.plot4[[mod]])
    dev.off()
  }
  
  
  
  ######################################################################################################
  # final plot with 6 columns: one for each recession model
  # plot4paper.final = plot_grid(t.plot4[[1]] , t.plot4[[2]] ,  t.plot4[[3]] ,  t.plot4[[4]] ,  t.plot4[[5]] ,  t.plot4[[6]] , 
  #                              ncol = 6, nrow =1,
  #                              rel_widths = c(1, 1,1,1,1,1))
  plot4paper.final = plot_grid(t.plot4[[1]] , t.plot4[[2]] ,  t.plot4[[3]] ,
                               ncol = 3, nrow =1,
                               rel_widths = c(1, 1, 1)) 
  pdf(paste0(dir.segm.recessions,"/Figure_article_recession_models_comparison.pdf"), 15,20 ,useDingbats=F)
  print(plot4paper.final)
  dev.off()
  ########################################################################################################
}






































####################################################################
prior.mcmc.bam = function(mcmc,  distrib, prior){
  ####################################################################
  if ((distrib == "'Uniform'")){
    #print("ok")
    priors.param =  dunif(mcmc, 
                          min = prior[1], 
                          max = prior[2], 
                          log = TRUE)
  } else if ((distrib == "'LogNormal'")){
    priors.param =   dlnorm(mcmc, 
                            meanlog = prior[1],  
                            sdlog   = prior[2],  
                            log     = TRUE)
  } else if ((distrib == "'Gaussian'")){
    priors.param =  dnorm(mcmc, 
                          mean = prior[1], 
                          sd   = prior[2], 
                          log  = TRUE)
  } else if (distrib =="'FlatPrior'"){
    priors.param = 0
  }
  return(priors.param)
}






























# fonction that computes the accuracy of the segmentation:
#######################################################################################
accuracy.rec = function(time, ts, nshift.found, nshift.true, t2, t97, tMAP , dir) {
####################################################################################### 
  TP =0; TN=1; FP=0; FN=0;      j =1; l = 1
  nobs.synthet = length(time)
  shift.time.TP=NULL; real.shift.time.TP=NULL; 
  tMAP.ok = vector(mode="list", length=nshift.true) 
  IND.MAP.ok = vector(mode="list", length=nshift.true) 
  TP=0; FP=0; FN=0; TN=0;  tMAP.TP=0; tTRUE.TP=0;
  type.point =rep("TN",length =nobs.synthet)
  type.MAP =rep("FP",length =nshift.found) 
  MAP.taken=rep(FALSE, length =nshift.found) 
  cross.taken = rep(FALSE, length =nshift.true) 
  index.ts.map = NULL;
  
  # find the True Positive points:
  ################################
  for (j in 1:nshift.true) {
    k =0; 
    for (s in 1:nshift.found) {
      if (((t2[s] <= ts[j]) & (t97[s] >= ts[j]))) {
        k= k+1
        tMAP.ok[[j]][k] = tMAP[s]
        IND.MAP.ok[[j]][k] = s
      }  
    }
    if (is.null(tMAP.ok[[j]]) == FALSE ) {
      TP = TP+1
      differ = abs(ts[j]-tMAP.ok[[j]])
      INDEX = which.min(differ)
      tMAP.TP[TP] = tMAP.ok[[j]][INDEX]
      tTRUE.TP[TP] = ts[j]
      type.MAP[IND.MAP.ok[[j]][INDEX]] = "TP"
      MAP.taken[IND.MAP.ok[[j]][INDEX]]  =TRUE
      differ.TP = abs(ts[j] - time)
      INDEX.TP = which.min(differ.TP)
      type.point[INDEX.TP] = "TP"
      index.ts.map[TP] = IND.MAP.ok[[j]][INDEX]
      #==============================> TP  number of true cp classified as cp !!
      cross.taken[j] = TRUE
    }
  }
  
  #find the False Negative points:
  ################################
  if (length(tMAP.TP) >1) {
    for (tt in 2:length(tMAP.TP)) {
      if (tMAP.TP[tt]== tMAP.TP[tt-1]) {
        # many true change points for one found cp!
        # find the true cp nearest to the found cp:
        aa = which.min(c(abs(tTRUE.TP[tt] - tMAP.TP[tt]), abs(tTRUE.TP[tt-1] - tMAP.TP[tt-1])))
        if (aa ==1) {
          next
        } else {
          type.MAP[tt] = "FN"
          MAP.taken[tt]  =FALSE
          differ.TP = abs(ts[tt] - time)
          INDEX.TP = which.min(differ.TP)
          type.point[INDEX.TP] = "FN"
          TP =TP-1
          FN= FN+1
        }
        
      }
    }
  }
  
  # find the False Positive points:
  #################################
  for (i in 2:nobs.synthet) {
    if ((type.point[i] != "TP") &  (type.point[i] != "FN")){
      for (s in 1:nshift.found) {
        if ((time[i] >= tMAP[s]) & (time[i-1] < tMAP[s])) {
          if ((type.MAP[s] != "TP") & (MAP.taken[s]==FALSE)) {
            FP = FP +1 
            type.point[i] = "FP" #====================> FP  number of true non-cp classified as cp !!
            break
          }
        }
      }
      for (j in 1:nshift.true) {
        if (((time[i] >= ts[j]) & (time[i-1] < ts[j])) & (cross.taken[j] == FALSE)) {
          FN = FN +1
          type.point[i] = "FN"  #=====================> FN number of true cp classified as non-cp !!
        }
      }
    }
    #print("ok")
  }
  
  # the True Negative points:
  ##########################
  TN= nobs.synthet - FP - FN - TP  #==================> TN number of true non cp classified as non-cp 
  if (tail(ts,1) > tail(time,1)) {
    FN=FN +1
  } 
  
  # Final results:
  #################################
  print(c("TN=",TN, "FP=", FP, "FN=", FN, "TP=",TP))
  points = data.frame(time = time , y= rep(0, nobs.synthet), type = type.point)
  plot.accuracy = ggplot()+
                  scale_x_continuous(name="Time [year]", expand = c(0.1,0.1)) +
                  scale_y_continuous(name=element_blank(), limits = c(-1,1), expand = c(0,0)) +         
                  geom_point(aes(x= ts, y=rep(0, nshift.true)), size = 4, shape =4, stroke=2, col ="black" )+
                  geom_vline(xintercept = tMAP, col="blue")+
                  geom_rect(mapping= aes(xmin= t2, xmax=t97 ,ymin=-Inf, ymax=Inf), fill="blue", alpha=0.1)+
                  annotate(geom="text", x=time, y=rep(-0.1, nobs.synthet), label=type.point, color="red", size=3)+       
                  geom_point(data = points, aes(x=time, y=y, color =points$type), size=3, shape =16)+
                  theme_bw(base_size = 25)
  ggsave(plot.accuracy, filename =paste0(dir,"/accuracy_points.png"),  device = "png", width = 24, height =8,
         dpi = 400, units = "in")
  # calculation of Accuracy and RMSE:
  accuracy = (TP+TN)/(TP+FP+FN+TN)
  print(c("accuracy=",accuracy))

  return(list(accuracy = accuracy, 
              tMAP.TP  = tMAP.TP, 
              tTRUE.TP = tTRUE.TP,
              point.type = data.frame(TN=TN, FP=FP, FN=FN, TP=TP) ))
}














































#####################################################################################################
plot.performance.model.comparison = function(model.names, 
                                             model.titles,
                                             chi.test,
                                             dir.case_study,
                                             priors, 
                                             which.recession,
                                             bt.from.gaugingsss,
                                             pdf.ts.gaugings,
                                             data.annotate.off,
                                             time.limits
                                            ){
#####################################################################################################
  #General directories:
  dir.regression = paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression")
  dir.rec.pool = paste0(dir.regression,"/Pooling")
  dir.segm.recessions = paste0(dir.case_study,"/Results/segmentation_recessions/segmentation")
  npar=NULL; df.mcmc.LL =NULL;  t.b.comparison =NULL; t.deltab.comparison =NULL;
  DIC.rec =NULL; RMSE.b1.rec =NULL; RMSE.b2.rec =NULL; Accuracy.rec =NULL;  accu.rec=NULL; rmse.b1.rec =NULL; rmse.b2.rec =NULL;
  t.plot = NULL; X1.new=NULL; X1 =NULL;
  shift.times.recessions =NULL;
  df.ts =NULL;
  shift.rmse.rec =NULL; shift.rootmse.rec =NULL;
  colllor =c("blue", NA, NA)
  shapppe = c(15,    19 ,  17)
  type.plot = c(NA, "red", "green")
  linetype.models  =  c("solid", "dashed" , "dotdash", "dotted", "longdash", 
                        "solid", "dashed" , "dotdash", "dotted", "longdash", 
                        "solid")
  color.models = c("red", "blue", "green", "black", "orange", "cyan", "brown", "green", "violet", "yellow", "gray")
  
  
  for (m in 1:length(model.names)){
    title <- paste0("M", m)
    t.plot[[m]]  = ggplot()  + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0,1)) + 
    ggtitle(title) + 
    xlab(limni.labels[1])+ 
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+  
    #theme_set(theme_gray(base_size = 20))+
    theme_classic(base_size = 20) +
    theme(
        axis.line.y = element_blank()
        ,axis.line.x = element_blank()
        ,axis.title=element_text(size= 20,face="bold")
        ,plot.title = element_text(face = "bold", size = 15, color= "black", margin=margin(0,0,0,0))
        ,panel.grid.major=element_blank()
        ,panel.grid.minor=element_blank()
        ,panel.background=element_blank()
        ,legend.text=element_text(size=20)
        ,legend.title=element_text(size=30)
        ,legend.key.size=unit(1.5, "cm")
        ,legend.position="none"
        ,plot.margin=unit(c(0, 1, 0, 0.5),"cm")
        # ,axis.text.x =element_blank()
        ,axis.text = element_blank()
        ,axis.ticks = element_blank())+
    scale_x_continuous(name=element_blank(), expand = c(0,0), limits =time.limits)
  }
  
  
  
  
  
  
  
  
  
  # do a specific computation for the results of each model and of each chi value:  
  ################################################################################
  for (chi.i in 1:length(chi.test)){
  ##################################
       print(paste0("CHI :  ", chi.test[chi.i], "  ##############################################################################################"))
    
       DIC.rec[[chi.i]]= list()
       RMSE.b1.rec[[chi.i]]= list()
       RMSE.b2.rec[[chi.i]]= list()
       Accuracy.rec[[chi.i]]= list()
       shift.rootmse.rec[[chi.i]]= list()
       accu.rec[[chi.i]]= list()
       rmse.b1.rec[[chi.i]]= list()
       rmse.b2.rec[[chi.i]]= list()
       shift.rmse.rec[[chi.i]]= list()
       X1.new[[chi.i]] = list()
       X1[[chi.i]] = list()
       df.ts[[chi.i]] = list()
       
       mu.results.df.rec = NULL; gamma.results.df.rec = NULL;  Q2.mu.rec  = NULL; mu.res.rec  = NULL; 
       Q97.mu.rec  = NULL;  Data.segm.rec = NULL;  mcmc.seg.rec = NULL;
       gamma_segm_recess = NULL;  tau.results.df.rec = NULL; df.shift.times.rec = NULL;
       df.shift.times.plus.rec = NULL; ts.res.before.rec = NULL; ts.res.rec = NULL; ts.res.plus.rec = NULL; Q2.ts.rec = NULL;
       Q97.ts.rec = NULL; nS.ok.rec = NULL; pdf.ts.rec = NULL;
       data.annotate.recess = NULL; data.annotate.recess.adjust = NULL;
       shift.times.recessions[[chi.i]] = list();
       Ysim = NULL;
       data.tmp.2 =NULL;
       
       
       for (mod in 1:length(model.names)) {
       #########################################
           print(paste0("compute DIC for model: ",model.names[mod],"... ",mod,"/",length(model.names)," ********************"))
           # specific directories for each test:
           print(paste0("chi ",chi.test[chi.i]))
           dir.rec.pool.model = paste0(dir.rec.pool,"/test_", model.names[mod])   
           dir.rec.pool.test = paste0(dir.rec.pool.model,"/chi_",chi.test[chi.i])
           dir.segm.recessions.model = paste0(dir.segm.recessions,"/", model.names[mod])
           dir.segm.recessions.model.param = paste0(dir.segm.recessions.model,"/chi_",chi.test[chi.i])
           data.rec.curve = read.table(paste0(dir.rec.pool.test,"/Curves_Data.txt"),header=TRUE,dec=".", sep="")
           ncurv =  tail(data.rec.curve$Period, 1) #length(which.recession) 
           # read data mcmc:  
           data.MCMC.cooked=as.matrix(read.table(paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"),header=TRUE,dec=".", sep=""))
           # extract the prior from the mcmc (this depends on the model):
           priors.gamma = 0
           len.mcmc = length(data.MCMC.cooked[,1])
           logpost = data.MCMC.cooked[,ncol(data.MCMC.cooked)]
           npar = ncol(data.MCMC.cooked)-1
           loglikelihood=0;
           singlelikelihoods = 0;
 
           
           if (model.names[[mod]] == "1expWithAsympt"){
           ############################################
                 all.parameters = c("a1(k)", "b1", "a2(k)", "g1", "g2")
                 par.var = c("a1", "a2")
                 priors.param = matrix(NA, nrow = length(data.MCMC.cooked[,1]), ncol = npar)
                 for (ii in 1:ncurv){
                   priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                      distrib = as.character(priors[[mod]][3]),
                                                      prior   = c(as.numeric(priors[[mod]][1]), 
                                                                  as.numeric(priors[[mod]][2]))
                   )
                 }
                 priors.param[,ncurv+1] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                         distrib = priors[[mod]][8],
                                                         prior   = c(as.numeric(priors[[mod]][6]),
                                                                     as.numeric(priors[[mod]][7]))
                                                         )
                 for (ii in (ncurv+2):(ncurv*2+1)){
                   priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                      distrib = priors[[mod]][13],
                                                      prior   = c(as.numeric(priors[[mod]][11]),
                                                                  as.numeric(priors[[mod]][12]))
                   )
                 }  
                 priors.param[, (npar-1)] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,(npar-1)],
                                                           distrib = "'Uniform'",
                                                           prior   = c(0, 100))
                 
                 priors.param[,npar] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,npar],
                                                      distrib = "'Uniform'",
                                                      prior   = c(0, 100))
             
                 

             
           } else if ((model.names[[mod]] =="expexp")|(model.names[[mod]] =="hyperb")|(model.names[[mod]] =="Coutagne")){  
           ################################################################################################################
                 parameters = c("a1(k)", "b1", "n1", "a2(k)", "g1", "g2")
                 par.var = c("a1", "a2")
                 priors.param = matrix(NA, nrow = length(data.MCMC.cooked[,1]), ncol = npar)
                 for (ii in 1:ncurv){
                   priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                      distrib = as.character(priors[[mod]][3]),
                                                      prior   = c(as.numeric(priors[[mod]][1]), 
                                                                  as.numeric(priors[[mod]][2]))
                   )
                 }
                 priors.param[,ncurv+1] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                         distrib = priors[[mod]][8],
                                                         prior   = c(as.numeric(priors[[mod]][6]),
                                                                     as.numeric(priors[[mod]][7]))
                 )
                 priors.param[,ncurv+2] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+2],
                                                         distrib = priors[[mod]][13],
                                                         prior   = c(as.numeric(priors[[mod]][11]),
                                                                     as.numeric(priors[[mod]][12]))
                 )
                 for (ii in (ncurv+3):(ncurv*2+2)){
                    priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                      distrib = priors[[mod]][18],
                                                      prior   = c(as.numeric(priors[[mod]][16]),
                                                                  as.numeric(priors[[mod]][17]))
                   )
                 }  
                 priors.param[, (npar-1)] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,(npar-1)],
                                                           distrib = "'Uniform'",
                                                           prior   = c(0, 100))
                 
                 priors.param[,npar] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,npar],
                                                      distrib = "'Uniform'",
                                                      prior   = c(0, 100))
             
             
             
           } else if (model.names[[mod]] =="2expWithAsympt_bis"){
           #####################################################################################
                 parameters = c("a1(k)", "b1", "a2(k)", "b2", "a3(k)", "g1", "g2")
                 par.var = c("a1", "a2", "a3")
                 priors.param = matrix(NA, nrow = length(data.MCMC.cooked[,1]), ncol = npar)
                 for (ii in 1:ncurv){
                   priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                      distrib = as.character(priors[[mod]][3]),
                                                      prior   = c(as.numeric(priors[[mod]][1]), 
                                                                  as.numeric(priors[[mod]][2]))
                   )
                 }
                 priors.param[,ncurv+1] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                         distrib = priors[[mod]][8],
                                                         prior   = c(as.numeric(priors[[mod]][6]),
                                                                     as.numeric(priors[[mod]][7]))
                 )
                 for (ii in (ncurv+2):(ncurv*2+1)){
                   priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                      distrib = as.character(priors[[mod]][13]),
                                                      prior   = c(as.numeric(priors[[mod]][11]), 
                                                                  as.numeric(priors[[mod]][12]))
                   )
                 }
                 priors.param[,2*ncurv+2] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+2],
                                                         distrib = priors[[mod]][18],
                                                         prior   = c(as.numeric(priors[[mod]][16]),
                                                                     as.numeric(priors[[mod]][17]))
                 )
                 for (ii in (2*ncurv+3):(ncurv*3+2)){
                   priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                      distrib = priors[[mod]][23],
                                                      prior   = c(as.numeric(priors[[mod]][21]),
                                                                  as.numeric(priors[[mod]][22]))
                   )
                 }  
                 priors.param[, (npar-1)] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,(npar-1)],
                                                           distrib = "'Uniform'",
                                                           prior   = c(0, 100))
                 
                 priors.param[,npar] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,npar],
                                                      distrib = "'Uniform'",
                                                      prior   = c(0, 100))
           
                 
                 
                 
           } else if ((model.names[[mod]] =="2expWithAsympt")|(model.names[[mod]] =="2expWithAsympt_rel")){
           #####################################################################################
             parameters = c("a1(k)", "b1", "a2", "b2", "a3(k)", "g1", "g2")
             par.var = c("a1", "a3")
             priors.param = matrix(NA, nrow = length(data.MCMC.cooked[,1]), ncol = npar)
             for (ii in 1:ncurv){
               priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                  distrib = as.character(priors[[mod]][3]),
                                                  prior   = c(as.numeric(priors[[mod]][1]), 
                                                              as.numeric(priors[[mod]][2]))
               )
             }
             priors.param[,ncurv+1] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][8],
                                                     prior   = c(as.numeric(priors[[mod]][6]),
                                                                 as.numeric(priors[[mod]][7]))
             )
             priors.param[,ncurv+2] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+2],
                                                     distrib = priors[[mod]][13],
                                                     prior   = c(as.numeric(priors[[mod]][11]),
                                                                 as.numeric(priors[[mod]][12]))
             )
             priors.param[,ncurv+3] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+3],
                                                     distrib = priors[[mod]][18],
                                                     prior   = c(as.numeric(priors[[mod]][16]),
                                                                 as.numeric(priors[[mod]][17]))
             )
             for (ii in (ncurv+4):(ncurv*2+3)){
               priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                  distrib = priors[[mod]][23],
                                                  prior   = c(as.numeric(priors[[mod]][21]),
                                                              as.numeric(priors[[mod]][22]))
               )
             }  
             priors.param[, (npar-1)] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,(npar-1)],
                                                       distrib = "'Uniform'",
                                                       prior   = c(0, 100))
             
             priors.param[,npar] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,npar],
                                                  distrib = "'Uniform'",
                                                  prior   = c(0, 100))
             
          
             
             
             
           } else if (model.names[[mod]] =="3expWithAsympt"){
           #####################################################################################
             parameters = c("a1(k)", "b1", "a2", "b2", "a3", "b3", "a4(k)" , "g1", "g2")
             par.var = c("a1", "a4")
             priors.param = matrix(NA, nrow = length(data.MCMC.cooked[,1]), ncol = npar)
             for (ii in 1:ncurv){
               priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                  distrib = as.character(priors[[mod]][3]),
                                                  prior   = c(as.numeric(priors[[mod]][1]), 
                                                              as.numeric(priors[[mod]][2]))
               )
             }
             priors.param[,ncurv+1] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][8],
                                                     prior   = c(as.numeric(priors[[mod]][6]),
                                                                 as.numeric(priors[[mod]][7]))
             )
             priors.param[,ncurv+2] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][13],
                                                     prior   = c(as.numeric(priors[[mod]][11]),
                                                                 as.numeric(priors[[mod]][12]))
             )
             priors.param[,ncurv+3] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][18],
                                                     prior   = c(as.numeric(priors[[mod]][16]),
                                                                 as.numeric(priors[[mod]][17]))
             )
             priors.param[,ncurv+4] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][23],
                                                     prior   = c(as.numeric(priors[[mod]][21]),
                                                                 as.numeric(priors[[mod]][22]))
             )
             priors.param[,ncurv+5] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][28],
                                                     prior   = c(as.numeric(priors[[mod]][26]),
                                                                 as.numeric(priors[[mod]][27]))
             )
             for (ii in (ncurv+6):(ncurv*2+5)){
               priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                  distrib = priors[[mod]][33],
                                                  prior   = c(as.numeric(priors[[mod]][31]),
                                                              as.numeric(priors[[mod]][32]))
               )
             }  
             priors.param[, (npar-1)] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,(npar-1)],
                                                       distrib = "'Uniform'",
                                                       prior   = c(0, 100))
             
             priors.param[,npar] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,npar],
                                                  distrib = "'Uniform'",
                                                  prior   = c(0, 100))
             
             
             
             
           } else if (model.names[[mod]] =="3expWithAsympt_bis"){
             #####################################################################################
             parameters = c("a1(k)", "b1", "a2(k)", "b2", "a3", "b3", "a4(k)" , "g1", "g2")
             par.var = c("a1", "a2", "a4")
             priors.param = matrix(NA, nrow = length(data.MCMC.cooked[,1]), ncol = npar)
             for (ii in 1:ncurv){
               priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                  distrib = as.character(priors[[mod]][3]),
                                                  prior   = c(as.numeric(priors[[mod]][1]), 
                                                              as.numeric(priors[[mod]][2]))
               )
             }
             priors.param[,ncurv+1] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][8],
                                                     prior   = c(as.numeric(priors[[mod]][6]),
                                                                 as.numeric(priors[[mod]][7]))
             )
             for (ii in (ncurv+2):(2*ncurv+1)){
               priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                  distrib = as.character(priors[[mod]][13]),
                                                  prior   = c(as.numeric(priors[[mod]][11]), 
                                                              as.numeric(priors[[mod]][12]))
               )
             }
             priors.param[,ncurv*2+2] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][18],
                                                     prior   = c(as.numeric(priors[[mod]][16]),
                                                                 as.numeric(priors[[mod]][17]))
             )
             priors.param[,ncurv*2+3] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][23],
                                                     prior   = c(as.numeric(priors[[mod]][21]),
                                                                 as.numeric(priors[[mod]][22]))
             )
             priors.param[,ncurv*2+4] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ncurv+1],
                                                     distrib = priors[[mod]][28],
                                                     prior   = c(as.numeric(priors[[mod]][26]),
                                                                 as.numeric(priors[[mod]][27]))
             )
             for (ii in (ncurv*2+5):(ncurv*3+4)){
               priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                                  distrib = priors[[mod]][33],
                                                  prior   = c(as.numeric(priors[[mod]][31]),
                                                              as.numeric(priors[[mod]][32]))
               )
             }  
             priors.param[, (npar-1)] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,(npar-1)],
                                                       distrib = "'Uniform'",
                                                       prior   = c(0, 100))
             
             priors.param[,npar] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,npar],
                                                  distrib = "'Uniform'",
                                                  prior   = c(0, 100))
             
             
             
        
          
       } else if ((model.names[[mod]] =="Coutagne_bis")|(model.names[[mod]] =="expexp_bis")|(model.names[[mod]] =="hyperb_bis")){
         #####################################################################################
         parameters = c("a1(k)", "b1(k)", "n1", "a2(k)", "g1", "g2")
         par.var = c("a1", "b1","a2")
         priors.param = matrix(NA, nrow = length(data.MCMC.cooked[,1]), ncol = npar)
         for (ii in 1:ncurv){
           priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                              distrib = as.character(priors[[mod]][3]),
                                              prior   = c(as.numeric(priors[[mod]][1]), 
                                                          as.numeric(priors[[mod]][2]))
           )
         }
         for (ii in (ncurv+1):(2*ncurv)){
           priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                              distrib = as.character(priors[[mod]][8]),
                                              prior   = c(as.numeric(priors[[mod]][6]), 
                                                          as.numeric(priors[[mod]][7]))
           )
         }
         priors.param[,2*ncurv+1] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,2*ncurv+1],
                                                 distrib = priors[[mod]][13],
                                                 prior   = c(as.numeric(priors[[mod]][11]),
                                                             as.numeric(priors[[mod]][12]))
         )
         for (ii in (2*ncurv+2):(3*ncurv + 1)){
           priors.param[,ii] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,ii],
                                              distrib = priors[[mod]][18],
                                              prior   = c(as.numeric(priors[[mod]][16]),
                                                          as.numeric(priors[[mod]][17]))
           )
         }  
         priors.param[, (npar-1)] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,(npar-1)],
                                                   distrib = "'Uniform'",
                                                   prior   = c(0, 100))
         
         priors.param[,npar] = prior.mcmc.bam(mcmc    = data.MCMC.cooked[,npar],
                                              distrib = "'Uniform'",
                                              prior   = c(0, 100))
       }
           
           # computation of the loglikelihood for the computation of DIC:
           logprior = 0
           for (ll in 1:len.mcmc){
             logprior[ll] = sum(priors.param[ll,])   
           }
           loglikelihood = logpost - logprior
           varLogLikelihood = var(loglikelihood)
           MeanDev = -2*mean(loglikelihood)
           ###
           # df.mcmc.LL[[mod]][[chi.i]] = data.frame(loglikelihood = loglikelihood, logprior=logprior, logposterior= logpost)
           
           DIC.rec[[chi.i]][[mod]] = MeanDev + 2*varLogLikelihood       # Gelman 2004 "Bayesian data analysis"
           print(paste0("DIC=", DIC.rec[[chi.i]][[mod]]))
           #############################################################
  
           
           # Accuracy:
           print(paste0("compute Accuracy ********************"))
           dir.rec.segm.test.param           = paste0(dir.segm.recessions.model.param, "/", tail(par.var,1))
           mu.results.df.rec[[mod]]          = read.table(file=paste0(dir.rec.segm.test.param, "/mu.results.df.txt"), header=TRUE)
           gamma.results.df.rec[[mod]]       = read.table(file=paste0(dir.rec.segm.test.param, "/gamma.results.df.txt"), header=TRUE)
           Q2.mu.rec[[mod]]                  = read.table(file=paste0(dir.rec.segm.test.param , "/Q2.mu.txt"), header=TRUE)
           mu.res.rec[[mod]]                 = read.table(file=paste0(dir.rec.segm.test.param , "/mu.res.txt"), header=TRUE)
           Q97.mu.rec[[mod]]                 = read.table(file=paste0(dir.rec.segm.test.param , "/Q97.mu.txt"), header=TRUE)
           Data.segm.rec[[mod]]              = read.table(file=paste0(dir.rec.segm.test.param , "/Data.segm.rec.txt"), header=TRUE)
           mcmc.seg.rec[[mod]]               = read.table(file=paste0(dir.rec.segm.test.param, "/mcmc_segmentation.txt"), header=TRUE)
           gamma_segm_recess[[mod]]          = c(mean(mcmc.seg.rec[[mod]]$Y1_gamma1), std(mcmc.seg.rec[[mod]]$Y1_gamma1))
           
           if (length(mu.results.df.rec[[mod]]$mu.MAP) >1) {
              tau.results.df.rec[[mod]]      = read.table(file=paste0(dir.rec.segm.test.param, "/tau.results.df.txt"), header=TRUE)
              df.shift.times.rec[[mod]]      = read.table(file=paste0(dir.rec.segm.test.param, "/df.shift.times.txt"), header=TRUE)
              df.shift.times.plus.rec[[mod]] = read.table(file=paste0(dir.rec.segm.test.param, "/df.shift.times.plus.txt"), header=TRUE)
              ts.res.before.rec[[mod]]       = read.table(file=paste0(dir.rec.segm.test.param, "/ts.res.before.txt"), header=TRUE)
              ts.res.rec[[mod]]              = read.table(file=paste0(dir.rec.segm.test.param, "/ts.res.txt"), header=TRUE)
              ts.res.plus.rec[[mod]]         = read.table(file=paste0(dir.rec.segm.test.param, "/ts.res.plus.txt"), header=TRUE)
              Q2.ts.rec[[mod]]               = read.table(file=paste0(dir.rec.segm.test.param, "/Q2.ts.txt"), header=TRUE)
              Q97.ts.rec[[mod]]              = read.table(file=paste0(dir.rec.segm.test.param, "/Q97.ts.txt"), header=TRUE)
              nS.ok.rec[[mod]]               = nrow(df.shift.times.rec[[mod]])+1
              pdf.ts.rec[[mod]]              = mcmc.seg.rec[[mod]][, (nS.ok.rec[[mod]]+1):(2*nS.ok.rec[[mod]]-1)]
              shift.times.recessions[[chi.i]][[mod]]  = read.table(file =paste0(dir.rec.segm.test.param,"/df.shift.times.txt"), header = TRUE)
              data.annotate.recess[[mod]]    = data.frame(t     = shift.times.recessions[[chi.i]][[mod]]$ts.res,
                                                          tstart= shift.times.recessions[[chi.i]][[mod]]$Q2.ts,
                                                          tend  = shift.times.recessions[[chi.i]][[mod]]$Q97.ts)
              #data.annotate.recess.adjust[[mod]][[param]] = data.frame(t.adj = shift.times.recessions[[mod]][[param]]$ts.real)
              data.annotate.recess.adjust[[mod]] = data.frame(t.adj = shift.times.recessions[[chi.i]][[mod]]$ts.res)
           } else {
              nS.ok.rec[[mod]]                   = 1
              shift.times.recessions[[chi.i]][[mod]]  = list()
              pdf.ts.rec[[mod]]                  = list()
              data.annotate.recess[[mod]]        = list()
              data.annotate.recess.adjust[[mod]] = list()
              tau.results.df.rec[[mod]] = list()
           } 
           shift.time.true  = bt.from.gaugingsss$t.shift.for.b$treal
           shift.time.found = shift.times.recessions[[chi.i]][[mod]]$ts.res
           nshift.true      = length(shift.time.true)
           nshift.found     = length(shift.time.found) 
           Nrecess          = length(Data.segm.rec[[mod]]$t)
           
           if (nshift.found > 0) {
              Accuracy.rec[[chi.i]][[mod]] <- accuracy.rec(time         =  Data.segm.rec[[mod]]$t, 
                                                           ts           =  shift.time.true , 
                                                           nshift.found =  nshift.found,
                                                           nshift.true  =  nshift.true,
                                                           t2           =  shift.times.recessions[[chi.i]][[mod]]$Q2.ts, 
                                                           t97          =  shift.times.recessions[[chi.i]][[mod]]$Q97.ts,
                                                           tMAP         =  shift.times.recessions[[chi.i]][[mod]]$ts.res,
                                                           dir          =  dir.rec.segm.test.param)
                                                       
              shift.rootmse.rec[[chi.i]][[mod]] = rmserr(x=Accuracy.rec[[chi.i]][[mod]]$tMAP.TP,
                                                         y=Accuracy.rec[[chi.i]][[mod]]$tTRUE.TP, 
                                                         summary = TRUE)
              Data.segm.rec[[mod]] = cbind(Data.segm.rec[[mod]], 
                                           b1.gaug = rep(0, length(Data.segm.rec[[mod]]$t)),
                                           b2.gaug = rep(0, length(Data.segm.rec[[mod]]$t))
                                           )
              
              for (bbb in 1:length(Data.segm.rec[[mod]]$t)) {
                if (Data.segm.rec[[mod]]$t[bbb]       < bt.from.gaugingsss$t.shift.for.b$treal[1]) {
                    Data.segm.rec[[mod]]$b1.gaug[bbb] = bt.from.gaugingsss$bt1.df$maxpost[1]
                    Data.segm.rec[[mod]]$b2.gaug[bbb] = bt.from.gaugingsss$bt2.df$maxpost[1]
                }
                for (ccc in 2:length(bt.from.gaugingsss$t.shift.for.b$treal)) {
                  if ((Data.segm.rec[[mod]]$t[bbb]   <  bt.from.gaugingsss$t.shift.for.b$treal[ccc]) &
                      (Data.segm.rec[[mod]]$t[bbb]   >= bt.from.gaugingsss$t.shift.for.b$treal[ccc-1])) {
                    Data.segm.rec[[mod]]$b1.gaug[bbb] = bt.from.gaugingsss$bt1.df$maxpost[ccc]
                    Data.segm.rec[[mod]]$b2.gaug[bbb] = bt.from.gaugingsss$bt2.df$maxpost[ccc]
                  }
                }
                if ((Data.segm.rec[[mod]]$t[bbb]   >= tail(bt.from.gaugingsss$t.shift.for.b$treal, 1))) {
                  Data.segm.rec[[mod]]$b1.gaug[bbb] = tail(bt.from.gaugingsss$bt1.df$maxpost,1)
                  Data.segm.rec[[mod]]$b2.gaug[bbb] = tail(bt.from.gaugingsss$bt2.df$maxpost,1)
                }
              }
              RMSE.b1.rec[[chi.i]][[mod]] = rmserr(x       = Data.segm.rec[[mod]]$b1.gaug,
                                                   y       = Data.segm.rec[[mod]]$maxpost/100, 
                                                   summary = TRUE)
              RMSE.b2.rec[[chi.i]][[mod]] = rmserr(x       = Data.segm.rec[[mod]]$b2.gaug,
                                                   y       = Data.segm.rec[[mod]]$maxpost/100, 
                                                   summary = TRUE)
              accu.rec[[chi.i]][[mod]] =  Accuracy.rec[[chi.i]][[mod]]$accuracy; 
              shift.rmse.rec[[chi.i]][[mod]] = shift.rootmse.rec[[chi.i]][[mod]]$rmse;
              rmse.b1.rec[[chi.i]][[mod]] =  RMSE.b1.rec[[chi.i]][[mod]]$rmse;
              rmse.b2.rec[[chi.i]][[mod]] =  RMSE.b2.rec[[chi.i]][[mod]]$rmse;
              
           } else {
              nshift.found = 0
              Accuracy.rec[[chi.i]][[mod]] = NULL
              RMSE.b1.rec[[chi.i]][[mod]]  = NULL
              RMSE.b2.rec[[chi.i]][[mod]]  = NULL
              shift.rootmse.rec[[chi.i]][[mod]] =NULL
              shift.rmse.rec[[chi.i]][[mod]] = NULL;
              accu.rec[[chi.i]][[mod]] = (length(Data.segm.rec[[mod]]$t) - nshift.true)/length(Data.segm.rec[[mod]]$t)
              rmse.b1.rec[[chi.i]][[mod]]  = 0
              rmse.b1.rec[[chi.i]][[mod]]  = 0
            }
           

           ##########################################################################################
           # PLot of the pdf of the shifts:
           #preparing datasets of pdf of shift times:
           if (!is.null(tau.results.df.rec[[mod]]$tau.MAP)) {
             X1[[chi.i]][[mod]]  =  pdf.ts.rec[[mod]]    
             if (length(tau.results.df.rec[[mod]]$tau.MAP) ==1) {
                 X1.new[[chi.i]][[mod]]  = X1[[chi.i]][[mod]] 
                 df.ts[[chi.i]][[mod]] = data.frame(xx = unlist(X1.new[[chi.i]][[mod]]))
                 
                 t.plot[[mod]]  = t.plot[[mod]]  + 
                                  geom_density(data = df.ts[[chi.i]][[mod]],
                                               aes(x= xx ,
                                               ..scaled.. ),
                                               fill = colllor[chi.i], 
                                               colour=type.plot[chi.i], 
                                               alpha=0.4)
               
             } else {
                 X1.new[[chi.i]][[mod]]  = data.frame(time= X1[[chi.i]][[mod]][,1], 
                                                    ord = rep(1, length(X1[[chi.i]][[mod]][,1])))
                 for (orderr in 1:length(tau.results.df.rec[[mod]]$tau.MAP)) {
                       X1.new[[chi.i]][[mod]]  = rbind(X1.new[[chi.i]][[mod]] , 
                                                 data.frame(time= X1[[chi.i]][[mod]][,orderr], 
                                                            ord = rep(orderr, length(X1[[chi.i]][[mod]][,1]))))
                 }
                 df.ts[[chi.i]][[mod]] = data.frame(xx = unlist(X1.new[[chi.i]][[mod]]$time), 
                                                    group = unlist(X1.new[[chi.i]][[mod]]$ord)
                                                    )
                 t.plot[[mod]]  = t.plot[[mod]]  +
                                  geom_density(data = df.ts[[chi.i]][[mod]], 
                                  aes(x =xx, 
                                      ..scaled.. , 
                                      group = group),  
                                      fill = colllor[chi.i], 
                                      colour=type.plot[chi.i], alpha=0.4)
                 
             }
           }
           
           # Quantiles for tot uncertainty for each parameter:
           Ysim[[mod]] = list();
           data.tmp.2[[mod]]  = list()
           for (jj in 1:nS.ok.rec[[mod]]) {
             Ysim[[mod]][[jj]] =0
             for (ii in 1: length(mcmc.seg.rec[[mod]][,jj])) {
               Ysim[[mod]][[jj]][ii] = mcmc.seg.rec[[mod]][ii,jj] + 
                        rnorm(1,  mean = 0, sd = mcmc.seg.rec[[mod]][, 2*nS.ok.rec[[mod]]])
             }
             data.tmp.2[[mod]][[jj]] = quantile(Ysim[[mod]][[jj]], probs=c(0.025,0.975) ) #, na.rm=TRUE)
           } 
           
       }  # end loop in models!
    
       
       ##################
       # beta comparison:
       ##################
       best.models = c(1, 2, 3, 4, 5, 9)
       color.models = c("red", "blue", "green", "violet", "cyan", "brown", "green", "lightblue", "orange", "yellow", "gray")
       bt.from.gaugingsss.modified = list(bt1.df       = bt.from.gaugingsss$bt1.df[-6,], 
                                          t.shift.for.b = bt.from.gaugingsss$t.shift.for.b[-6,]) 
       
       t.b.comparison[[chi.i]] <- ggplot()+
         scale_x_continuous(expand = c(0,0), 
                            limits =time.limits) +
         scale_y_continuous(expand = c(0,0), 
                            limits = c(-1, 0.5), 
                            breaks=seq(-1, 0.5, 0.5)) +
         ylab(TeX("$\\hat{b_1} $")) +
         xlab("Time [day]") +
         coord_cartesian(clip = 'off') +
         #geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray60", size = 0.2) + # limni
         theme_light(base_size=20)+
         theme(axis.text=element_text(size=20)
               ,axis.title=element_text(size=30, face="bold")
               ,panel.grid.major=element_blank()
               ,panel.grid.minor=element_blank()
               ,legend.text=element_text(size=20)
               ,legend.title=element_text(size=30)
               ,legend.key.size=unit(1.5, "cm")
               ,legend.position="none"
               ,plot.margin=unit(c(2, 0.5, 4.5, 1),"cm")
               ,axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=35, face ="bold")) +
         geom_segment(mapping= aes(x = c(time.limits[1], 
                                         bt.from.gaugingsss.modified$t.shift.for.b$treal), 
                                   y = bt.from.gaugingsss.modified$bt1.df$maxpost, 
                                   xend = c(bt.from.gaugingsss.modified$t.shift.for.b$treal, 
                                            time.limits[2]),  
                                   yend = bt.from.gaugingsss.modified$bt1.df$maxpost), 
                      color = "black", linetype ="solid",
                      size = 3)+
         annotate("text",
                  x = 780,
                  y = bt.from.gaugingsss.modified$bt1.df$maxpost[1] + 0.002, 
                  label = 'atop(bold("From gaugings"))', 
                  color = "black", parse=TRUE, fontface= 2,
                  size=6)
       
       
       for (mmmmm in best.models){ #length(model.names)) {
         if ( nS.ok.rec[[mmmmm]] > 1) {
           for (jj in 1:nS.ok.rec[[mmmmm]]){
             t.b.comparison[[chi.i]]  = t.b.comparison[[chi.i]]   
               # annotate("rect", 
               #          xmin = unlist(ts.res.before.rec[[mmmmm]])[[jj]], 
               #          xmax = unlist(ts.res.plus.rec[[mmmmm]])[[jj]] , 
               #          ymin = data.tmp.2[[mmmmm]][[jj]][1]/100,
               #          ymax = data.tmp.2[[mmmmm]][[jj]][2]/100,
               #          fill = "red", alpha = 0.2)
           }
           df.segment.beta = data.frame(x        = unlist(ts.res.before.rec[[mmmmm]]), 
                                        y        = unlist(mu.res.rec[[mmmmm]])/100, 
                                        xend     = unlist(ts.res.plus.rec[[mmmmm]]), 
                                        yend     = unlist(mu.res.rec[[mmmmm]])/100)
           
           t.b.comparison[[chi.i]]  =  t.b.comparison[[chi.i]] + 
                                       geom_segment(data = df.segment.beta,
                                                    mapping=aes(x        = x, 
                                                                y        = y, 
                                                                xend     = xend, 
                                                                yend     = yend),
                                                    col      = color.models[mmmmm], 
                                                    linetype = "11", #linetype.models[mmmmm],
                                                    size = 2)+
           annotate("text", 
                    x = df.segment.beta$x[1] + 110 ,
                    y = df.segment.beta$y[1] + 0.03, 
                    label = paste0("M", mmmmm), 
                    color = color.models[mmmmm],  fontface= 2,
                    size=6)
         } else {
           t.b.comparison[[chi.i]]  = t.b.comparison[[chi.i]]  
           #                            annotate("rect",  
           #                                              xmin = Data.segm.rec[[mmmmm]][[length(nS.ok.rec[[mmmmm]])]]$t[1], 
           #                                              xmax = tail(Data.segm.rec[[mmmmm]][[length(nS.ok.rec[[mmmmm]])]]$t,1),
           #                                              ymin = data.tmp.2[[length(nS.ok.rec[[mmmmm]])]][[1]][1]/100,
           #                                              ymax = data.tmp.2[[length(nS.ok.rec[[mmmmm]])]][[1]][2]/100,
           #                                              fill = "red", alpha=0.2)
         }
       }

       ########################
       # Delta beta comparison:
       ########################
       deltab.gaug=c(); deltab.mod = NULL; deltab.mod.df =NULL;
       shape.models = c(21, 22, 23, 24, 25, 21, 22, 23, 7, 25, 24);
       size.models =  rep(8, 11)
       #best.models = c(1, 2, 3, 4, 5, 9)
       color.mod.db = c("black","red", "blue", "green", "violet", "cyan", "orange")
       for (iii in 2:length(bt.from.gaugingsss.modified$bt1.df$maxpost)) {
         deltab.gaug = c(deltab.gaug, 
                         bt.from.gaugingsss.modified$bt1.df$maxpost[iii] - 
                         bt.from.gaugingsss.modified$bt1.df$maxpost[iii-1])
       }
       deltab.gaug.df = data.frame(deltab =deltab.gaug, 
                                   tshift = bt.from.gaugingsss.modified$t.shift.for.b$treal,
                                   color.mod = "black",
                                   shape.mod = 15, size.mod = 10)
       deltab.tot.df =  deltab.gaug.df
       for (mmmmm in best.models){ #length(model.names)) {
         if ( nS.ok.rec[[mmmmm]] > 1) {
           deltab.mod[[mmmmm]] = list()
           for (jj in 2:nS.ok.rec[[mmmmm]]){
             deltab.mod[[mmmmm]][[jj-1]]  =   
               unlist(mu.res.rec[[mmmmm]])[[jj]]/100 - 
               unlist(mu.res.rec[[mmmmm]])[[jj-1]]/100 
           }
           deltab.mod.df[[mmmmm]] = data.frame(deltab =  unlist(deltab.mod[[mmmmm]]), 
                                               tshift =  unlist(ts.res.rec[[mmmmm]]), 
                                               color.mod = color.models[mmmmm],
                                               shape.mod = shape.models[mmmmm],
                                               size.mod  = size.models[mmmmm]
                                               )
           deltab.tot.df = rbind(deltab.tot.df, deltab.mod.df[[mmmmm]] )
         }
       }

       ######
       t.deltab.comparison[[chi.i]] <- ggplot(data = deltab.tot.df)+
         scale_x_continuous(expand = c(0,0), 
                            limits = time.limits) +
         scale_y_continuous(expand = c(0,0), 
                            limits = c(-0.6, 0.1), 
                            breaks = c(seq(-0.6, 0, 0.1), 0.1)) +
         geom_hline(aes(yintercept = 0), linetype ="dashed", color="gray")+
         ylab(TeX("$\\Delta \\hat{b_1} $")) +
         xlab("Time [day]") +
         coord_cartesian(clip = 'off') +
         #geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "gray60", size = 0.2) + # limni
         theme_light(base_size=20)+
         theme(axis.text=element_text(size=20)
               ,axis.title=element_text(size=30, face="bold")
               ,panel.grid.major=element_blank()
               ,panel.grid.minor=element_blank()
               ,legend.text=element_text(size=30)
               ,legend.key.size=unit(1.5, "cm")
               ,legend.position="bottom"
               ,legend.direction= "horizontal"
               ,plot.margin=unit(c(2, 0.5, 0, 1),"cm")
               ,legend.title = element_blank()
               ,axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=35, 
                                            face ="bold")) +
         geom_point(
                    mapping = aes(y = deltab, 
                                  x= tshift,
                                  fill = color.mod),
                    shape = deltab.tot.df$shape.mod, 
                    size = deltab.tot.df$size.mod,
                    color = deltab.tot.df$color.mod) +
         scale_fill_manual(values =c("black","red", "blue", "green", "violet", "cyan", "orange"),
                           labels =c("From gaugings", "M1",  "M2" , "M3"  ,"M4" , "M5"  ,"M9"),
                           guide = guide_legend(override.aes = list(
                             linetype = c("blank","blank","blank","blank","blank","blank","blank"),
                             shape = c(15, 21, 22, 23, 24, 25, 7),
                             #size = c(8,6,6,6,6,6,6),
                             size = c(10,8,8,8,8,8,8),
                             color = c("black","red", "blue", "green", "violet", "cyan", "orange"))))
         #geom_vline(data= deltab.gaug.df, aes(xintercept = tshift), linetype="dashed", color="gray")

       # t.deltab.comparison[[chi.i]]  = t.deltab.comparison[[chi.i]] +
       # scale_colour_manual(name = element_blank(), 
       #                     values = colss.mod,
       #                     breaks = c("M1","M2","M3","M4","M5","M9", "colgaug"),
       #                     labels = c("M1","M2","M3","M4","M5","M9", "From Gaug"))
                           # guide = guide_legend(override.aes = list(
                           #                              #linetype = c("blank","blank","blank","blank","blank","blank"),
                           #                              shape = c(21 ,22 ,23 ,24 ,25  ,7))))
       # plot chi30:
       
       t.b.comparison.tot = plot_grid( t.b.comparison[[chi.i]] ,  t.deltab.comparison[[chi.i]],
                                       labels = c("a)", "b)"), label_size = 35, label_fontface = "bold",
                                       ncol =2, nrow = 1, rel_widths = c(1, 1))
       pdf(paste0(dir.regression,"/deltab_b_comparison.pdf"), 20, 12 ,useDingbats=F)
       print(t.b.comparison.tot)
       dev.off()
  } # end loop in chi values!
  #################################################################################### end of iterations

  
  
  
  
  
  
  
  
  
  
  
  
  X.gaug  =  pdf.ts.gaugings    
    if (ncol( pdf.ts.gaugings)==1) {
      X.gaug.new = X.gaug
    } else {
      X.gaug.new = data.frame(time= X.gaug[,1], 
                              ord = rep(1, length(X.gaug[,1])))
      for (orderr in 1:(ncol( pdf.ts.gaugings))) {
        X.gaug.new  = rbind(X.gaug.new, 
                            data.frame(time = X.gaug[,orderr], 
                                       ord = rep(orderr, length(X.gaug[,1]))))
      }
    }
    t.plot.gaug  = ggplot()  + 
      scale_y_continuous(name = "", expand = c(0,0), limits = c(0,1)) + 
      ggtitle("From gaugings") + 
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+  
      #theme_set(theme_gray(base_size = 20))+
      #theme_light(base_size = 20) +
      theme_classic(base_size = 20) +
      theme(
            axis.line.y = element_blank()
            ,axis.line.x = element_blank()
            ,axis.title=element_text(size= 20,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "black", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,panel.background=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0, 1, 0, 0.5),"cm")
            # ,axis.text.x =element_blank()
            ,axis.text = element_blank()
            ,axis.ticks = element_blank())+
      scale_x_continuous(name=element_blank(), expand = c(0,0), limits =time.limits)+
      geom_density(aes(x = X.gaug.new$time, 
                       ..scaled.. , 
                       group = X.gaug.new$ord),  
                   fill = unlist("black"), 
                   colour=NA, alpha=0.4)

    #----------------------
    #PLOT official dates
    #----------------------
    t.plot.off <- ggplot() +
      scale_y_continuous(name = "", expand = c(0,0), limits = c(0,1))+
      ggtitle("Official dates of RC update")+
      coord_cartesian(clip = 'off')+
      #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
      #theme_set(theme_grey(base_size = 20))+
      theme_classic(base_size = 20) +
      theme( axis.text=element_text(size= 15)
            ,axis.line.y = element_blank()
            ,axis.line.x = element_line(color="gray10", size= 0.5)
            ,axis.title=element_text(size= 20,face="bold")
            ,plot.title = element_text(face = "bold", size = 15, color= "black", margin=margin(0,0,0,0))
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            #,panel.background=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin=unit(c(0, 1, 0, 0.5),"cm")
            # ,axis.text.x =element_blank()
            ,axis.text.y = element_blank()
            ,axis.ticks.y = element_blank()
            ,axis.ticks.x = element_line(color="gray10", size= 0.5)) +
      scale_x_continuous(name="Time [days]", expand = c(0,0), 
                         limits =time.limits)+
      geom_point(data = data.annotate.off, aes(x = xeffect, y = rep(0, length(data.annotate.off$xeffect))), 
                                               color= "black", size = 4, shape =4, stroke=2 )
    

  ts.shift.total.plot = plot_grid( t.plot[[1]] ,  t.plot[[2]] , 
                                   t.plot[[3]] ,  t.plot[[4]] , 
                                   t.plot[[5]] ,  t.plot[[6]] , 
                                   t.plot[[7]] ,  t.plot[[8]] , 
                                   t.plot[[9]] ,  t.plot[[10]] , 
                                   t.plot[[11]] , 
                                   #t.plot[[12]] , 
                                   t.plot.gaug,
                                   t.plot.off,
                                   ncol =1, nrow = 14)
  
  pdf(paste0(dir.regression,"/tshift_comparison.pdf"), 10, 10 ,useDingbats=F)
  print(ts.shift.total.plot)
  dev.off()
  

  
  
  
  
  
  
  
  ##############################################################################################################
  # tshift_comparison
  ###################
  ddd=NULL;   gglegend = NULL;  
  ts.shift.total.plot2=NULL; ts.shift.total.plot3 =NULL; 
  color.m = c("red", "blue", "green"); 
  typelin = c("solid", "solid", "solid");
  chi.1 <- c("chi = 10 cm");
  chi.2 = c("chi = 30 cm"); 
  chi.3 = c("chi = 50 cm")
  chi2 <- c(TeX("$\\chi = 10 \\; cm$"), TeX("$\\chi = 30 \\; cm$"), TeX("$\\chi = 50 \\; cm$"))
  fromgaug <- "Rating shift times obtained from gaugings"
  off.times <- "Official dates of RC update"
  colss <- c("chi = 10 cm" = "red", 
             "chi = 30 cm" = "blue",
             "chi = 50 cm" = "green",
             "Rating shift times obtained from gaugings" ="black",
             "Official dates of RC update" = "black")
  
  oders.models.chi10 = c(1,2,4,5,11,3,9,6,8,10,7)   #put the models in the plot inorder of their performance
  
  # plot chi10:
    ts.shift.total.plot2[[1]]   = ggplot() +
      scale_y_continuous(name = "Different stage recession models", 
                         expand = c(0,0), 
                         limits = c(0,length(model.titles)+2),  
                         breaks =seq(0.5, length(model.titles) + 1.5, 1), 
                         label = c("Official dates", "From Gaugings" , model.titles[oders.models.chi10])) +
      geom_hline(yintercept = seq(1, length(model.titles)+1, 1), color="gray70", size= 0.2, linetype ="dashed") +
      #geom_hline(yintercept = 2, color="black", size= 0.2) +
      coord_cartesian(clip = 'off')+
      theme_bw(base_size = 20) +
      theme(    axis.text=element_text(size= 20)
                ,axis.text.y=element_text(size= 20, face ="bold")
                ,axis.title=element_text(size= 30,face="bold")
                ,panel.grid.major=element_blank()
                ,panel.grid.minor=element_blank()
                #,panel.background=element_blank()
                ,legend.text=element_text(size=30)
                ,legend.title=element_text(size=30)
                ,legend.key.size=unit(1.5, "cm")
                ,legend.direction = "horizontal"
                ,legend.position ="bottom"
                ,axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))
                ,plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
      scale_x_continuous(name="Time [days]", expand = c(0,0), 
                         limits =time.limits) + 
      geom_point(data = data.annotate.off, 
                 aes(x = xeffect, y = rep(0, length(data.annotate.off$xeffect))), 
                 color= "black", size = 3, shape =4, stroke= 3) +
      geom_segment(aes(x = bt.from.gaugingsss$t.shift.for.b$treal, 
                       xend = bt.from.gaugingsss$t.shift.for.b$treal, 
                       y =1, yend =2 ), color= "black", size= 1.1) +
      geom_vline(aes(xintercept=-1000, color = chi.1),
                 show.legend = T, size= 1.5, linetype = "solid")
    
    ddd[[1]]  = list()  
    for (mmm in 1:length(model.names)){
      if (!is.null(unlist(shift.times.recessions[[1]][[oders.models.chi10[mmm]]]$ts.res))) {
        ddd[[1]][[mmm]] = data.frame(xx =unlist(shift.times.recessions[[1]][[oders.models.chi10[mmm]]]$ts.res), 
                                     yy = mmm+1)
        ts.shift.total.plot2[[1]] = ts.shift.total.plot2[[1]] +
          geom_segment(data = ddd[[1]][[mmm]], 
                       aes(x = xx ,
                           xend = xx,
                           y =yy, 
                           yend = yy+1), 
                       color = color.m[1], 
                       size  = 1.1,
                       linetype = "solid")  #typelin[ccc])+
      }
    }
    ts.shift.total.plot2[[1]] = ts.shift.total.plot2[[1]] +
                                  scale_colour_manual(name = element_blank(),
                                                      values = colss[1],
                                                      breaks = chi.1,
                                                      labels = unname(chi2[1]),
                                                      guide = guide_legend(override.aes = list(
                                                                            linetype = c("solid"),
                                                                             shape = c(NA))))
    # plot chi30:
    oders.models.chi30 = c(1,2,5,3,9,4,7,6,10,8,11)   #put the models in the plot inorder of their performance
    ts.shift.total.plot2[[2]]   = ggplot() +
      scale_y_continuous(name = "Different stage recession models", 
                         expand = c(0,0), 
                         limits = c(0,length(model.titles)+2),  
                         breaks =seq(0.5, length(model.titles) + 1.5, 1), 
                         label = c("Official dates", "From Gaugings" , model.titles[oders.models.chi30])) +
      geom_hline(yintercept = seq(1, length(model.titles)+1, 1), color="gray70", size= 0.2, linetype ="dashed") +
      #geom_hline(yintercept = 2, color="black", size= 0.2) +
      coord_cartesian(clip = 'off')+
      theme_bw(base_size = 20) +
      theme(    axis.text=element_text(size= 20)
                ,axis.text.y=element_text(size= 20, face ="bold")
                ,axis.title=element_text(size= 30,face="bold")
                ,panel.grid.major=element_blank()
                ,panel.grid.minor=element_blank()
                #,panel.background=element_blank()
                ,legend.text=element_text(size=30)
                ,legend.title=element_text(size=30)
                ,legend.key.size=unit(1.5, "cm")
                ,legend.direction = "horizontal"
                ,legend.position ="bottom"
                ,axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))
                ,axis.title.y = element_blank()
                ,plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
      scale_x_continuous(name="Time [days]", expand = c(0,0), 
                         limits =time.limits) + 
      geom_point(data = data.annotate.off, 
                 aes(x = xeffect, y = rep(0, length(data.annotate.off$xeffect))), 
                 color= "black", size = 3, shape =4, stroke= 3) +
      geom_segment(aes(x    = bt.from.gaugingsss$t.shift.for.b$treal, 
                       xend = bt.from.gaugingsss$t.shift.for.b$treal, 
                       y    = 1, yend =2 ), color= "black", size= 1.1) +
      geom_vline(aes(xintercept=-1000, color = chi.2),
                 show.legend = T, size= 1.5, linetype = "solid")
    
    ddd[[2]]  = list()  
    for (mmm in 1:length(model.names)){
      if (!is.null(unlist(shift.times.recessions[[2]][[oders.models.chi30[mmm]]]$ts.res))) {
        ddd[[2]][[mmm]] = data.frame(xx =unlist(shift.times.recessions[[2]][[oders.models.chi30[mmm]]]$ts.res), 
                                     yy = mmm+1)
        ts.shift.total.plot2[[2]] = ts.shift.total.plot2[[2]] +
          geom_segment(data = ddd[[2]][[mmm]], 
                       aes(x = xx ,
                           xend = xx,
                           y =yy, 
                           yend = yy+1), 
                       color = color.m[2], 
                       size  = 1.1,
                       linetype = "solid")  #typelin[ccc])+
      }
    }
    ts.shift.total.plot2[[2]] = ts.shift.total.plot2[[2]] +
      scale_colour_manual(name = element_blank(),
                          values = colss[2],
                          breaks = chi.2,
                          labels = unname(chi2[2]),
                          guide = guide_legend(override.aes = list(
                            linetype = c("solid"),
                            shape = c(NA))))
    
    
    # plot chi50:
    oders.models.chi50 = c(2,1,4,9,11,3,5,6,8,10,7)  #put the models in the plot inorder of their performance
    ts.shift.total.plot2[[3]]   = ggplot() +
      scale_y_continuous(name = "Different stage recession models", 
                         expand = c(0,0), 
                         limits = c(0,length(model.titles)+2),  
                         breaks =seq(0.5, length(model.titles) + 1.5, 1), 
                         label = c("Official dates", "From Gaugings" , model.titles[oders.models.chi50])) +
      geom_hline(yintercept = seq(1, length(model.titles)+1, 1), color="gray70", size= 0.2, linetype ="dashed") +
      #geom_hline(yintercept = 2, color="black", size= 0.2) +
      coord_cartesian(clip = 'off')+
      theme_bw(base_size = 20) +
      theme(    axis.text=element_text(size= 20)
                ,axis.text.y=element_text(size= 20, face ="bold")
                ,axis.title=element_text(size= 30,face="bold")
                ,panel.grid.major=element_blank()
                ,panel.grid.minor=element_blank()
                #,panel.background=element_blank()
                ,legend.text=element_text(size=30)
                ,legend.title=element_text(size=30)
                ,legend.key.size=unit(1.5, "cm")
                ,legend.direction = "horizontal"
                ,legend.position ="bottom"
                ,axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))
                ,axis.title.y = element_blank()
                ,plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
      scale_x_continuous(name="Time [days]", expand = c(0,0), 
                         limits =time.limits) + 
      geom_point(data = data.annotate.off, 
                 aes(x = xeffect, y = rep(0, length(data.annotate.off$xeffect))), 
                 color= "black", size = 3, shape =4, stroke= 3) +
      geom_segment(aes(x = bt.from.gaugingsss$t.shift.for.b$treal, 
                       xend = bt.from.gaugingsss$t.shift.for.b$treal, 
                       y =1, yend =2 ), color= "black", size= 1.1) +
      geom_vline(aes(xintercept=-1000, color = chi.3),
                 show.legend = T, size= 1.5, linetype = "solid")
    
    ddd[[3]]  = list()  
    for (mmm in 1:length(model.names)){
      if (!is.null(unlist(shift.times.recessions[[3]][[oders.models.chi50[mmm]]]$ts.res))) {
        ddd[[3]][[mmm]] = data.frame(xx =unlist(shift.times.recessions[[3]][[oders.models.chi50[mmm]]]$ts.res), 
                                     yy = mmm+1)
        ts.shift.total.plot2[[3]] = ts.shift.total.plot2[[3]] +
          geom_segment(data = ddd[[3]][[mmm]], 
                       aes(x = xx ,
                           xend = xx,
                           y =yy, 
                           yend = yy+1), 
                       color = color.m[3], 
                       size  = 1.1,
                       linetype = "solid")  #typelin[ccc])+
      }
    }
    ts.shift.total.plot2[[3]] = ts.shift.total.plot2[[3]] +
      scale_colour_manual(name = element_blank(),
                          values = c(colss[3]),
                          breaks = c(chi.3),
                          labels = unname(chi2[3]),
                          guide = guide_legend(override.aes = list(
                            linetype = c("solid"),
                            shape = c(NA))))  

  ts.shift.total.plot.3 = plot_grid(            ts.shift.total.plot2[[1]], 
                                                ts.shift.total.plot2[[2]],
                                                ts.shift.total.plot2[[3]],
                                                ncol =3, nrow = 1)
  pdf(paste0(dir.regression,"/tshift_comparison_tot.pdf"), 20, 12 ,useDingbats=F)
  print(ts.shift.total.plot.3)
  dev.off()
  
  
  
  

  
  
  
  

  
  
  
  ##############################################################################################################
  # Deltab comparison:
  ####################
  ddd=NULL;   gglegend = NULL;  
  ts.shift.total.plot2=NULL; ts.shift.total.plot3 =NULL; 
  color.m = c("red", "blue", "green"); 
  typelin = c("solid", "solid", "solid");
  chi.1 <- c("chi = 10 cm");
  chi.2 = c("chi = 30 cm"); 
  chi.3 = c("chi = 50 cm")
  chi2 <- c(TeX("$\\chi = 10 \\; cm$"), TeX("$\\chi = 30 \\; cm$"), TeX("$\\chi = 50 \\; cm$"))
  fromgaug <- "Rating shift times obtained from gaugings"
  off.times <- "Official dates of RC update"
  colss <- c("chi = 10 cm" = "red", 
             "chi = 30 cm" = "blue",
             "chi = 50 cm" = "green",
             "Rating shift times obtained from gaugings" ="black",
             "Official dates of RC update" = "black")
  
  oders.models.chi10 = c(1,2,4,5,11,3,9,6,8,10,7)   #put the models in the plot inorder of their performance
  
  # plot chi10:
  ts.shift.total.plot2[[1]]   = ggplot() +
    scale_y_continuous(name = "Different stage recession models", 
                       expand = c(0,0), 
                       limits = c(0,length(model.titles)+2),  
                       breaks =seq(0.5, length(model.titles) + 1.5, 1), 
                       label = c("Official dates", "From Gaugings" , model.titles[oders.models.chi10])) +
    geom_hline(yintercept = seq(1, length(model.titles)+1, 1), color="gray70", size= 0.2, linetype ="dashed") +
    #geom_hline(yintercept = 2, color="black", size= 0.2) +
    coord_cartesian(clip = 'off')+
    theme_bw(base_size = 20) +
    theme(    axis.text=element_text(size= 20)
              ,axis.text.y=element_text(size= 20, face ="bold")
              ,axis.title=element_text(size= 30,face="bold")
              ,panel.grid.major=element_blank()
              ,panel.grid.minor=element_blank()
              #,panel.background=element_blank()
              ,legend.text=element_text(size=20)
              ,legend.title=element_text(size=30)
              ,legend.key.size=unit(1.5, "cm")
              ,legend.direction = "horizontal"
              ,legend.position ="bottom"
              ,axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))
              ,plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
    scale_x_continuous(name="Time [days]", expand = c(0,0), 
                       limits =time.limits) + 
    geom_point(data = data.annotate.off, 
               aes(x = xeffect, y = rep(0, length(data.annotate.off$xeffect))), 
               color= "black", size = 3, shape =4, stroke= 3) +
    geom_segment(aes(x = bt.from.gaugingsss$t.shift.for.b$treal, 
                     xend = bt.from.gaugingsss$t.shift.for.b$treal, 
                     y =1, yend =2 ), color= "black", size= 1.1) +
    geom_vline(aes(xintercept=-1000, color = chi.1),
               show.legend = T, size= 1.5, linetype = "solid")
  
  ddd[[1]]  = list()  
  for (mmm in 1:length(model.names)){
    if (!is.null(unlist(shift.times.recessions[[1]][[oders.models.chi10[mmm]]]$ts.res))) {
      ddd[[1]][[mmm]] = data.frame(xx =unlist(shift.times.recessions[[1]][[oders.models.chi10[mmm]]]$ts.res), 
                                   yy = mmm+1)
      ts.shift.total.plot2[[1]] = ts.shift.total.plot2[[1]] +
        geom_segment(data = ddd[[1]][[mmm]], 
                     aes(x = xx ,
                         xend = xx,
                         y =yy, 
                         yend = yy+1), 
                     color = color.m[1], 
                     size  = 1.1,
                     linetype = "solid")  #typelin[ccc])+
    }
  }
  ts.shift.total.plot2[[1]] = ts.shift.total.plot2[[1]] +
    scale_colour_manual(name = element_blank(),
                        values = colss[1],
                        breaks = chi.1,
                        labels = unname(chi2[1]),
                        guide = guide_legend(override.aes = list(
                          linetype = c("solid"),
                          shape = c(NA))))
  # plot chi30:
  oders.models.chi30 = c(1,2,5,3,9,4,7,6,10,8,11)   #put the models in the plot inorder of their performance
  ts.shift.total.plot2[[2]]   = ggplot() +
    scale_y_continuous(name = "Different stage recession models", 
                       expand = c(0,0), 
                       limits = c(0,length(model.titles)+2),  
                       breaks =seq(0.5, length(model.titles) + 1.5, 1), 
                       label = c("Official dates", "From Gaugings" , model.titles[oders.models.chi30])) +
    geom_hline(yintercept = seq(1, length(model.titles)+1, 1), color="gray70", size= 0.2, linetype ="dashed") +
    #geom_hline(yintercept = 2, color="black", size= 0.2) +
    coord_cartesian(clip = 'off')+
    theme_bw(base_size = 20) +
    theme(    axis.text=element_text(size= 20)
              ,axis.text.y=element_text(size= 20, face ="bold")
              ,axis.title=element_text(size= 30,face="bold")
              ,panel.grid.major=element_blank()
              ,panel.grid.minor=element_blank()
              #,panel.background=element_blank()
              ,legend.text=element_text(size=20)
              ,legend.title=element_text(size=30)
              ,legend.key.size=unit(1.5, "cm")
              ,legend.direction = "horizontal"
              ,legend.position ="bottom"
              ,axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))
              ,axis.title.y = element_blank()
              ,plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
    scale_x_continuous(name="Time [days]", expand = c(0,0), 
                       limits =time.limits) + 
    geom_point(data = data.annotate.off, 
               aes(x = xeffect, y = rep(0, length(data.annotate.off$xeffect))), 
               color= "black", size = 3, shape =4, stroke= 3) +
    geom_segment(aes(x    = bt.from.gaugingsss$t.shift.for.b$treal, 
                     xend = bt.from.gaugingsss$t.shift.for.b$treal, 
                     y    = 1, yend =2 ), color= "black", size= 1.1) +
    geom_vline(aes(xintercept=-1000, color = chi.2),
               show.legend = T, size= 1.5, linetype = "solid")
  
  ddd[[2]]  = list()  
  for (mmm in 1:length(model.names)){
    if (!is.null(unlist(shift.times.recessions[[2]][[oders.models.chi30[mmm]]]$ts.res))) {
      ddd[[2]][[mmm]] = data.frame(xx =unlist(shift.times.recessions[[2]][[oders.models.chi30[mmm]]]$ts.res), 
                                   yy = mmm+1)
      ts.shift.total.plot2[[2]] = ts.shift.total.plot2[[2]] +
        geom_segment(data = ddd[[2]][[mmm]], 
                     aes(x = xx ,
                         xend = xx,
                         y =yy, 
                         yend = yy+1), 
                     color = color.m[2], 
                     size  = 1.1,
                     linetype = "solid")  #typelin[ccc])+
    }
  }
  ts.shift.total.plot2[[2]] = ts.shift.total.plot2[[2]] +
    scale_colour_manual(name = element_blank(),
                        values = colss[2],
                        breaks = chi.2,
                        labels = unname(chi2[2]),
                        guide = guide_legend(override.aes = list(
                          linetype = c("solid"),
                          shape = c(NA))))
  
  
  # plot chi50:
  oders.models.chi50 = c(2,1,4,9,11,3,5,6,8,10,7)  #put the models in the plot inorder of their performance
  ts.shift.total.plot2[[3]]   = ggplot() +
    scale_y_continuous(name = "Different stage recession models", 
                       expand = c(0,0), 
                       limits = c(0,length(model.titles)+2),  
                       breaks =seq(0.5, length(model.titles) + 1.5, 1), 
                       label = c("Official dates", "From Gaugings" , model.titles[oders.models.chi50])) +
    geom_hline(yintercept = seq(1, length(model.titles)+1, 1), color="gray70", size= 0.2, linetype ="dashed") +
    #geom_hline(yintercept = 2, color="black", size= 0.2) +
    coord_cartesian(clip = 'off')+
    theme_bw(base_size = 20) +
    theme(    axis.text=element_text(size= 20)
              ,axis.text.y=element_text(size= 20, face ="bold")
              ,axis.title=element_text(size= 30,face="bold")
              ,panel.grid.major=element_blank()
              ,panel.grid.minor=element_blank()
              #,panel.background=element_blank()
              ,legend.text=element_text(size=20)
              ,legend.title=element_text(size=30)
              ,legend.key.size=unit(1.5, "cm")
              ,legend.direction = "horizontal"
              ,legend.position ="bottom"
              ,axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))
              ,axis.title.y = element_blank()
              ,plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
    scale_x_continuous(name="Time [days]", expand = c(0,0), 
                       limits =time.limits) + 
    geom_point(data = data.annotate.off, 
               aes(x = xeffect, y = rep(0, length(data.annotate.off$xeffect))), 
               color= "black", size = 3, shape =4, stroke= 3) +
    geom_segment(aes(x = bt.from.gaugingsss$t.shift.for.b$treal, 
                     xend = bt.from.gaugingsss$t.shift.for.b$treal, 
                     y =1, yend =2 ), color= "black", size= 1.1) +
    geom_vline(aes(xintercept=-1000, color = chi.3),
               show.legend = T, size= 1.5, linetype = "solid")
  
  ddd[[3]]  = list()  
  for (mmm in 1:length(model.names)){
    if (!is.null(unlist(shift.times.recessions[[3]][[oders.models.chi50[mmm]]]$ts.res))) {
      ddd[[3]][[mmm]] = data.frame(xx =unlist(shift.times.recessions[[3]][[oders.models.chi50[mmm]]]$ts.res), 
                                   yy = mmm+1)
      ts.shift.total.plot2[[3]] = ts.shift.total.plot2[[3]] +
        geom_segment(data = ddd[[3]][[mmm]], 
                     aes(x = xx ,
                         xend = xx,
                         y =yy, 
                         yend = yy+1), 
                     color = color.m[3], 
                     size  = 1.1,
                     linetype = "solid")  #typelin[ccc])+
    }
  }
  ts.shift.total.plot2[[3]] = ts.shift.total.plot2[[3]] +
    scale_colour_manual(name = element_blank(),
                        values = c(colss[3]),
                        breaks = c(chi.3),
                        labels = unname(chi2[3]),
                        guide = guide_legend(override.aes = list(
                          linetype = c("solid"),
                          shape = c(NA))))  
  
  ts.shift.total.plot.3 = plot_grid(            ts.shift.total.plot2[[1]], 
                                                ts.shift.total.plot2[[2]],
                                                ts.shift.total.plot2[[3]],
                                                ncol =3, nrow = 1)
  pdf(paste0(dir.regression,"/tshift_comparison_tot.pdf"), 20, 12 ,useDingbats=F)
  print(ts.shift.total.plot.3)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  ################################################################################################
  # DIC COMPARISON
  #####################
  # Initialise the plot:
  DIC.plot <- ggplot() + 
    theme_light(base_size = 20)+
    scale_x_continuous(name = "Recession curve model", expand = c(0,0), 
                       limits =c(0, length(model.titles)+1), breaks = seq(1,length(model.titles)), labels =model.titles ) + 
    scale_y_continuous(name = "Deviance Information Criterion") +
    coord_cartesian(clip = 'off') +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), 
          plot.background = element_rect(fill ="transparent", color = NA),
          panel.grid.major=element_line(size=0.6, linetype = "dashed"), 
          panel.grid.minor=element_blank(),
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill ="transparent"), 
          #axis.line = element_line(colour = "black"),
          axis.text.x =   element_text(face ="bold"),
          axis.ticks = element_line(colour = "black"),
          axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
          plot.margin=unit(c(1,0.5,1.5,1.05),"cm"),
          text = element_text(size=20),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.background = element_rect(colour = "transparent", fill = "transparent"),
          #      axis.line = element_line(colour = "black"),
          #legend.title = element_text(colour="blue", size=16, face="bold"),
          legend.position="bottom",
          legend.direction = "horizontal")+
    geom_point(aes(x= seq(1:length(model.names)), y = unlist(DIC.rec[[1]]), color = chi.1),
               size= 6,  pch = shapppe[1], na.rm = TRUE)+
    geom_line(aes(x= seq(1:length(model.names)), y = unlist(DIC.rec[[1]])), size= 1, 
               color = color.m[1], na.rm = TRUE, linetype = "dashed")+
    geom_point(aes(x= seq(1:length(model.names)), y = unlist(DIC.rec[[2]]), color = chi.2), 
               size= 6, pch = shapppe[2], na.rm = TRUE) +
    geom_line(aes(x= seq(1:length(model.names)), y = unlist(DIC.rec[[2]])), size= 1, 
              color = color.m[2], na.rm = TRUE, linetype = "dashed") +
    geom_point(aes(x= seq(1:length(model.names)), y = unlist(DIC.rec[[3]]), color = chi.3), 
               size= 6,  pch = shapppe[3], na.rm = TRUE) +
    geom_line(aes(x= seq(1:length(model.names)), y = unlist(DIC.rec[[3]])), size= 1, 
              color = color.m[3], na.rm = TRUE, linetype = "dashed") +
    scale_colour_manual(name = element_blank(),
                        values = c(colss[1], colss[2], colss[3]),
                        breaks = c(chi.1, chi.2, chi.3),
                        labels = unname(chi2),
                        guide = guide_legend(override.aes = list(
                          linetype = c("blank", "blank", "blank"),
                          shape = c(shapppe[1], shapppe[2], shapppe[3]))))
# geom_point(aes(x= seq(1:length(model.names)), y = unlist(DIC.rec[[3]])), size= 5, 
    #            color = colllor[chi.i],  pch = shapppe[chi.i], na.rm = TRUE)+
  pdf(paste0(dir.regression,"/DIC_model_comparison.pdf"), 10, 7 ,useDingbats=F)
  print(DIC.plot)
  dev.off()
  
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ############################################################################################################################
  #                                          RMSE b1 COMPARISON
  ############################################################################################################################
  RMSE.b1.plot <- ggplot() + 
    theme_light(base_size = 20)+
    scale_x_continuous(name = "Recession curve model", expand = c(0,0), 
                       limits =c(0, length(model.titles)+1), breaks = seq(1,length(model.titles)), labels =model.titles ) + 
    scale_y_continuous(name = TeX("RMSE (between $\\beta$ and $b_1$)")) +
    coord_cartesian(clip = 'off') +
    theme(plot.title = element_text(hjust = 0.5), 
          plot.background = element_rect(fill ="transparent", color = NA),
          panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
          panel.grid.minor=element_blank(),
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill ="transparent"), 
          #axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.margin=unit(c(1,0.5,1.5,1.05),"cm"),
          text = element_text(size=20),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.background = element_rect(colour = "transparent", fill = "transparent"),
          #      axis.line = element_line(colour = "black"),
          #legend.title = element_text(colour="blue", size=16, face="bold"),
          legend.position="none") +
    geom_point(aes(x= seq(1:length(model.names)), y = unlist(rmse.b1.rec[[1]])), size= 5,
               color = colllor[1],  pch = shapppe[1], na.rm = TRUE)+
    geom_line(aes(x= seq(1:length(model.names)), y = unlist(rmse.b1.rec[[1]])), size= 1,
               color = colllor[1],  linetype = "dashed", na.rm = TRUE)+
    geom_point(aes(x= seq(1:length(model.names)), y = unlist(rmse.b1.rec[[2]])), size= 5,
               color = colllor[2],  pch = shapppe[2], na.rm = TRUE)+
    geom_line(aes(x= seq(1:length(model.names)), y = unlist(rmse.b1.rec[[2]])), size= 1,
              color = colllor[2],  linetype = "dashed", na.rm = TRUE)+
    theme(text = element_text(size=14),
          axis.text = element_text(size=10))
    # geom_point(aes(x= seq(1:length(model.names)), y = unlist(rmse.b1.rec[[3]][[mod]])), size= 5,
    #            color = colllor[3],  pch = shapppe[3], na.rm = TRUE)
    pdf(paste0(dir.regression,"/RMSE_b1_model_comparison.pdf"), 10, 10 ,useDingbats=F)
    print(RMSE.b1.plot)
    dev.off()
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # RMSE.b2.plot <- ggplot() + 
    #   theme_light(base_size = 20)+
    #   scale_x_continuous(name = "Recession curve model", expand = c(0,0), 
    #                      limits =c(0, length(model.titles)+1), breaks = seq(1,length(model.titles)), labels =model.titles ) + 
    #   scale_y_continuous(name = TeX("RMSE (between $\\beta$ and $b_2$)")) +
    #   coord_cartesian(clip = 'off') +
    #   theme(plot.title = element_text(hjust = 0.5), 
    #         plot.background = element_rect(fill ="transparent", color = NA),
    #         panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
    #         panel.grid.minor=element_blank(),
    #         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #         panel.background = element_rect(fill ="transparent"), 
    #         #axis.line = element_line(colour = "black"),
    #         axis.ticks = element_line(colour = "black"),
    #         plot.margin=unit(c(1,0.5,1.5,1.05),"cm"),
    #         text = element_text(size=20),
    #         legend.key = element_rect(colour = "transparent", fill = "transparent"),
    #         legend.background = element_rect(colour = "transparent", fill = "transparent"),
    #         #      axis.line = element_line(colour = "black"),
    #         #legend.title = element_text(colour="blue", size=16, face="bold"),
    #         legend.position="none")
    # 
    # RMSE.shift.plot <- ggplot() + 
    #   theme_light(base_size = 20)+
    #   scale_x_continuous(name = "Recession curve model", expand = c(0,0), 
    #                      limits =c(0, length(model.titles)+1), breaks = seq(1,length(model.titles)), labels =model.titles ) + 
    #   scale_y_continuous(name = "RMSE of shift times") +
    #   coord_cartesian(clip = 'off') +
    #   theme(plot.title = element_text(hjust = 0.5), 
    #         plot.background = element_rect(fill ="transparent", color = NA),
    #         panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
    #         panel.grid.minor=element_blank(),
    #         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #         panel.background = element_rect(fill ="transparent"), 
    #         #axis.line = element_line(colour = "black"),
    #         axis.ticks = element_line(colour = "black"),
    #         plot.margin=unit(c(1,0.5,1.5,1.05),"cm"),
    #         text = element_text(size=20),
    #         legend.key = element_rect(colour = "transparent", fill = "transparent"),
    #         legend.background = element_rect(colour = "transparent", fill = "transparent"),
    #         #      axis.line = element_line(colour = "black"),
    #         #legend.title = element_text(colour="blue", size=16, face="bold"),
    #         legend.position="none")
    # 
    # Accuracy.plot <- ggplot() + 
    #   theme_light(base_size = 20)+
    #   scale_x_continuous(name = "Recession curve model", expand = c(0,0), 
    #                      limits =c(0, length(model.titles)+1), breaks = seq(1,length(model.titles)), labels =model.titles ) + 
    #   scale_y_continuous(name = "Accuracy of the segmentation") +
    #   coord_cartesian(clip = 'off') +
    #   theme(plot.title = element_text(hjust = 0.5), 
    #         plot.background = element_rect(fill ="transparent", color = NA),
    #         panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
    #         panel.grid.minor=element_blank(),
    #         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #         panel.background = element_rect(fill ="transparent"), 
    #         #axis.line = element_line(colour = "black"),
    #         axis.ticks = element_line(colour = "black"),
    #         plot.margin=unit(c(1,0.5,1.5,1.05),"cm"),
    #         text = element_text(size=20),
    #         legend.key = element_rect(colour = "transparent", fill = "transparent"),
    #         legend.background = element_rect(colour = "transparent", fill = "transparent"),
    #         legend.direction = "horizontal",
    #         legend.text = element_text(size = 20))
    # 
  # RMSE.b2.plot = RMSE.b2.plot +
  #   geom_point(aes(x= seq(1:length(model.names)), y = unlist(rmse.b2.rec[[1]])), size= 5,
  #              color = colllor[1],  pch = shapppe[1], na.rm = TRUE)+
  #   geom_line(aes(x= seq(1:length(model.names)), y = unlist(rmse.b2.rec[[1]])), size= 1,
  #              color = colllor[1],   linetype = "dashed", na.rm = TRUE)+
  #   geom_point(aes(x= seq(1:length(model.names)), y = unlist( rmse.b2.rec[[2]])), size= 5,
  #              color = colllor[2],  pch = shapppe[2], na.rm = TRUE)+
  #   geom_line(aes(x= seq(1:length(model.names)), y = unlist(rmse.b2.rec[[2]])), size= 1,
  #              color = colllor[2],   linetype = "dashed", na.rm = TRUE)+
  #   theme(text = element_text(size=14),
  #         axis.text = element_text(size=10))
  #   # geom_point(aes(x= seq(1:length(model.names)), y = unlist(rmse.b2.rec[[3]][[mod]])), size= 5,
  #   #            color = colllor[3],  pch = shapppe[3], na.rm = TRUE)
  # 
  # 
  # 
  # 
  # RMSE.shift.plot = RMSE.shift.plot +
  #   geom_point(aes(x= seq(1:length(model.names)), y = unlist(shift.rmse.rec[[1]])), size= 5,
  #              color = colllor[1],  pch = shapppe[1], na.rm = TRUE)+
  #   geom_line(aes(x= seq(1:length(model.names)), y = unlist(shift.rmse.rec[[1]])), size= 1,
  #             color = colllor[1],   linetype = "dashed", na.rm = TRUE)+
  #   geom_point(aes(x= seq(1:length(model.names)), y = unlist( shift.rmse.rec[[2]])), size= 5,
  #              color = colllor[2],  pch = shapppe[2], na.rm = TRUE)+
  #   geom_line(aes(x= seq(1:length(model.names)), y = unlist(shift.rmse.rec[[2]])), size= 1,
  #             color = colllor[2],   linetype = "dashed", na.rm = TRUE)+
  #   theme(text = element_text(size=14),
  #         axis.text = element_text(size=10))
  # # geom_point(aes(x= seq(1:length(model.names)), y = unlist(rmse.b2.rec[[3]][[mod]])), size= 5,
  # #            color = colllor[3],  pch = shapppe[3], na.rm = TRUE)
  # 
  # # legend:
  # chi1 <- "chi = 10 cm"
  # chi2 <- "chi = 30 cm"
  # cols <- c("chi = 10 cm" = "red", 
  #           "chi = 30 cm" = "green")
  # # colllor =c("red", "green", "blue")
  # # shapppe = c(15,    19 ,     17)
  # Accuracy.plot = Accuracy.plot +
  #   geom_point(aes(x= seq(1:length(model.names)), y = unlist(accu.rec[[1]]), color = chi1),
  #              size= 5, pch = shapppe[1], na.rm = TRUE)+
  #   geom_line(aes(x= seq(1:length(model.names)), y = unlist(accu.rec[[1]])), size= 1,
  #              color = colllor[1],   linetype = "dashed",  na.rm = TRUE)+
  #   geom_point(aes(x= seq(1:length(model.names)), y = unlist(accu.rec[[2]]), color = chi2), 
  #              size= 5, pch = shapppe[2], na.rm = TRUE)  +
  #   geom_line(aes(x= seq(1:length(model.names)), y = unlist(accu.rec[[2]])), size= 1,
  #             color = colllor[2],   linetype = "dashed",  na.rm = TRUE) + 
  #   scale_colour_manual(name = element_blank(), 
  #                       values=c(cols),
  #                       breaks=c(chi1,
  #                                chi2
  #                                ),
  #                       labels = c(chi1 = TeX("$\\chi = 10 cm$"), 
  #                                  chi2 = TeX("$\\chi = 30 cm$")),
  #                       #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
  #                       guide = guide_legend(override.aes = list(
  #                         linetype = c("blank", "blank"),
  #                         shape = c(15, 19))))+
  #   theme(text = element_text(size=14),
  #         axis.text = element_text(size=10))
  # 
  # Accuracy.plot.without.legend = Accuracy.plot + theme(legend.position="none")
  #   # geom_point(aes(x= seq(1:length(model.names)), y = unlist(accu.rec[[3]])), size= 5,
  #   #            color = colllor[3],  pch = shapppe[3], na.rm = TRUE)
  # 
  # 
  # pdf(paste0(dir.regression,"/RMSE_b2_model_comparison.pdf"), 10, 10 ,useDingbats=F)
  # print(RMSE.b2.plot)
  # dev.off()
  # 
  # pdf(paste0(dir.regression,"/Accuracy_model_comparison.pdf"), 10, 10 ,useDingbats=F)
  # print(Accuracy.plot)
  # dev.off()
  # 
  # pdf(paste0(dir.regression,"/rmse_shift_times_model_comparison.pdf"), 10, 10 ,useDingbats=F)
  # print(RMSE.shift.plot)
  # dev.off()
  # 
  # performance.plot = plot_grid(DIC.plot, Accuracy.plot.without.legend, RMSE.b1.plot, RMSE.b2.plot, #RMSE.shift.plot,
  #                              ncol =2, nrow = 2, labels = c("a)", "b)", "c)", "d)"), label_size = 20)
  #
  # # + the legend:
  # #########################################################################
  # Legend.title <- ggdraw() + 
  #                 draw_label("Legend",
  #                 fontface = 'bold',
  #                 x = 0, hjust = 0,
  #                 size = 20) +
  #                 theme(# add margin on the left of the drawing canvas,
  #                 # so title is aligned with left edge of first plot
  #                 plot.margin = margin(10, 0, 10, 0))
  # glegend <- cowplot::get_legend( Accuracy.plot)
  # Legends =  plot_grid(Legend.title, 
  #                      glegend,
  #                      ncol = 1,
  #                      nrow = 2,
  #                      # rel_heights values control vertical title margins
  #                      rel_heights = c(0.3, 1))
  # 
  # performance.plot.and.legend =  plot_grid(performance.plot, 
  #                                          Legends,
  #                                          ncol =1,
  #                                          nrow = 2,
  #                                          rel_heights = c(1,0.2))
  # pdf(paste0(dir.regression,"/performance_rec_model_comparison.pdf"), 15, 10 ,useDingbats=F)
  # print(performance.plot.and.legend)
  # dev.off()
  
}




















































############################################################################################################
plot.segmentation.recession.one.model = function(dir.rec.pool.test,  
                                                 dir.rec.segm.test,
                                                 rec.mod,
                                                 stage.limits,
                                                 stage.scale.shift,
                                                 asymptote.limits,
                                                 limits.Y.lambda,
                                                 limits.Y.alpha,
                                                 limits.x.recess,
                                                 x.name,
                                                 plot.b.from.gaugings,
                                                 plot.recession.uncert,
                                                 plot.gamma.uncertainty,
                                                 plot.dates,
                                                 date_origin,
                                                 initial.time,
                                                 gaugings,
                                                 limni.labels,
                                                 grid_limni.ylim,
                                                 df.limni) {
  ##############################################################################################################
  # Recession models studied :
  #***************************
  # 1) 1exp =>                 h(t) = alpha1 *exp(-lambda1* t) + beta
  # 2) 2exp =>                 h(t) = alpha1 *exp(-lambda1* t) + alpha1 *exp(-lambda1* t) + beta
  # 3) 3exp =>                 h(t) = alpha1 *exp(-lambda1* t) + alpha2 *exp(-lambda2* t) + alpha3*exp(-lambda3* t) + beta
  # 4) double exp (horton) =>  h(t) = alpha1 *exp(-lambda1* t^n1) + beta
  # 5) hyperb   =>             h(t) = alpha1 /(1+lambda1*t)^2
  # 6) Boussinesq (1903) =>    h(t) = alpha1 *(1 + lambda1*t)^(-2)
  # 7) Coutagne
  # model.names = c("expexp", "hyperb", "expexp", "hyperb", "expexp", "hyperb")
  

  # preparing  the dataframe for riverbed estimated through gaugings :
  # if (plot.b.from.gaugings == TRUE){
  #   dir.sed.transp = paste0(dir.case_study,"/Results/segmentation_sed_transp")
  #   dir.sed.transp.SPD = paste0(dir.sed.transp,"/SPD")
  #   bt1.df = read.table(paste0(dir.sed.transp.SPD,"/bt1_df.txt"), header =TRUE)
  #   bt2.df = read.table(paste0(dir.sed.transp.SPD,"/bt2_df.txt"), header =TRUE)
  #   ts.gaug =  read.table(paste0(dir.sed.transp.SPD,"/STEP1_shift_times.txt"), header =TRUE)
  #   ts.before.gaug =  sort(c(0, ts$t.adj))
  #   ts.plus.gaug = sort(c(ts$t.adj, tail(df.limni$t_limni,1)))
  # }
  
  # recess model with tot uncertainty:
  #*******************************************
  Recession.model = function(theta, model, t){ 
  #*******************************************
    h=0*t
    if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[8],theta[9]))
    } else if (model =="2expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="2expWithAsympt_bis"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="2expWithAsympt_rel"){
      h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="1expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[4],theta[5]))
    } else if ((model =="expexp")|(model =="expexp_bis")){
      h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } else if ((model =="hyperb")|(model =="hyperb_bis")){
      h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } 
    return(res.h)
  }
  
  # recess model for maxpost:
  #************************************************
  Recess.Maxpost.model = function(theta, model, t){
    #************************************************
    h=0*t
    if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
    }  else if (model =="2expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
    }  else if (model =="2expWithAsympt_bis"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
    }  else if (model =="2expWithAsympt_rel"){
      h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
        theta[5] 
    } else if ((model =="1expWithAsympt")|(model =="1expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) + theta[3] 
    } else if ((model =="expexp")|(model =="expexp_bis")){
      h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
    } else if ((model =="hyperb")|(model =="hyperb_bis")){
      h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
    } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
    }
    return(h)
  }
  
  #initialisation of the lists of plot objects:
  reg.pool.plot=NULL; reg.pool.plot2=NULL;  seg.rec.plot=NULL; t.plot=NULL; t.plot2 =NULL; t.plot3 =NULL
  t.plot4 = NULL; title.model=NULL; dir.rec.segm.test.param=NULL;
  tau.results.df.rec =NULL; mu.results.df.rec=NULL; gamma.results.df.rec=NULL; df.shift.times.rec=NULL
  df.shift.times.plus.rec=NULL; ts.res.before.rec=NULL; ts.res.rec=NULL; ts.res.plus.rec=NULL; Q2.ts.rec=NULL
  Q97.ts.rec=NULL; Q2.mu.rec=NULL; mu.res.rec=NULL; Q97.mu.rec=NULL; Data.segm.rec=NULL; nS.ok.rec=NULL;
  mcmc.seg.rec=NULL; pdf.ts.rec=NULL; gamma_segm_recess=NULL; gaugings.df.recess=NULL;
  shift.times.recessions=NULL; data.annotate.recess=NULL; data.annotate.recess.adjust=NULL; 
  X=NULL; X1=NULL; data.tmp= NULL;data.tmp.2= NULL; Ysim=NULL; time.adjust.before =NULL; 
  time.adjust.plus=NULL; parameters =NULL; parameters.names =NULL;  model.title =0
  
  # do a specific subfigure column for the results of one selected model only: 
  #####
  model.names = rec.mod
    # read data mcmc:  
    data.MCMC.cooked  = as.matrix(read.table(paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"), header=TRUE,dec=".", sep=""))
    data.MCMC.MaxPost = as.numeric(read.table(paste0(dir.rec.pool.test,"/Results_Summary.txt"), row.names=1,dec=".",sep="", skip = 16))
    colfunc           = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))
    summary.rec       = read.table(file=paste0(dir.rec.pool.test,"/Results_Summary.txt"), header=TRUE)
    mcmc.rec          = read.table(file=paste0(dir.rec.pool.test,"/Results_MCMC_cooked.txt"), header=TRUE)
    residuals.rec     = read.table(file=paste0(dir.rec.pool.test,"/Results_Residuals.txt"), header=TRUE)
    curves.data.rec   = read.table(file=paste0(dir.rec.pool.test,"/Curves_Data.txt"), header=TRUE)
    
    ###########################################################################################
    # Recessions regression plot:
    ###########################################################################################
    print("plot 1")
    nsample      = length(data.MCMC.cooked[,1])
    tgrid        = seq(limits.x.recess[1],  
                       limits.x.recess[2],  0.5) 
    Ncurves.pool = tail(curves.data.rec$Period,1)
    
    #Initialisation:
    if (model.names =="3expWithAsympt"){
      ########################################
      model.title = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2 e^{-\\lambda_2 t} +  \\alpha_3 e^{-\\lambda_3 t} + \\beta^{(k)}$")
      b1.var =FALSE
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                                  #a1, b1,
                     Ncurves.pool +2, Ncurves.pool + 3,                    #a2, b2
                     Ncurves.pool +4, Ncurves.pool + 5,                      #a3, b3,
                     Ncurves.pool +5 + i,                            #a4 (asymptotic level parameter)
                     Ncurves.pool*2 + 5+1, Ncurves.pool*2 + 5+2)  #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
    } else if (model.names =="3expWithAsympt_bis"){
      #####################################################
      model.title = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2^{(k)} e^{-\\lambda_2 t} + \\alpha_3 e^{-\\lambda_3 t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                               #a1(var), b1,
                     Ncurves.pool +1 + i, 2*Ncurves.pool + 2,         #a2(var), b2
                     2*Ncurves.pool + 3, 2*Ncurves.pool + 4,          #a3, b3,
                     2*Ncurves.pool + 4 + i,                          #a4(var) (asymptotic level parameter)
                     3*Ncurves.pool + 4 + 1, 3*Ncurves.pool + 6)      #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }  
      
      
    }  else if (model.names =="2expWithAsympt_bis"){
      #####################################################
      model.title = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2^{(k)} e^{-\\lambda_2 t}  + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=8) # 7 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=8)      # 7 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                     Ncurves.pool + 1 + i, Ncurves.pool*2 +1 + 1,   #a2, b2
                     Ncurves.pool*2 +1+1+i,                       #a3,
                     Ncurves.pool*3 + 2 + 1, Ncurves.pool*3 +2+2)  #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
      
    }  else if ((model.names=="2expWithAsympt")|(model.names =="2expWithAsympt_rel")){
      ##############################################################################################
      model.title = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2 e^{-\\lambda_2  t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=8) # 7 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=8)      # 7 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                             #a1, b1,
                     Ncurves.pool +2, Ncurves.pool +3,              #a2, b2
                     Ncurves.pool + 3 + i,                           #a3,
                     Ncurves.pool*2 + 3 + 1, Ncurves.pool*2 +3 +2)  #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
      
    } else if (model.names =="1expWithAsympt"){
      ################################################
      model.title = TeX("$h(t,k) = \\alpha^{(k)} \\; e^{-\\lambda t} + \\beta^{(k)}$")
      MCMC.save    = matrix(NA, nrow=Ncurves.pool*nsample, ncol=6) # 5 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=6)      # 5 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                     Ncurves.pool +1+ i,                      #a2
                     Ncurves.pool*2 +1+1, Ncurves.pool*2 +1 + 2)  #gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
      
    } else if ((model.names=="expexp")|(model.names =="hyperb")|(model.names =="Coutagne")) {
      ##########################################################################################################
      model.title = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t}^{\\eta t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool + 1,                          # a1(k), b1,
                     Ncurves.pool + 2,                             # n1
                     Ncurves.pool + 2 + i,                         # a2(k) (asymptotic level parameter)
                     Ncurves.pool*2 + 2 +1, Ncurves.pool*2 +2 +2)  # gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      }
      
      
    } else if ((model.names =="expexp_bis")|(model.names =="hyperb_bis")|(model.names =="Coutagne_bis")){
      ####################################################################################################################
      model.title = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t}^{\\eta^{(k)} t} + \\beta^{(k)}$")
      MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
      MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period 
      for(i in 1:Ncurves.pool){
        col.num = c( i, Ncurves.pool + i,                     # a1(k), b1(k),
                     Ncurves.pool*2 + 1,                      # n1
                     Ncurves.pool*2 + 1 + i,                  # a2 (asymptotic level parameter)
                     Ncurves.pool*3 + 2, Ncurves.pool*3 +3)   # gamma1, gamma2
        MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
        MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
      } 
    }
    
    #######################################################################################################
    # Apply recession model (total uncertainty and maxpost):
    message("- Plotting all Recession curves !!!  Wait ... "); flush.console()
    # Rec.Post    = apply(MCMC.save,    MARGIN=1, Recession.model,      model=model.names[mod],  t=tgrid) #add structural error:
    # Rec.MaxPost = apply(MaxPost.save, MARGIN=1, Recess.Maxpost.model, model=model.names[mod],  t=tgrid) # Maximum posterior 
    # # Quantiles: 
    # List.Rec.quants = list(NULL)
    # for(i in 1:Ncurves.pool){
    #   data.tmp = apply(Rec.Post[,(nsample*(i-1)+1):(nsample*i)], 
    #                    MARGIN=1, quantile, probs=c(0.025, 0.975),  na.rm=TRUE)
    #   List.Rec.quants[[i]] = data.frame(cbind(tgrid, 
    #                                           t(data.tmp - stage.scale.shift),
    #                                           Rec.MaxPost[,i] - stage.scale.shift))  #!!!!!!!!!!!!!! 
    #   colnames(List.Rec.quants[[i]]) = c("t", "inf", "sup", "maxpost")
    # }
    # 
    # #prepare plot:
    inter.per     = seq(1,Ncurves.pool,1) 
    nRec          = length(inter.per)
    palette.per   = colfunc(Ncurves.pool)
    palette.per   = palette.per[1:Ncurves.pool]
    data.rec.obs  = read.table(paste0(dir.rec.pool.test,"/Curves_Data.txt"), header=TRUE,dec=".",sep="") # Gaugings loading
    ylim.wind     = c(stage.limits[1],    stage.limits[2])
    xlim.wind     = c(limits.x.recess[1], limits.x.recess[2])
    pos.num = function(x.int){
               inter = which(x.int==data.rec.obs$period);
               return(inter)}
    inter.rec.obs       = unlist(sapply(inter.per, pos.num), recursive = TRUE, use.names = TRUE)
    data.rec.obs$Period = as.factor(data.rec.obs$Period)
    data.Rec            = read.table(file = paste0(dir.rec.pool.test,"/Rec_SPD_env.txt"), header=TRUE)
    data.Rec$Period     = factor(data.Rec$Period)
    PltData.rec         = data.Rec[data.Rec$inf > -200,]
    
    
    ###########################################
    reg.pool.plot = ggplot(PltData.rec) +
      # geom_smooth(aes(x=t,
      #                 y=maxpost,
      #                 ymax=sup, 
      #                 ymin=inf, 
      #                 group=Period, 
      #                 fill=Period), 
      #             size=0.1, stat='identity', alpha=0.1)  + #alpha=0.1
      geom_path(aes(x=t,   
                    y=maxpost, 
                    group=Period, 
                    colour=Period), size=0.1) +
      ### recession curve obs:
      geom_linerange(aes(x      = time, 
                         ymax   = (h + 2*uh) - stage.scale.shift, 
                         ymin   = (h-2*uh)   - stage.scale.shift, 
                         colour = Period),
                     data=data.rec.obs, size=0.2)+
      geom_point(aes(x=time, 
                     y=h-stage.scale.shift,
                     colour=Period), 
                 data=data.rec.obs, shape=16, size=2)+
      ### Labels
      xlab("Recession time [days]")+
      ylab("Stage h [cm]") +
      labs(colour = "Period") +
      scale_colour_manual(values = palette.per) +
      # scale_fill_manual(values = palette.per)+
      #scale_y_continuous(breaks=breaks.wind,labels=labels.wind)+
      # coord_cartesian(ylim = ylim.wind, 
      #                 xlim = xlim.wind) +
      scale_x_continuous(expand = c(0,0))+  #, limits = xlim.wind) +
      scale_y_continuous(expand = c(0,0)) +
      ### Theme
      theme_light(base_size=20)+
      theme(axis.text         = element_text(size=15)
            ,axis.title       = element_text(size=20,face="plain")
            ,plot.title       = element_text(size=25, face="bold", hjust = 0.5)
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,legend.text      = element_text(size=20)
            ,legend.title     = element_text(size=30)
            ,legend.key.size  = unit(1.5, "cm")
            ,legend.position  = "none"
            ,plot.margin      = unit(c(2, 0.5, 0, 1),"cm")
            ,axis.title.y     = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0)))

    
    
    
    
    
    
    
    
    #####################################################################################
    # Segmentation for the specific model mod and all parameters:
    #####################################################################################
    print("plot 2")
    message("- Plotting results of segmentation applied to the time series of beta (asymptote)")
    # initialise:
    dir.rec.segm.test.param  = list() 
    tau.results.df.rec       = list()
    mu.results.df.rec        = list() 
    gamma.results.df.rec     = list()  
    df.shift.times.rec       = list() 
    df.shift.times.plus.rec  = list()
    ts.res.before.rec        = list() 
    ts.res.rec               = list()
    ts.res.plus.rec          = list()  
    Q2.ts.rec                = list() 
    Q97.ts.rec               = list()    
    Q2.mu.rec                = list()
    mu.res.rec               = list() 
    Q97.mu.rec               = list()
    Data.segm.rec            = list() 
    nS.ok.rec                = list() 
    mcmc.seg.rec             = list()
    pdf.ts.rec               = list() 
    gamma_segm_recess        = list() 
    gaugings.df.recess       = list()   
    shift.times.recessions       = list() 
    data.annotate.recess         = list() 
    data.annotate.recess.adjust  = list() 
    X                            = list()
    X1                           = list() 
    data.tmp.2                   = list() 
    Ysim                         = list()
    time.adjust.before           = list()
    time.adjust.plus             = list()
    seg.rec.plot                 = list()

    
    if (model.names == "3expWithAsympt"){
      # parameters = c("a1", "a2", "a3", "a4")
      # parameters.names = c(TeX("$\\alpha_{1}$"), TeX("$\\alpha_{2}$"), TeX("$\\alpha_{3}$"), TeX("$\\beta$"))
      parameters = c("a1", "a4")
      parameters.names = c(TeX("$\\alpha_{1}$"), TeX("$\\beta \\left[ cm \\right]$"))
      
    } else if (model.names =="hyperb"){  #Drogue model
      parameters = c("a1", "a2")
      parameters.names = c(TeX("$\\alpha$"),  TeX("$\\beta [cm]$"))
      
    } else if (model.names =="expexp"){  #Horton model
      parameters = c("a1", "a2")
      parameters.names  = c(TeX("$\\alpha$"), TeX("$\\beta [cm]$"))
      
    } else if (model.names =="2expWithAsympt_bis"){
      parameters = c("a1", "a2", "a3")
      parameters.names = c(TeX("$\\alpha_{1}$"), TeX("$\\alpha_{2}$"), TeX("$\\beta [cm]$"))
      
    } else if (model.names  =="2expWithAsympt"){
      parameters  = c("a1", "a3")
      parameters.names  = c(TeX("$\\alpha_{1}$"), TeX("$\\beta [cm]$"))
      
    } else if ((model.names  =="1expWithAsympt")){
      parameters  = c("a1", "a2")
      parameters.names  = c(TeX("$\\alpha$"), TeX("$\\beta [cm]$"))
      
    } else if ((model.names  =="2expWithAsympt_rel")){
      parameters  = c("a1", "a3")
      parameters.names  = c(TeX("$\\alpha_{1}$"), TeX("$\\beta [cm]$"))
      
    } else if (model.names  =="Coutagne"){
      parameters  = c("a1", "a2")
      parameters.names  = c(TeX("$\\alpha$"),  TeX("$\\beta [cm]$"))
      
    } else if (model.names  == "3expWithAsympt_bis"){
      parameters  = c("a1", "a2", "a4")
      parameters.names  = c(TeX("$\\alpha_{1}$"),TeX("$\\alpha_{2}$"), TeX("$\\beta [cm]$"))
      
    } else if ((model.names  == "Coutagne_bis")|(model.names  == "expexp_bis")|(model.names  == "hyperb_bis")){
      parameters  = c("a1", "b1", "a2")
      parameters.names  = c(TeX("$\\alpha$"), TeX("$\\lambda$"),  TeX("$\\beta [cm]$"))
    }
    
    
    ##########################################################################
    param  = length(parameters)   # plot only last parameter: asymptotic stage
    ##########################################################################
      dir.rec.segm.test.param[[param]]   = paste0(dir.rec.segm.test,"/", parameters[param])
      mu.results.df.rec[[param]]         = read.table(file = paste0(dir.rec.segm.test.param[[param]], "/mu.results.df.txt"), header=TRUE)
      gamma.results.df.rec[[param]]      = read.table(file = paste0(dir.rec.segm.test.param[[param]], "/gamma.results.df.txt"), header=TRUE)
      Q2.mu.rec[[param]]                 = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/Q2.mu.txt"), header=TRUE)
      mu.res.rec[[param]]                = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/mu.res.txt"), header=TRUE)
      Q97.mu.rec[[param]]                = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/Q97.mu.txt"), header=TRUE)
      Data.segm.rec[[param]]             = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/Data.segm.rec.txt"), header=TRUE)
      mcmc.seg.rec[[param]]              = read.table(file=paste0(dir.rec.segm.test.param[[param]],   "/mcmc_segmentation.txt"), header=TRUE)
      gamma_segm_recess[[param]]         = c(mean(mcmc.seg.rec[[param]]$Y1_gamma1), std(mcmc.seg.rec[[param]]$Y1_gamma1))
      #
      #
      if (length(mu.results.df.rec[[param]]$mu.MAP) >1) {
        tau.results.df.rec[[param]]      = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/tau.results.df.txt"), header=TRUE)
        df.shift.times.rec[[param]]      = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/df.shift.times.txt"), header=TRUE)
        df.shift.times.plus.rec[[param]] = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/df.shift.times.plus.txt"), header=TRUE)
        ts.res.before.rec[[param]]       = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/ts.res.before.txt"), header=TRUE)
        ts.res.rec[[param]]              = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/ts.res.txt"), header=TRUE)
        ts.res.plus.rec[[param]]         = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/ts.res.plus.txt"), header=TRUE)
        Q2.ts.rec[[param]]               = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/Q2.ts.txt"), header=TRUE)
        Q97.ts.rec[[param]]              = read.table(file=paste0(dir.rec.segm.test.param[[param]] ,  "/Q97.ts.txt"), header=TRUE)
        nS.ok.rec[[param]]               = nrow(df.shift.times.rec[[param]])+1
        pdf.ts.rec[[param]]              = mcmc.seg.rec[[param]][, (nS.ok.rec[[param]]+1):(2*nS.ok.rec[[param]]-1)]
        shift.times.recessions[[param]]  = read.table(file =paste0(dir.rec.segm.test.param[[param]] , "/df.shift.times.txt"), header = TRUE)
        data.annotate.recess[[param]]    = data.frame(t      = shift.times.recessions[[param]]$ts.res,
                                                      tstart = shift.times.recessions[[param]]$Q2.ts,
                                                      tend   = shift.times.recessions[[param]]$Q97.ts,
                                                      t.adj  = shift.times.recessions[[param]]$ts.res,
                                                      t.real = shift.times.recessions[[param]]$ts.real)
        # plot dates of shifts:
        if (plot.dates == TRUE){
          dates.recess = 0
          for (dat in 1:length(data.annotate.recess[[param]]$t.adj)) {
            # dates.recess[dat]   =  data.annotate.recess[[param]]$t.adj[dat] +  initial.time
            dates.recess[dat]     =  data.annotate.recess[[param]]$t.real[dat]  +  initial.time
          }
          RealDates.recess        = as.Date(dates.recess,  origin = date_origin)
          RealDates.recess.new    = strptime(as.character(RealDates.recess), "%Y-%m-%d" )
          RealDates.recess.newnew = format( RealDates.recess.new, "%d/%m/%Y")
        } else {
          RealDates.recess.newnew = NULL
        }
      #########
      } else {
      #########
        nS.ok.rec[[param]]                   = 1
        shift.times.recessions[[param]]      = list()
        pdf.ts.rec[[param]]                  = list()
        data.annotate.recess[[param]]        = list()
        data.annotate.recess.adjust[[param]] = list()
        RealDates.recess.newnew              = NULL
      } 
      
      
      

      # gaugings to plot in the stage record with the periods derived from segmentation:
      gaugings.df.recess[[param]] = read.table(file =paste0(dir.rec.segm.test.param[[param]],  "/data_with_periods.txt"), header = TRUE)
      # Quantiles for tot uncertainty for each parameter:
      Ysim[[param]]       = list();
      data.tmp.2[[param]] = list()
      for (jj in 1:nS.ok.rec[[param]]) {
        Ysim[[param]][[jj]] =0
        for (ii in 1: length(mcmc.seg.rec[[param]][,jj])) {
          Ysim[[param]][[jj]][ii] = mcmc.seg.rec[[param]][ii,jj] + 
            rnorm(1, 
                  mean = 0, 
                  sd = mcmc.seg.rec[[param]][, 2*nS.ok.rec[[param]]])
        }
        if (param == length(parameters)) {
          data.tmp.2[[param]][[jj]] = quantile(Ysim[[param]][[jj]], probs=c(0.025,0.975) ) #, na.rm=TRUE)
        } else {
          delete.mcmc = which(Ysim[[param]][[jj]] <0)
          data.tmp.2[[param]][[jj]] = quantile(Ysim[[param]][[jj]][-delete.mcmc], probs=c(0.025,0.975) ) #, na.rm=TRUE)
        }
      } 
      # Preparing vector with pdf of shift times :
      if (length(pdf.ts.rec[[param]]) > 0) {
        X1[[param]] = pdf.ts.rec[[param]]
        if (nS.ok.rec[[param]]==2) {
          X[[param]] = X1[[param]]
        } else {
          X[[param]] = data.frame(time= X1[[param]][,1], 
                                         ord = rep(1, length(X1[[param]][,1])))
          for (orderr in 1:(ncol(pdf.ts.rec[[param]]))) {
            X[[param]] = rbind(X[[param]], 
                                      data.frame(time= X1[[param]][,orderr],
                                                 ord = rep(orderr, length(X1[[param]][,1]))))
          }
        }
      }  
      
      # Second  plots: time series of parameters of the regression model with the shift times and b(t)
      #print(param)
      if (param == length(parameters)) {
        inddd = which((Data.segm.rec[[param]]$X97.5. < stage.limits[2]) & 
                        (Data.segm.rec[[param]]$X2.5. > stage.limits[1]))
        
      } else if (parameters[param] == "b1") {
        inddd = which((Data.segm.rec[[param]]$X97.5. < limits.Y.lambda[2]) & 
                        (Data.segm.rec[[param]]$X2.5. > limits.Y.lambda[1]))
      } else {
        inddd = which((Data.segm.rec[[param]]$X97.5. < limits.Y.alpha[2]) & 
                        (Data.segm.rec[[param]]$X2.5. > limits.Y.alpha[1]))
      } 
      dddf = Data.segm.rec[[param]][inddd,]
      
      dddf = Data.segm.rec[[param]]
      
      obs.uncertainty.y = TRUE
      limni.time.limits = c(df.limni$t_limni[1], tail(df.limni$t_limni,1))
      
      ##################################################################
      if ( nS.ok.rec[[param]] > 1) {
      ##################################################################
        seg.rec.plot[[param]]  <- ggplot()
        obs.uncertainty.y = TRUE
        if (obs.uncertainty.y == TRUE) {
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] +  
            # geom_errorbar(data = Data.segm.rec[[param]],
            #               aes(x=t, 
            #                   ymin= (mean - 2*stdev), 
            #                   ymax = (mean+2*stdev)),
            #               size = 0.1, width=20, col= "black")
            geom_errorbar(data = dddf,
                          aes(x    = t, 
                              ymin = X2.5., 
                              ymax = X97.5.), col= "gray30")
        }
        # time.adjust.before[[param]] = c(0, data.annotate.recess[[param]]$t.adj)
        # time.adjust.plus[[param]]   = c(data.annotate.recess[[param]]$t.adj, 
        #                                        limni.time.limits[2])
        
        time.adjust.before[[param]] = c(0, data.annotate.recess[[param]]$t.real)
        time.adjust.plus[[param]]   = c(data.annotate.recess[[param]]$t.real, limni.time.limits[2])
        
        
        if (plot.gamma.uncertainty ==TRUE){
          for (jj in 1:nS.ok.rec[[param]]){
            seg.rec.plot[[param]] =  seg.rec.plot[[param]] + 
              annotate("rect", 
                       xmin = time.adjust.before[[param]][jj], 
                       xmax = time.adjust.plus[[param]][jj], 
                       ymin = max(data.tmp.2[[param]][[jj]][1]), #, stage.limits[1]),
                       ymax = min(data.tmp.2[[param]][[jj]][2]), #, stage.limits[2]),
                       fill = "pink", alpha = 0.3)
          } 
        }
        
        round.choose = function(x, roundTo, dir=1){
          if (dir ==1){
            if (x >0){
              x+(roundTo - x %% roundTo)
            } else {
              x-(x%% roundTo)
            }
          } else if (dir==0){
            if (x >0){
               x-(x%% roundTo)
            } else {
              x+(roundTo - x %% roundTo)
            }
          }
        }
        
        
        seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
          # annotate("rect",
          #      xmin = time.adjust.before[[param]],
          #      xmax = time.adjust.plus[[param]],
          #      ymin = Q2.mu.rec[[param]]$x,
          #      ymax = Q97.mu.rec[[param]]$x, 
          #      fill ="red", alpha=0.2) +
          
          # geom_point(data = Data.segm.rec[[param]], aes(x = t, y = mean), size = 2, pch=1) +
          # geom_line(data = Data.segm.rec[[param]], aes(x = t, y = mean), size = 0.1, color = "blue")+
          geom_point(data = dddf, aes(x = t, y = mean), size = 3, pch=21, fill="gray30", color= "gray30")
        
        if (param == length(parameters)) {
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
            geom_segment(aes(
              x    = time.adjust.before[[param]],
              xend = time.adjust.plus[[param]],
              y    = mu.res.rec[[param]]$x,
              yend = mu.res.rec[[param]]$x),
              color ="red", size =1.5)
          
        } else if (parameters[param] == "b1") {
          geom_segment(aes(
            x    = time.adjust.before[[param]],
            xend = time.adjust.plus[[param]],
            y    = mu.res.rec[[param]]$x,
            yend = mu.res.rec[[param]]$x),
            color ="red", size =1.5)
        } else {
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
            geom_segment(aes(
              x    = time.adjust.before[[param]],
              xend = time.adjust.plus[[param]],
              y    = mu.res.rec[[param]]$x,
              yend = mu.res.rec[[param]]$x),
              color ="red", size =1.5)
        } 
        seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
          theme_light(base_size = 20) +
          theme(plot.title = element_text(hjust = 0.5))+
          scale_x_continuous(name = x.name, expand = c(0,0), limits = limni.time.limits) +
          # geom_segment(x = time.adjust.before[[param]],
          #              xend = time.adjust.plus[[param]],
          #              y = mu.res.rec[[param]]$x,
          #              yend = mu.res.rec[[param]]$x, 
          #              col="red") +
          #geom_vline(xintercept = read.res.rec$ts.res, col="black", lwd =0.5, linetype= "dashed")+
          # geom_vline(aes(xintercept = data.annotate.recess[[param]]$t.adj), linetype= "solid", lwd = 1.5, color="blue") +
          geom_vline(aes(xintercept = data.annotate.recess[[param]]$t.real), linetype= "solid", lwd = 1.5, color="blue") +
          coord_cartesian(clip = "off") +
          theme_light(base_size=20)+
          theme(axis.text         = element_text(size=15)
                ,axis.title       = element_text(size=20, face="plain")
                ,panel.grid.major = element_blank() #element_line(size=1.2)
                ,panel.grid.minor = element_blank()
                ,legend.text      = element_text(size=20)
                ,legend.title     = element_text(size=30)
                #,plot.margin     = unit(c(0.5,0.5,0.5,0.5),"cm")
                ,legend.key.size  = unit(1.5, "cm")
                ,legend.position  = "none"
                ,plot.margin      = unit(c(2, 0.5, 0, 1),"cm")
                ,axis.title.y     = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=20, face="bold"))
        
      ##############
      } else {
      #############
        seg.rec.plot[[param]] = ggplot()
        if (obs.uncertainty.y == TRUE) {
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] +  
            # geom_errorbar(data = Data.segm.rec[[param]],
            #               aes(x=t, 
            #                   ymin= (mean - 2*stdev), 
            #                   ymax = (mean+2*stdev)),
            #               size = 0.3, width=40, col= "black")
            geom_errorbar(data = dddf,
                          aes(x    = t, 
                              ymin = X2.5., 
                              ymax = X97.5.),
                          #size = 0.2, width=80,
                          col= "gray30")
        }
        seg.rec.plot[[param]] =  seg.rec.plot[[param]] +  
          #geom_line(data = Data.segm.rec[[param]], aes(x = t, y = mean), size = 0.1, color = "blue")+
          geom_point(data = Data.segm.rec[[param]], aes(x = t, y = mean),  size = 3, pch=21, fill="gray30", color= "gray30") +
          theme_light(base_size = 20) +  
          coord_cartesian(clip = "off") +
          theme(plot.title = element_text(hjust = 0.5))
        
        if (plot.gamma.uncertainty ==TRUE){
          seg.rec.plot[[param]] =  seg.rec.plot[[param]]+
            annotate("rect",  
                     xmin = Data.segm.rec[[param]]$t[1], 
                     xmax = tail(Data.segm.rec[[param]]$t,1),
                     ymin = data.tmp.2[[param]][[1]][1],  #max(data.tmp.2[[param]][[1]][1],  stage.limits[1]),
                     ymax = data.tmp.2[[param]][[1]][2],  #min(data.tmp.2[[param]][[1]][2],  stage.limits[2]),
                     fill = "pink", alpha=0.3)
        }
        # annotate("rect", 
        #          xmin= Data.segm.rec[[param]]$t[1], 
        #          xmax= tail(Data.segm.rec[[param]]$t,1),
        #          ymin= Q2.mu.rec[[param]]$x,
        #          ymax= Q97.mu.rec[[param]]$x, 
        #          fill="red", alpha=0.2) 
        
        if (param == length(parameters)) {
        # plot segments means:
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
              geom_segment(aes(
                x    = Data.segm.rec[[param]]$t[1],
                xend = tail(Data.segm.rec[[param]]$t,1),
                y    = mu.res.rec[[param]]$x,
                yend = mu.res.rec[[param]]$x),
                color ="red", size =1.5)
          
        } else if (parameters[param] == "b1") {
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
            geom_segment(aes(
              x    = Data.segm.rec[[param]]$t[1],
              xend = tail(Data.segm.rec[[param]]$t,1),
              y    = mu.res.rec[[param]]$x,
              yend = mu.res.rec[[param]]$x),
              color ="red", size =1.5)
        } else {
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
            geom_segment(aes(
                     x    = Data.segm.rec[[param]]$t[1],
                     xend = tail(Data.segm.rec[[param]]$t,1),
                     y    = mu.res.rec[[param]]$x,
                     yend = mu.res.rec[[param]]$x),
                     color ="red", size =1.5)
        } 
        # geom_segment(x     = Data.segm.rec[[param]]$t[1], 
        #              xend  = tail(Data.segm.rec[[param]]$t,1),
        #              y     = mu.res.rec[[param]]$x, 
        #              yend  = mu.res.rec[[param]]$x, 
        #              color = "red") +
        seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
        scale_x_continuous(name = x.name, expand = c(0,0), limits =limni.time.limits)+
        theme_light(base_size=20) +
        theme(axis.text         = element_text(size=15)
              ,axis.title       = element_text(size=20, face="plain")
              ,panel.grid.major = element_blank() #element_line(size=1.2)
              ,panel.grid.minor = element_blank()
              ,legend.text      = element_text(size=20)
              ,legend.title     = element_text(size=30)
              #,plot.margin     = unit(c(0.5,0.5,0.5,0.5),"cm")
              ,legend.key.size  = unit(1.5, "cm")
              ,legend.position  = "none"
              ,plot.margin      = unit(c(2, 0.5, 0, 1),"cm")
              ,axis.title.y     = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=20, face="bold"))
      }
      
      ####################################################################
      if (param == length(parameters)) {
        if (plot.b.from.gaugings == TRUE) {
          # uncomment this part if you want to plot the segments of "b1" obtained from gaugings!!!!!!
          seg.rec.plot[[param]] =  seg.rec.plot[[param]] + 
            # geom_segment(mapping= aes(x = c(limits.X[1], bt.from.gaugingsss$t.shift.for.b$treal), 
            #                           y = bt.from.gaugingsss$bt1.df$maxpost*100, 
            #                           xend = c(bt.from.gaugingsss$t.shift.for.b$treal, limits.X[2]),  
            #                           yend = bt.from.gaugingsss$bt1.df$maxpost*100), 
            #              color = "black", 
            #              size = 0.5,
            #              linetype = "dashed") +
            # geom_segment(mapping= aes(x = c(limits.X[1], bt.from.gaugingsss$t.shift.for.b$treal), 
            #                           y = bt.from.gaugingsss$bt2.df$maxpost*100, 
            #                           xend = c(bt.from.gaugingsss$t.shift.for.b$treal, limits.X[2]),  
            #                           yend = bt.from.gaugingsss$bt2.df$maxpost*100), 
          #              color = "blue", 
          #              size = 0.5,
          #              linetype = "dashed")+
          theme(axis.text         = element_text(size=15)
                ,axis.title       = element_text(size=20, face="plain")
                ,panel.grid.major = element_blank() #element_line(size=1.2)
                ,panel.grid.minor = element_blank()
                ,legend.text      = element_text(size=20)
                ,legend.title     = element_text(size=30)
                #,plot.margin     = unit(c(0.5,0.5,0.5,0.5),"cm")
                ,legend.key.size  = unit(1.5, "cm")
                ,legend.position  = "none"
                ,plot.margin      = unit(c(2, 0.5, 0, 1),"cm")
                ,axis.title.y     = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=20))
          
          # geom_rect(mapping = aes(xmin = ts.before.gaug, 
          #                         xmax = ts.plus.gaug, 
          #                         ymin = bt2.df$X2.5.,
          #                         ymax = bt2.df$X97.5.), 
          #           fill="red", alpha=0.3)+
          # geom_vline(aes(xintercept = ts.gaug$t.adj), color = "black",
          #            lwd =0.7, linetype = "dotdash")
        }
        
        # fix the axis scale y for the asymptotic level:
        seg.rec.plot[[param]]   =  seg.rec.plot[[param]] + 
          scale_y_continuous(name =  parameters.names[param], 
                           #limits =c(asymptote.limits[1], asymptote.limits[2]), 
                           expand = c(0.01, 0.01))
      ###########
      }  else if (parameters[param] == "b1") { 
      ###########
        seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
          scale_y_continuous(name =  parameters.names[param],
                             limits =c(limits.Y.lambda[1], limits.Y.lambda[2]),  ##############    change this !!!!!!!!!!!!!!!!!!!!!!!!
                             expand = c(0.01,0.01))
        
      #########  
      }  else {
      #########
        # automatic y scale axis:
        seg.rec.plot[[param]] =  seg.rec.plot[[param]] +
          scale_y_continuous(name =  parameters.names[param],
                             limits =c(limits.Y.alpha[1], limits.Y.alpha[2]),  ##############    change this !!!!!!!!!!!!!!!!!!!!!!!!
                             expand = c(0.01,0.01))
      }
      # TO VERIFY::
      # seg.rec.plot = seg.rec.plot+
      #   coord_cartesian(clip = 'off')+
      #   theme_bw(base_size=20)+
      #   theme(axis.text=element_text(size=15)
      #         ,axis.title=element_text(size=20,face="plain")
      #         ,panel.grid.major=element_blank()
      #         ,panel.grid.minor=element_blank()
      #         ,legend.text=element_text(size=20)
      #         ,legend.title=element_text(size=30)
      #         ,legend.key.size=unit(1.5, "cm")
      #         ,legend.position="none"
      #         ,plot.margin=unit(c(2, 0.5, 2, 2),"cm")
      #         ,axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l =0 ))
      #         ,axis.title.x = element_blank())
      # pdf(paste0(dir.segm.recessions,"/Figure_segment_model_",model.names[mod],"_chi",chi.test,"par_",parameters[param],".pdf"), 15,10 ,useDingbats=F)
      # print(seg.rec.plot[[param]])
      # dev.off()
      
    
    # combine the plots with the parameters time series:
    if ((any(model.names =="2expWithAsympt_bis"))|(any(model.names =="3expWithAsympt_bis"))|(any(model.names =="expexp_bis"))|(any(model.names =="hyperb_bis"))|(any(model.names =="Coutagne_bis"))){
      if (length(parameters)==2){
        # t.plot3 = plot_grid(seg.rec.plot[[1]], NULL, NULL, seg.rec.plot[[2]],
        #                         ncol = 1, nrow = 4, rel_heights = c(1, 1, 1, 1))
        t.plot3 = seg.rec.plot[[2]]
        
      } else if (length(parameters)==3){
        t.plot3 = seg.rec.plot[[3]]
      }
    } else { 
      t.plot3 = seg.rec.plot[[2]]
      
    } 
    
    
    
    
    
    
    
    
    
    ###########################################################################################
    #PLot of stage record with the shifts:
    print("plot 3")
    message("- Plotting stage record with gaugings and the river bed estimates with uncertainty")
    # filter the time series of stage record removing the long periods with missing data (putting a NA instead):
    dt_limni         = 0 
    new_NA_limni     = 0   
    t_limni_filtered = df.limni$t_limni
    h_limni_filtered = df.limni$h_limni
    
    for (tt in 2:length(df.limni$t_limni)){
        dt_limni[tt] =  df.limni$t_limni[tt] - df.limni$t_limni[tt-1]
    }
    
    for (tt in 2:length(df.limni$t_limni)){
        if (dt_limni[tt] > 10*mean(dt_limni)){
          t_limni_filtered[(tt+ new_NA_limni): (length(t_limni_filtered)+1)] <- c(t_limni_filtered[tt+ new_NA_limni] -10, 
                                                                                  t_limni_filtered[(tt+ new_NA_limni) : length(t_limni_filtered)])
          h_limni_filtered[(tt+ new_NA_limni): (length(h_limni_filtered)+1)] <- c(NA, 
                                                                                  h_limni_filtered[(tt+ new_NA_limni) : length(h_limni_filtered)])
          new_NA_limni = new_NA_limni +1 
        }
    } 
    df.limni_filtered = data.frame(t_limni = t_limni_filtered, h_limni = h_limni_filtered)
      
    # stage record plot:
    t.plot <- ggplot()
      # if (is.null(df.limni)==FALSE) {
      #   t.plot= t.plot + 
      #     geom_line(data = df.limni_filtered, aes(x = t_limni, y = h_limni), color = "gray70",size = 0.2)+
      #     scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      # } else {
      #   t.plot= t.plot + 
      #     scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0,tail(g.s$X.tP.,1)))
      # }  
      # 
      
    t.plot <- ggplot()+
      scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limni.time.limits) +
      scale_y_continuous(name=limni.labels[2], expand = c(0,0)) +
                        # limits = c(grid_limni.ylim[1],   grid_limni.ylim[2]), 
                        # breaks = seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3])) +
      ylab(limni.labels[2]) +
      coord_cartesian(clip = 'off') +
      geom_line(data = df.limni_filtered, aes(x = t_limni, y = h_limni), color = "gray60",size = 0.2) + # LIMNI !
      geom_point(data=gaugings.df.recess[[length(parameters)]],                 # gaugings
                 aes(x = t , y= h), size = 3, color = gaugings.df.recess[[length(parameters)]]$color) +
      theme_light(base_size=20) +
      theme(axis.text.x.top   = element_blank()
            ,axis.title       = element_text(size=20,face="plain")
            ,axis.title.x     = element_blank()  
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,legend.text      = element_text(size=20)
            ,legend.title     = element_text(size=30)
            ,legend.key.size  = unit(1.5, "cm")
            ,legend.position  = "none"
            ,plot.margin      = unit(c(2, 0.5, 0, 1),"cm")
            ,axis.title.y     = element_text(size=20, margin = margin(t = 0, r = 55, b = 0, l = 0))
            ,axis.text.x      = element_blank()
            ,axis.line.x      = element_blank())
    
    ###################
    if ( nS.ok.rec[[length(nS.ok.rec)]] >1) {
      if (plot.gamma.uncertainty == TRUE){
        for (jj in 1:nS.ok.rec[[length(nS.ok.rec)]]){
          t.plot  = t.plot  + 
            annotate("rect", 
                     xmin = time.adjust.before[[length(nS.ok.rec)]][jj], 
                     xmax = time.adjust.plus[[length(nS.ok.rec)]][jj], 
                     ymin = max(data.tmp.2[[length(nS.ok.rec)]][[jj]][1]/100), #stage.limits[1]/100),
                     ymax = min(data.tmp.2[[length(nS.ok.rec)]][[jj]][2]/100), #stage.limits[2]/100),
                     fill = "pink", alpha = 0.3)
        }
      }
      t.plot  = t.plot  + 
        # annotate("rect", 
        #          xmin = time.adjust.before[[length(nS.ok.rec)]], 
        #          xmax = time.adjust.plus[[length(nS.ok.rec)]], 
        #          ymin = unlist(Q2.mu.rec[[length(nS.ok.rec)]])/100,
        #          ymax = unlist(Q97.mu.rec[[length(nS.ok.rec)]])/100, 
        #          fill ="red", alpha=0.2)+
        #geom_vline(aes(xintercept = data.annotate.recess.adjust$t.adj), color = "blue", lwd =0.3, linetype = "dashed")+
        geom_segment(mapping=aes(x    = time.adjust.before[[length(parameters)]] , 
                                 y    = unlist(mu.res.rec[[length(parameters)]])/100, 
                                 xend = time.adjust.plus[[length(parameters)]], 
                                 yend = unlist(mu.res.rec[[length(parameters)]])/100),
                     col="red", size=1.5)
      #################
    } else {
      #################  
      if (plot.gamma.uncertainty ==TRUE){
        t.plot  = t.plot + annotate("rect",  
                                    xmin = Data.segm.rec[[length(nS.ok.rec)]]$t[1], 
                                    xmax = tail(Data.segm.rec[[length(nS.ok.rec)]]$t,1),
                                    ymin = max(data.tmp.2[[length(nS.ok.rec)]][[1]][1]/100), #stage.limits[1]/100),
                                    ymax = min(data.tmp.2[[length(nS.ok.rec)]][[1]][2]/100), #stage.limits[2]/100),
                                    fill = "pink", alpha=0.3)
      }
      # t.plot  = t.plot  + 
      #   annotate("rect", 
      #            xmin= Data.segm.rec[[length(nS.ok.rec)]]$t[1], 
      #            xmax= tail(Data.segm.rec[[length(nS.ok.rec)]]$t,1),
      #            ymin= unlist(Q2.mu.rec[[length(nS.ok.rec)]])/100,
      #            ymax= unlist(Q97.mu.rec[[length(nS.ok.rec)]])/100, fill="red", alpha=0.1)
    }
    if (plot.dates == TRUE){
      if (!is.null(RealDates.recess.newnew)){
        t.plot  <- t.plot  +
          annotate("text", 
                   x     = data.annotate.recess[[param]]$t.adj,
                   y     = grid_limni.ylim[1], vjust = 0 , 
                   label = RealDates.recess.newnew, color = "blue", size=4)
      }
    }
    # ##########################################################################################
    # # PLot of the pdf of the shifts:
    # t.plot2  <- ggplot() +
    #   scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,0.5))+
    #   xlab(limni.labels[1])+ 
    #   coord_cartesian(clip = 'off') + 
    #   geom_hline(yintercept = c( -0.5,0, 0.5), color="darkgray", linetype="dashed", size = 0.5)+
    #   theme_light(base_size=15)+
    #   theme(axis.text=element_text(size=10)
    #         ,axis.title=element_text(size=15,face="plain")
    #         ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
    #         ,legend.text=element_text(size=20),legend.title=element_text(size=30)
    #         ,legend.key.size=unit(1.5, "cm"),legend.position="none"
    #         ,plot.margin=unit(c(0, 0.5, 0, 1),"cm")
    #         ,axis.title.y = element_text(margin = margin(t = 0, r = 43, b = 0, l = 0))
    #         ,axis.text.y=element_blank()
    #         ,axis.ticks.y = element_blank()) +
    #   scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limits.X)
    # 
    # if (is.null(data.annotate.off)==FALSE) {
    #   t.plot2  <- t.plot2  +
    #     geom_point(data = data.annotate.off, 
    #                aes(x = xeffect, y = -1), color= "black", size = 3, shape =4, stroke=1.5 )
    # }
    # if (!is.null(pdf.ts.rec[[length(parameters)]])) {
    #   if (nS.ok.rec[[length(nS.ok.rec)]] ==2) {
    #     t.plot2  = t.plot2  + 
    #       geom_density(aes(x= X[[length(parameters)]], ..scaled.. /2),
    #                    fill="blue", 
    #                    colour=NA, alpha=0.4)
    #   } else {
    #     t.plot2  = t.plot2  +
    #       geom_density(aes(x= X[[length(parameters)]]$time, 
    #                        ..scaled.. /2, 
    #                        group =X[[length(parameters)]]$ord),  
    #                    fill= "blue", 
    #                    colour=NA, alpha=0.4)
    #   }
    #   t.plot2  = t.plot2  + 
    #     geom_segment(data = data.annotate.recess.adjust[[length(parameters)]], 
    #                  aes(x = t.adj, y = -0.5, yend =0, xend=t.adj ),
    #                  size = 0.9, color ="blue", linetype = "solid")
    # }
    # 
    
    ##########################################################################################
    print("plot 4")
    message("- Plotting pdf of the detected shift times, if any")
    # PLot of the pdf of the shifts:
    t.plot2  <- ggplot() +
      scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 0.5))+
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off') + 
      geom_hline(yintercept = c( 0.5), color="darkgray", linetype="dashed", size = 0.5)+
      theme_light(base_size=20)+
      theme(axis.text         = element_text(size=15)
            ,axis.title       = element_text(size=20,face="plain")
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,legend.text      = element_text(size=20)
            ,legend.title     = element_text(size=30)
            ,legend.key.size  = unit(1.5, "cm")
            ,legend.position  = "none"
            ,plot.margin      = unit(c(0, 0.5, 0, 1),"cm")
            ,axis.title.y     = element_text(margin = margin(t = 0, r = 68, b = 0, l = 0))
            ,axis.text.y      = element_blank()
            ,axis.ticks.y     = element_blank()) +
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =limni.time.limits)
    if ( nS.ok.rec[[length(nS.ok.rec)]] >1) { 
      if (nS.ok.rec[[length(nS.ok.rec)]] ==2) {
        t.plot2  = t.plot2  + 
          geom_density(aes(x= X[[length(parameters)]], ..scaled.. /2),
                       fill="blue", 
                       colour=NA, alpha=0.2)
      } else {
        t.plot2  = t.plot2  +
          geom_density(aes(x= X[[length(parameters)]]$time, 
                           ..scaled.. /2, 
                           group =X[[length(parameters)]]$ord),  
                       fill= "blue", 
                       colour=NA, alpha=0.2)
      }
    } else {
      t.plot2  = t.plot2  +
      annotate("text",
               x = limni.time.limits[2]/2,  
               y = 0.2, 
               label= "No shift detected", color = "gray", size=6)
    }
    
    
    ####################################################################################################################
    # ADD GAUGINGS SEGMENTATION RESULTS
    # if (plot_ts_gaugings==TRUE){
    #   if (mod==1){
    #     t.plot4 = plot_grid(reg.pool.plot, t.plot3, t.plot , t.plot2 , 
    #                                ncol = 1, nrow = 4, rel_heights = c(1.1, 1, 0.9, 0.15),  labels = c('a) ', 
    #                                                                                                    'b)' ,
    #                                                                                                    'c)', 
    #                                                                                                    ''), 
    #                                label_size = 25, label_fontface = "bold")
    #   } else {
    #     t.plot4 = plot_grid(reg.pool.plot, t.plot3, t.plot , t.plot2 , 
    #                                ncol = 1, nrow = 4, rel_heights = c(1.1, 1, 0.9, 0.15),  labels = c('', '' ,'', ''), 
    #                                label_size = 25, label_fontface = "plain")
    #   }
    #   
    # }

    
  ### COMBINING THREE PLOTS: a) recession estimates b) segmentation of asymptotic stage series c) stage record
  path.final.plot = paste0(dir.rec.segm.test,"/Recession_analysis_results.pdf")
  message("- Combining plots in unique plot + saving to:")
  print(path.final.plot)

  t.plot4 = plot_grid(reg.pool.plot, 
                             t.plot3,
                             t.plot,
                             t.plot2, 
                             ncol = 1, nrow = 4, 
                             rel_heights = c(1, 1, 0.9, 0.2),
                             labels = c('a) ',  'b)' ,   'c)' , ''),
                             label_size = 25, 
                             label_fontface = "bold")
  t.plot4
  pdf(path.final.plot, 15,16 ,useDingbats=F)
  print(t.plot4)
  dev.off()
  
  ########################################################################################################
  # # plot of the time series with detected shift times.
  # plot.time.shifts.gaugings(dir                       = dir.rec.pool.test, 
  #                           g.s                       = gaugings.df.recess[[length(parameters)]], 
  #                           data.annotate.off         = data.annotate.off, 
  #                           data.annotate.gaug        = data.annotate.gaug.1,
  #                           t_Gaug                    = t_Gaug, 
  #                           h_Gaug                    = h_Gaug, 
  #                           df.limni                  = stage.record, 
  #                           limni.labels              = limni.labels,
  #                           grid_limni.ylim           = grid_limni.ylim,
  #                           plot.shift.times.on.limni = FALSE,
  #                           pdf.ts                    = pdf.ts.results.1,
  #                           dates                     = shift.times.date) 
  ########################################################################################################
}





































############################################################################################################
plot.all.recessions = function(dir.rec.pool.test,  
                               rec.model,
                               stage.limits,
                               limits.x.recess,
                               stage.scale.shift,
                               plot.recession.uncert){
#############################################################################################################
  # Recession models studied :
  #***************************
  # 1) 1exp =>                 h(t) = alpha1 *exp(-lambda1* t) + beta
  # 2) 2exp =>                 h(t) = alpha1 *exp(-lambda1* t) + alpha1 *exp(-lambda1* t) + beta
  # 3) 3exp =>                 h(t) = alpha1 *exp(-lambda1* t) + alpha2 *exp(-lambda2* t) + alpha3*exp(-lambda3* t) + beta
  # 4) double exp (horton) =>  h(t) = alpha1 *exp(-lambda1* t^n1) + beta
  # 5) hyperb   =>             h(t) = alpha1 /(1+lambda1*t)^2
  # 6) Boussinesq (1903) =>    h(t) = alpha1 *(1 + lambda1*t)^(-2)
  # 7) Coutagne
  # model.names = c("expexp", "hyperb", "expexp", "hyperb", "expexp", "hyperb")
  
  
  # recess model with tot uncertainty:
  #********************************************************************************
  Recession.model = function(theta, model, t){ 
    #******************************************************************************
    h=0*t
    if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[8],theta[9]))
    } else if (model =="2expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="2expWithAsympt_bis"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="2expWithAsympt_rel"){
      h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
        theta[5] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[6],theta[7]))
    } else if (model =="1expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[4],theta[5]))
    } else if ((model =="expexp")|(model =="expexp_bis")){
      h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } else if ((model =="hyperb")|(model =="hyperb_bis")){
      h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
      res.h = sapply(h, 
                     function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
                     theta=c(theta[5],theta[6]))
    } 
    return(res.h)
  }
  
  # recess model for maxpost:
  #*************************************************************************************
  Recess.Maxpost.model = function(theta, model, t){
  #*************************************************************************************
    h=0*t
    if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5]*exp(-theta[6]*t) +
        theta[7] 
    }  else if (model =="2expWithAsympt"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
    }  else if (model =="2expWithAsympt_bis"){
      h = theta[1]*exp(-theta[2]*t) +
        theta[3]*exp(-theta[4]*t) + 
        theta[5] 
    }  else if (model =="2expWithAsympt_rel"){
      h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
        theta[5] 
    } else if ((model =="1expWithAsympt")|(model =="1expWithAsympt_bis")){
      h = theta[1]*exp(-theta[2]*t) + theta[3] 
    } else if ((model =="expexp")|(model =="expexp_bis")){
      h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
    } else if ((model =="hyperb")|(model =="hyperb_bis")){
      h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
    } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
    }
    return(h)
  }
  
  #initialisation of the lists of plot objects:
  reg.pool.plot=NULL; reg.pool.plot2=NULL;  seg.rec.plot=NULL; t.plot=NULL; t.plot2 =NULL; t.plot3 =NULL
  t.plot4 = NULL; title.model=NULL; dir.rec.segm.test.param=NULL;
  tau.results.df.rec =NULL; mu.results.df.rec=NULL; gamma.results.df.rec=NULL; df.shift.times.rec=NULL
  df.shift.times.plus.rec=NULL; ts.res.before.rec=NULL; ts.res.rec=NULL; ts.res.plus.rec=NULL; Q2.ts.rec=NULL
  Q97.ts.rec=NULL; Q2.mu.rec=NULL; mu.res.rec=NULL; Q97.mu.rec=NULL; Data.segm.rec=NULL; nS.ok.rec=NULL;
  mcmc.seg.rec=NULL; pdf.ts.rec=NULL; gamma_segm_recess=NULL; gaugings.df.recess=NULL;
  shift.times.recessions=NULL; data.annotate.recess=NULL; data.annotate.recess.adjust=NULL; 
  X=NULL; X1=NULL; data.tmp= NULL;data.tmp.2= NULL; Ysim=NULL; time.adjust.before =NULL; 
  time.adjust.plus=NULL; parameters =NULL; parameters.names =NULL;  model.title =0
  #one selected model only: 
  #####
  mod=1; model.names=0
  model.names[mod] = rec.model[mod]
  # read results of cooked mcmc:  
  data.MCMC.cooked  = as.matrix(read.table(paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"), header=TRUE,dec=".", sep=""))
  data.MCMC.MaxPost = as.numeric(read.table(paste0(dir.rec.pool.test,"/Results_Summary.txt"), row.names=1,dec=".",sep="", skip = 16))
  colfunc           = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))
  summary.rec       = read.table(file=paste0(dir.rec.pool.test,"/Results_Summary.txt"), header=TRUE)
  mcmc.rec          = read.table(file=paste0(dir.rec.pool.test,"/Results_MCMC_cooked.txt"), header=TRUE)
  residuals.rec     = read.table(file=paste0(dir.rec.pool.test,"/Results_Residuals.txt"), header=TRUE)
  curves.data.rec   = read.table(file=paste0(dir.rec.pool.test,"/Curves_Data.txt"), header=TRUE)
  
  ######################################################################################################
  # Recessions regression plot:
  ######################################################################################################
  nsample      = length(data.MCMC.cooked[,1])
  tgrid        = seq(0, round(max(residuals.rec[,1]), 1), 0.5)    #seq(limits.x.recess[1],  limits.x.recess[2],  0.5) 
  Ncurves.pool = tail(curves.data.rec$Period,1)
  #Initialisation:
  #########################################
  if (model.names[mod] =="3expWithAsympt"){
  #########################################
    model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2 e^{-\\lambda_2 t} +  \\alpha_3 e^{-\\lambda_3 t} + \\beta^{(k)}$")
    b1.var =FALSE
    MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period 
    for(i in 1:Ncurves.pool){
      col.num = c( i, Ncurves.pool+1,                                  #a1, b1,
                   Ncurves.pool +2, Ncurves.pool + 3,                    #a2, b2
                   Ncurves.pool +4, Ncurves.pool + 5,                      #a3, b3,
                   Ncurves.pool +5 + i,                            #a4 (asymptotic level parameter)
                   Ncurves.pool*2 + 5+1, Ncurves.pool*2 + 5+2)  #gamma1, gamma2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    }
    
  } else if (model.names[mod] =="3expWithAsympt_bis"){
    #####################################################
    model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2^{(k)} e^{-\\lambda_2 t} + \\alpha_3 e^{-\\lambda_3 t} + \\beta^{(k)}$")
    MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period 
    for(i in 1:Ncurves.pool){
      col.num = c( i, Ncurves.pool+1,                               #a1(var), b1,
                   Ncurves.pool +1 + i, 2*Ncurves.pool + 2,         #a2(var), b2
                   2*Ncurves.pool + 3, 2*Ncurves.pool + 4,          #a3, b3,
                   2*Ncurves.pool + 4 + i,                          #a4(var) (asymptotic level parameter)
                   3*Ncurves.pool + 4 + 1, 3*Ncurves.pool + 6)      #gamma1, gamma2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    }  
    
    
  }  else if (model.names[mod] =="2expWithAsympt_bis"){
    #####################################################
    model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2^{(k)} e^{-\\lambda_2 t}  + \\beta^{(k)}$")
    MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=8) # 7 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=8)      # 7 param + # of period 
    for(i in 1:Ncurves.pool){
      col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                   Ncurves.pool + 1 + i, Ncurves.pool*2 +1 + 1,   #a2, b2
                   Ncurves.pool*2 +1+1+i,                       #a3,
                   Ncurves.pool*3 + 2 + 1, Ncurves.pool*3 +2+2)  #gamma1, gamma2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    }
    
    
  }  else if ((model.names[mod] =="2expWithAsympt")|(model.names[mod] =="2expWithAsympt_rel")){
    ##############################################################################################
    model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t} + \\alpha_2 e^{-\\lambda_2  t} + \\beta^{(k)}$")
    MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=8) # 7 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=8)      # 7 param + # of period 
    for(i in 1:Ncurves.pool){
      col.num = c( i, Ncurves.pool+1,                             #a1, b1,
                   Ncurves.pool +2, Ncurves.pool +3,              #a2, b2
                   Ncurves.pool + 3 + i,                           #a3,
                   Ncurves.pool*2 + 3 + 1, Ncurves.pool*2 +3 +2)  #gamma1, gamma2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    }
    
    
  } else if (model.names[mod] =="1expWithAsympt"){
    ################################################
    model.title[mod] = TeX("$h(t,k) = \\alpha^{(k)} \\; e^{-\\lambda t} + \\beta^{(k)}$")
    MCMC.save    = matrix(NA, nrow=Ncurves.pool*nsample, ncol=6) # 5 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=6)      # 5 param + # of period 
    for(i in 1:Ncurves.pool){
      col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                   Ncurves.pool +1+ i,                      #a2
                   Ncurves.pool*2 +1+1, Ncurves.pool*2 +1 + 2)  #gamma1, gamma2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    }
    
    
  } else if ((model.names[mod] =="expexp")|(model.names[mod] =="hyperb")|(model.names[mod] =="Coutagne")) {
    ##########################################################################################################
    model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t}^{\\eta t} + \\beta^{(k)}$")
    MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period 
    for(i in 1:Ncurves.pool){
      col.num = c( i, Ncurves.pool + 1,                          # a1(k), b1,
                   Ncurves.pool + 2,                             # n1
                   Ncurves.pool + 2 + i,                         # a2(k) (asymptotic level parameter)
                   Ncurves.pool*2 + 2 +1, Ncurves.pool*2 +2 +2)  # gamma1, gamma2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    }
    
    
  } else if ((model.names[mod] =="expexp_bis")|(model.names[mod] =="hyperb_bis")|(model.names[mod] =="Coutagne_bis")){
    ####################################################################################################################
    model.title[mod] = TeX("$h(t,k) = \\alpha_1^{(k)} e^{-\\lambda_1 t}^{\\eta^{(k)} t} + \\beta^{(k)}$")
    MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
    MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period 
    for(i in 1:Ncurves.pool){
      col.num = c( i, Ncurves.pool + i,                     # a1(k), b1(k),
                   Ncurves.pool*2 + 1,                      # n1
                   Ncurves.pool*2 + 1 + i,                  # a2 (asymptotic level parameter)
                   Ncurves.pool*3 + 2, Ncurves.pool*3 +3)   # gamma1, gamma2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
    } 
  }
  
  #######################################################################################################
  # Apply recession model (total uncertainty and maxpost):
  message("- Plotting all Recession curves !!!  Wait ... "); flush.console()
  Rec.Post    = apply(MCMC.save,    MARGIN=1, Recession.model,      model=model.names[mod],  t=tgrid) #add structural error:
  Rec.MaxPost = apply(MaxPost.save, MARGIN=1, Recess.Maxpost.model, model=model.names[mod],  t=tgrid) # Maximum posterior 
  
  # Quantiles: 
  List.Rec.quants = list(NULL)
  for(i in 1:Ncurves.pool){
    data.tmp = apply(Rec.Post[,(nsample*(i-1)+1):(nsample*i)], 
                     MARGIN=1, quantile, probs=c(0.025, 0.975),  na.rm=TRUE)
    List.Rec.quants[[i]] = data.frame(cbind(tgrid, 
                                            t(data.tmp - stage.scale.shift),
                                            Rec.MaxPost[,i] - stage.scale.shift))  #!!!!!!!!!!!!!! 
    colnames(List.Rec.quants[[i]]) = c("t", "inf", "sup", "maxpost")
  }
  
  #prepare plot:
  inter.per=seq(1,Ncurves.pool,1) 
  nRec = length(inter.per)
  palette.per = colfunc(Ncurves.pool)
  palette.per = palette.per[1:Ncurves.pool]
  data.plot.Rec = data.frame(do.call("rbind", 
                                     List.Rec.quants[inter.per]), 
                             Period=rep(inter.per,each=length(tgrid)))
  inter.null=which(data.plot.Rec$maxpost==0)
  data.rec.obs = read.table(paste0(dir.rec.pool.test,"/Curves_Data.txt"),
                            header=TRUE,dec=".",sep="") # Gaugings loading
  
  ylim.wind = c(stage.limits[1],    stage.limits[2])
  xlim.wind = c(limits.x.recess[1], limits.x.recess[2])
  
  pos.num =function(x.int){
    inter=which(x.int==data.rec.obs$period);
    return(inter)}
  inter.rec.obs =unlist(sapply(inter.per, pos.num), recursive = TRUE, use.names = TRUE)
  data.Rec = data.plot.Rec #[-inter.null,]
  data.Rec$Period = factor(data.Rec$Period)
  data.rec.obs$Period = as.factor(data.rec.obs$Period)
  write.table( data.Rec, paste0(dir.rec.pool.test,"/Rec_SPD_env.txt"), sep ="\t", row.names=FALSE)
  PltData.rec <- data.Rec[data.Rec$inf > -200,] ########## !!!!!!!!!!!!!!!!!!!!!!!!!!!!! update this !!!!!
  
  
  ##############################################
  reg.pool.plot[[mod]] = ggplot(PltData.rec)
  if (plot.recession.uncert == TRUE) {
    reg.pool.plot[[mod]]  = reg.pool.plot[[mod]] +
    geom_smooth(aes(x=t,
                    y=maxpost,
                    ymax=sup,
                    ymin=inf,
                    group=Period,
                    fill=Period),
                size=0.1, stat='identity', alpha=0.1)
  }
  reg.pool.plot[[mod]] = reg.pool.plot[[mod]] +
    geom_path(aes(x = t,   
                  y = maxpost,
                  group  = Period, 
                  colour = Period), size=0.3) +
    ### recession curve obs:
    geom_linerange(aes(x = time, 
                       ymax  = (h + 2*uh)  - stage.scale.shift, 
                       ymin  = (h-2*uh)    - stage.scale.shift, 
                       colour= Period), data = data.rec.obs, size=0.5)+
    geom_point(aes(x=time, 
                   y=h-stage.scale.shift,
                   colour=Period), 
               data=data.rec.obs, shape=16, size=2)+
    ### Labels
    xlab("Recession time [days]") +
    ylab("Stage h [cm]") +
    labs(colour = "Period") +
    scale_colour_manual(values = palette.per)+
    scale_fill_manual(values = palette.per)+
    #scale_y_continuous(breaks=breaks.wind,labels=labels.wind)+
    #coord_cartesian(ylim=ylim.wind, xlim=xlim.wind)+
    scale_x_continuous(expand=c(0,0))+ #, limits = xlim.wind) +
    scale_y_continuous(expand=c(0,0)) +
    ### Theme
    theme_light(base_size=25)+
    theme(axis.text         = element_text(size=10)
          ,axis.title       = element_text(size=15, face="plain")
          ,panel.grid.major = element_blank() #element_line(size=1.2)
          ,panel.grid.minor = element_blank()
          ,legend.text      = element_text(size=20)
          ,legend.title     = element_text(size=30)
          #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "none"
          ,plot.margin      = unit(c(2, 0.5, 0, 1),"cm")
          ,axis.title.y     = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
  
  pdf(paste0(dir.rec.pool.test,"/Recession_estimation.pdf"), 10, 8 ,useDingbats=F)
  print( reg.pool.plot[[mod]] )
  dev.off()
  
  print(reg.pool.plot[[mod]])
  
########################################################################################################
}


