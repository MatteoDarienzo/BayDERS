
###################################################
# Module for the MANUAL recostruction of hydrograph
###################################################
# Author: Matteo Darienzo
# 10/01/2022

  # charge general settings:
  source(file.options.general)
  
  
  ### Lucite Dombe:
  #################
  # ticks_RC.y.log = c(0.1, 1, 10, 100, 1000, 10000)
  # grid_RC.ylim.log = c(0.1, 10000)
  # Hmin = -0.8
  # Hmax = 12
  # Hgrid  = seq(Hmin, Hmax, 0.1)
  # 
  # #treal                = c(0,      3071.82,   5509,     5997,     6118,    6697,    7274,   10000,   18575,  19572,   21095,   21854,  22136,   22626,   22975)
  # treal                = c(0,      3071.82,   5509,     5997,    6697,     6878,      7274,   10000,   18575,  19572,   21095,   21854,  22136,   22626,   22975)  
  #  
  # 
  # c.prior.corrected    = c(1.589,   1.648,    1.538,    1.57,     1.3,      1.3,    1.3,    1.3,     1.3,     1.3,     1.3,    1.3,     1.3,     1.3)
  # st_c.prior.corrected = c(0.01,    0.01,     0.01,     0.01,     0.01,     0.01,   0.01,   0.01,    0.01,    0.01,    0.01,   0.01,    0.01,    0.01)
  # 
  # a.prior.corrected    = c(38.96,  32.84,    39.69,    45.43,    50.04,    50.04,   50.04,  50.04,   50.04,  50.04,    50.04,   50.04,   50.04,   50.04)
  # st_a.prior.corrected = c(2,       2,        2,        2,          2  ,    2,       2,      2,      2,       2,       2,      2,      2,        2)
  # 
  # b.prior.corrected    = c(1,      0.756,   0.536,     0.485,   0.576,   0.765,  0.658,   1.2,      0.51,    0.3,     0.5,     0.3,     0.5,     0.8)
  # st_b.prior.corrected = c( 0.2,    0.2,    0.1,      0.1,      0.1,     0.1,    0.1,    0.2,    0.1,    0.1,     0.1,     0.1,     0.1,     0.2)
  
  
  
  
  
  ### Estaquinhq:
  ###############
  ticks_RC.y.log = c(0.1, 1, 10, 100, 1000, 10000)
  grid_RC.ylim.log = c(0.1, 10000)
  Hmin = -0.8
  Hmax = 12
  Hgrid  = seq(Hmin, Hmax, 0.1)
  
  #treal              = c(0,       2000,      5307,     5667,     5672,       6765,               8699,   18081,   21000,   22279)
  treal                = c(0,       2000,      3226,     5601,     5814,       5967,     6970,     8699,   18081,    21000,   22279)
  
  c.prior.corrected    = c(2.573,   2.573,    2.573,     2.48,     2.652,      2.698,    2.664,    2.664,   2.664,   2.664)
  st_c.prior.corrected = c(0.01,    0.01,      0.01,     0.01,     0.01,       0.01,     0.01,     0.01,   0.01,     0.01)
  
  a.prior.corrected    = c(13.4,   13.415,    13.415,    18.259,    10.918,    10.1,     10.811,   10.8,   10.8,      10.8)
  st_a.prior.corrected = c(2,       1,          1,        1,         1,          1  ,     1,       2,      2,         2  )
  
  b.prior.corrected    = c( 0.8,   0.023,    0.023,      0.192,     0.124,      0.08,    0.081,    0.0,   -0.5,     -0.8)
  st_b.prior.corrected = c( 0.2,    0.2,      0.1,       0.1,       0.1,         0.1,     0.1,      0.2,   0.2,      0.2)
  
  
  
  
  

  
  
  
  
  
  ############################################################################################################
  #                                               COMPUTATION
  ############################################################################################################
  source(paste0(dir.modules,"/module_BaRatin.R"))
  limni.all            =  cbind(stage.record, date = format(as.Date(limni$Date, format= "%d/%m/%Y"), "%Y-%m-%d"))
  limni.NA             =  na.omit(limni.all)
  nlimni               =  length(limni.NA$t_limni)
  
  g1.prior             = c(0, 0.1, 0.001)
  g2.prior             = c(0, 0.1, 0.001)
  
  tin.period    =  treal
  t_shifts_date =  as.Date(floor(tin.period + limni$Datenum[1]), origin = date_origin)
  data.annotate.off$date =  as.Date(floor(data.annotate.off$xeffect + limni$Datenum[1]), origin = date_origin) 
  stage.record.date = as.Date(floor(stage.record$t_limni + limni$Datenum[1]), origin = date_origin) 
  # treal_after   =  c(treal[-1])
  # treal_before  =  c(treal[-length(treal)])
  treal_after   =  c(t_shifts_date[-1])
  treal_before  =  c(t_shifts_date[-length(t_shifts_date)])
  
  #b.prior.corrected    = c( 0.8,   0.8, 0.023,   0.192,      0.124,      0.08,    0.081,    0.0,   -0.5,     -0.8)
  #st_b.prior.corrected = c( 0.2,    0.2,  0.2,    0.1,       0.1,        0.2,     0.2,      0.2,   0.2,      0.2) 
  
  limni.shifts.plot <- ggplot() +
    geom_line(aes(x = stage.record.date, y=stage.record$h_limni), color = "lightblue",    size =0.2) +    
    # geom_line(aes(x=stage.record$t_limni, y=stage.record$h_limni), color = "lightblue",    size =0.2)+
    # scale_x_continuous(expand=c(0,0))+
    scale_x_date(expand=c(0,0), breaks= function(x) seq.Date(from = min(x), 
                                                             to = max(x),
                                                             by = "5 years")) +
    labs()+
    xlab(limni.labels[1]) +
    ylab(limni.labels[2])+
    # geom_point(aes(x = t_Gaug,  y = h_Gaug) , pch=21,    fill= "blue",   size = 2) + 
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text         = element_text(size=10)
          ,axis.title       = element_text(size=15)
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,legend.text      = element_text(size=20)
          ,legend.title     = element_text(size=30)
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "none"
          ,plot.margin      = unit(c(0.5,0.5,0.2, 1),"cm"))+
    # geom_vline(xintercept   = treal, color="blue", size=1, linetype="solid")+
    # geom_vline(xintercept = data.annotate.off$xeffect, color="red", size=1, linetype="solid")
    
    geom_vline(xintercept   = t_shifts_date, color="blue", size=1, linetype="solid")+
    geom_vline(xintercept   = data.annotate.off$date, color ="red", size=1, linetype="solid")+
    geom_rect(aes(xmin = treal_before,
                  xmax = treal_after,
                  ymin = b.prior.corrected - 2*st_b.prior.corrected,
                  ymax = b.prior.corrected+ 2*st_b.prior.corrected) , fill ="red", alpha=0.2)+
    geom_segment(aes(x= treal_before, xend= treal_after,  y= b.prior.corrected,  yend = b.prior.corrected) , color ="red", size=1)
  
  dirplotsQt = paste0(dir.case_study,"/Results")
  ggsave(limni.shifts.plot, filename =paste0(dirplotsQt,"/all_shift_times.png"), bg = "transparent", width = 14, height =6, dpi = 200)
  plot(limni.shifts.plot)
  
  
  
  #########################################################################################################
  ngrid  = length(Hgrid)
  dir.config = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2")
  write.table(Hgrid, file =paste0(dir.config,"/Hgrid.txt"), col.names = FALSE, row.names = FALSE)
  #write.table(limni.NA, file =paste0(dir.config,"/limni.txt"), col.names = FALSE, row.names = FALSE)
  # initialise:
  limn = hgrid = Qt.env = Qt.spag = RC.env = RC.spag = RC.summary =RC.mcmc = q2prior = q97prior = q50prior = NULL
  activat.stage.1 = b1.init = NULL
  Qt.2prior = Qt.50prior = Qt.97prior =NULL
  RC = Qt = Qt.cum = NULL
  
  ########################################################################################################
  for (t in 2:length(tin.period)){
    ######################################################################################################
    # for each period run BaRatin app.
    # extract the limni for the current period:
    #limni (stage record) of the current period "P":
    
    t_limni.P = h_limni.P  = index_limni = c()
    t_limni_date.P = c()
    j =0
    for (i in 1:nlimni){
      if ((limni.NA$t_limni[i] >= tin.period[t-1]) & (limni.NA$t_limni[i] < tin.period[t]))  {
        j = j+1
        t_limni.P[j]      = limni.NA$t_limni[i]
        t_limni_date.P[j] = format(as.Date(limni.NA$date[i], format="%Y-%m-%d"), "%Y-%m-%d")
        h_limni.P[j]      = limni.NA$h_limni[i]
        index_limni[j]    = i 
      }
    }   
    
    
    if (!is.null(t_limni.P)) {
    ##########################
      write(h_limni.P, file =paste0(dir.config,"/limni.txt"), ncolumns = 1, sep = "\n")
      nobs.limni.P  = length(h_limni.P)
      write.table(data4BaRatin[, 1:4], file =paste0(dir.config,"/Gaugings_data.txt"),  sep="\t",row.names=FALSE, col.names = TRUE)
      message("************************************************************************")
      message(paste0("Period ==> ", t-1 , "  from ", t_limni_date.P[1] , " (", t_limni.P[1],
                     ") to ", tail(t_limni_date.P,1) , " (", tail(t_limni.P, 1), ")") )
      message("************************************************************************")
      # Configure files for BaRatin.
      BaRatin_config(dir               = dir.config,
                     nsim              = 10000,
                     propagat          = FALSE,
                     b.distr           = b.distr,
                     a.distr           = "Gaussian",
                     c.distr           = c.distr,
                     a.prior           = a.prior.corrected[t-1],
                     st_a.prior        = st_a.prior.corrected[t-1],
                     c.prior           = c.prior.corrected[t-1],
                     st_c.prior        = st_c.prior.corrected[t-1],
                     b.prior           = b.prior.corrected[t-1],
                     st_b.prior        = st_b.prior.corrected[t-1],    #priors
                     Bw.prior          = Bw.prior,
                     Cr.prior          = Cr.prior,
                     g.prior           = g.prior,
                     Bc.prior          = Bc.prior,
                     KS.prior          = KS.prior,
                     S0.prior          = S0.prior,
                     st_Bw.prior       = st_Bw.prior,
                     st_Cr.prior       = st_Cr.prior,
                     st_g.prior        = st_g.prior,
                     st_Bc.prior       = st_Bc.prior,
                     st_KS.prior       = st_KS.prior,
                     st_S0.prior       = st_S0.prior,
                     ncontrol          = ncontrols,
                     M                 = M,
                     nobs              = 0,
                     Ncycles           = 100,
                     ngrid             = ngrid,
                     nlimni            = nobs.limni.P,
                     predictionRC      = FALSE,
                     predictionQt      = TRUE,
                     predictionPrior   = TRUE,
                     simMCMC           = FALSE,
                     mcmc.prior        = 1000,
                     remnant.err.model = remnant.err.model,
                     g1.prior          = g1.prior,
                     g2.prior          = g2.prior,
                     g1.distr.type     = g1.distr.type,
                     g2.distr.type     = g2.distr.type)
      
      message("***************************************************************"); flush.console()
      message(c("Applying BaRatin to period!!!  Wait ... ")); flush.console()
      message("***************************************************************"); flush.console()
      setwd(paste0(dir_code,"/BaM_exe"))
      system2(paste0(dir_code,"/BaM_exe/BaM_Segmentation2.exe"))
      
      # Read results:
      ###############
      limn[[t]]                    = read.table(paste0(dir.config,"/limni.txt"),               sep="\t",  header=FALSE)
      hgrid[[t]]                   = read.table(paste0(dir.config,"/Hgrid.txt"),               sep="\t",  header=FALSE)
      
      ########################
      # print RC :
      ########################
      RC.env[[t]]                  = read.table(paste0(dir.config,"/Qrc_Prior.env"),           sep="",    header=TRUE) 
      RC.spag[[t]]                 = read.table(paste0(dir.config,"/Qrc_Prior.spag"),          sep="",    header=FALSE)
      RC.summary[[t]]              = read.table(paste0(dir.config,"/Results_Summary.txt"),     sep="",    header=TRUE)
      RC.mcmc[[t]]                 = read.table(paste0(dir.config,"/Results_MCMC_Cooked.txt"), sep="",    header=TRUE)
      RC.spag[[t]][RC.spag[[t]]    == -666.666] = NA
      activat.stage.1[[t]]         = RC.summary[[t]][16,1] - 2*RC.summary[[t]][11,1]  
      b1.init[[t]]                 = RC.summary[[t]][16,1]
      q2prior[[t]]                 = 0; 
      q97prior[[t]]                = 0; 
      q50prior[[t]]                = 0;
      for (i in 1:length(RC.spag[[t]]$V1)) {
        q2prior[[t]][i]  = quantile(RC.spag[[t]][i,], probs = c(0.025), na.rm = TRUE)
        q97prior[[t]][i] = quantile(RC.spag[[t]][i,], probs = c(0.975), na.rm = TRUE)
        q50prior[[t]][i] = quantile(RC.spag[[t]][i,], probs = c(0.50) , na.rm = TRUE)
      }
      
      RC[[t]] = data.frame(h        = hgrid[[t]], 
                           RC.MAP   = q50prior[[t]],
                           q2total  = q2prior[[t]],
                           q97total = q97prior[[t]])
      
      if (length(which(RC[[t]]$V1 < activat.stage.1[[t]])) > 0){
        RC[[t]][1: tail(which(RC[[t]]$V1 < activat.stage.1[[t]]),1), ] = NA
      }
      RC[[t]][,2:4][ RC[[t]][,2:4] <= grid_RC.ylim.log[1] ] <- grid_RC.ylim.log[1]  #there are 4 columns.
      
      RC.plot <- ggplot() +
        scale_y_continuous(name = expression("Discharge Q [m3/s]"))+
        scale_x_continuous(name = expression("Stage h [m]"), limits = c(grid_RC.xlim[1], grid_RC.xlim[2])) +
        labs(x = "Stage [m]", y = "Discharge [m3/s]") +
        theme_bw(base_size=15)+
        theme(axis.text=element_text(size=10)
              ,axis.title=element_text(size=15, face="bold")
              #,axis.title.x = element_text(size=15, face="bold")
              ,panel.grid.major=element_blank()
              ,panel.grid.minor=element_blank()
              ,legend.text=element_text(size=20)
              ,legend.title=element_text(size=30)
              ,legend.key.size=unit(1.5, "cm")
              ,legend.position="none"
              ,plot.margin= unit(c(2, 1, 0, 1),"cm")
              ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
        geom_ribbon(data = RC[[t]], aes(x= V1, ymin = q2total, ymax=q97total), fill ="red",  alpha=0.4)+
        geom_line(data   = RC[[t]], aes(x= V1, y    = RC.MAP), size=0.3, color="black") 
      
      RC.plot.log = RC.plot  +
        scale_y_log10(name = expression("Discharge Q [m3/s]"),
                      breaks=ticks_RC.y.log, labels=ticks_RC.y.log, limits=grid_RC.ylim.log, expand = c(0,0)) +
        annotation_logticks(base = 10, sides = "l", scaled = TRUE,colour = "black",     size = 0.8, linetype = 1)
      
      ggsave(RC.plot, filename = paste0(dirplotsQt,"/RC_it",t,".png"),  width = 16, height =8, dpi = 200)
      ggsave(RC.plot.log, filename = paste0(dirplotsQt,"/RClog_it",t,".png"),  width = 16, height =8, dpi = 200)
      
      
      
      # Hydrograph:
      Qt.env[[t]]                  = read.table(paste0(dir.config,"/Qt_prior.env"),        sep="",    header=TRUE)
      Qt.spag[[t]]                 = read.table(paste0(dir.config,"/Qt_prior.spag"),       sep="",    header=FALSE)
      Qt.spag[[t]][Qt.spag[[t]]==-666.666] <- NA 
      Qt.2prior[[t]]               = 0
      Qt.97prior[[t]]              = 0  
      Qt.50prior[[t]]              = 0
      for (i in 1:length(Qt.spag[[t]]$V1)) {
        Qt.2prior[[t]][i]  = unlist(quantile(Qt.spag[[t]][i,], probs = c(0.025) , na.rm = TRUE))
        Qt.97prior[[t]][i] = unlist(quantile(Qt.spag[[t]][i,], probs = c(0.975) , na.rm = TRUE))
        Qt.50prior[[t]][i] = unlist(quantile(Qt.spag[[t]][i,], probs = c(0.50)  ,  na.rm = TRUE))
      }
      
      tt=2
      while(tt <=length(t_limni.P)){
        if ((t_limni.P[tt] - t_limni.P[tt-1]) > 2 ){
          ghost = tt
          print(tt)
          print("add NA")
          t_limni_date.P  = c(t_limni_date.P[1: (ghost-1)],   format(mean(c(as.Date(t_limni_date.P[tt],     format = "%Y-%m-%d"),
                                                                            as.Date(t_limni_date.P[tt-1],   format = "%Y-%m-%d"))), "%Y-%m-%d"),
                              t_limni_date.P[ghost:length(t_limni_date.P)])
          t_limni.P       = c(t_limni.P[1: (ghost-1)],       (t_limni.P[tt] +t_limni.P[tt-1])/2 ,           t_limni.P[ghost:length(t_limni.P)])
          h_limni.P       = c(h_limni.P[1: (ghost-1)],        NA,                                           h_limni.P[ghost:length(h_limni.P)])
          Qt.2prior[[t]]  = c(Qt.2prior[[t]][1: (ghost-1)],   NA,                                           Qt.2prior[[t]][ghost:length( Qt.2prior[[t]])])
          Qt.50prior[[t]] = c(Qt.50prior[[t]][1: (ghost-1)],  NA,                                           Qt.50prior[[t]][ghost:length( Qt.50prior[[t]])])
          Qt.97prior[[t]] = c(Qt.97prior[[t]][1: (ghost-1)],  NA,                                           Qt.97prior[[t]][ghost:length( Qt.97prior[[t]])])
          
          tt = tt+1
        }
        tt = tt+1
      }
      
      Qt[[t]] = data.frame( t        =  t_limni.P,
                            tdate    =  t_limni_date.P,    #format(as.Date(t_limni_date.P,  format= "%d/%m/%Y"), "%Y-%m-%d"),
                            h        =  h_limni.P,
                            Qt.MAP   =  unlist(Qt.50prior[[t]]),
                            q2total  =  unlist(Qt.2prior[[t]]),
                            q97total =  unlist(Qt.97prior[[t]]))
      
      
      Qt[[t]][,4:6][ Qt[[t]][,4:6] <= grid_RC.ylim.log[1]] <- grid_RC.ylim.log[1]
      
      if (t > 2) { 
        ############
        if ((t_limni.P[1] - tail(Qt.cum[[t-1]]$t,1)) > 3) {
          print("add NA")
          Qt[[t]] = rbind(data.frame(t        = (tail(Qt.cum[[t-1]]$t,1)  +  t_limni.P[1])/2 ,
                                     tdate    = format(mean(c(as.Date(tail(Qt.cum[[t-1]]$tdate,1), format ="%Y-%m-%d"), 
                                                              as.Date(t_limni_date.P[1],  format = "%Y-%m-%d"))), "%Y-%m-%d"),
                                     h        = NA, 
                                     Qt.MAP   = NA, 
                                     q2total  = NA, 
                                     q97total = NA),  
                          Qt[[t]])
          #Qt.cum$tdate =  format(as.Date( Qt.cum$tdate, format ="%Y-%m-%m"), "%d/%m/%Y")
          
        } else {
          #Qt.cum$tdate =  format(as.Date( Qt.cum$tdate, format ="%Y-%m-%m"), "%d/%m/%Y")
        }
        Qt.cum[[t]]  = rbind(Qt.cum[[t-1]],  Qt[[t]])
        #########    
      }  else {
        #########
        # only for first period:
        Qt.cum[[t]] =   Qt[[t]]
      }
      
      
      index_last_t_limni =  tail(index_limni,1)
      
      
      
      #############
    } else {   # no data in the period:
      #############
      Qt.cum[[t]] = rbind(Qt.cum[[t-1]], 
                          data.frame(t        = (limni.NA$t_limni[index_last_t_limni]  +  limni.NA$t_limni[index_last_t_limni+1])/2,
                                     tdate    = format(mean(c(as.Date(limni.NA$date[index_last_t_limni],      format = "%Y-%m-%d"), 
                                                              as.Date(limni.NA$date[index_last_t_limni+1],    format = "%Y-%m-%d"))), "%Y-%m-%d"),
                                     h        = NA, 
                                     Qt.MAP   = NA, 
                                     q2total  = NA, 
                                     q97total = NA))
      
      index_last_t_limni = index_last_t_limni
    }
    
    
    
    
    
    Qt.cum[[t]]$tdate = as.Date(Qt.cum[[t]]$tdate, format ="%Y-%m-%d")
    
    ########################
    # print HYDROGRAPH Q(t):
    ########################
    Qt.plot <- ggplot() +
      scale_x_date(name = "Time [year]",  expand=c(0,0))+
      #scale_x_continuous(name = expression("Time [day]"), expand=c(0,0))+    
      scale_y_continuous(name = expression("Discharge Q [m3/s]"), expand=c(0,0))+
      labs(x = "Time [day]", y = "Discharge [m3/s]") +
      theme_bw(base_size=15)+
      theme(axis.text=element_text(size=10)
            ,axis.title=element_text(size=15, face="bold")
            #,axis.title.x = element_text(size=15, face="bold")
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin= unit(c(2, 1, 0, 1),"cm")
            ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      geom_ribbon(data = Qt.cum[[t]], aes(x= tdate, ymin = q2total, ymax= q97total), fill="red", alpha =0.4)+
      geom_line(data = Qt.cum[[t]], aes(x= tdate, y = Qt.MAP),   size=0.5, color ="black")
    #geom_vline(xintercept = tin.period[-1], color="blue", linetype = "dashed", size =1)
    
    
    Qt.plot.log = Qt.plot  +
      scale_y_log10(name = expression("Discharge Q [m3/s]"),
                    breaks=ticks_RC.y.log, labels=ticks_RC.y.log, limits=grid_RC.ylim.log, expand = c(0,0)) +
      annotation_logticks(base = 10, sides = "l", scaled = TRUE,colour = "black",     size = 0.8, linetype = 1)
    
    ggsave(Qt.plot, filename = paste0(dirplotsQt,"/Qt_it",t,".png"),  width = 16, height =6, dpi = 200)
    ggsave(Qt.plot.log, filename = paste0(dirplotsQt,"/Qtlog_it",t,".png"), width = 16, height =6, dpi = 200)
    
    
    Qt[[t]]$tdate = as.Date(Qt[[t]]$tdate, format ="%Y-%m-%d")
    
    Qt.plot.period <- ggplot() +
      scale_x_date(name = "Time [year]",  expand=c(0,0))+
      #scale_x_continuous(name = expression("Time [day]"), expand=c(0,0))+    
      scale_y_continuous(name = expression("Discharge Q [m3/s]"), expand=c(0,0))+
      labs(x = "Time [day]", y = "Discharge [m3/s]") +
      theme_bw(base_size=15)+
      theme(axis.text=element_text(size=10)
            ,axis.title=element_text(size=15, face="bold")
            #,axis.title.x = element_text(size=15, face="bold")
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin= unit(c(2, 1, 0, 1),"cm")
            ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      geom_ribbon(data = Qt[[t]], aes(x= tdate, ymin = q2total, ymax= q97total), fill="red", alpha =0.4)+
      geom_line(data = Qt[[t]], aes(x= tdate, y = Qt.MAP),   size=0.5, color ="black")
    
    ggsave(Qt.plot.period, filename = paste0(dirplotsQt,"/Qt_period",t,".png"),  width = 16, height =6, dpi = 200)
    
    
    
    
    
    hydrograph <- read.csv2(paste0(dir.case_study,"/",file_gaugings), fileEncoding="UTF-8", quote="", sep=";",dec=".",header=TRUE, na.strings=c(";","NA"))
    hydrograph$Date = as.Date(hydrograph$Date, format ="%d/%m/%Y")
    
    Qt.plot.off <- ggplot() +
      scale_x_date(name = "Time [year]",  expand=c(0,0))+
      #scale_x_continuous(name = expression("Time [day]"), expand=c(0,0))+    
      scale_y_continuous(name = expression("Discharge Q [m3/s]"), expand=c(0,0))+
      labs(x = "Time [day]", y = "Discharge [m3/s]") +
      theme_bw(base_size=15)+
      theme(axis.text=element_text(size=10)
            ,axis.title=element_text(size=15, face="bold")
            #,axis.title.x = element_text(size=15, face="bold")
            ,panel.grid.major=element_blank()
            ,panel.grid.minor=element_blank()
            ,legend.text=element_text(size=20)
            ,legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm")
            ,legend.position="none"
            ,plot.margin= unit(c(2, 1, 0, 1),"cm")
            ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
      geom_ribbon(data = Qt.cum[[t]], aes(x= tdate, ymin = q2total, ymax= q97total), fill="red", alpha =0.5)+
      geom_line(data = Qt.cum[[t]], aes(x= tdate, y = Qt.MAP),   size=0.5, color ="black") +
      geom_line(data = hydrograph, aes(x= Date,  y = Q),   size=0.5, color ="lightblue")
    
    
    
    ggsave(Qt.plot.off, filename = paste0(dirplotsQt,"/Qt_new_vs_old.png"), width = 16, height =6, dpi = 200)      
  }
  ################
  # end of loop.
  ################
  
  
  # save hydrograph to csv file:
  #Qt.cum$tdate = format(as.Date(Qt.cum$tdate, format = "%Y-%m-%d"), "%d/%m/%Y")
  write.table(Qt.cum[[length(Qt.cum)]][,2:6], file = paste0(dir.config,"/hydrograph.csv"), row.names = FALSE, dec = ".", sep = ";", quote = FALSE)
  
  
  message("
############################################
#              All done !                  #  
############################################
")
  


  
  
  




# #save results of BaRatin: 
# list.of.files <- c(
#   paste0(dir.BaRatin.config,"/limni.txt"),
#   paste0(dir.BaRatin.config,"/Hgrid.txt"),
#   paste0(dir.BaRatin.config,"/Results_Summary.txt"),
#   paste0(dir.BaRatin.config,"/Config_Model.txt"),
#   paste0(dir.BaRatin.config,"/Qrc_Prior.spag"),
#   paste0(dir.BaRatin.config,"/Qrc_Prior.env"), 
#   paste0(dir.BaRatin.config,"/Gaugings_data.txt"),
#   paste0(dir.BaRatin.config,"/Qt_Prior.spag"),
#   paste0(dir.BaRatin.config,"/Qt_Prior.env")
# )
# for (ll in 1:length(list.of.files)) {
#   file.copy(list.of.files[ll], paste0(dir.rt.it[[rt]],"/BaRatin/from_ST_priors"), overwrite = TRUE)
# } 
####################################################################################################


