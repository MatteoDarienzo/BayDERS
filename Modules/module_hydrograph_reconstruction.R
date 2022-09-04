
###################################################
# Module for the MANUAL recostruction of hydrograph
###################################################
# Author: Matteo Darienzo
# 10/01/2022

  # charge general settings:
  source(file.options.general)
  







### Rurrenabaqye Bolivia:
#########################
ticks_RC.y.log = c(100, 1000, 10000, 100000)
grid_RC.ylim.log = c(100, 100000)
Hmin = -2
Hmax = 7
Hgrid  = seq(Hmin, Hmax, 0.01)

############### Using new segmentation :
Nperiods = 19
treal                = c(0,  400, 659.731, 1045.41, 1383.09, 2464.74, 2773.1, 3147.35, 3442.81, 3561.86, 3792.02, 3905.77, 4249.33, 5034.28, 5676.13, 6565.11, 6892.89, 7031.92, 7167.04, 8089.958)
c.prior.corrected    = rep(list(c(1.697270, 1.650040)), Nperiods)
st_c.prior.corrected = rep(list(c(0.02, 0.026458)), Nperiods)

a.prior.corrected    = rep(list(c(359.831, 432.682000)), Nperiods)
st_a.prior.corrected = rep(list(c(40, 61)), Nperiods)

b2.prior.corrected    = 1.153920
st_b2.prior.corrected = 0.303400

b1.prior.corrected    = c(-0.83, -0.43, -0.83, -0.65, -0.7354, -0.7477, -0.75942, -0.907, -0.916, -0.984, -1.2, -1.21, -1.38, -1.64, -1.34, -2.03, -1.64, -0.87, -1.55)
st_b1.prior.corrected = c(0.25,   0.1,   0.15, 0.143,  0.112,   0.178,    0.145,   0.1602, 0.2,    0.19,   0.2,   0.19,  0.2,  0.22,  0.1,   0.234, 0.248, 0.276, 0.25)
b.prior.corrected = st_b.prior.corrected = NULL
for (i in 1:length(b1.prior.corrected)){
   b.prior.corrected[[i]]    = c(b1.prior.corrected[i],    b2.prior.corrected) 
   st_b.prior.corrected[[i]] = c(st_b1.prior.corrected[i], st_b2.prior.corrected)   
}
g1.distr.type          = 'Uniform'
g2.distr.type          = 'Uniform'
g1.prior               = c(0,   97,   90)
g2.prior               = c(0, 0.008, 0.007)
dirplotsQt = paste0(dir.case_study,"/Results/Reconstruction")

# c.prior.corrected    = rep(list(c(1.697270)), Nperiods)
# st_c.prior.corrected = rep(list(c(0.01)), Nperiods)
# 
# a.prior.corrected    = rep(list(c(359.831)), Nperiods)
# st_a.prior.corrected = rep(list(c(5)), Nperiods)
# b.prior.corrected    = list(-1.008380, -0.646025,	-0.735455,	-0.747683, -0.759419,	-0.906989, -0.916295, -0.983358,	-1.198930, 
#                             -1.207530, -1.382520,	-1.589190, -2.030680, -1.641240,	-0.869200,	-1.548870, -1.548870, -1.548870)
# st_b.prior.corrected = list(0.143089,	0.142877,	0.111940, 0.177593,	0.145189, 0.160194,	0.199864,	0.189920,	0.204962, 
#                             0.191067,	0.200128,	0.211174, 0.233940,	0.247046,	0.275630, 0.230131, 0.230131, 0.230131)




# ########## official 2014
# Nperiods              = 1
# treal                 = c(0, 8089.958)
# c.prior.corrected     = rep(list(c(1.59, 1.23)),  Nperiods)
# st_c.prior.corrected  = rep(list(c(0.001, 0.001)), Nperiods)
# 
# a.prior.corrected     = list(c(587, 1144))
# st_a.prior.corrected  = rep(list(c(120, 140)),   Nperiods)
# 
# b2.prior.corrected    = 2.25
# st_b2.prior.corrected = 0.4
# 
# b1.prior.corrected    = c(-0.66)
# st_b1.prior.corrected = c(0.3)
# b.prior.corrected     = st_b.prior.corrected = NULL
# for (i in 1:length(b1.prior.corrected)){
#   b.prior.corrected[[i]]    = c(b1.prior.corrected[i],    b2.prior.corrected) 
#   st_b.prior.corrected[[i]] = c(st_b1.prior.corrected[i], st_b2.prior.corrected)   
# }
# g1.distr.type          = 'Uniform'
# g2.distr.type          = 'Uniform'
# g1.prior               = c(0,   1,   0.5)
# g2.prior               = c(0, 0.00001, 0.000005)
# dir.create(paste0(dir.case_study,  "/Results/segmentation_gaugings/Official_segmentation/SPD_test1/2014"))
# dirplotsQt = paste0(dir.case_study,"/Results/segmentation_gaugings/Official_segmentation/SPD_test1/2014")
# 
# 
# 
# 
# ########## official  2012 
# Nperiods              = 1
# treal                 = c(0,  8089.958)
# c.prior.corrected     = rep(list(c(1.88, 1.88)),  Nperiods)
# st_c.prior.corrected  = rep(list(c(0.001, 0.001)), Nperiods)
# 
# a.prior.corrected     = list(c(423, 464))
# st_a.prior.corrected  = rep(list(c(10, 10)),   Nperiods)
# 
# b2.prior.corrected    = 1.9
# st_b2.prior.corrected = 0.4
# 
# b1.prior.corrected    = c(-1.35)
# st_b1.prior.corrected = c(0.2)
# b.prior.corrected     = st_b.prior.corrected = NULL
# for (i in 1:length(b1.prior.corrected)){
#   b.prior.corrected[[i]]    = c(b1.prior.corrected[i],    b2.prior.corrected) 
#   st_b.prior.corrected[[i]] = c(st_b1.prior.corrected[i], st_b2.prior.corrected)   
# }
# g1.distr.type          = 'Uniform'
# g2.distr.type          = 'Uniform'
# g1.prior               = c(0,   1,   0.5)
# g2.prior               = c(0, 0.00001, 0.000005)
# dir.create(paste0(dir.case_study,  "/Results/segmentation_gaugings/Official_segmentation/SPD_test1/2012"))
# dirplotsQt = paste0(dir.case_study,"/Results/segmentation_gaugings/Official_segmentation/SPD_test1/2012")
# 


		




   
  # ### Lucite Dombe:
  # #################
  # ticks_RC.y.log = c(0.1, 1, 10, 100, 1000, 10000)
  # grid_RC.ylim.log = c(0.1, 10000)
  # Hmin = -0.8
  # Hmax = 12
  # Hgrid  = seq(Hmin, Hmax, 0.1)
  # 
  # #treal                = c(0,      3071.82,   5509,     5997,     6118,    6697,    7274,   10000,   18575,  19572,   21095,   21854,  22136,   22626,   22975)
  # treal                = c(0,                 5509,     5997,    6697,     6878,      7274,   10000,   18575,  19572,   21095,   21854,  22136,   22626,   22975)
  # 
  # 
  # c.prior.corrected    = c(1.589,     1.589,   1.648,    1.538,    1.57,     1.3,   1.3,     1.3,     1.3,     1.3,     1.3,     1.3,     1.3)
  # st_c.prior.corrected = c(0.01,      0.01,     0.01,     0.01,     0.01,   0.01,   0.01,    0.01,    0.01,    0.01,   0.01,    0.01,    0.01)
  # 
  # a.prior.corrected    = c(39,        38.96,    32.84,    39.69,    45.43,    50.04,  50.04,  50.04,   50.04,   50.04,  50.04,   50.04,   50.04)
  # st_a.prior.corrected = c(2,         1,        1,          1  ,    1,        1,       2,     2,       2,       2,       2,      2,        2)
  # 
  # b.prior.corrected    = c(1,         0.536,     0.485,   0.576,    0.765,    0.658,   1.2,    0.51,    0.3,     0.5,     0.3,     0.5,     0.8)
  # st_b.prior.corrected = c( 0.2,      0.1,       0.1,      0.1,      0.1,      0.1,    0.2,    0.1,     0.1,     0.1,     0.1,     0.1,     0.2)
  # 
  # g1.prior             = c(0, 0.0001, 0.00001)
  # g2.prior             = c(0, 0.0000001, 0.00000001)
  # dirplotsQt = paste0(dir.case_study,"/Results")


  
  
  
  ### Estaquinhq:
  ###############
  # ticks_RC.y.log = c(0.1, 1, 10, 100, 1000, 10000)
  # grid_RC.ylim.log = c(0.1, 10000)
  # Hmin = -0.8
  # Hmax = 12
  # Hgrid  = seq(Hmin, Hmax, 0.1)
  # 
  # #treal              = c(0,       2000,      5307,     5667,     5672,       6765,               8699,   18081,   21000,   22279)
  # treal                = c(0,       2000,      3226,     5601,     5814,       5967,     6970,     8699,   18081,    21000,   22279)
  # 
  # c.prior.corrected    = c(2.573,   2.573,    2.573,     2.48,     2.652,      2.698,    2.664,    2.664,   2.664,   2.664)
  # st_c.prior.corrected = c(0.01,    0.01,      0.01,     0.01,     0.01,       0.01,     0.01,     0.01,   0.01,     0.01)
  # 
  # a.prior.corrected    = c(13.4,   13.415,    13.415,    18.259,    10.918,    10.1,     10.811,   10.8,   10.8,      10.8)
  # st_a.prior.corrected = c(2,       1,          1,        1,         1,          1  ,     1,       2,      2,         2  )
  # 
  # b.prior.corrected    = c( 0.8,   0.023,    0.023,      0.192,     0.124,      0.08,    0.081,    0.0,   -0.5,     -0.8)
  # st_b.prior.corrected = c( 0.2,    0.2,      0.1,       0.1,       0.1,         0.1,     0.1,      0.2,   0.2,      0.2)
  #
  # g1.prior             = c(0, 0.0001, 0.00001)
  # g2.prior             = c(0, 0.0000001, 0.00000001)
  # dirplotsQt = paste0(dir.case_study,"/Results")

  

  
  
  










  
  
  ############################################################################################################
  #                                               COMPUTATION
  ############################################################################################################
  source(paste0(dir.modules,"/module_BaRatin.R"))


  stage.record           = df.limni
  limni.all              = cbind(stage.record, date = format(as.Date(stage.record$t_limni.true,  origin = date_origin), "%Y-%m-%d"))
  limni.NA               = na.omit(limni.all)
  nlimni                 = length(limni.NA$t_limni)
  tin.period             = treal
  t_shifts_date          = as.Date(floor(tin.period + stage.record$t_limni.true[1]), origin = date_origin, format = "%Y-%m-%d")
  data.annotate.off$date = as.Date(floor(data.annotate.off$xeffect + stage.record$t_limni.true[1]), origin = date_origin) 
  stage.record.date      = as.Date(floor(stage.record$t_limni + stage.record$t_limni.true[1]), origin = date_origin) 
  write.table(t_shifts_date, file = paste0(dirplotsQt,"/RC_shift_dates.csv"), row.names = FALSE, dec = ".", sep = ";", quote = FALSE)
  # treal_after          = c(treal[-1])
  # treal_before         = c(treal[-length(treal)])
  treal_after            = c(t_shifts_date[-1])
  treal_before           = c(t_shifts_date[-length(t_shifts_date)])
  b1.prior.final         = sapply(b.prior.corrected,"[[",1)
  st_b1.prior.final      = sapply(st_b.prior.corrected,"[[",1)
  
  # Plot stage record with offset b over time:
  limni.shifts.plot <- ggplot() +
    geom_line(aes(x = stage.record.date, y=stage.record$h_limni), color = "lightblue",    size =0.2) +    
    # geom_line(aes(x=stage.record$t_limni, y=stage.record$h_limni), color = "lightblue",    size =0.2)+
    # scale_x_continuous(expand=c(0,0))+
    scale_x_date(expand=c(0,0), breaks= function(x) seq.Date(from = min(x), to = max(x), by = "5 years")) +
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
                  ymin = b1.prior.final - 2*st_b1.prior.final,
                  ymax = b1.prior.final + 2*st_b1.prior.final) , fill ="red", alpha=0.2)+
    geom_segment(aes(x= treal_before, xend= treal_after,  y= b1.prior.final,  yend = b1.prior.final) , color ="red", size=1)
  

  ggsave(limni.shifts.plot, filename =paste0(dirplotsQt,"/all_shift_times.png"), bg = "transparent", width = 14, height =6, dpi = 200)
  plot(limni.shifts.plot)
  
  
  
  #########################################################################################################
  ngrid  = length(Hgrid)
  dir.config = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2")
  write.table(Hgrid, file =paste0(dir.config,"/Hgrid.txt"), col.names = FALSE, row.names = FALSE)
  #write.table(limni.NA, file =paste0(dir.config,"/limni.txt"), col.names = FALSE, row.names = FALSE)
  # Initialise arrayys and lists:
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
      nobs_gaug = length(data4BaRatin$h)  # not used
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
                     a.prior           = a.prior.corrected[[t-1]],
                     st_a.prior        = st_a.prior.corrected[[t-1]],
                     c.prior           = c.prior.corrected[[t-1]],
                     st_c.prior        = st_c.prior.corrected[[t-1]],
                     b.prior           = b.prior.corrected[[t-1]],
                     st_b.prior        = st_b.prior.corrected[[t-1]],    #priors
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
                     nobs              = nobs_gaug,
                     Ncycles           = 100,
                     ngrid             = ngrid,
                     nlimni            = nobs.limni.P,
                     predictionRC      = FALSE,
                     predictionQt      = FALSE,
                     predictionPrior   = TRUE,
                     simMCMC           = TRUE,
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
        theme_bw(base_size=20)+
        theme(axis.text=element_text(size=10)
              ,axis.title=element_text(size=15, face="bold")
              #,axis.title.x = element_text(size=15, face="bold")
              #,panel.grid.major=element_blank()
              #,panel.grid.minor=element_blank()
              ,legend.text=element_text(size=20)
              ,legend.title=element_text(size=30)
              ,legend.key.size=unit(1.5, "cm")
              ,legend.position="none"
              ,plot.margin= unit(c(2, 1, 0, 1),"cm")
              ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
        geom_ribbon(data = RC[[t]], aes(x= V1, ymin = q2total, ymax = q97total), fill ="red",  alpha=0.2)+
        geom_line(data   = RC[[t]], aes(x= V1, y    = RC.MAP), size=1, color="black") 
      
      RC.plot.log = RC.plot  +
          scale_y_log10(name = expression("Discharge Q [m3/s]"),
                        breaks=ticks_RC.y.log, labels=ticks_RC.y.log, limits=grid_RC.ylim.log, expand = c(0,0)) +
          annotation_logticks(base = 10, sides = "l", scaled = TRUE,colour = "black",     size = 0.8, linetype = 1)
      
      ggsave(RC.plot,     filename = paste0(dirplotsQt,"/RC_it",t,".png"),  width = 12, height =8, dpi = 200)
      ggsave(RC.plot.log, filename = paste0(dirplotsQt,"/RClog_it",t,".png"),  width = 12, height =8, dpi = 200)
      
      
      
      
      
      
      
    
      ####################
      # print HYDROGRAPH :
      ####################
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
    
    # plot Q(t):
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
    
    ggsave(Qt.plot,     filename = paste0(dirplotsQt,"/Qt_it",t,".png"),  width = 16, height =6, dpi = 200)
    ggsave(Qt.plot.log, filename = paste0(dirplotsQt,"/Qtlog_it",t,".png"), width = 16, height =6, dpi = 200)
    
    Qt[[t]]$tdate = as.Date(Qt[[t]]$tdate, format ="%Y-%m-%d")
    # Q(t) for the current period only:
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
  }
  ################
  # end of loop.
  ################
  
  
  
  
  
  
  
  dateQ = format(as.Date(stage.record$t_limni.true,  origin = date_origin), "%Y-%m-%d %H:%S")
  Qt.fin = t_limni.date
  
  # Save final new hydrograph to csv file:
  ########################################
  # Qt.cum$tdate = format(as.Date(Qt.cum$tdate, format = "%Y-%m-%d"), "%d/%m/%Y")
  write.table(Qt.cum[[length(Qt.cum)]][,2:6], file = paste0(dirplotsQt,"/hydrograph.csv"), row.names = FALSE, dec = ".", sep = ";", quote = FALSE)
  
  
  
  
  
  # Plot all RCs together:
  ########################
  colorr = colo
  RC.all.plot <- ggplot() +
    scale_y_continuous(name = expression("Discharge Q [m3/s]"))+
    scale_x_continuous(name = expression("Stage h [m]"))+ #, limits = c(grid_RC.xlim[1], grid_RC.xlim[2])) +
    labs(x = "Stage [m]", y = "Discharge [m3/s]") +
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=10)
          ,axis.title=element_text(size=15, face="bold")
          #,axis.title.x = element_text(size=15, face="bold")
          #,panel.grid.major=element_blank()
          #,panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,plot.margin= unit(c(2, 1, 0, 1),"cm")
          ,axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
  
  for( i in 2:length(RC)){
    RC.all.plot = RC.all.plot+
      geom_ribbon(data = RC[[i]], aes(x= V1, ymin = q2total, ymax=q97total), fill ="red",  alpha=0.1)+
      geom_line(data   = RC[[i]], aes(x= V1, y    = RC.MAP), size=1, color= "black") #colorr[i]) 
    
  }
  ggsave(RC.all.plot, filename = paste0(dirplotsQt,"/RC_all.png"),  width = 11, height =8, dpi = 200)
  RC.all.plot = RC.all.plot  +
    scale_y_log10(name = expression("Discharge Q [m3/s]"),
                  breaks=ticks_RC.y.log, labels=ticks_RC.y.log, limits=grid_RC.ylim.log, expand = c(0,0)) +
    annotation_logticks(base = 10, sides = "l", scaled = TRUE,colour = "black",     size = 0.8, linetype = 1)
  ggsave(RC.all.plot, filename = paste0(dirplotsQt,"/RClog_all.png"),  width = 11, height =8, dpi = 200)
  
  
  
  
  
  
  
  
  
  
  
  
  
  # COMPARISON WITH OLD HYDROGRAPH:
  ##################################
  # hydrograph      = read.csv2(paste0(dir.case_study,"/", "E188_Rio_Buzi_at_Estaquinha_discharge.csv"), fileEncoding="UTF-8", quote="", sep=",",dec=".", header=TRUE, na.strings=c(";","NA"))
  # hydrograph$Date = as.Date(hydrograph$X, format ="%m/%d/%Y")

  # hydrograph      = read.csv2(paste0(dir.case_study,"/", "E246_Rio_Lucite_at_Dombe_discharge.csv"), fileEncoding="UTF-8", quote="", sep=";",dec=".", header=TRUE, na.strings=c(";","NA"))
  # hydrograph$Date = as.Date(hydrograph$Date, format ="%d/%m/%Y")

  
  #########################################################################################################################################
  Qt.plot.off = ggplot() +
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
                geom_line(data = Qt.cum[[t]], aes(x= tdate, y = Qt.MAP),   size=0.5, color ="black")+
                geom_line(data = hydrograph, aes(x= Date,  y = Discharge),   size=0.5, color ="lightblue")
                ggsave(Qt.plot.off, filename = paste0(dirplotsQt,"/Qt_new_vs_old.png"), width = 16, height =6, dpi = 200)
  
  
                
                
  # Zoom on some periods:
  ########################
   Qt.plot.off_1 = ggplot() +
                  scale_x_date(name = "Time [year]",  expand=c(0,0), limits = c(as.Date("1957-01-01"), as.Date("1961-01-01")))+
                  #scale_x_continuous(name = expression("Time [day]"), expand=c(0,0))+
                  scale_y_continuous(name = expression("Discharge Q [m3/s]"), expand=c(0,0), limits = c(0, 4000))+
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
                  geom_line(data = Qt.cum[[t]], aes(x= tdate, y = Qt.MAP),   size=1, color ="black")+
                  geom_line(data = hydrograph, aes(x= Date,  y = Discharge),   size=1, color ="lightblue")
                ggsave(Qt.plot.off_1, filename = paste0(dirplotsQt,"/Qt_new_vs_old_Period1.png"), width = 16, height =6, dpi = 200)
                
                
                Qt.plot.off_2 = ggplot() +
                  scale_x_date(name = "Time [year]",  expand=c(0,0), limits = c(as.Date("1965-01-01"), as.Date("2000-01-01")))+
                  #scale_x_continuous(name = expression("Time [day]"), expand=c(0,0))+
                  scale_y_continuous(name = expression("Discharge Q [m3/s]"), expand=c(0,0), limits = c(0, 4000))+
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
                  geom_line(data = Qt.cum[[t]], aes(x= tdate, y = Qt.MAP),   size=1, color ="black")+
                  geom_line(data = hydrograph, aes(x= Date,  y = Discharge),   size=1, color ="lightblue")
                ggsave(Qt.plot.off_2, filename = paste0(dirplotsQt,"/Qt_new_vs_old_Period2.png"), width = 16, height =6, dpi = 200)
                
                
                Qt.plot.off_3 = ggplot() +
                  scale_x_date(name = "Time [year]",  expand=c(0,0), limits = c(as.Date("2005-01-01"), as.Date("2017-01-01")))+
                  #scale_x_continuous(name = expression("Time [day]"), expand=c(0,0))+
                  scale_y_continuous(name = expression("Discharge Q [m3/s]"), expand=c(0,0), limits = c(0, 4000))+
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
                  geom_line(data = Qt.cum[[t]], aes(x= tdate, y = Qt.MAP),   size=1, color ="black")+
                  geom_line(data = hydrograph, aes(x= Date,  y = Discharge),   size=1, color ="lightblue")
                ggsave(Qt.plot.off_3, filename = paste0(dirplotsQt,"/Qt_new_vs_old_Period3.png"), width = 16, height =6, dpi = 200)
                
                
                
                
               
                
                
          # Comparison with two official series:
          ######################################
          hydrograph_new = read.csv2(paste0(dir.case_study,
                                     "/Results/Reconstruction/hydrograph.csv"), 
                                     fileEncoding="UTF-8", quote="", sep=";",dec=".", header=TRUE, na.strings=c(";","NA"))
          hydrograph_new$Date = as.Date(hydrograph_new$tdate,  format ="%Y-%m-%d")
          hydrograph1  = read.csv2(paste0(dir.case_study,
                                   "/Results/segmentation_gaugings/Official_segmentation/SPD_test1/2014/hydrograph.csv"), 
                                   fileEncoding="UTF-8", quote="", sep=";",dec=".", header=TRUE, na.strings=c(";","NA"))
          hydrograph1$Date = as.Date(hydrograph1$tdate,  format ="%Y-%m-%d")
          hydrograph2  = read.csv2(paste0(dir.case_study,
                                   "/Results/segmentation_gaugings/Official_segmentation/SPD_test1/2012/hydrograph.csv"), 
                                   fileEncoding="UTF-8", quote="", sep=";",dec=".", header=TRUE, na.strings=c(";","NA"))    
          hydrograph2$Date = as.Date(hydrograph2$tdate,  format ="%Y-%m-%d")
                
          Qt.plot.event = ggplot() +
            scale_x_date(name = "Date",  expand=c(0,0), limits = c(as.Date("2013-12-23"), as.Date("2014-03-01")), labels=date_format("%d/%m/%Y"))+
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
            geom_ribbon(data = hydrograph_new, aes(x= Date, ymin = q2total, ymax= q97total), fill="red", alpha =0.2)+
            geom_line(data   = hydrograph_new, aes(x= Date, y = Qt.MAP),   size=1, color ="black")+
            
            geom_line(data   = hydrograph1,  aes(x= Date,  y = Qt.MAP),   size=1, color ="black", linetype = "dotted")+
            geom_line(data   = hydrograph2,  aes(x= Date,  y = Qt.MAP),   size=1, color ="black", linetype = "dashed") 
            Qt.plot.event
            ggsave(Qt.plot.event, filename = paste0(dir.case_study, "/Results/Reconstruction/Qt_new_vs_old.png"), width = 16, height =6, dpi = 200)    
                
                
                
                 

  
  
    
  
  message("
############################################
#              All done !                  #  
############################################
")
####################################################################################################
  
  
  


  
  
  



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

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
