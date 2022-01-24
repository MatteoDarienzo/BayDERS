############################################################################################
#                              Input Data Reading ...
############################################################################################
# author: Matteo Darienzo, INRAE
# purpose: 
# initialise the project, directories and read all input data for computation, for instance:
# - stage record, 
# - stage-discharge gaugings,
# - official dates of RC update (if available)




#********************************************************************************************
#Read the directories:
#********************************************************************************************
dir.case_study =  paste0(dir_code,"/Case_studies/", case_study_name) #creation of the folder with the case study. 
dir.exe        =  paste0(dir_code,"/BaM_exe")     #read directory of BaM.exe.
message("
#####################################
#         Welcome to BayDERS        #  
#####################################
# This program includes methods for #
# detecting stage-discharge rating  #
# shifts by using gaugings and/or   #
# stage record.                     #
#                                   #
# Developed by:                     #
# Matteo Darienzo (INRAE Lyon).     #
# Under supervision of:             #
# B. Renard, J. Le Coz, M. Lang.    # 
# Version v1. 2021                  #
# Funded by EDF, CNR, SCHAPI        #   
#####################################
# References:                       #
# Le Coz et al., 2014)              # 
# Mansanarez et al., 2019)          #
# Darienzo et al., 2021             #
#####################################
# Case study:                        ")



#********************************************************************************************
#                        Include all modules for computation:                               *
#********************************************************************************************
source(paste0(dir.case_study,"/Options_General.R"))    
source(paste0(dir.modules,"/module_gaugings_segmentation.R"))
source(paste0(dir.modules,"/module_gaugings_segmentation_plots.R"))
source(paste0(dir.modules,"/module_BaRatin.R"))
source(paste0(dir.modules,"/module_BaRatinSPD.R"))
source(paste0(dir.modules,"/module_BaM_results.R"))
source(paste0(dir.modules,"/module_shift_generator.R"))
source(paste0(dir.modules,"/module_prior_propagation.R"))
source(paste0(dir.modules,"/module_recession_analysis.R"))
source(paste0(dir.modules,"/module_recessions_analysis_plots.R"))
source(paste0(dir.modules,"/module_sediment_transport.R"))
source(paste0(dir.modules,"/module_sediment_transport_plots.R"))
source(paste0(dir.modules,"/module_retrospective_analysis.R"))
source(paste0(dir.modules,"/module_retrospective_analysis_plots.R"))
# source(paste0(dir.modules,"/Real_time_plots.R"))





print(station.name)
message("Period of study:           ")
print(data.period)
message("
- Reading input data and general options. 
  Please, wait ..."); flush.console()







if (file_gaugings        != "synthet_gaugings.csv") {
   dir.create(paste0(dir.case_study,"/Results"))   #creation of the folder with Results within the case study folder.
}













##############################################################################
#              Limnigraph (stage record, water level): read stage data h(t):
##############################################################################
if (file_limni != FALSE) {
  
  #limni <- read.csv2(paste(dir.case_study,"/",file_limni,sep=""),  header= TRUE)
  
  limni <- read.csv2(paste0(dir.case_study,"/",file_limni),fileEncoding="UTF-8", 
                     quote="", sep=";", dec=".", header= TRUE, na.strings=c(",", " ", "NA") )
  
  
  limni = na.omit(limni)
  
  limni = limni[seq(1, length(limni[,1]), limni_filter), ] 

  if (tLimni.type =="date") {
  #**************************
  t_limni          = limni[,tLimni.col]
  origin.numeric   = as.numeric(as.POSIXct(as.POSIXlt(0, origin = date_origin)))/86400
  t_limni.numeric  = 0; t_limni.date=0;
  message("- Converting stage record dates (if <1970).") 
  message("  Please wait ...")
  pb <- txtProgressBar(min = 0,               # Minimum value of the progress bar
                       max = length(t_limni), # Maximum value of the progress bar
                       style = 3,             # Progress bar style (also available style = 1 and style = 2)
                       width = 50,            # Progress bar width. Defaults to getOption("width")
                       char = "=")            # Character used to create the bar
  
  for (i in 1:length(t_limni)){
    t_limni.date[i]    = as.character(as.POSIXct(as.POSIXlt(t_limni[i], origin = date_origin, 
                                                            tz="Europe/Paris", tryFormats = c("%Y-%m-%d %H:%M:%S",
                                                                                              "%Y-%m-%d %H:%M",
                                                                                              "%Y-%m-%d",
                                                                                              
                                                                                              "%d-%m-%Y %H:%M:%S",
                                                                                              "%d-%m-%Y %H:%M",
                                                                                              "%d-%m-%Y",
                                                                                              
                                                                                              "%d/%m/%Y %H:%M:%S",
                                                                                              "%d/%m/%Y %H:%M",
                                                                                              "%d/%m/%Y",
                                                                                              
                                                                                              "%Y/%m/%d %H:%M:%S",
                                                                                              "%Y/%m/%d %H:%M",
                                                                                              "%Y/%m/%d"
                                                                                              )),
                                                 tryFormats = c("%Y-%m-%d",
                                                                "%d/%m/%Y %H:%M"))) 
    t_limni.numeric[i] = as.numeric(as.POSIXct(as.POSIXlt(t_limni[i], origin = date_origin, 
                                                          tz="Europe/Paris", tryFormats = c("%Y-%m-%d %H:%M:%S",
                                                                                            "%Y-%m-%d %H:%M",
                                                                                            "%Y-%m-%d",
                                                                                            
                                                                                            "%d-%m-%Y %H:%M:%S",
                                                                                            "%d-%m-%Y %H:%M",
                                                                                            "%d-%m-%Y",
                                                                                            
                                                                                            "%d/%m/%Y %H:%M:%S",
                                                                                            "%d/%m/%Y %H:%M",
                                                                                            "%d/%m/%Y",
                                                                                            
                                                                                            "%Y/%m/%d %H:%M:%S",
                                                                                            "%Y/%m/%d %H:%M",
                                                                                            "%Y/%m/%d"
                                                                                            )),
                                               tryFormats = c("%Y-%m-%d",
                                                              "%d/%m/%Y")))/86400 
    setTxtProgressBar(pb, i)
  }  
  close(pb)
  t_limni.numeric2 = t_limni.numeric  - origin.numeric
  
  dt_increase = 0.1     # if limni data have equal timing we can shift the timings. change this if you want to increase the delay [day] 
  while ((length(which(duplicated(t_limni.numeric2)==T))) > 0) {
    t_limni.numeric2[which(duplicated(t_limni.numeric2) == T)] =  t_limni.numeric2[which(duplicated(t_limni.numeric2) == T)] + dt_increase 
  }
  
  
  t_limni          = t_limni.numeric2
  t_limni.true     = t_limni.numeric2

  
  
  } else if (tLimni.type == "numeric"){
  #************************************ 
    dt_increase = 0.1     # if limni data have equal timing we can shift the timings. change this if you want to increase the delay [day] 
    while ((length(which(duplicated(limni[,tLimni.col])==T))) > 0) {
      limni[,tLimni.col][which(duplicated(limni[,tLimni.col]) == T)] =  limni[,tLimni.col][which(duplicated(limni[,tLimni.col]) == T)] + dt_increase 
    }
    
    t_limni      = limni[,tLimni.col] - limni[1,tLimni.col]
    t_limni.true = limni[,tLimni.col]
    
    
    # transforming times in dates:
    t_limni.date = as.Date(floor(t_limni.true), origin = date_origin)
    t_limni.date = as.POSIXct(as.POSIXlt(t_limni.date))
    t_limni.date = t_limni.date  + (t_limni.true - floor(t_limni.true))* 86400      # it is in seconds

  
  } else {
    message("***** Error: tLimni.type not available! *****")
  }
  
  

  #conversion of the stage record to stage in "m":
  if (u.m.limni == "cm") {
    h_limni = limni[,hLimni.col]/100
  } else if (u.m.limni == "mm") {
    h_limni = limni[,hLimni.col]/1000
  } else if (u.m.limni == "m") {
    h_limni = limni[,hLimni.col]
  } else if (u.m.limni == "ft") {
    h_limni = limni[,hLimni.col]*0.3048
  }
  df.limni   = data.frame(t_limni, h_limni,  t_limni.true =t_limni.true)  #dataframe of limni
  nobs.limni = length(h_limni)
  print(paste0("1) Stage record correctly loaded (",file_limni,")"))
  ##################################
  
  
} else {
  # no stage record availability!
  print("1) Stage record not available !") 
  print("   recession analysis cannot be performed !") 
  print("   adjustment of shift times cannot be done !")
  limni        = NULL;
  h_limni      = NULL; 
  t_limni      = NULL;
  df.limni     = NULL;
  t_limni.date = NULL;
  
}
















###################################################################
#                      Gaugings:
###################################################################
if (file_gaugings != FALSE) {
  if (file_gaugings        != "synthet_gaugings.csv") {
    Gaugings.raw <- read.csv2(paste0(dir.case_study,"/",file_gaugings),
                              fileEncoding="UTF-8", quote="",
                              sep=";",dec=".",header=TRUE, na.strings=c(";","NA"))
    #Gaugings.raw = na.omit(Gaugings.raw)
    
  } else {
    Gaugings.raw <- read.csv2(paste0(dir.single.dataset,"/",file_gaugings),
                              fileEncoding="UTF-8", quote="",
                              sep=";",dec=".",header=TRUE, na.strings=c(";","NA"))
    #Gaugings.raw = na.omit(Gaugings.raw)
  }
  Gaugings.raw = Gaugings.raw[seq(1, nrow(Gaugings.raw), gaugings_filter), ]   # uncomment this if you want to select only specific rows

  first.gaug = 1      # to be changed in future versions.
  last.gaug  = "last" # to be changed in future versions.
  if (last.gaug == "last") {
    last = length(Gaugings.raw[,1])
  } else {
    last = last.gaug
  }
  Gaugings = Gaugings.raw[first.gaug:last,]
  
  
  
  
  if (tGaug.type== "date"){
  #************************  
    message("- Converting gaugings dates (if < 1970).") 
    message("  Please wait ...")
    t_gaug           = Gaugings[, tGaug.col]
    origin.numeric   = as.numeric(as.POSIXct(as.POSIXlt(0, origin = date_origin)))/86400
    t_gaug.numeric   = 0; t_gaug.date=c();
    pb <- txtProgressBar(min = 0,               # Minimum value of the progress bar
                         max = length(t_gaug), # Maximum value of the progress bar
                         style = 3,             # Progress bar style (also available style = 1 and style = 2)
                         width = 50,            # Progress bar width. Defaults to getOption("width")
                         char = "=")            # Character used to create the bar.
    for (i in 1:length(t_gaug)){
      t_gaug.date[i]    = as.character(as.POSIXct(as.POSIXlt(t_gaug[i], origin = date_origin, 
                                                             tz="Europe/Paris",
                                                             tryFormats = c("%Y-%m-%d %H:%M:%S",
                                                                            "%Y-%m-%d %H:%M",
                                                                            "%Y-%m-%d",
                                                                            
                                                                            "%d-%m-%Y %H:%M:%S",
                                                                            "%d-%m-%Y %H:%M",
                                                                            "%d-%m-%Y",
                                                                            
                                                                            "%d/%m/%Y %H:%M:%S",
                                                                            "%d/%m/%Y %H:%M",
                                                                            "%d/%m/%Y",
                                                                            
                                                                            "%Y/%m/%d %H:%M:%S",
                                                                            "%Y/%m/%d %H:%M",
                                                                            "%Y/%m/%d"
                                                                            )), tz="Europe/Paris",
                                                  tryFormats = c("%Y-%m-%d",
                                                                 "%d/%m/%Y %H:%M")))
      t_gaug.numeric[i] = as.numeric(as.POSIXct(as.POSIXlt(t_gaug[i], origin = date_origin, 
                                                           tz="Europe/Paris", tryFormats = c("%Y-%m-%d %H:%M:%S",
                                                                                             "%Y-%m-%d %H:%M",
                                                                                             "%Y-%m-%d",
                                                                                             
                                                                                             "%d-%m-%Y %H:%M:%S",
                                                                                             "%d-%m-%Y %H:%M",
                                                                                             "%d-%m-%Y",
                                                                                             
                                                                                             "%d/%m/%Y %H:%M:%S",
                                                                                             "%d/%m/%Y %H:%M",
                                                                                             "%d/%m/%Y",
                                                                                             
                                                                                             "%Y/%m/%d %H:%M:%S",
                                                                                             "%Y/%m/%d %H:%M",
                                                                                             "%Y/%m/%d"
                                                                                             )),
                                                tryFormats = c("%Y-%m-%d",
                                                               "%d/%m/%Y %H:%M")))/86400 
      setTxtProgressBar(pb, i)
    }  
    close(pb)
    t_gaug.numeric2 = t_gaug.numeric  - origin.numeric
    
    
    dt_increase = 0.1     # cif gaugings data have equal timing we can shift the timings. change this if you want to increase the delay [day] 
    while ((length(which(duplicated(t_gaug.numeric2)==T))) > 0) {
      t_gaug.numeric2[which(duplicated(t_gaug.numeric2) == T)] =  t_gaug.numeric2[which(duplicated(t_gaug.numeric2) == T)] + dt_increase 
    }

    t_gaug          = t_gaug.numeric2
    t_gaug.true     = t_gaug.numeric2
    
    
    if (file_limni != FALSE) {
      if (t_gaug.numeric2[1] <= t_limni.numeric2[1]) {
        t_gaug       = t_gaug.numeric2  - t_gaug.numeric2[1]
        t_limni      = t_limni.numeric2 - t_gaug.numeric2[1]
        initial.time = t_gaug.numeric2[1]
        
      } else {
        t_Gaug       = t_gaug.numeric2  - t_limni.numeric2[1]
        t_limni      = t_limni.numeric2 - t_limni.numeric2[1]
        initial.time = t_limni.numeric2[1]
      }
      df.limni   = data.frame(t_limni, h_limni,  t_limni.true = t_limni.true)  #dataframe of limni
    } else if (file_limni == FALSE) {
      t_Gaug         = t_gaug.numeric2  - t_gaug.numeric2[1]
      initial.time   = t_gaug.numeric2[1]
    }
    gaug.date = t_gaug.date
    
    
    
  } else if (tGaug.type =="numeric"){
  ###################################
    # t_gaug        = t_limni.numeric2 - t_limni.numeric2[1]
    t_gaug.true     = Gaugings[,tGaug.col]
    Gaug.time       = Gaugings[,tGaug.col]
    
    if (file_limni != FALSE) {
      if (Gaug.time[1] <= t_limni.true[1]) {
        t_Gaug       = Gaug.time - Gaug.time[1]
        initial.time =  Gaug.time[1]
      } else {
        t_Gaug       = Gaug.time - t_limni.true[1]
        initial.time = t_limni[1]
      }
    } else if (file_limni == FALSE) {
      t_Gaug       = Gaug.time - Gaug.time[1]
      initial.time = Gaug.time[1]
    }
    gaug.date = as.Date(floor(Gaug.time), origin = date_origin)
    gaug.date = as.POSIXct(as.POSIXlt(gaug.date))
    gaug.date = gaug.date  + (Gaug.time - floor(Gaug.time)) * 24 * 3600

    
    
    
  } else {
  ########
    message("***** Fatal Error: 'tGaug.type' not available! *****")
  }
  
  #print(paste0("Initial time of simulation = ", initial.time))
  
  
  Q_Gaug    <- Gaugings[,QGaug.col]
  Q_BaRatin <- matrix(-9999,length(Q_Gaug),1)
  Q_BaRatin[1:length(Q_Gaug)] = Q_Gaug[1:length(Q_Gaug)] 
  if (uQ.absolute == FALSE) {
    uQ_Gaug <- (Gaugings[,uQGaug.col]/100*Q_Gaug)/2
  } else {
    uQ_Gaug <- Gaugings[,uQGaug.col]
  }
  
  if (u.m.Hgaug == "cm") {
    h_Gaug <- Gaugings[,hGaug.col]/100
  } else if (u.m.Hgaug == "mm") {
    h_Gaug <- Gaugings[,hGaug.col]/1000
  } else if (u.m.Hgaug == "m") {
    h_Gaug <- Gaugings[,hGaug.col]
  } else if (u.m.Hgaug == "ft") {
    h_Gaug <- Gaugings[,hGaug.col]*0.3048
  }
  if (u.m.Qgaug == "m^3.s-1") {
    Q_Gaug <- Gaugings[,QGaug.col]
  } else if (u.m.Qgaug == "ft^3.s-1") {
    Q_Gaug <- Gaugings[,QGaug.col]*0.028317
  }
  gaug.df <- data.frame(t_Gaug,h_Gaug,Q_Gaug,uQ_Gaug)   # gaugings dataframe
  
  
  
  data4BaRatin <- data.frame("h"      = h_Gaug,
                             "Q"      = Q_BaRatin,
                             "uQ"     = uQ_Gaug,
                             "Period" = 1 , 
                             "t"      = t_Gaug,
                             "t.true" = t_gaug.true)
  nobs.gaug = length(Q_Gaug)
  print(paste0("2) Gaugings correctly loaded (",file_gaugings,")"))
  
  
  
  
  
  #The stage-discharge plot (discharge in log scale):
  #**************************************************
  gaugings.save  =  "gaugings" 
  print("   Plotting stage-discharge gaugings.")
  gaugings.plot = ggplot(data= data4BaRatin) + 
    coord_trans(xlim = grid_RC.xlim, 
                ylim = grid_RC.ylim.log) +
    scale_y_log10(na.value = -10,  
                  breaks   = ticks_RC.y.log,   
                  labels   = ticks_RC.y.log) +
    annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black", size = 0.5, linetype = 1)+
    scale_x_continuous(breaks = seq(grid_RC.xlim[1], grid_RC.xlim[2],  grid_RC.xstep)) +  
    #theme_light(base_size = 15)+
    ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
    xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]")))+
    geom_point(aes(x = h , y= Q), fill ="blue", pch=21, size = 1.5)+
    geom_errorbar(aes(x = h,  ymin =Q-2*uQ , ymax =Q +2*uQ), color= "blue", width=.05, size = 0.3) +
    theme_bw(base_size=10)+
    theme(  axis.text        = element_text(size=15)
            ,axis.title       = element_text(size=15, face="bold")
            ,panel.grid.major = element_line(size=1)
            ,legend.text      = element_text(size=20)
            ,legend.title     = element_text(size=20)
            ,legend.key.size  = unit(1.5, "cm")
            ,legend.position  = "none")
  #scale_colour_gradientn(colors=rainbow(7))
  
  ggsave(gaugings.plot, filename =paste0(dir.case_study,"/", gaugings.save, ".png"), bg = "transparent", width = 10, height =6, dpi = 200)
  plot(gaugings.plot)
  
} else {
  # no gaugings availability!
  print("2) Gaugings not available !")
  print("   segmentation of gaugings cannot be performed !")
  data4BaRatin   = NULL;
  h_Gaug         = NULL; 
  Q_Gaug         = NULL;
  uQ_Gaug        = NULL;
  t_Gaug         = NULL;
  Gaugings       = NULL;
}





















###########################################################################################
#                       Official dates of RC update:
###########################################################################################
if (official.shift.times != FALSE) {
  if (file_gaugings        != "synthet_gaugings.csv") {
    
    officialShifts <- read.csv2(paste0(dir.case_study, "/",  official.shift.times),
                                fileEncoding="UTF-8-BOM", quote="", sep= ";", dec=".", header=TRUE)
    
    if (tOfficial.type== "date"){
      #************************  
      t_official       = officialShifts[,tOfficial.col]
      origin.numeric   = as.numeric(as.POSIXct(as.POSIXlt(0, origin = date_origin)))/86400
      t_limni.numeric  = 0; t_limni.date=0;
      message("- Converting gaugings dates (if < 1970).") 
      message("  Please wait ...")
      t_gaug           = Gaugings$Date
      origin.numeric   = as.numeric(as.POSIXct(as.POSIXlt(0, origin = date_origin)))/86400
      t_gaug.numeric   = 0; t_gaug.date=c();
      # pb <- txtProgressBar(min = 0,               # Minimum value of the progress bar
      #                      max = length(t_gaug), # Maximum value of the progress bar
      #                      style = 3,             # Progress bar style (also available style = 1 and style = 2)
      #                      width = 50,            # Progress bar width. Defaults to getOption("width")
      #                      char = "=")            # Character used to create the bar.
      # for (i in 1:length(t_gaug)){
      #   t_gaug.date[i]    = as.character(as.POSIXct(as.POSIXlt(t_gaug[i], origin = date_origin, 
      #                                                          tz="Europe/Paris", tryFormats = c("%Y-%m-%d %H:%M:%S",
      #                                                                                            "%Y-%m-%d %H:%M",
      #                                                                                            "%Y-%m-%d",
      #                                                                                            
      #                                                                                            "%d-%m-%Y %H:%M:%S",
      #                                                                                            "%d-%m-%Y %H:%M",
      #                                                                                            "%d-%m-%Y",
      #                                                                                            
      #                                                                                            "%d/%m/%Y %H:%M:%S",
      #                                                                                            "%d/%m/%Y %H:%M",
      #                                                                                            "%d/%m/%Y",
      #                                                                                            
      #                                                                                            "%Y/%m/%d %H:%M:%S",
      #                                                                                            "%Y/%m/%d %H:%M",
      #                                                                                            "%Y/%m/%d"
      #                                                          )),
      #                                               tz="Europe/Paris",
      #                                               tryFormats = c("%Y-%m-%d",
      #                                                              "%d/%m/%Y")))
      #   t_gaug.numeric[i] = as.numeric(as.POSIXct(as.POSIXlt(t_gaug[i], origin = date_origin, 
      #                                                        tz="Europe/Paris", tryFormats = c("%Y-%m-%d %H:%M:%S",
      #                                                                                          "%Y-%m-%d %H:%M",
      #                                                                                          "%Y-%m-%d",
      #                                                                                          
      #                                                                                          "%d-%m-%Y %H:%M:%S",
      #                                                                                          "%d-%m-%Y %H:%M",
      #                                                                                          "%d-%m-%Y",
      #                                                                                          
      #                                                                                          "%d/%m/%Y %H:%M:%S",
      #                                                                                          "%d/%m/%Y %H:%M",
      #                                                                                          "%d/%m/%Y",
      #                                                                                          
      #                                                                                          "%Y/%m/%d %H:%M:%S",
      #                                                                                          "%Y/%m/%d %H:%M",
      #                                                                                          "%Y/%m/%d"
      #                                                        )),
      #                                             tryFormats = c("%Y-%m-%d",
      #                                                            "%d/%m/%Y")))/86400 
      #   setTxtProgressBar(pb, i)
      # }  
      # close(pb)
    }
    
    
    
    
    if (!is.null(Gaugings)) {
      if (!is.null(t_limni)){
      if (Gaugings[1,tGaug.col] <= t_limni[1]) {
        officialShiftsTime <- officialShifts[,1] - Gaugings[1,tGaug.col]
      } else {
        officialShiftsTime <- officialShifts[,1] - limni[1,tLimni.col]
      }
      } else {
        officialShiftsTime <- officialShifts[,1] - Gaugings[1,tGaug.col]
      }
    } else {
      officialShiftsTime <- officialShifts[,1] - limni[1,tLimni.col]
    }
    
  } else {
    officialShifts <- read.csv2(paste0(dir.single.dataset ,"/",  official.shift.times),
                                fileEncoding="UTF-8", quote="", sep=";", dec=".", header=TRUE)
    officialShiftsTime <- officialShifts[,1] - Gaugings[1,tGaug.col]
  }
  
  
  
  
  
  if (!is.null(officialShiftsTime)) {
    data.annotate.off =  data.frame(xeffect = officialShiftsTime, 
                                    xpotent = officialShiftsTime)
  } else {
    data.annotate.off = NULL
  }
  
  
  
  
  print(paste0("3) Official shift times correctly loaded (",official.shift.times,")"))
  ###################  
} else {
  ###################
  print("3) Official shift times not available!")
  print("   Performance evaluation not possible !!!")
  officialShiftsTime = NULL
  data.annotate.off  = NULL
}
















################################################################################################
# Plot the stage record with the gauged stage:
################################################################################################
limni.labels            =  c("Time [days]", "Stage h [m]") 
limni.save              =  "Stage_record" # Name of the study Limni (for plots)

if (!is.null(df.limni)) {
  print("   Plotting stage record with gaugings.")
  limni.plot <- ggplot() +
    geom_line(aes(x=t_limni, y=df.limni$h_limni), color = "lightblue",    size =0.2)+
    #scale_x_continuous(expand=c(0,0))+
    #scale_x_date(expand=c(0,0))+
    labs()+
    xlab(limni.labels[1]) +
    ylab(limni.labels[2])+
    geom_point(aes(x = t_Gaug,  y = h_Gaug) , pch=21,    fill= "blue",   size = 2) + 
    # geom_point(aes(x = gaug.date[1:index.lastgauging.retro], 
    #                y= h_Gaug[1:index.lastgauging.retro]) ,
    #            pch=21, fill= "blue", size = 3)+
    # geom_point(aes(x = gaug.date[index.lastgauging.retro+1:length(t_Gaug)],
    #                y = h_Gaug[index.lastgauging.retro+1:length(t_Gaug)]) ,
    #            pch=21, fill= "green", size = 3)+
    # annotate("text", x= c(limni.date[index.start.realtime/2],
    #                       limni.date[index.start.realtime+(length(limni.date) - index.start.realtime+1)/2],
    #                       limni.date[index.start.realtime +1]),
    #          y= c(grid_limni.ylim[2], grid_limni.ylim[2], grid_limni.ylim[2]), 
    #          label= c("Retrospective \n Analysis", "Real Time", "Initialisation \n t0"),
  #          color = c("blue", "green", "red"), size=5) +
  coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text         = element_text(size=15)
          ,axis.title       = element_text(size=15)
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,legend.text      = element_text(size=20)
          ,legend.title     = element_text(size=30)
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "none"
          ,plot.margin      = unit(c(0.5,0.5,0.2, 1),"cm"))
  
  if (!is.null(data.annotate.off)) {
    limni.plot = limni.plot+ 
      geom_vline(xintercept = data.annotate.off$xeffect, color="red", size=0.5, linetype="dashed")
    
  }
  # limni.plot = limni.plot+
  #   geom_vline(xintercept = limni.date[index.start.realtime], 
  #              color="red", size=1.5, linetype="dashed")
  
  ggsave(limni.plot, filename =paste0(dir.case_study,"/",limni.save,".png"),
         bg = "transparent", width = 14, height =6, dpi = 200)
  plot(limni.plot)
}     

























########################################################################
#                   grid limits for plots :
########################################################################
print("4) Preparing grids for plots. Please wait ...")
grid.RC  = "Manual"  

# grid for x axis:
if (is.null(h_limni) == FALSE) {
  Hmin = min(c(min(h_limni, na.rm= TRUE), min(h_Gaug, na.rm= TRUE)),  na.rm= TRUE)
  Hmax = max(c(max(h_limni, na.rm= TRUE), max(h_Gaug, na.rm= TRUE)),  na.rm= TRUE)
} else {
  Hmin = min(h_Gaug,  na.rm= TRUE)
  Hmax = max(h_Gaug,  na.rm= TRUE)
}
generic.grid = seq(-1000000,1000000, 0.25)

for (i in 1:length(generic.grid)) {
  if((Hmin >= generic.grid[i]) & (Hmin <= generic.grid[i+1])) {
    if (grid.RC == "Manual") {
      Hmin_grid = grid_RC.xlim[1]
    } else { 
      Hmin_grid = generic.grid[i]
    }
  }
}
for (i in 1:length(generic.grid)) {
  if((Hmax >= generic.grid[i]) & (Hmax <= generic.grid[i+1])) {
    if (grid.RC == "Manual") {
      Hmax_grid = grid_RC.xlim[2]
    } else { 
      Hmax_grid = generic.grid[i+1]
    }
  }
}  




# grid for y axis:
if (!is.null(Q_Gaug)) {
  Qmin = min(Q_Gaug,  na.rm= TRUE)
  Qmax = max(Q_Gaug,  na.rm= TRUE) + max(uQ_Gaug,  na.rm= TRUE)
  generic.grid.Qlog = c(0.0001, 0.001, 0.01,0.1, 1, 10,100, 1000, 10000, 100000)
  generic.grid.Q = seq(0,100000, 100)
  for (i in 1:length(generic.grid.Qlog)) {
    if((Qmin >= generic.grid.Qlog[i]) & (Qmin <= generic.grid.Qlog[i+1])) {
      Qmin_grid.log = generic.grid.Qlog[i]
    }
  }
  for (i in 1:length(generic.grid.Qlog)) {
    if((Qmax >= generic.grid.Qlog[i]) & (Qmax <= generic.grid.Qlog[i+1])) {
      Qmax_grid.log = generic.grid.Qlog[i+1]
    }
  }
  for (i in 1:length(generic.grid.Q)) {
    if((Qmax >= generic.grid.Q[i]) & (Qmax <= generic.grid.Q[i+1])) {
      if (grid.RC == "Manual") {
        Qmax_grid = grid_RC.ylim[2]
      } else { 
        Qmax_grid = generic.grid.Q[i+1]
      }
    }
  }
  ticks.log <- generic.grid.Qlog[which(generic.grid.Qlog >= Qmin_grid.log  & generic.grid.Qlog <= Qmax_grid.log)]
  
  #time grid:
  time.grid = seq(0, 10000000, 1000)
  for (i in 1:length(time.grid)) {
    if((time.grid[i+1] > tail(t_Gaug,1) & (time.grid[i] < tail(t_Gaug,1)))) {
      tmax_grid = time.grid[i+1]
    }
  }
  
  
} else {
  # no gaugings:
  Qmin_grid         = grid_RC.ylim[1]
  Qmax_grid         = grid_RC.ylim[2]
  Qmin_grid.log     = grid_RC.ylim.log[1]
  Qmax_grid.log     = grid_RC.ylim.log[2]
  ticks.log         = ticks_RC.y.log       = c(0.1, 1, 10, 100, 1000, 10000)    # grid for RC-log plots
  time.grid = seq(0, 10000000, 1000)
  for (i in 1:length(time.grid)) {
    if((time.grid[i+1] > tail(df.limni$t_limni,1) & (time.grid[i] < tail(df.limni$t_limni,1)))) {
      tmax_grid = time.grid[i+1]
    }
  }
}




















#########
# COLORS:
#########
# generation of colors for the plots: 
#each color will represent a distinct period of RC stability.
################################
color.generation = function(n) {
  ################################
  # Generate the vector "colo" with the colors names for the different periods.
  # Notice that the colors need a high contrast !!
  colo.max = rep(c("red","blue","chartreuse", "orange","darkviolet",  "forestgreen",  
                   "cyan", "saddlebrown","pink", "black",  "yellow2",
                   "navyblue", "magenta", "lightsteelblue4", "khaki", "darkolivegreen", 
                   "lightsalmon", "mediumvioletred", "orange3", "indianred4", "brown4",
                   "black" , "green", "lightblue3", "hotpink4", "lightgoldenrod4"), 100)  # max 100x26 colors !!
  colo=colo.max[1:n]
  return(colo)
}
max.number.periods = 500   # maximum number of colors for the periods and plots !!!!
colo               = color.generation(n = max.number.periods)  






############
# RC PRIORS:
############
# Reading priors for BaRatin: RC estimation
message("You have selected the following RC priors:")
message("- Hydraulic configuration:")
print(c(paste0("  ", ncontrols, " controls defined: "), control.type))
print("  Matrix of controls:") 
print(M)
Sys.sleep(1.5)


#Parameter a:
message("- RC parameter 'a':")
if (a.distr== "LogNormal"){
  adistr = "LN"
} else if (a.distr== "Gaussian"){
  adistr = "N"
} else if (a.distr== "Uniform"){
  adistr = "U"
} else {
  adistr = "pdf"
} 

if (propagat == TRUE){
  cat(paste0("  Propagate the priors of geometric/physical properties: \n",
             "  (e.g. Channel width, longitudinal slope, Roughness coeff.)\n"))
  cat(paste0("  Prior pdf of parameter 'a' = ", a.distr, "\n  **************************************\n"))
  
  for (c in 1:ncontrols){
    cat(paste0("  Priors of control ",c,":\n"))
    if (control.type[c] == "rect.channel") {
      if (a.distr== "LogNormal"){
        # Transform priors from gaussian to lognormal:
        # function which takes gaussian mean and stdev and gives lognormal mean and stdev
        st_Bc.prior[c] = round(Transf_Gauss_lognorm(E = Bc.prior[c], stdev = st_Bc.prior[c])$sd, digits = 2)
        st_KS.prior[c] = round(Transf_Gauss_lognorm(E = KS.prior[c], stdev = st_KS.prior[c])$sd, digits = 2)
        st_S0.prior[c] = round(Transf_Gauss_lognorm(E = S0.prior[c], stdev = st_S0.prior[c])$sd, digits = 2)
      }
      cat(paste0("  Bc = ",adistr, "(",  Bc.prior[c], ", ", st_Bc.prior[c], ")\n"))
      cat(paste0("  Ks = ",adistr, "(",  KS.prior[c], ", ", st_KS.prior[c], ")\n"))
      cat(paste0("  S0 = ",adistr, "(",  S0.prior[c], ", ", st_S0.prior[c], ")\n\n"))       
    } else if (control.type[c] == "rect.weir") {
      if (a.distr== "LogNormal"){
        # Transform priors from gaussian to lognormal:
        # function which takes gaussian mean and stdev and gives lognormal mean and stdev
        st_Cr.prior[c] = round(Transf_Gauss_lognorm(E = Cr.prior[c], stdev = st_Cr.prior[c])$sd, digits = 2)
        st_g.prior[c]  = round(Transf_Gauss_lognorm(E = g.prior[c],  stdev = st_g.prior[c])$sd,  digits = 2)
        st_Bw.prior[c] = round(Transf_Gauss_lognorm(E = Bw.prior[c], stdev = st_Bw.prior[c])$sd, digits = 2)
      }
      cat(paste0("  Cr = ",adistr, "(",  Cr.prior[c], ", ", st_Cr.prior[c], ")\n"))
      cat(paste0("  g  = ",adistr, "(",  g.prior[c],  ", ", st_g.prior[c],  ")\n"))
      cat(paste0("  Bw = ",adistr, "(",  Bw.prior[c], ", ", st_Bw.prior[c], ")\n\n"))       
    } else {
      message("**** Fatal Input error: you have selected a wrong control type!")
      message("control types available: 'rect.channel' and 'rect.weir'.")
      message("Please, check again!")
    }
  }
  
} else {
  cat(paste0("  Prior pdf of parameter 'a' = ", a.distr, "\n  **************************************\n"))
  for (c in 1:ncontrols){
    cat(paste0("  a",c," = "))
    if (a.distr== "LogNormal"){
      # Transform priors from gaussian to lognormal:
      # function which takes gaussian mean and stdev and gives lognormal mean and stdev
      st_a.prior[c] = round(Transf_Gauss_lognorm(E = a.prior[c], stdev = st_a.prior[c])$sd, digits = 4)
    }
    cat(paste0(adistr, "(",  a.prior[c], ", ", st_a.prior[c], ")\n"))
  }
} 

Sys.sleep(1.5)





#Parameter b:
message("- RC parameter 'b':")
if (b.distr== "LogNormal"){
  bdistr = "LN"
} else if (b.distr== "Gaussian"){
  bdistr = "N"
} else if (b.distr== "Uniform"){
  bdistr = "U"
} else {
  bdistr = "?pdf"
} 
cat(paste0("  Prior pdf of parameter 'b' = ", b.distr, "\n  **************************************\n"))
for (c in 1:ncontrols){
  cat(paste0("  b",c," = "))
  if (b.distr== "LogNormal"){
    # Transform priors from gaussian to lognormal:
    # function which takes gaussian mean and stdev and gives lognormal mean and stdev
    st_b.prior[c] = round(Transf_Gauss_lognorm(E = b.prior[c], stdev = st_b.prior[c])$sd, digits = 4)
  }
  cat(paste0(bdistr, "(",  b.prior[c], ", ", st_b.prior[c], ")\n"))
}
Sys.sleep(1.5)





#Parameter c:
message("- RC parameter 'c':")
if (c.distr== "LogNormal"){
  cdistr = "LN"
} else if (c.distr== "Gaussian"){
  cdistr = "N"
} else if (c.distr== "Uniform"){
  cdistr = "U"
} else {
  cdistr = "?pdf"
} 
cat(paste0("  Prior pdf of parameter 'c' = ", c.distr, "\n  **************************************\n"))
for (c in 1:ncontrols){
  cat(paste0("  c",c," = "))
  if (c.distr== "LogNormal"){
    # Transform priors from gaussian to lognormal:
    # function which takes gaussian mean and stdev and gives lognormal mean and stdev
    st_c.prior[c] = round(Transf_Gauss_lognorm(E = c.prior[c], stdev = st_c.prior[c])$sd, digits = 4)
  }
  cat(paste0(cdistr, "(",  c.prior[c], ", ", st_c.prior[c], ")\n"))
}

Sys.sleep(1)



# #Remnant error model:
# #--------------------
# #Parameters gamma1 and gamma2 (structural error model):
# message("- RC parameters 'gamma1 and gamma2':")
# if (g1.distr.type== "LogNormal"){
#   g1distr = "LN"
# } else if (g1.distr.type== "Gaussian"){
#   g1distr = "N"
# } else if (g1.distr.type== "Uniform"){
#   g1distr = "U"
# } else {
#   g1distr = "?pdf"
# } 
# if (g2.distr.type== "LogNormal"){
#   g2distr = "LN"
# } else if (g2.distr.type== "Gaussian"){
#   g2distr = "N"
# } else if (g2.distr.type== "Uniform"){
#   g2distr = "U"
# } else {
#   g2distr = "?pdf"
# } 
# 
# cat(paste0("  Prior pdf of parameter 'gamma1' = ", g1.distr.type, "\n  ******************************************\n"))
# cat(paste0("  gamma1"," = "))
# cat(paste0(g1distr, "(",  g1.prior[1], ", ", g1.prior[2], ")\n"))
# 
# 
# 
# cat(paste0("  gamma2"," = "))
# cat(paste0(g1distr, "(",  g2.prior[1], ", ", g2.prior[2], ")\n"))  
# 


message("- Creating a new folder 'Results' within the case study")
message("  folder that will contain all results of the analysis.")
dir.create(paste0(dir.case_study,"/Results"))




message("
############################################
#              All done !                  #  
############################################
# Check the two plots just created.        #
# If everything looks fine, you can proceed# 
# with the application of the detection    #
# methods:                                 #
# 1) Segmentation of gaugings              #
# 2) Stage-recession analysis              #
# 3) sediment transport proxy analysis.    #
############################################

")



