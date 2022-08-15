#############################################################################################################
# Functions for the retrospective analysis:
#############################################################################################################
# Author : MAtteo Darienzo ,
# Organisation : INRAE
# Date of last update : 08/10/2020
#
# Module for recession analysis. It contains the functions necessary to extract and estimate 
# the stage recessions from the stage record.




#######################################################
# FINDING THE ALL THE LoCAL MINIMUM (when h(i)<h(i-1)):
#######################################################
localmin <- function(t, h, toll) {
###########################
  hmin = 0; 
  tmin = 0  ;
  hmin[1] = h[1]; 
  tmin[1] = t[1]; 
  tmax = length(t) ; 
  i = 1; k = 1;
  while (i < tmax) {
    i = i+1
    if (h[i] <= h[i-1] + toll) {
      k = k + 1 
      hmin[k] = h[i]
      tmin[k] = t[i]
    }
  }
  min_loc2 = data.frame(tmin, hmin)
  return(min_loc2)
}





#############################################################################
# FINDING THE ALL THE MINIMUM of THE MINIMUM (when hmin(i)< hmin(i-1)):  hrec
rec <- function(tmin, hmin, chi, delta.t.max,  delta.t.min, toll) {
#############################################################################
  # initialize:
  hrec       = trec =  trecmax = hrecmax = 0 
  hrec[1]    = hmin[1] 
  trec[1]    = tmin[1]
  hrecmax[1] = hmin[1]
  trecmax[1] = tmin[1]
  j          = 2 
  m          = 1 
  k          = 1
  min_loc3   = data.frame(tmin,  hmin)
  elim.index = c()
  
  
  for (ccc in 2:length( min_loc3$tmin)) {
    if (( min_loc3$tmin[ccc] -  min_loc3$tmin[ccc-1]) < delta.t.min) {
      elim.index = c(elim.index, ccc)
    }
  }
  if (length(elim.index) > 1) {
    min_loc4                = min_loc3[-elim.index,]
  } else {
    min_loc4                = min_loc3
  }
  
  j_max = length(min_loc4$tmin) 
  
  # loop in time series of "tmin":
  while (j <= j_max) {
    # print(paste0('timestep = ', tmin[j] - trec[m]))
    # print(paste0('jump =     ', hmin[j] - hrec[m]))
    
    if (min_loc4$hmin[j] <= hrec[m] + toll) {
      # ok, recession point
      m = m + 1 
      hrec[m] = min_loc4$hmin[j]
      trec[m] = min_loc4$tmin[j]
      
    } else {

      if (((abs(min_loc4$hmin[j] - hrec[m])) >= chi ) | ((min_loc4$tmin[j]-trec[m]) >= delta.t.max)) { # stop the recession
        #print("step > deltamax or > chi")
        m = m + 1
        k = k + 1
        hrec[m] = hmin[j]
        trec[m] = tmin[j]
        hrecmax[k] = hrec[m]   #not used yet
        trecmax[k] = trec[m]   #not used yet
      }
    }
    j = j+1
  }
  min_loc5 = data.frame(trec,  hrec)

  return(min_loc5)
}








###############################################################
# Extract all recession curves:
extract_curve <- function(trec, hrec, chi, delta.t.max, toll) {
###############################################################
  # max number of recessions:
  tmaxrec = length(trec) 
  
  # initialize:
  RCi = 0 ; hcurve=0 ; tcurve =0 ; 
  p = 1 ; l = 1 ; n = 1 ; RCi = 1 ; 
  hfin = 0 ; tfin = 0;
  nc = 0 ; ts = 0 ; s = 1 ; 
  hpeak = 0  ; tpeak = 0 ; t.real =0;
  hcurve[1] = hrec[1];
  tcurve[1] = trec[1];
  tpeak[1]  = trec[1]; 
  hpeak[1]  = hrec[1];
  tfinal=0; hfinal=0; tf.real=0;
  #gradient.max = -0.5
  
  
  #loop on trec:
  while (l < tmaxrec) {
    l = l+1  # index of recession dataset
    
    if ((hrec[l] <= hrec[l-1] + toll)  & (abs(hrec[l] - hrec[l-1]) < chi) &  (trec[l] - trec[l-1] < delta.t.max)) {
        # ongoing recession:
        n = n +1
        nc = nc +1
        hcurve[l]  =  hrec[l]
        tcurve[l]  =  trec[l] - trec[l-1] + tcurve[l-1]
        t.real[l]  =  trec[l]
        RCi[l]     =  RCi[l-1]
        tpeak[l]   =  tpeak[l-1]
        hpeak[l]   =  hpeak[l-1]
        tfinal[l]  =  0
        hfinal[l]  =  0
        tf.real[l] =  0
      
    } else {
      # stop ongoing recession ==> new recession !!!
      p  = p+1   
      nc = 0
      
      # # remove all values with angular coeff high:
      # eliminate.index2 =0;
      # gradient.max = -0.5 # max angular coefficient coefficient
      # for (ccc in 1:(length( min_loc4$trec)-5)) {
      #   if ((min_loc4$hrec[ccc] - mean(min_loc4$hrec[ccc:(ccc+5)]))/(min_loc4$trec[ccc+5] - min_loc4$trec[ccc]) > 
      #       gradient.max) {
      #     eliminate.index2 = c(eliminate.index2, ccc)
      #   }
      # }
      # if (length(eliminate.index2) > 1) {
      #   min_loc5                = min_loc4[-eliminate.index2,]
      # } else {
      #   min_loc5 = min_loc4
      # }
      # 
      hcurve[l]    = hrec[l]
      tcurve[l]    = 0
      t.real[l]    = trec[l]
      hpeak[l]     = hrec[l]
      tpeak[l]     = trec[l]
      tfinal[l-1]  = tcurve[l-1]
      hfinal[l-1]  = hrec[l-1]
      tf.real[l-1] = hrec[l-1]
      tfinal[l]    = 0
      hfinal[l]    = 0
      tf.real[l]   = 0
      RCi[l]       = p
    }
  }
  curves = data.frame(tcurve, 
                      hcurve, 
                      RCi, 
                      t.real,
                      tpeak, 
                      hpeak,
                      tfinal,
                      hfinal,
                      tf.real)
  return(curves)
}


















#########################################################################################################
recession.selection <- function(  dir.exe,
                                  dir.segment.rec,
                                  file.options.general,
                                  file.options.recess,
                                  stage.record
                                  ) {
##########################################################################################################
  
  if (is.null(stage.record)){
    message("***************************************************************************")
    message("***** ERROR: stage record (with times and values), 'df.limni' is required !!! ")
    message("***************************************************************************")
  }
  
  # Read inputs and options for computation:
  source(file.options.general)
  source(file.options.recess)
  ############################
  # other hardcoded settings:
  index.period           = c(1:length(stage.record$t_limni))   # indexes of stage record to analyse
  stage.limni.u.m.       = "m"                                 # "m" or "cm"  unity of the stage record (df.limni$h_limni)
  tmax.rec               = 150                                 # max recession time (for plot) NOT USED !
  hmax.rec               = 100                                 # max recession stage (for plot) NOT USED !
  #limits.x.recess        = c(0, 150, 30)                       # [c(min, max, step)] limits for the recession period in days. by default = c(0, 150, 30).
  stage.scale.shift      = 1000                               # is in [cm]: parameter used to shift stage values and avoid negative values: h = h + stage.scale.shift
  BayesianOption         = 2
  toll                   = 0.1   # [cm] tolerance for the recession extraction. 
  #############################
  
  # directories:
  dir.segment.rec.test1  = paste0(dir.segment.rec,"/",name.folder.results.recession)
  dir.extraction         = paste0(dir.segment.rec.test1,"/1_curves_extraction")
  dir.BaM.rec            = paste0(dir.exe,"/Recession_h") 
  dir.BaM.rec.pool       = paste0(dir.exe,"/Recession_h_pooling") 
  dir.create(paste0(dir.segment.rec,"/", name.folder.results.recession))
  dir.create(paste0(dir.segment.rec.test1,"/1_curves_extraction"))
  
  # saving file with the selected options for the current test:
  file.options = paste0(dir.segment.rec.test1,"/options.txt")
  cat(paste0("stage.limni.u.m. =" , stage.limni.u.m.,  "        # unity of the stage record (df.limni$h_limni)"), file = file.options,  sep="\n")   
  cat(paste0("stage.scale.shift =", stage.scale.shift, "      # parameter used to shift stage values and avoid negative values [cm]:"), file = file.options, append = TRUE, sep="\n")
  cat(paste0("uh.rec =",            uh.rec,            "             # stdev of the stage obs error [in cm]"),         file = file.options, append = TRUE, sep="\n")  
  cat(paste0("tburn.rec =",         tburn.rec,         "          # discard the first n days of the studied recession curve"),    file = file.options, append = TRUE, sep="\n")        
  cat(paste0("Nmin.rec =",          Nmin.rec,          "          # min number of data in a recession curve"),      file = file.options, append = TRUE, sep="\n")     
  cat(paste0("tgood =",             tgood,             "             # min length of the recession in days"),          file = file.options, append = TRUE, sep="\n")         
  cat(paste0("delta.t.min =",       delta.t.min,       "          # min days between two recess data"),             file = file.options, append = TRUE, sep="\n")        
  cat(paste0("delta.t.max =",       delta.t.max,       "          # max days between two recess data"),             file = file.options, append = TRUE, sep="\n")          
  cat(paste0("chi =",               chi,               "                    # max stage rise between two recess data"),       file = file.options, append = TRUE, sep="\n")     
  cat(paste0("index.period =",      paste(index.period[1], tail(index.period,1), sep=";"),  "     # indexes of stage record to analyse"), file = file.options, append = TRUE, sep="\n")   
  cat(paste0("BayesianOption =",    BayesianOption,    "        # Type of Bayesian model (only 2 is available!)"), file = file.options, append = TRUE, sep="\n")   
  # cat("********************************************",  file = file.options, append = TRUE, sep="\n")
  # cat(paste0("recession model =",   rec.model),        file = file.options, append = TRUE, sep="\n") # recession model: "exp", "expexp", "hyperb", etc.
  # cat(paste0("jump.pos/neg.rec=",   paste(jump.pos.rec, jump.neg.rec), sep=";"), file = file.options, append = TRUE, sep="\n") 
  # cat(paste0("Ncycles.mcmc.rec=",   Ncycles.mcmc.rec), file = file.options, append = TRUE, sep="\n") 
  # cat(paste0("nmcmc.rec=",          nmcmc.rec),        file = file.options, append = TRUE, sep="\n") 
  # cat(paste0("nslim.rec=",          nslim.rec),  file = file.options, append = TRUE, sep="\n") 
  # cat(paste0("prior.gamma.rec =",   paste(prior.gamma.rec, sep=";")),  file = file.options, append = TRUE, sep="\n") 
  # cat(paste0("prior.param.rec=",    prior.param.rec),  file = file.options, append = TRUE, sep="\n") 

  
  
  
  #*********************************************************************************************************
  # ! stage.scale.shift is in cm !
  # we need to be sure that stage record is always positive, because the remnant error model associated 
  # to the recession model (in BaM) is linearly dependent on stage: err = g1 + g2*h
  # for this reason we shift upward the stage record of a quantity "stage.scale.shift" = 10000 cm.
  
  # now shift and convert to cm:
  if (stage.limni.u.m. == "cm"){
      stage.rec.df.raw = data.frame(t = stage.record$t_limni[index.period],
                                    h = stage.record$h_limni[index.period] +  stage.scale.shift)
  } else if (stage.limni.u.m. == "m") {
      stage.rec.df.raw = data.frame(t = stage.record$t_limni[index.period],
                                    h = stage.record$h_limni[index.period]*100 +  stage.scale.shift)
  } else if (stage.limni.u.m. == "mm") {
      stage.rec.df.raw = data.frame(t = stage.record$t_limni[index.period],
                                    h = stage.record$h_limni[index.period]/10 +  stage.scale.shift)
  } else {
      print("You have inserted a wrong  u.m. for the stage record, put 'm' or 'cm' or 'mm'.")
  }
  
  # check for na:
  if (any(is.na(stage.rec.df.raw$h))){
    stage.rec.df = stage.rec.df.raw[-which(is.na(stage.rec.df.raw$h)),]
  } else {
    stage.rec.df = stage.rec.df.raw
  }
  
  
  # # test to study recession events:
  # #plot(stage.rec.df$t, stage.rec.df$h)
  period2study = c(0 , 25500)
  indexess = seq(which.min(abs(stage.rec.df$t - period2study[1])), which.min(abs(stage.rec.df$t - period2study[2])), 1)
  test.rec = stage.rec.df[indexess,]
  #gtest=ggplot() + geom_point(aes(x=test.rec$t, y=test.rec$h)) + geom_line(aes(x=test.rec$t, y=test.rec$h))

  
  
  
  # ALL CURVES EXTRACTION START:
  #********************************************************************************************************************
  message('- Total number of stage data:')
  print(length(stage.rec.df$t))
  
  # find all local decreasing points:
  message('- Find all stage data decreasing: h(t) < h(t-1)')
  min_loc.h   = localmin( t = stage.rec.df$t, 
                          h = stage.rec.df$h,
                          toll =toll)
  indexessx   = seq(which.min(abs(min_loc.h$t - period2study[1])), 
                    which.min(abs(min_loc.h$t - period2study[2])), 1)
  #gtest       = gtest + geom_point(aes(x=min_loc.h$t[indexessx], y=min_loc.h$h[indexessx]), color ="red", size= 2)  

  
  # find all those points in recession phase:
  message('- Among all decreasing stage data find all recession data.')
  recess.h    = rec(tmin        = min_loc.h$tmin, 
                    hmin        = min_loc.h$hmin, 
                    chi         = chi,
                    delta.t.max = delta.t.max, 
                    delta.t.min = delta.t.min, 
                    toll        = toll)
  indexessxx  = seq(which.min(abs(recess.h$trec - period2study[1])), 
                    which.min(abs(recess.h$trec - period2study[2])), 1)
  #gtest       = gtest + geom_point(aes(x=recess.h$trec[indexessxx], y=recess.h$hrec[indexessxx]), color ="green", size= 2)  

  
  # separate and extract recessions:
  message('- Separate recessions using parameters "chi" and deltatmax".')
  curves.h    = extract_curve(trec         = recess.h$trec,  
                              hrec         = recess.h$hrec,
                              chi          = chi, 
                              delta.t.max  = delta.t.max,
                              toll         = toll)

  indexessxxx = seq(which.min(abs(curves.h$t.real - period2study[1])),
                    which.min(abs(curves.h$t.real - period2study[2])),1)
  #gtest       = gtest+geom_point(aes(x=curves.h$t.real[indexessxxx], y=curves.h$hcurve[indexessxxx]), color = "blue", size= 2)
  
  
  # GRID:
  #######
  hmin      = min(curves.h$hcurve); 
  hmax      = max(curves.h$hcurve); 
  hmin_grid = hmax_grid = 0
  hgrid     = seq(-10000,100000,100)
  tgridd    = seq(0,10000,10)
  tmax      = max(curves.h$tcurve); 
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
  
  # all recession stage values:
  Data_h = data.frame(  RCi      = curves.h$RCi, 
                        hcurve   = curves.h$hcurve, 
                        tcurve   = curves.h$tcurve, 
                        treal    = curves.h$t.real, 
                        tpeak    = curves.h$tpeak,
                        hpeak    = curves.h$hpeak, 
                        tfinal   = curves.h$tfinal,
                        hfinal   = curves.h$hfinal)
  message('- Number of stage-recession data:')
  print(length(curves.h$hcurve))
    
    
    
  # detect all time final of the recession curves:
    #***********************************************
    tf =  hf =  index.h =  icurve.h =  hpeak.h =  tpeak.h =0
    for (ii in 1:length(curves.h$hcurve)) {
      if (curves.h$hfinal[ii] != 0) {
           icurve.h          = icurve.h + 1 
           tf[icurve.h]      = curves.h$tfinal[ii]   # Final time t of a recession
           hf[icurve.h]      = curves.h$hfinal[ii]   # Final stage h of a recession
           hpeak.h[icurve.h] = curves.h$hpeak[ii]    # peak of a stage-recession
           tpeak.h[icurve.h] = curves.h$tpeak[ii]    # time of the peak of a recession
           index.h[icurve.h] = ii                    # recession curve index
      } else {
           icurve.h = icurve.h
      }
    }
    
    
    message('- Number of potential stage-recessions:')
    print(icurve.h)
    
    # Saving results of recession curves extraction (in stage h) :
    message('- Saving intermediate results of recession curves extraction in file "Extract_rec_curves.csv".')
    write.table(Data_h, file=paste0(dir.extraction,"/Extract_rec_curves.csv"), sep=";")
    
    # Select some specific recessions based on user defined options:
    #***************************************************************
    #initialisation:
    asymptote.h  = curve_good.h = hpeakgood.h =  t.real.good.h = index.good.h = 0
    iiii = Nrec_longer_100 = results.regression = 0
    asym.df.temp = d.h =  d.h.selected=  d.h.selected.with.true.time = Resultss = quantiles.rec = NULL
    asymptote.df = data.frame(NULL)
    gggg         = ggplot()
    ncurves.h    = icurve.h
    message('\n\n\n**************************')
    message('Extracted stage-recessions')
    message('**************************')
    message('Loop on all recessions. Select only recessions with:')
    message('- recession length >= tgood')
    message('- number of recession points >= Nmin')
    message('- gradient <= gradient.max')
    
    
    
    #####################################################################################################
    for (i in 2:ncurves.h) {     # LOOP on all available curves:     
    #####################################################################################################
        d.h[[i]] = data.frame( t  = curves.h$tcurve[(index.h[i-1]+1):(index.h[i]-1)],
                               h  = curves.h$hcurve[(index.h[i-1]+1):(index.h[i]-1)],
                               uh = uh.rec) #abs(0/100*curves.h$hcurve[(index.h[i-1]+1):(index.h[i]-1)]))
        #*******************************************
        # NOW, REMOVE VALUES WITH TOO HIGH GRADIENT:
        #*******************************************
        # initialise:
        eliminate.index = c()  # first value always removed !!! 
        coeff_ang       = 0
        gradient.max    = gradient.max
        moving_average  = 4
        stop            = 1
        
        if (length(d.h[[i]]$t) > moving_average){ # if at least "moving_average" points !!
            # compute the angular coefficient:
            for (ccc in ((moving_average+1):(length( d.h[[i]]$t )))) {
               #linn = lm(formula = d.h[[i]]$h[(ccc-moving_average):ccc] ~  d.h[[i]]$t[(ccc-moving_average):ccc])
               linn = lm(formula = d.h[[i]]$h[1:ccc] ~  d.h[[i]]$t[1:ccc])
               coeff_ang[ccc]  = linn$coefficients[[2]]
               # coeff_ang[ccc] = (mean(d.h[[i]]$h[(ccc - moving_average):ccc]) - mean(d.h[[i]]$h[ccc])) /
               #                  (d.h[[i]]$t[ccc - moving_average]  - d.h[[i]]$t[ccc])
               if (coeff_ang[ccc]  <=  gradient.max) {
                  stop = ccc
               }
            }
            eliminate.index = c(seq(1:stop))
        }
        
        
        # Remove some data:
        ###################
        if (length(eliminate.index) > 1) {
          # deltan = nobs.h - length(eliminate.index)
          # if (deltan < Nmin.rec){
          #   eliminate.index = eliminate.index[ (Nmin.rec- deltan +1) : length(eliminate.index)]
          # }
          init_time_rec =  d.h[[i]]$t[1]
          d.h[[i]]      =  d.h[[i]][- eliminate.index,]
          # now shift the recession times earlier:
          d.h[[i]]$t    =  d.h[[i]]$t - d.h[[i]]$t[1] # + init_time_rec
        }

        
        # Select only recessions with length > tgood and with a number of point > Nmin :
        #********************************************************************************************************************
        if (length(d.h[[i]]$t) > Nmin.rec) {
        #********************************************************************************************************************
           if (tail(d.h[[i]]$t,1) >= tgood) {  
           #***************************************************************************************************************** 
                curve_good.h                  = curve_good.h + 1
                hpeakgood.h[curve_good.h]     = hpeak.h[i]
                t.real.good.h[curve_good.h]   = tpeak.h[i]
                index.good.h[curve_good.h]    = index.h[i]
                iiii[curve_good.h]            = i
                Nburn.rec                     = max(tail(which(curves.h$tcurve[(index.h[i-1]+1):(index.h[i]-1)] < tburn.rec), 1), 0)
                d.h.selected[[curve_good.h]]  = data.frame(t  = curves.h$tcurve[(index.h[i-1] +  1 + Nburn.rec):(index.h[i]-1)],
                                                           h  = curves.h$hcurve[(index.h[i-1]+ 1 + Nburn.rec):(index.h[i]-1)],
                                                           uh = uh.rec) 
                d.h.selected.with.true.time[[curve_good.h]]  = data.frame(t       = curves.h$tcurve[(index.h[i-1] +1 + Nburn.rec):(index.h[i]-1)],
                                                                          h       = curves.h$hcurve[(index.h[i-1] +1 + Nburn.rec):(index.h[i]-1)],
                                                                          uh      = uh.rec,
                                                                          tpeak   = t.real.good.h[curve_good.h],
                                                                          treal   = curves.h$tcurve[(index.h[i-1] +1 + Nburn.rec):(index.h[i]-1)] + t.real.good.h[curve_good.h],
                                                                          ind.rec = rep(curve_good.h, length(curves.h$tcurve[(index.h[i-1] + 1+Nburn.rec):(index.h[i]-1)])))
                nobs.h       = length(d.h.selected[[curve_good.h]]$t)
                print(paste0("Recession # ", curve_good.h, ' (', i, ') with ', nobs.h, ' points' ))
           
                # } else if ((rec.model == "1expWithAsymptNorm") |(rec.model == "2expWithAsymptNorm")) {
                #   d.h.selected[[curve_good.h]]  = data.frame(curves.h$tcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)],
                #                                              curves.h$Qcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)]/curves.h$Qcurve[(index.h[i-1]+Nburn.rec)],
                #                                              uh.rec/100*curves.h$Qcurve[(index.h[i-1]+Nburn.rec):(index.h[i]-1)])
                # }
                # eliminate.index =0
                # for (ccc in 2:length( d.h.selected[[curve_good.h]]$t)) {
                #   if ((d.h.selected[[curve_good.h]]$t[ccc] - d.h.selected[[curve_good.h]]$t[ccc-1]) < delta.t.min) {
                #        eliminate.index = c(eliminate.index, ccc)
                #   }
                # }
                # if (eliminate.index > 0) {
                #   d.h.selected[[curve_good.h]]                 = d.h.selected[[curve_good.h]][-eliminate.index,]
                #   d.h.selected.with.true.time[[curve_good.h]]  = d.h.selected.with.true.time[[curve_good.h]][-eliminate.index,]
                # }
           
                # SAVE VALUES:
                #*************
                curve_data.h = paste0(dir.BaM.rec.pool,"/Curves_Data.txt")
                write.table(d.h.selected[[curve_good.h]] , 
                            file = curve_data.h, append = FALSE, sep = "\t", eol = "\n", 
                            na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh"))
           
                # find out the longest recessions:
                if (tail(d.h.selected.with.true.time[[curve_good.h]]$t, 1) >= 80) {       ####### Change this eventually !!!!!!
                   Nrec_longer_100 = Nrec_longer_100 + 1
                } 
           
                # if ((tail(d.h[[i]]$t,1) >= 170)){  # JUST AS A TEST!!
                #   #cat(paste0("curve_good.h = ", curve_good.h, '\n', "tin          = ", t.real.good.h[curve_good.h]))
                #   #print(d.h.selected.with.true.time[[curve_good.h]])
                #   #gggg = gggg + geom_point(data = d.h.selected.with.true.time[[curve_good.h]], aes(x= t, y =h))
                # }
           }
        }
    }
    ################################################
    # write recessions to file:
    df.curves   = bind_rows( d.h.selected, .id = "column_label")
    df.curves.t = bind_rows( d.h.selected.with.true.time, .id = "column_label")
    if (is.null(d.h.selected)) {
       err1 = print("******* ERROR: extracted recessions dataframe is NULL ! please change settings.")
       return(err1)
    } else {
       write.table(df.curves,file=paste0(dir.BaM.rec.pool, "/Curves_Data.txt"), sep = "\t", eol = "\n",
                na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh", "Period"))

    
    
    ############################
    # plot extracted recessions:
    ############################
    message("************************************")
    message("Results of Extracted recessions")
    message("***********************************")    
    Nrec          = length(df.curves.t$t)
    message("- Number of extracted stage-recession data:")
    print(Nrec)
    Nk            = tail(df.curves.t$ind.rec,1)
    message("- Number of extracted stage-recession periods:")
    print(Nk)
    message("- Plotting extracted results to:")
    path.plot.extracted.rec = paste0(dir.extraction,"/Figure_recession_selection_chi",chi,".pdf")
    print(path.plot.extracted.rec)
    rec.plot.test = NULL
    rec.plot.test = plot.extracted.recessions.paper(Data_rec          = df.curves.t,
                                                    stage.record      = stage.record,
                                                    hmin_grid         = hmin_grid, 
                                                    hmax_grid         = hmax_grid, 
                                                    dir.extraction    = dir.extraction, 
                                                    chi               = chi,
                                                    tmax.rec          = tmax.rec,
                                                    hmax.rec          = hmax.rec,
                                                    nobs              = Nrec, 
                                                    Nperiods          = Nk, 
                                                    Nrec_longer_100   = Nrec_longer_100,
                                                    stage.scale.shift = stage.scale.shift)
    
    pdf(path.plot.extracted.rec, 14, 6 ,useDingbats=F)
    print(rec.plot.test$plot4paper.extraction)
    dev.off()
    Ncurves = length(d.h.selected)
    
    
    min_rec = data.frame(time_stage_min = double(),  stage_min = double())
    for (rec in 1:length(d.h.selected.with.true.time)) {
      min_rec = rbind(min_rec, data.frame(time_stage_min = tail(d.h.selected.with.true.time[[rec]]$treal, 1),
                                          stage_min =      tail(d.h.selected.with.true.time[[rec]]$h, 1)))
    }
    # Plot extracted curves against recession time t*:
    #*************************************************
    rec.minimum.plot <- ggplot(data =min_rec, aes(x = time_stage_min, y = stage_min - stage.scale.shift)) + 
                        #geom_line(data = stage.record, aes(t_limni, h_limni*100), size =0.1, color = "Gray80") +
                        geom_point(size = 2.5) +
                        coord_cartesian(clip = "off")+
                        scale_x_continuous(name = "Time [days]",  expand = c(0,0)) +
                        scale_y_continuous(name = "Minimum Recession Stage H [cm]", expand = c(0,0)) +
                        theme_bw(base_size = 15) +
                        theme( axis.text             = element_text(size=15)
                              ,axis.title           = element_text(size=15)
                              ,plot.margin          = unit(c(0.5,0.5, 0.5, 2),"cm"))   
    path.plot.minimum.rec = paste0(dir.extraction,"/Figure_recession_minimum_stages.pdf")
    pdf(path.plot.minimum.rec, 10, 5 ,useDingbats=F)
    print(rec.minimum.plot)
    dev.off()
    
    print("****************")
    print("   All done!    ")
    print("****************")
    return(list(d.h.selected  = d.h.selected, 
              t.real.good.h = t.real.good.h,
              curve.h       = curves.h, 
              data.curve    = Data_h, 
              index.curve   = icurve.h, 
              hpeak         = hpeak.h, 
              tpeak         = tpeak.h, 
              index         = index.h, 
              tfinal        = tf, 
              hfinal        = hf,
              hmin          = hmin_grid, 
              hmax          = hmax_grid,
              tmax          = tmax_grid,
              Ncurves       = Ncurves))
    }
}







































#########################################################################################################
recession.regression   <-  function(dir.exe,
                                    dir.segment.rec,
                                    file.options.general,
                                    file.options.recess,
                                    data.recess,
                                    initial.time.rec,
                                    which.recession){
#########################################################################################################
  # read inputs and options for computation:
  source(file.options.general)
  source(file.options.recess)
  
  
  
  ####################
  # Hardcoded options
  BayesianOption         = 2                             # Type of Bayesian model (1 = simple, 2 = pooling). only 2 is available!
  limits.Y.alpha         = c(0, 1000, 200)               # Limits for the alpha (initial stage) plot.
  limits.Y.lambda        = c(0, 40, 10)                  # Limits for the lambda (recession rate) plot.
  plot.b.from.gaugings   = FALSE                         # plot the river bed estimation obtained by using gaugings
  x.name                 = "Time [days]"                 # x-axis label for plot of time series
  y.name                 = "Asymptotic stage [m]"        # y-axis label for plot of time series of asymptotic stage
  save.all.results       = TRUE                          # TRUE = save all segmentation computations
  plot.gamma.uncertainty = TRUE                          # plot structural uncertainty on the segmentation.
  segm_all_parameters    = FALSE
  stage.scale.shift      = 1000                          # is in [cm]: parameter used to shift stage values and avoid negative values: h = h + stage.scale.shift
  stage.limits           = c(-150, 50, 50)               # limits for recession stage is in "cm" !!!
  limits.y               = stage.limits
  limits.x.recess        = c(0, 150, 30)                 # [c(min, max, step)] limits for the recession period in days. by default = c(0, 150, 30).
  limits.x               = limits.x.recess
  plot.recession.uncert  = TRUE
  limni.time.limits      = NULL
  asymptote.limits       = stage.limits
  #####################
  
  

  
  # Directories :
  ###############
  dir.segment.rec.test1   = paste0(dir.segment.rec,"/",name.folder.results.recession)
  dir.BaM.recession.pool  = paste0(dir.exe,"/Recession_h_pooling") 
  dir.estim.recessions    = paste0(dir.segment.rec.test1,"/2_curves_estimation")
  dir.segm.recessions     = paste0(dir.segment.rec.test1,"/3_curves_segmentation")
  
  dir.create( paste0(dir.segment.rec,"/", name.folder.results.recession))
  dir.create( paste0(dir.segment.rec.test1,"/2_curves_estimation"))
  dir.create( paste0(dir.segment.rec.test1,"/3_curves_segmentation"))
  dir.create( paste0(dir.estim.recessions,"/Pooling"))

  message("Recessions exponential regression - using BaM !!!"); flush.console()
  Ncurves       = length(data.recess[[1]])
  colfunc       = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))   
  output_file_h = paste0(dir.estim.recessions,"/Param_rec_h.csv")
  dir.rec.pool  = paste0(dir.estim.recessions,"/Pooling")
  Ncurves.pool  = length(which.recession)
  d.h           = data.frame(cbind(data.recess[[1]][[which.recession[1]]], 
                                   Period = rep(which.recession[1], length(data.recess[[1]][[which.recession[1]]]$t))))
  for (icurve in which.recession[-1]) { 
    if (is.null(data.recess[[1]][[icurve]])){
         print("error in input dataset: please check the period column")
    }
    d.h = rbind(d.h,  cbind(data.recess[[1]][[icurve]],  Period = rep(icurve, length(data.recess[[1]][[icurve]]$t))))
  }   
  d.h$Period    = d.h$Period - which.recession[1] + 1 
  nobs.pooling  = length(d.h$t)
  message(paste0("Total number of recessions to estimate = ", Ncurves.pool, "... Wait ... ")); flush.console();
  message(paste0("Total number of mcmc = ", nmcmc.rec*Ncycles.mcmc.rec)); flush.console()
  message(paste0("Total number of stage-recession data = ", nobs.pooling)); flush.console()
  
  
  
  
  
  
  # Start Bayesian inference of stage-recession model 
  # for each selected model:
  # Be sure that the recession options file is entirely filled, 
  # in particular the options for the recession estimation:
  # priors for recession model parameters
  read.reg.rec = list()
  ###################################
  for (irec in 1:length(rec.model)) {
  ###################################
            rec.mod    = rec.model[[irec]]
            prior.rec  = prior.param.rec[[irec]]
            #shift the pprior of asymptotic parameter:
            prior.rec[length(prior.rec)-4] =  as.numeric(prior.rec[length(prior.rec)-4])  + stage.scale.shift
            prior.rec[length(prior.rec)-3] =  as.numeric(prior.rec[length(prior.rec)-3])  + stage.scale.shift  
            prior.rec[length(prior.rec)-1] =  as.numeric(prior.rec[length(prior.rec)-1])  + stage.scale.shift
            dir.rec.pool.model = paste0(dir.rec.pool,"/model_",rec.mod)   
            dir.rec.pool.test  = paste0(dir.rec.pool.model,"/chi_",chi)
            dir.create(paste0(dir.rec.pool,"/model_",rec.mod))
            dir.create(paste0(dir.rec.pool.model,"/chi_",chi))
            message("###########################################################################################"); flush.console()
            message(paste0("Stage-recession model =", rec.mod)); flush.console()
            message("###########################################################################################"); flush.console()            

            
            #######################################
            if (estim.plot.results.only == TRUE) {
            #######################################
               if ((file.exists(paste0(dir.rec.pool.test,"/Results_Summary.txt")) == FALSE)|
                  (file.exists(paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt")) == FALSE)|
                  (file.exists(paste0(dir.rec.pool.test,"/Results_Residuals.txt")) == FALSE)){
                    plot.results.new= as.logical(dlgInput(c("You have selected the option to read results only",
                                       "(estim.plot.results.only =  TRUE) in the Options_Recession file", 
                                       "But there are no results to read and process.", " ",
                                       "--> If you want to perform also the recession estimation press 'F' or 'FALSE'",
                                       "--> Otherwise press T, check settings file and run again!"), 
                                     Sys.info()[" "])$res)
                    if (plot.results.new == FALSE){
                         estim.plot.results.only = FALSE
                    } else {
                         err1 = print("******* ERROR: here are no results to process, please put estim.plot.results.only = FALSE")
                         return(err1)
                    }
               }
            }
  
            
            #######################################
            if (estim.plot.results.only == FALSE) {
            #######################################
            #initialisation:
            asymptote.h    = 0; curve_good.h= 0;  hpeakgood.h= 0; t.real.good.h =0; index.good.h = 0; 
            asym.h.maxpost = 0; asym.h.stdev= 0;  asym.h.Q10 = 0; asym.h.Q90 = 0; asym.h.mean = 0;
            asym.h.Q2      = 0; asym.h.Q95  = 0;  results.regression = 0;
            theta1.maxpost = 0; theta1.stdev= 0;  theta1.Q10 = 0; theta1.Q90 = 0; theta1.mean = 0;
            theta2.maxpost = 0; theta2.stdev= 0;  theta2.Q10 = 0; theta2.Q90 = 0; theta2.mean = 0;
            theta3.maxpost = 0; theta3.stdev= 0;  theta3.Q10 = 0; theta3.Q90 = 0; theta3.mean = 0;
            theta4.maxpost = 0; theta4.stdev= 0;  theta4.Q10 = 0; theta4.Q90 = 0; theta4.mean = 0;
            theta5.maxpost = 0; theta5.stdev= 0;  theta5.Q10 = 0; theta5.Q90 = 0; theta5.mean = 0;
            theta6.maxpost = 0; theta6.stdev= 0;  theta6.Q10 = 0; theta6.Q90 = 0; theta6.mean = 0;
            asymptote.df   = data.frame(NULL)
            asym.df.temp   = NULL  
            Resultss       = NULL
            quantiles.rec  = NULL
            setwd(dir.exe)
            start_time = Sys.time()
            curve_data = "Recession_h_pooling/Curves_Data.txt"
            write.table(d.h, file = curve_data, append = FALSE, sep = "\t", eol = "\n",
                        na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh", "Period"))
            if (prior.gamma.rec[7] == "'Uniform'"){
                prior.gamma.rec[5] =  as.numeric(prior.gamma.rec[5])/stage.scale.shift
                prior.gamma.rec[6] =  as.numeric(prior.gamma.rec[6])/stage.scale.shift
                prior.gamma.rec[8] =  as.numeric(prior.gamma.rec[8])/stage.scale.shift
            }
            # Launch BaM application in Bayesian pooling:
            message("***************************************************************"); flush.console()
            message("Applying BaM Regression by pooling method !!!  Please, wait ... "); flush.console()
            message("***************************************************************"); flush.console()
            # Configurate the files for BaM (Folder "/BaM_exe/Recession_h_pooling")
            BaM_config.pooling(dir.exe     = dir.exe,
                               model       = rec.mod, 
                               nobs        = nobs.pooling, 
                               ncycles     = Ncycles.mcmc.rec,
                               nmcmc       = nmcmc.rec, 
                               ncurves     = Ncurves.pool, 
                               jump.pos    = jump.pos.rec, 
                               jump.neg    = jump.neg.rec, 
                               slim        = nslim.rec, 
                               burn        = nburn.rec,
                               prior       = prior.rec,
                               gamma.model = gamma.model.rec,
                               prior.gamma = prior.gamma.rec)
            # Launch BaM.exe (Benjamin Renard)
            system2("BaM_recession_multi_model_final.exe")
            ##############################################
            # Read results from BaM :
            list.files.pool <- c( paste0(dir.BaM.recession.pool,"/Results_MCMC_Cooked.txt"),
                                  paste0(dir.BaM.recession.pool,"/Results_Residuals.txt"),
                                  paste0(dir.BaM.recession.pool,"/Results_Summary.txt"),
                                  paste0(dir.BaM.recession.pool,"/Config_Model.txt"),
                                  paste0(dir.BaM.recession.pool,"/Curves_Data.txt"))
            for (i in 1:length(list.files.pool)) {
                  file.copy(list.files.pool[i], dir.rec.pool.test ,overwrite = TRUE)
            }
            end_time <- Sys.time()
            print(c("computat. time for Bayesian regression of all recessions by pooling =", end_time - start_time))
            
            
            ########
            } else {
            ########
                 message("*********************"); flush.console()
                 message("READ RESULTS ONLY !!!"); flush.console()
                 message("*********************"); flush.console()
            }

            
            
            
            
            
            
            
            ########################
            # PLOT RESULTS:
            ########################  
            setwd(dir.rec.pool.test)
            Results_summary.pool    = read.table(paste0(dir.rec.pool.test,"/Results_Summary.txt"),     header =TRUE)
            Results_mcmc.pool       = read.table(paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"), header =TRUE)
            Results_residuals.pool  = read.table(paste0(dir.rec.pool.test,"/Results_Residuals.txt"),   header =TRUE)
            nparam                  = ncol(Results_mcmc.pool)
            vertical.length         = nparam/3
            # plot mcmc traceplots and density plots to check convergence and representiveness:
            message("- Plotting mcmc trace and density plots for visual check ..."); flush.console()
            mcmc = MCMCplot(doLogPost = T,
                            doPar     = T,
                            doDPar    = F, 
                            MCMCfile  = paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"), 
                            type      = "trace",  #="trace", # "histogram","density","scatterplot"
                            xlab      = '',
                            ylab      = '',
                            ncol      = 6, 
                            prior     = NULL,
                            burn      = 0, 
                            slim      = 1,
                            theme_size= 15)
            mcmc2 = MCMCplot(doLogPost = T,
                             doPar     = T,
                             doDPar    = F, 
                             MCMCfile  = paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"), 
                             type      = 'density', 
                             prior     = NULL,
                             xlab      = '',
                             ylab      = '',
                             ncol      = 6, 
                             burn      = 0, 
                             slim      = 1,
                             theme_size= 15)  
            pdf(paste0(dir.rec.pool.test,"/mcmc_recessions.pdf"), 40, vertical.length, useDingbats=F)
            print(plot_grid(mcmc, mcmc2,
                            nrow=1, ncol = 2, 
                            labels = c("Trace plots", "Density plots"),
                            label_size = 25))
            dev.off()
            # test convergence:
            ###################
            #conv = Convergence.test(dir.seg  = dir.rec.pool.test , 
            #                        npar     = nparam, 
            #                        dir.plot = dir.rec.pool.test)
            
            # Residuals:
            ############
            message("- Plotting stage-recessions residuals ..."); flush.console()
            leg.obs  = "Observed recessions data"
            leg.sim  = "Simulated recessions data"
            leg.obs.sim  = c("Observed recessions data" = "black", "Simulated recessions data" = "orange")
            residuals.plot = ggplot(Results_residuals.pool) +
                             geom_point(aes(x = X1_obs,  y = Y1_sim, fill = leg.sim), pch =21, size = 1.8) +
                             geom_point(aes(x = X1_obs,  y = Y1_obs, fill = leg.obs), pch =21, size = 1) +
                             xlab("Recession time [day]") + ylab("Stage [m]") + 
                             scale_fill_manual(name=element_blank(), values=leg.obs.sim, breaks=c(leg.obs,leg.sim),labels=c(leg.obs,leg.sim)) +
                             theme_bw() + theme(legend.position = "bottom")
                             ggsave(residuals.plot, filename =paste0(dir.rec.pool.test, "/Residuals.png"),
                                    bg = "transparent", width = 8, height =4, dpi = 200)
            
            
                             
            # save the asymptotic stage parameter:
            ###############################################################################   
            if (rec.mod == "1expWithAsympt"){   #a1(k), b1, a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 1 + which.recession]
                #write.table(a1.mcmc, file = "a1.mcmc.txt", append = FALSE, row.names = FALSE)
            } else if (rec.mod == "2expWithAsympt"){ # a1(k), b1, a2, b2, a3(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                a2.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                b2.mcmc   = Results_mcmc.pool[, Ncurves.pool + 3]
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 3 + which.recession]
            } else if (rec.mod == "2expWithAsympt_bis"){ # a1(k), b1, a2(k), b2, a3(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                a2.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1+ which.recession]
                b2.mcmc   = Results_mcmc.pool[, 2*Ncurves.pool + 2]
                asym.mcmc = Results_mcmc.pool[, 2*Ncurves.pool + 2 + which.recession]
            } else if (rec.mod == "2expWithAsympt_rel"){ # a1(k), b1, a2, b2, a3(k)
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 3 + which.recession]
            } else if (rec.mod == "3expWithAsympt"){   # a1(k), b1, a2, b2, a3, b3, a4(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                a2.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                b2.mcmc   = Results_mcmc.pool[, Ncurves.pool + 3]
                a3.mcmc   = Results_mcmc.pool[, Ncurves.pool + 4]
                b3.mcmc   = Results_mcmc.pool[, Ncurves.pool + 5]
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 5 + which.recession]
            } else if (rec.mod == "3expWithAsympt_bis"){   # a1(k), b1, a2, b2, a3, b3, a4(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                a2.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1+ which.recession]
                b2.mcmc   = Results_mcmc.pool[, 2*Ncurves.pool + 2]
                a3.mcmc   = Results_mcmc.pool[, 2*Ncurves.pool + 3]
                b3.mcmc   = Results_mcmc.pool[, 2*Ncurves.pool + 4]
                asym.mcmc = Results_mcmc.pool[, 2*Ncurves.pool + 4 + which.recession]
            } else if (rec.mod == "expexp"){    # a1(k), b1, n1, a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                n1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 2 + which.recession]
            } else if (rec.mod == "expexp_bis"){    # a1(k), b1(k), n1, a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                n1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                asym.mcmc = Results_mcmc.pool[, 2*Ncurves.pool + 1 + which.recession]
            } else if (rec.mod == "hyperb"){  # a1(k), b1, n1, a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                n1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 2 + which.recession]
            } else if (rec.mod == "hyperb_bis"){  # a1(k), b1, n1, a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                n1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                asym.mcmc = Results_mcmc.pool[, 2*Ncurves.pool + 1 + which.recession]
            } else if (rec.mod == "Boussinesq"){  # a1(k), b1,  a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 1 + which.recession]
            } else if (rec.mod == "Coutagne"){   # a1(k), b1, n1, a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                n1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                asym.mcmc = Results_mcmc.pool[, Ncurves.pool + 2 + which.recession]
            } else if (rec.mod == "Coutagne_bis"){   # a1(k), b1, n1, a2(k)
                a1.mcmc   = Results_mcmc.pool[, which.recession]
                b1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 1]
                n1.mcmc   = Results_mcmc.pool[, Ncurves.pool + 2]
                asym.mcmc = Results_mcmc.pool[, 2*Ncurves.pool + 1 + which.recession]
            }
                             
            # Asymptote:
            asym.df.pool        = data.frame(t(quantile(x = asym.mcmc[,1], p = c(0.025, 0.5, 0.975)))); 
            names(asym.df.pool) = c("q2", "mean", "q97");
            asym.df.pool        = cbind(asym.df.pool , stdev= std(asym.mcmc[,1]), t = data.recess[[2]][1])
            for (aa in 2:Ncurves.pool) {
                 asym.df.temp        = data.frame(t(quantile(x = asym.mcmc[,aa], p = c(0.025, 0.5, 0.975))))
                 names(asym.df.temp) = c("q2", "mean", "q97")
                 asym.df.pool        = rbind(asym.df.pool, cbind(  asym.df.temp,  stdev= std(asym.mcmc[,aa]), t = data.recess[[2]][aa]))
            }
            write.table(asym.df.pool, 
                        file = paste0(dir.rec.pool.test, "/df.asymptote.txt"), append = FALSE, sep = "\t", eol = "\n", 
                        na = "NA", dec = ".", row.names = FALSE, col.names=c("q2", "mean", "q97", "stdev", "time"))  
            
            # plot the time series of the asymptotic stage parameter:
            error.max = 0
            error.min = 0
            for (iii in 1:length(asym.df.pool$t)) {
                 error.min[iii] = asym.df.pool$q2[iii]  - stage.scale.shift   # )),  limits.y[1])
                 error.max[iii] = asym.df.pool$q97[iii] - stage.scale.shift   # )),  limits.y[2])
            }
            asym.df.pool$error.max = error.max
            asym.df.pool$error.min = error.min
            message("- Plotting asymptotic stage time series ..."); flush.console()
            asympt.plot =  ggplot()+
                           geom_point(data = asym.df.pool, aes(x = t, y = mean - stage.scale.shift), colour = "black", size = 3) +
                           geom_line(data  = asym.df.pool, aes(x = t, y = mean - stage.scale.shift), colour = "gray90", size = 0.2) +
                           geom_errorbar(data = asym.df.pool, aes(x = t, ymin = error.min, ymax = error.max), width=50, size = 0.3)+
                           scale_y_continuous( expand = c(0,0)) + #, limits= c(limits.y[1], limits.y[2])) +
                           theme_bw(base_size = 15)+
                           xlab("Time [day]") + ylab("Asymptotic level [cm]") +
                           theme( panel.grid.major     = element_blank() 
                                 ,panel.grid.minor     = element_blank())
            
                           pdf(paste0(dir.rec.pool.test,"/Figure_asymptote_time_series.pdf"), 12, 6 ,useDingbats=F)
                           print(asympt.plot)
                           dev.off()
                           
           # saving results of the Bayesian recession estimation:
           #****************************************************************************************************
           read.reg.rec[[irec]]  =  read.results.regression.rec(    dir.recess          = dir.rec.pool.test,
                                                                    dir.segm.recessions = dir.segm.recessions,
                                                                    dir.BaM.rec.pool    = dir.BaM.recession.pool,
                                                                    rec.model           = rec.mod,
                                                                    BayesianOption      = BayesianOption,
                                                                    which.recession     = which.recession,
                                                                    time.rec            = initial.time.rec,
                                                                    is.b1.var           = FALSE,
                                                                    time.limits         = limni.time.limits,
                                                                    asymptote.limits    = asymptote.limits,
                                                                    stage.scale.shift   = stage.scale.shift,
                                                                    chi                 = chi)
           #****************************************************************************************************
           plot.all.recessions(dir.rec.pool.test     = dir.rec.pool.test,
                               rec.model             = rec.mod,
                               rec.model.gamma       = gamma.model.rec,
                               stage.limits          = stage.limits,
                               limits.x.recess       = limits.x.recess,
                               stage.scale.shift     = stage.scale.shift,
                               plot.recession.uncert = plot.recession.uncert)
           
           message("- Recession estimation correctly performed, and saved on folder:"); flush.console()
           message(paste0("  ", dir.rec.pool.test)); flush.console()
  }
  #*************************************************************************************************************
  print("****************")
  print("   All done!    ")
  print("****************")
  return(read.reg.rec)
}





























# individual regression model
######################################################################################
BaM_config1.h <- function(tgrid, nobs, 
                          ncycles, nmcmc, nslim, jump.pos,  jump.neg, nburn,
                          model, pred, asympt.before, remnant,
                          prior, prior.gamma) {                            
######################################################################################
  write.table(tgrid, file =paste0(dir_code,"/BaM_exe/Recession_h/tgrid.txt"), 
              col.names = FALSE, row.names = FALSE)
  ngrid = length(tgrid)
  #-----------------------------------------------------------------
  #creation of Config_BaM.txt
  file.bam = paste0(dir_code,"/BaM_exe/Config_BaM.txt")
  cat('"Recession_h/"', file =file.bam , sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"', file = file.bam , sep="\n", append = TRUE)    
  cat('"Config_Model.txt"', file = file.bam , sep="\n", append = TRUE)
  cat('"Config_xtra.txt"', file = file.bam , sep="\n", append = TRUE)
  cat('"Config_Data.txt"', file = file.bam, sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file.bam , sep="\n", append = TRUE)                                       
  cat('"Config_MCMC.txt"', file = file.bam , sep="\n", append = TRUE)                                            
  cat('"Config_Cooking.txt"', file = file.bam , sep="\n", append = TRUE)
  cat('"Config_Summary.txt"', file = file.bam , sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"', file = file.bam, sep="\n", append = TRUE)
  if(pred==TRUE) {
    cat('"Config_Pred_Master.txt"', file = file.bam , sep="\n", append = TRUE)
  } else if(pred==FALSE) {
    cat('""', file = file.bam , sep="\n", append = TRUE)
  }
  #----------------------------------------------------------------------
  file.name = paste(dir_code,"/BaM_exe/Recession_h/Config_Model.txt",sep="")
  cat('"Recession_h"', file =file.name,sep="\n")
  cat(1, file = file.name, append = TRUE,sep="\n")
  cat(1, file =file.name, append = TRUE,sep="\n")
  
  
  if (model == "1expWithAsympt") { 
  ###############################################################################################
    cat(3, file =file.name, append = TRUE,sep="\n")       #Number of paramaters  5
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    cat(prior[3],  file =file.name, append = TRUE,sep="\n")
    cat(prior[1],  file =file.name, append = TRUE,sep=",")
    cat(",",       file =file.name, append = TRUE, sep=",")
    cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    #-------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    cat(prior[8], file =file.name, append = TRUE,sep="\n")
    cat(prior[6], file =file.name, append = TRUE, sep=",")
    cat(",",      file =file.name, append = TRUE, sep=",")
    cat(prior[7], file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"a2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    cat(prior[13],   file =file.name, append = TRUE,sep="\n")
    cat(prior[11],   file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    # file extra:
    file.xtra = paste(dir_code,"/BaM_exe/Recession_h/Config_xtra.txt",sep="")
    cat("exp", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    
    
  } else if ((model == "2expWithAsympt")|(model == "2expWithAsympt_bis")) { # =========================================> 2exp
    ############################################################################
    cat(5, file =file.name, append = TRUE,sep="\n")     #Number of paramaters  5
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    cat(prior[3],  file =file.name, append = TRUE,sep="\n")
    cat(prior[1],  file =file.name, append = TRUE,sep=",")
    cat(",",       file =file.name, append = TRUE, sep=",")
    cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    #------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    cat(prior[8], file =file.name, append = TRUE,sep="\n")
    cat(prior[6], file =file.name, append = TRUE, sep=",")
    cat(",",      file =file.name, append = TRUE, sep=",")
    cat(prior[7], file =file.name, append = TRUE, sep="\n")
    #------------------------------------------------------
    cat('"a2"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    cat(prior[13],   file =file.name, append = TRUE,sep="\n")
    cat(prior[11],    file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"b2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    cat(prior[18],   file =file.name, append = TRUE,sep="\n")
    cat(prior[16] ,  file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"a3"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[24],   file =file.name, append = TRUE,sep="\n")
    cat(prior[23],   file =file.name, append = TRUE,sep="\n")
    cat(prior[21],   file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[22],   file =file.name, append = TRUE, sep="\n")
    # file extra:
    file.xtra = paste0(dir_code,"/BaM_exe/Recession_h/Config_xtra.txt")
    cat("exp", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(2, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    

  } else if ((model == "2expWithAsympt_rel")) { # =============================================> 2exp terms (whose one is relative)
    ###########################################################
    cat(5, file =file.name, append = TRUE,sep="\n")       #Number of paramaters  5
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    cat(prior[3],  file =file.name, append = TRUE,sep="\n")
    cat(prior[1],  file =file.name, append = TRUE,sep=",")
    cat(",",       file =file.name, append = TRUE, sep=",")
    cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    #-------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    cat(prior[8], file =file.name, append = TRUE,sep="\n")
    cat(prior[6], file =file.name, append = TRUE, sep=",")
    cat(",",      file =file.name, append = TRUE, sep=",")
    cat(prior[7], file =file.name, append = TRUE, sep="\n")
    #------------------------------------------------------
    cat('"a2"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    cat(prior[13],   file =file.name, append = TRUE,sep="\n")
    cat(prior[11],    file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"b2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    cat(prior[18],   file =file.name, append = TRUE,sep="\n")
    cat(prior[16] ,  file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"a3"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[24],   file =file.name, append = TRUE,sep="\n")
    cat(prior[23],   file =file.name, append = TRUE,sep="\n")
    cat(prior[21],   file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[22],   file =file.name, append = TRUE, sep="\n")
    # file extra:
    file.xtra = paste(dir_code,"/BaM_exe/Recession_h/Config_xtra.txt",sep="")
    cat("exp_rel", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    
    #######################################################################
  } else if ((model == "3expWithAsympt")|(model == "3expWithAsympt_bis")) { #============ ==> 3exp
    ######################################################################
    cat(7, file =file.name, append = TRUE,sep="\n")       #Number of paramaters  7
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    cat(prior[3],  file =file.name, append = TRUE,sep="\n")
    cat(prior[1],  file =file.name, append = TRUE,sep=",")
    cat(",",       file =file.name, append = TRUE, sep=",")
    cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    #------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    cat(prior[8], file =file.name, append = TRUE,sep="\n")
    cat(prior[6], file =file.name, append = TRUE, sep=",")
    cat(",",      file =file.name, append = TRUE, sep=",")
    cat(prior[7], file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"a2"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    cat(prior[13],   file =file.name, append = TRUE,sep="\n")
    cat(prior[11],    file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"b2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    cat(prior[18],   file =file.name, append = TRUE,sep="\n")
    cat(prior[16] ,  file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"a3"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[24],   file =file.name, append = TRUE,sep="\n")
    cat(prior[23],   file =file.name, append = TRUE,sep="\n")
    cat(prior[21],   file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[22],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"b3"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[29],   file =file.name, append = TRUE,sep="\n")
    cat(prior[28],   file =file.name, append = TRUE,sep="\n")
    cat(prior[26],   file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[27],   file =file.name, append = TRUE, sep="\n")
    #--------------------------------------------------------
    cat('"a4"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[34],   file =file.name, append = TRUE,sep="\n")
    cat(prior[33],   file =file.name, append = TRUE,sep="\n")
    cat(prior[31],   file =file.name, append = TRUE, sep=",")
    cat(",",         file =file.name, append = TRUE, sep=",")
    cat(prior[32],   file =file.name, append = TRUE, sep="\n")
    # file extra:
    file.xtra = paste0(dir_code,"/BaM_exe/Recession_h/Config_xtra.txt")
    cat("exp", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(3, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    
    
    
    
  } else if ((model == "expexp")|(model == "expexp_bis")) { # ======================================================> 1 double exp
  ###########################################################
    cat(4, file =file.name, append = TRUE,sep="\n")       #Number of paramaters  5
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[3],  file =file.name, append = TRUE,sep="\n")
      cat(prior[1],  file =file.name, append = TRUE,sep=",")
      cat(",",       file =file.name, append = TRUE, sep=",")
      cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE,sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"n1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_n1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE,sep="\n")
      cat(prior[11],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"a2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[18],   file =file.name, append = TRUE,sep="\n")
      cat(prior[16],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    }
    # file extra:
    file.xtra = paste(dir_code,"/BaM_exe/Recession_h/Config_xtra.txt",sep="")
    cat("expexp", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    
    
  } else if ((model == "hyperb")|(model == "hyperb_bis")) { # ==============================================> 1 double exp
  ###########################################################
    cat(4, file =file.name, append = TRUE,sep="\n")    #Number of paramaters  5
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[3],  file =file.name, append = TRUE,sep="\n")
      cat(prior[1],  file =file.name, append = TRUE,sep=",")
      cat(",",       file =file.name, append = TRUE, sep=",")
      cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    }
    #------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE,sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"n1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_n1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE,sep="\n")
      cat(prior[11],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #---------------------------------------------------------
    cat('"a2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[18],   file =file.name, append = TRUE,sep="\n")
      cat(prior[16],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    }
    # file extra:
    file.xtra = paste(dir_code,"/BaM_exe/Recession_h/Config_xtra.txt",sep="")
    cat("hyperb", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    
    
    
  } else if ((model == "Coutagne")|(model == "Coutagne_bis")) { #======================================================> Coutagne
    cat(4, file =file.name, append = TRUE,sep="\n")     #Number of paramaters  5
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[3],  file =file.name, append = TRUE,sep="\n")
      cat(prior[1],  file =file.name, append = TRUE,sep=",")
      cat(",",       file =file.name, append = TRUE, sep=",")
      cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE,sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"n1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_n1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE,sep="\n")
      cat(prior[11],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #---------------------------------------------------------
    cat('"a2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[18],   file =file.name, append = TRUE,sep="\n")
      cat(prior[16],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    }
    # file extra:
    file.xtra = paste0(dir_code,"/BaM_exe/Recession_h/Config_xtra.txt")
    cat("Coutagne", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
  }
  
  
  #----------------------------------------------------------------------  
  file.name2 = paste(dir_code,"/BaM_exe/Recession_h/Config_Data.txt",sep="")
  name = "Recession_h\\Curves_Data.txt"
  cat(name, file =file.name2,  sep="\n") 	               	# path to data file
  cat(1, file =file.name2,sep="\n",append = TRUE)      	 	# number of header lines
  cat(nobs, file =file.name2,sep="\n",append = TRUE)       # Nobs, number of rows in data file (excluding header lines)
  cat(3, file =file.name2,sep="\n",append = TRUE)      		# number of columns in the data file
  cat(1, file =file.name2,sep="\n",append = TRUE)        	# columns for X (observed inputs) in data file - comma-separated if several (order: t, h, T, T_smooth)
  cat(0, file =file.name2,sep="\n",append = TRUE)          # 8,9,10,11,12 !!! columns for Xu (random uncertainty in X, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="\n",append = TRUE)      	  # columns for Xb (systematic uncertainty in X, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="\n",append = TRUE)          # columns for Xb_indx (index of systematic errors in X - use 0 for a no-error assumption)
  cat(2, file =file.name2,sep="\n",append = TRUE)          # 16,18,20,21 columns for Y (observed outputs) in data file - comma-separated if several (order: Q, g0, cos(pi*g0), sin(pi*g0))
  cat(3, file =file.name2,sep="\n",append = TRUE)          # 17,19,22,23 columns for Yu (uncertainty in Y, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="\n",append = TRUE)      		# columns for Yb (systematic uncertainty in Y, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="",append = TRUE)         		# columns for Yb_indx (index of systematic errors in Y - use 0 for a no-error assumption)
  #----------------------------------------------------------------------
  file.mcmc = paste(dir_code,"/BaM_exe/Recession_h/Config_MCMC.txt",sep="")
  cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
  cat(nmcmc, file = file.mcmc, append = TRUE,sep="\n")   #Nadapt
  cat(ncycles, file = file.mcmc, append = TRUE,sep="\n")  #Ncycles
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #minMoveRate
  cat(0.5, file =file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
  cat(jump.neg, file =file.mcmc, append = TRUE,sep="\n")    #DownMult
  cat(jump.pos, file =file.mcmc, append = TRUE,sep="\n")    #UpMult
  cat(0, file =file.mcmc, append = TRUE,sep="\n")      #mode for init jump distr
  cat("****", file =file.mcmc, append = TRUE,sep="\n") 
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #MultFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     #RC MultiFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")      #Remnant MultiFact
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  #---------------------------------------------------------------------- 
  file.Pred1 = paste(dir_code,"/BaM_exe/Recession_h/Config_Pred_Master.txt",sep="")
  cat('3', file =file.Pred1,sep="\n")
  cat("'Config_Pred_Maxpost.txt'", file = file.Pred1, append = TRUE,sep="\n")
  cat("'Config_Pred_ParamU.txt'", file =file.Pred1, append = TRUE,sep="\n")
  cat("'Config_Pred_TotalU.txt'", file =file.Pred1, append = TRUE,sep="\n")
  #---------------------------------------------------------------------- 
  file.Pred3 = paste(dir_code,"/BaM_exe/Recession_h/Config_Pred_Maxpost.txt",sep="")
  cat("'Recession_h\\tgrid.txt'", file =file.Pred3, sep="\n")
  cat(ngrid, file =file.Pred3,sep="\n", append = TRUE)
  cat("1", file = file.Pred3, append = TRUE,sep="\n")   #n of spaghetti
  cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
  cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1", file = file.Pred3, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'ht_Maxpost.spag'", file = file.Pred3, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
  cat("'ht_Maxpost.env'", file = file.Pred3, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred3, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  #---------------------------------------------------------------------- 
  file.Pred4 = paste(dir_code,"/BaM_exe/Recession_h/Config_Pred_ParamU.txt",sep="")
  cat("'Recession_h\\tgrid.txt'", file =file.Pred4, sep="\n")
  cat(ngrid, file =file.Pred4,sep="\n", append = TRUE)
  cat("1", file = file.Pred4, append = TRUE,sep="\n")   #n of spaghetti
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
  cat(".false.", file = file.Pred4, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1", file = file.Pred4, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'ht_ParamU.spag'", file = file.Pred4, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
  cat("'ht_ParamU.env'", file = file.Pred4, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred4, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  #---------------------------------------------------------------------- 
  file.Pred5 = paste0(dir_code,"/BaM_exe/Recession_h/Config_Pred_TotalU.txt")
  cat("'Recession_h\\tgrid.txt'", file =file.Pred5,sep="\n")
  cat(ngrid,              file = file.Pred5, append = TRUE, sep="\n")                   #!!! n of grid
  cat("1",                file = file.Pred5, append = TRUE,sep="\n")                    #!!! n of spaghetti
  cat(".true.",           file = file.Pred5, append = TRUE,sep="\n")                    #!!! Propagate parametric uncertainty?
  cat(".true.",           file = file.Pred5, append = TRUE,sep="\n")                    #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1",               file = file.Pred5, append = TRUE,sep="\n")                    #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'ht_TotalU.spag'", file = file.Pred5, append = TRUE,sep="\n")                    #!!! Files containing spaghettis for each output variable (size nY)
  cat(".true.",           file = file.Pred5, append = TRUE,sep="\n")                    #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.",           file = file.Pred5, append = TRUE,sep="\n")                    #!!! Post-processing: create envelops? (size nY)
  cat("'ht_TotalU.env'",  file = file.Pred5, append = TRUE,sep="\n")                    #!!! Post-processing: name of envelop files (size nY)
  cat(".true.",           file = file.Pred5, append = TRUE,sep="\n")                    #!!! Print progress in console during computations?
  cat(".false." ,         file = file.Pred5, append = TRUE,sep="\n")                    #!!! Do state prediction? (size nState)
  #-------------------------------------------------------
  file.remnant = paste0(dir_code,"/BaM_exe/Recession_h/Config_RemnantSigma.txt")
  cat("'Linear'",     file = file.remnant, sep="\n")     #! Function f used in sdev=f(Qrc) 
  cat(2,              file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
  cat("gamma1",       file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
  cat(prior.gamma[4], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
  cat(prior.gamma[3], file = file.remnant, append = TRUE, sep="\n")            #! Prior distribution
  cat(prior.gamma[1], file = file.remnant, append = TRUE, sep=",")
  cat(",",            file = file.remnant, append = TRUE, sep=",")
  cat(prior.gamma[2], file = file.remnant, append = TRUE, sep="\n")
  #------------------------------------------------------------------
  cat("gamma2",       file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
  cat(prior.gamma[8], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
  cat(prior.gamma[7], file = file.remnant, append = TRUE, sep="\n")         #! Initial Guess
  cat(prior.gamma[5], file = file.remnant, append = TRUE, sep=",")
  cat(",",            file = file.remnant, append = TRUE, sep=",")
  cat(prior.gamma[6], file = file.remnant, append = TRUE, sep="\n")
  # cat(gamma[2], file =file.remnant, append = TRUE,sep="\n")
  # cat('"VAR"', file =file.remnant, append = TRUE,sep="\n")
  # cat('"Config_g2_VAR.txt"', file =file.remnant, append = TRUE,sep="\n")
  #-------------------------------------------------------------------
  
  #*********************************************************************************
  file.Pred6 = paste(dir_code,"/BaM_exe/Recession_h/Config_RunOptions.txt",sep="")
  cat(".true.", file = file.Pred6, sep="\n")  #Do MCMC?
  cat(".true.", file = file.Pred6, append = TRUE, sep="\n")  #Do MCMC summary?
  cat(".true.", file = file.Pred6, append = TRUE,sep="\n")  #Do Residual diagnostics?
  if (pred==TRUE) {
    cat(".true.", file = file.Pred6, append = TRUE,sep="\n")  #Do Predictions?
  } else {
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")  #Do Predictions?
  }
  ###################################################################   COOKING CONFIG
  file.cooking = paste(dir_code,"/BaM_exe/Recession_h/Config_Cooking.txt",sep="")
  cat("Results_MCMC_Cooked.txt" , file =file.cooking ,sep="\n")    #Result file
  cat(nburn, file =file.cooking, append = TRUE, sep="\n")            #Burn factor
  cat(nslim, file =file.cooking, append = TRUE, sep="\n")             #Nslim
  ###################################################################   RESIDUALS CONFIG
  file.residuals = paste(dir_code,"/BaM_exe/Recession_h//Config_Residuals.txt",sep="")
  cat("Results_Residuals.txt" , file =file.residuals ,sep="\n")    #Result file
  ###################################################################   SUMMARY CONFIG
  file.summary = paste(dir_code,"/BaM_exe/Recession_h//Config_Summary.txt",sep="")
  cat("Results_Summary.txt" , file =file.summary ,sep="\n")    #Result file
  
}

















######################################################################################
BaM_config.pooling <- function( dir.exe ,
                                model, 
                                nobs, 
                                ncycles,
                                nmcmc,
                                ncurves, 
                                jump.pos , 
                                jump.neg , 
                                slim, 
                                burn,
                                prior,
                                gamma.model,
                                prior.gamma) {  
######################################################################################
# This function writes the configuration files for BaM.exe for the estimation of
# stage-recession model through Bayesian pooling approach.
  tgrid     =  seq(0,60,0.1)
  ngrid     =  length(tgrid)
  theta2    =  -log(2)  + log(log(2))
  St_theta2 =  1
  theta4    =  -log(30) + log(log(2))
  St_theta4 =  1
  
  write.table(tgrid, file = paste0(dir.exe,"/Recession_h_pooling/tgrid.txt"), 
              col.names = FALSE, row.names = FALSE)

  
  #creation of Config_BaM.txt
  file.bam = paste(dir.exe,"/Config_BaM.txt",sep="")
  cat('"Recession_h_pooling/"',    file = file.bam, sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"',   file = file.bam, sep="\n", append = TRUE)    
  cat('"Config_Model.txt"',        file = file.bam, sep="\n", append = TRUE)
  cat('"Config_xtra.txt"',         file = file.bam, sep="\n", append = TRUE)
  cat('"Config_Data.txt"',         file = file.bam, sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file.bam, sep="\n", append = TRUE)                                       
  cat('"Config_MCMC.txt"',         file = file.bam, sep="\n", append = TRUE)                                            
  cat('"Config_Cooking.txt"',      file = file.bam, sep="\n", append = TRUE)
  cat('"Config_Summary.txt"',      file = file.bam, sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"',    file = file.bam, sep="\n", append = TRUE)
  cat('""',                        file = file.bam, sep="\n", append = TRUE)
  #----------------------------------------------------------------------
  file.name = paste0(dir.exe,"/Recession_h_pooling/Config_Model.txt")
  cat('"Recession_h"', file = file.name, sep="\n")
  cat(1,               file = file.name, append = TRUE, sep="\n")
  cat(1,               file = file.name, append = TRUE, sep="\n")

  
  if (model == "1expWithAsympt") { # ===================================================================> 1exp
    cat(3,             file =file.name, append = TRUE, sep="\n")  #Number of paramaters  5
    cat('"a1"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[4],      file =file.name, append = TRUE, sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[3],    file =file.name, append = TRUE, sep="\n")
      cat(prior[1],    file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[2],    file =file.name, append = TRUE, sep="\n")
    }
    #-----------------------------------------------------
    cat('"b1"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[9],      file =file.name, append = TRUE, sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[8],    file =file.name, append = TRUE, sep="\n")
      cat(prior[6],    file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[7],    file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------
    cat('"a2"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[14],     file =file.name, append = TRUE, sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE, sep="\n")
      cat(prior[11],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------------------------
    if (prior[5] == "var"){
      file.name.var1 = paste0(dir.exe,"/Recession_h_pooling/Config_a1_VAR.txt")
      cat(ncurves,     file =file.name.var1, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var1, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a1_',    file =file.name.var1, append = TRUE, sep="")
        cat(i,         file =file.name.var1, append = TRUE, sep="")  
        cat('"',       file =file.name.var1, append = TRUE, sep="\n")  
        cat(prior[4],  file =file.name.var1, append = TRUE, sep="\n")
        cat(prior[3],  file =file.name.var1, append = TRUE, sep="\n")
        cat(prior[1],  file =file.name.var1, append = TRUE, sep=",")
        cat(",",       file =file.name.var1, append = TRUE, sep=",")
        cat(prior[2],  file =file.name.var1, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[10] == "var"){
      file.name.var2 = paste0(dir.exe,"/Recession_h_pooling/Config_b1_VAR.txt")
      cat(ncurves,     file =file.name.var2, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var2, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b1_',    file =file.name.var2, append = TRUE, sep="")
        cat(i,         file =file.name.var2, append = TRUE, sep="")  
        cat('"',       file =file.name.var2, append = TRUE, sep="\n")  
        cat(prior[9],  file =file.name.var2, append = TRUE, sep="\n")
        cat(prior[8],  file =file.name.var2, append = TRUE, sep="\n")
        cat(prior[6],  file =file.name.var2, append = TRUE, sep=",")
        cat(",",       file =file.name.var2, append = TRUE, sep=",")
        cat(prior[7],  file =file.name.var2, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[15] == "var"){
      file.name.var3 = paste0(dir.exe,"/Recession_h_pooling/Config_a2_VAR.txt")
      cat(ncurves,     file =file.name.var3, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var3, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a2_',    file =file.name.var3, append = TRUE, sep="")
        cat(i,         file =file.name.var3, append = TRUE, sep="")  
        cat('"',       file =file.name.var3, append = TRUE, sep="\n")  
        cat(prior[14], file =file.name.var3, append = TRUE, sep="\n")
        cat(prior[13], file =file.name.var3, append = TRUE, sep="\n")
        cat(prior[11], file =file.name.var3, append = TRUE, sep=",")
        cat(",",       file =file.name.var3, append = TRUE, sep=",")
        cat(prior[12], file =file.name.var3, append = TRUE, sep="\n")
      }
    }
    #---------------------------------------------------------------
    # file extra:
    file.xtra = paste0(dir.exe,"/Recession_h_pooling/Config_xtra.txt")
    cat("exp",         file =file.xtra,  sep="\n")   # recession model: e.g. "exp", "expexp", "Boussinesq", "hyperb"
    cat(1,             file = file.xtra, append = TRUE, sep="\n")  # number of superposed terms 
    
    
  
    
  } else if ((model == "2expWithAsympt")|(model == "2expWithAsympt_bis")) { # ==================================> 2exp
    cat(5,             file =file.name, append = TRUE, sep="\n") # Number of paramaters  5
    cat('"a1"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[4],      file =file.name, append = TRUE, sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[3],    file =file.name, append = TRUE, sep="\n")
      cat(prior[1],    file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[2] ,   file =file.name, append = TRUE, sep="\n")
    }
    #-----------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE,sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"a2"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE,sep="\n")
      cat(prior[11],    file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------------
    cat('"b2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_b2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[18],   file =file.name, append = TRUE,sep="\n")
      cat(prior[16] ,  file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"a3"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[24],   file =file.name, append = TRUE,sep="\n")
    if (prior[25] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a3_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[23],   file =file.name, append = TRUE,sep="\n")
      cat(prior[21],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[22],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------------------------
    if (prior[5] == "var"){
      file.name.var1 = paste0(dir.exe,"/Recession_h_pooling/Config_a1_VAR.txt")
      cat(ncurves, file =file.name.var1, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var1, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a1_', file =file.name.var1, append = TRUE,sep="")
        cat(i, file =file.name.var1, append = TRUE,sep="")  
        cat('"', file =file.name.var1, append = TRUE,sep="\n")  
        cat(prior[4] , file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[3], file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[1], file =file.name.var1, append = TRUE,sep=",")
        cat(",",file =file.name.var1, append = TRUE, sep=",")
        cat(prior[2] , file =file.name.var1, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[10] == "var"){
      file.name.var2 = paste0(dir.exe,"/Recession_h_pooling/Config_b1_VAR.txt")
      cat(ncurves, file =file.name.var2, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var2, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b1_', file =file.name.var2, append = TRUE,sep="")
        cat(i, file =file.name.var2, append = TRUE,sep="")  
        cat('"', file =file.name.var2, append = TRUE,sep="\n")  
        cat(prior[9] , file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[8], file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[6], file =file.name.var2, append = TRUE,sep=",")
        cat(",",file =file.name.var2, append = TRUE, sep=",")
        cat(prior[7] , file =file.name.var2, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[15] == "var"){
      file.name.var3 = paste0(dir.exe,"/Recession_h_pooling/Config_a2_VAR.txt")
      cat(ncurves,     file =file.name.var3, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var3, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a2_',    file =file.name.var3, append = TRUE,sep="")
        cat(i,         file =file.name.var3, append = TRUE,sep="")  
        cat('"',       file =file.name.var3, append = TRUE,sep="\n")  
        cat(prior[14], file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[13], file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[11], file =file.name.var3, append = TRUE,sep=",")
        cat(",",       file =file.name.var3, append = TRUE, sep=",")
        cat(prior[12], file =file.name.var3, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[20] == "var"){
      file.name.var4 = paste0(dir.exe,"/Recession_h_pooling/Config_b2_VAR.txt")
      cat(ncurves,     file =file.name.var4, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var4, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b2_',    file =file.name.var4, append = TRUE,sep="")
        cat(i,         file =file.name.var4, append = TRUE,sep="")  
        cat('"',       file =file.name.var4, append = TRUE,sep="\n")  
        cat(prior[19], file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[18], file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[16], file =file.name.var4, append = TRUE,sep=",")
        cat(",",       file =file.name.var4, append = TRUE, sep=",")
        cat(prior[17], file =file.name.var4, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[25] == "var"){
      file.name.var5 = paste0(dir.exe,"/Recession_h_pooling/Config_a3_VAR.txt")
      cat(ncurves,     file =file.name.var5, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var5, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a3_',    file =file.name.var5, append = TRUE,sep="")
        cat(i,         file =file.name.var5, append = TRUE,sep="")  
        cat('"',       file =file.name.var5, append = TRUE,sep="\n")  
        cat(prior[24], file =file.name.var5, append = TRUE,sep="\n")
        cat(prior[23], file =file.name.var5, append = TRUE,sep="\n")
        cat(prior[21], file =file.name.var5, append = TRUE,sep=",")
        cat(",",       file =file.name.var5, append = TRUE,sep=",")
        cat(prior[22], file =file.name.var5, append = TRUE,sep="\n")
      }
    }
    #----------------------------------------------------------------
    # file extra:
    file.xtra = paste0(dir.exe,"/Recession_h_pooling/Config_xtra.txt")
    cat("exp", file = file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(2,     file = file.xtra, append = TRUE, sep="\n")  # number of superposed terms 
    
    
  } else if (model == "2expWithAsympt_rel") { # =========================================================> 2exp relatives
    cat(5,          file =file.name, append = TRUE,sep="\n")  #Number of paramaters  5
    cat('"a1"',     file =file.name, append = TRUE,sep="\n")
    cat(prior[4],   file =file.name, append = TRUE,sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"',  file =file.name, append = TRUE,sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[3], file =file.name, append = TRUE,sep="\n")
      cat(prior[1], file =file.name, append = TRUE,sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[2], file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------
    cat('"b1"',     file =file.name, append = TRUE,sep="\n")
    cat(prior[9],   file =file.name, append = TRUE,sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"',  file =file.name, append = TRUE,sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE,sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"a2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"',   file =file.name, append = TRUE,sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[13], file =file.name, append = TRUE,sep="\n")
      cat(prior[11], file =file.name, append = TRUE, sep=",")
      cat(",",       file =file.name, append = TRUE, sep=",")
      cat(prior[12], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"b2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[19],   file =file.name, append = TRUE,sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"',   file =file.name, append = TRUE,sep="\n")
      cat('"Config_b2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[18], file =file.name, append = TRUE,sep="\n")
      cat(prior[16], file =file.name, append = TRUE, sep=",")
      cat(",",       file =file.name, append = TRUE, sep=",")
      cat(prior[17], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"a3"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[24],   file =file.name, append = TRUE,sep="\n")
    if (prior[25] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a3_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[23],   file =file.name, append = TRUE,sep="\n")
      cat(prior[21],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[22],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------------------------
    if (prior[5] == "var"){
      file.name.var1 = paste0(dir.exe,"/Recession_h_pooling/Config_a1_VAR.txt")
      cat(ncurves, file =file.name.var1, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var1, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a1_', file =file.name.var1, append = TRUE,sep="")
        cat(i, file =file.name.var1, append = TRUE,sep="")  
        cat('"', file =file.name.var1, append = TRUE,sep="\n")  
        cat(prior[4] , file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[3], file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[1], file =file.name.var1, append = TRUE,sep=",")
        cat(",",file =file.name.var1, append = TRUE, sep=",")
        cat(prior[2] , file =file.name.var1, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[10] == "var"){
      file.name.var2 = paste0(dir.exe,"/Recession_h_pooling/Config_b1_VAR.txt")
      cat(ncurves, file =file.name.var2, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var2, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b1_', file =file.name.var2, append = TRUE,sep="")
        cat(i, file =file.name.var2, append = TRUE,sep="")  
        cat('"', file =file.name.var2, append = TRUE,sep="\n")  
        cat(prior[9] , file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[8], file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[6], file =file.name.var2, append = TRUE,sep=",")
        cat(",",file =file.name.var2, append = TRUE, sep=",")
        cat(prior[7] , file =file.name.var2, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[15] == "var"){
      file.name.var3 = paste0(dir.exe,"/Recession_h_pooling/Config_a2_VAR.txt")
      cat(ncurves, file =file.name.var3, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var3, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a2_', file =file.name.var3, append = TRUE,sep="")
        cat(i, file =file.name.var3, append = TRUE,sep="")  
        cat('"', file =file.name.var3, append = TRUE,sep="\n")  
        cat(prior[14] , file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[13], file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[11], file =file.name.var3, append = TRUE,sep=",")
        cat(",",file =file.name.var3, append = TRUE, sep=",")
        cat(prior[12] , file =file.name.var3, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[20] == "var"){
      file.name.var4 = paste0(dir.exe,"/Recession_h_pooling/Config_b2_VAR.txt")
      cat(ncurves, file =file.name.var4, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var4, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b2_', file =file.name.var4, append = TRUE,sep="")
        cat(i, file =file.name.var4, append = TRUE,sep="")  
        cat('"', file =file.name.var4, append = TRUE,sep="\n")  
        cat(prior[19] , file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[18], file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[16], file =file.name.var4, append = TRUE,sep=",")
        cat(",",file =file.name.var4, append = TRUE, sep=",")
        cat(prior[17] , file =file.name.var4, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[25] == "var"){
      file.name.var5 = paste0(dir.exe,"/Recession_h_pooling/Config_a3_VAR.txt")
      cat(ncurves, file =file.name.var5, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var5, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a3_', file =file.name.var5, append = TRUE,sep="")
        cat(i, file =file.name.var5, append = TRUE,sep="")  
        cat('"', file =file.name.var5, append = TRUE,sep="\n")  
        cat(prior[24] , file =file.name.var5, append = TRUE,sep="\n")
        cat(prior[23], file =file.name.var5, append = TRUE,sep="\n")
        cat(prior[21], file =file.name.var5, append = TRUE,sep=",")
        cat(",",file =file.name.var5, append = TRUE, sep=",")
        cat(prior[22] , file =file.name.var5, append = TRUE, sep="\n")
      }
    }
    #-----------------------------------------------------------------
    # file extra:
    file.xtra = paste0(dir.exe,"/Recession_h_pooling/Config_xtra.txt")
    cat("exp_rel", file =file.xtra,  sep="\n")    # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(2,         file = file.xtra, append = TRUE, sep="\n")  # number of superposed terms 
    
    

  } else if ((model == "3expWithAsympt")|(model == "3expWithAsympt_bis")) { #=====================================> 3exp
    cat(7,          file =file.name, append = TRUE,sep="\n")   #Number of paramaters  7
    cat('"a1"',     file =file.name, append = TRUE,sep="\n")
    cat(prior[4],   file =file.name, append = TRUE,sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"',  file =file.name, append = TRUE,sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[3], file =file.name, append = TRUE,sep="\n")
      cat(prior[1], file =file.name, append = TRUE,sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[2], file =file.name, append = TRUE, sep="\n")
    }
    #------------------------------------------------------
    cat('"b1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[9], file =file.name, append = TRUE,sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE,sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------
    cat('"a2"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[14],   file =file.name, append = TRUE,sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"',   file =file.name, append = TRUE,sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[13], file =file.name, append = TRUE,sep="\n")
      cat(prior[11], file =file.name, append = TRUE, sep=",")
      cat(",",       file =file.name, append = TRUE, sep=",")
      cat(prior[12], file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"b2"',       file =file.name, append = TRUE,sep="\n")
    cat(prior[19],    file =file.name, append = TRUE,sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"',    file =file.name, append = TRUE,sep="\n")
      cat('"Config_b2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[18],  file =file.name, append = TRUE,sep="\n")
      cat(prior[16],  file =file.name, append = TRUE, sep=",")
      cat(",",        file =file.name, append = TRUE, sep=",")
      cat(prior[17],  file =file.name, append = TRUE, sep="\n")
    }
    #---------------------------------------------------------
    cat('"a3"',       file =file.name, append = TRUE,sep="\n")
    cat(prior[24],    file =file.name, append = TRUE,sep="\n")
    if (prior[25] == "var"){
      cat('"VAR"',    file =file.name, append = TRUE,sep="\n")
      cat('"Config_a3_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[23],  file =file.name, append = TRUE,sep="\n")
      cat(prior[21],  file =file.name, append = TRUE, sep=",")
      cat(",",        file =file.name, append = TRUE, sep=",")
      cat(prior[22],  file =file.name, append = TRUE, sep="\n")
    }
    #----------------------------------------------------------
    cat('"b3"',      file =file.name, append = TRUE,sep="\n")
    cat(prior[29],   file =file.name, append = TRUE,sep="\n")
    if (prior[30] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_b3_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[28],   file =file.name, append = TRUE,sep="\n")
      cat(prior[26],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[27],   file =file.name, append = TRUE, sep="\n")
    }
    #---------------------------------------------------------
    cat('"a4"',       file =file.name, append = TRUE,sep="\n")
    cat(prior[34],    file =file.name, append = TRUE,sep="\n")
    if (prior[35] == "var"){
      cat('"VAR"',    file =file.name, append = TRUE,sep="\n")
      cat('"Config_a4_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[33],  file =file.name, append = TRUE,sep="\n")
      cat(prior[31],  file =file.name, append = TRUE, sep=",")
      cat(",",        file =file.name, append = TRUE, sep=",")
      cat(prior[32],  file =file.name, append = TRUE, sep="\n")
    }
    #---------------------------------------------------------
    if (prior[5] == "var"){
      file.name.var1 = paste0(dir.exe,"/Recession_h_pooling/Config_a1_VAR.txt")
      cat(ncurves,    file =file.name.var1, sep="\n")  #Number of paramaters
      cat(4,          file =file.name.var1, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a1_',   file =file.name.var1, append = TRUE,sep="")
        cat(i,        file =file.name.var1, append = TRUE,sep="")  
        cat('"',      file =file.name.var1, append = TRUE,sep="\n")  
        cat(prior[4], file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[3], file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[1], file =file.name.var1, append = TRUE,sep=",")
        cat(",",      file =file.name.var1, append = TRUE, sep=",")
        cat(prior[2], file =file.name.var1, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[10] == "var"){
      file.name.var2 = paste0(dir.exe,"/Recession_h_pooling/Config_b1_VAR.txt")
      cat(ncurves, file =file.name.var2, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var2, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b1_', file =file.name.var2, append = TRUE,sep="")
        cat(i, file =file.name.var2, append = TRUE,sep="")  
        cat('"', file =file.name.var2, append = TRUE,sep="\n")  
        cat(prior[9] , file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[8], file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[6], file =file.name.var2, append = TRUE,sep=",")
        cat(",",file =file.name.var2, append = TRUE, sep=",")
        cat(prior[7] , file =file.name.var2, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[15] == "var"){
      file.name.var3 = paste0(dir.exe,"/Recession_h_pooling/Config_a2_VAR.txt")
      cat(ncurves, file =file.name.var3, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var3, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a2_', file =file.name.var3, append = TRUE,sep="")
        cat(i, file =file.name.var3, append = TRUE,sep="")  
        cat('"', file =file.name.var3, append = TRUE,sep="\n")  
        cat(prior[14] , file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[13], file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[11], file =file.name.var3, append = TRUE,sep=",")
        cat(",",file =file.name.var3, append = TRUE, sep=",")
        cat(prior[12] , file =file.name.var3, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[20] == "var"){
      file.name.var4 = paste0(dir.exe,"/Recession_h_pooling/Config_b2_VAR.txt")
      cat(ncurves, file =file.name.var4, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var4, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b2_', file =file.name.var4, append = TRUE,sep="")
        cat(i, file =file.name.var4, append = TRUE,sep="")  
        cat('"', file =file.name.var4, append = TRUE,sep="\n")  
        cat(prior[19] , file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[18], file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[16], file =file.name.var4, append = TRUE,sep=",")
        cat(",",file =file.name.var4, append = TRUE, sep=",")
        cat(prior[17] , file =file.name.var4, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[25] == "var"){
      file.name.var5 = paste0(dir.exe,"/Recession_h_pooling/Config_a3_VAR.txt")
      cat(ncurves, file =file.name.var5, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var5, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a3_', file =file.name.var5, append = TRUE,sep="")
        cat(i, file =file.name.var5, append = TRUE,sep="")  
        cat('"', file =file.name.var5, append = TRUE,sep="\n")  
        cat(prior[24] , file =file.name.var5, append = TRUE,sep="\n")
        cat(prior[23], file =file.name.var5, append = TRUE,sep="\n")
        cat(prior[21], file =file.name.var5, append = TRUE,sep=",")
        cat(",",file =file.name.var5, append = TRUE, sep=",")
        cat(prior[22] , file =file.name.var5, append = TRUE, sep="\n")
      }
    }
    #-------------------------------------------------------------------------
    if (prior[30] == "var"){
      file.name.var6 = paste0(dir.exe,"/Recession_h_pooling/Config_b3_VAR.txt")
      cat(ncurves,      file =file.name.var6, sep="\n")  #Number of paramaters
      cat(4,            file =file.name.var6, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a3_',     file =file.name.var6, append = TRUE, sep="")
        cat(i,          file =file.name.var6, append = TRUE, sep="")  
        cat('"',        file =file.name.var6, append = TRUE, sep="\n")  
        cat(prior[29],  file =file.name.var6, append = TRUE, sep="\n")
        cat(prior[28],  file =file.name.var6, append = TRUE, sep="\n")
        cat(prior[26],  file =file.name.var6, append = TRUE, sep=",")
        cat(",",        file =file.name.var6, append = TRUE, sep=",")
        cat(prior[27],  file =file.name.var6, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[35] == "var"){
      file.name.var7 = paste0(dir.exe,"/Recession_h_pooling/Config_a4_VAR.txt")
      cat(ncurves,      file =file.name.var7, sep="\n")  #Number of paramaters
      cat(4,            file =file.name.var7, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a4_',     file =file.name.var7, append = TRUE, sep="")
        cat(i,          file =file.name.var7, append = TRUE, sep="")  
        cat('"',        file =file.name.var7, append = TRUE, sep="\n")  
        cat(prior[34],  file =file.name.var7, append = TRUE, sep="\n")
        cat(prior[33],  file =file.name.var7, append = TRUE, sep="\n")
        cat(prior[31],  file =file.name.var7, append = TRUE, sep=",")
        cat(",",        file =file.name.var7, append = TRUE, sep=",")
        cat(prior[32],  file =file.name.var7, append = TRUE, sep="\n")
      }
    }
    #-----------------------------------------------------------------
    # file extra:
    file.xtra = paste0(dir.exe,"/Recession_h_pooling/Config_xtra.txt")
    cat("exp", file = file.xtra,  sep="\n")             # recession model
    cat(3,     file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    
    
  } else if ((model == "expexp")|(model == "expexp_bis")) { # ======================================================> 1 double exp
    cat(4,             file =file.name, append = TRUE, sep="\n")      #Number of paramaters  5
    cat('"a1"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[4],      file =file.name, append = TRUE, sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[3],    file =file.name, append = TRUE, sep="\n")
      cat(prior[1],    file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[2] ,   file =file.name, append = TRUE, sep="\n")
    }
    #-----------------------------------------------------
    cat('"b1"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[9],      file =file.name, append = TRUE, sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[8],    file =file.name, append = TRUE, sep="\n")
      cat(prior[6],    file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[7],    file =file.name, append = TRUE, sep="\n")
    }
    #------------------------------------------------------
    cat('"n1"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[14],     file =file.name, append = TRUE, sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_n1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE, sep="\n")
      cat(prior[11],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"a2"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[19],     file =file.name, append = TRUE, sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[18],   file =file.name, append = TRUE, sep="\n")
      cat(prior[16],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------------------------
    if (prior[5] == "var"){
      file.name.var1 = paste0(dir.exe,"/Recession_h_pooling/Config_a1_VAR.txt")
      cat(ncurves,     file =file.name.var1, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var1, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves){
        cat('"a1_',    file =file.name.var1, append = TRUE, sep="")
        cat(i,         file =file.name.var1, append = TRUE, sep="")  
        cat('"',       file =file.name.var1, append = TRUE, sep="\n")  
        cat(prior[4],  file =file.name.var1, append = TRUE, sep="\n")
        cat(prior[3],  file =file.name.var1, append = TRUE, sep="\n")
        cat(prior[1],  file =file.name.var1, append = TRUE, sep=",")
        cat(",",       file =file.name.var1, append = TRUE, sep=",")
        cat(prior[2],  file =file.name.var1, append = TRUE, sep="\n")
      }
    }
    #-------------------------------------------------------------------------
    if (prior[10] == "var"){
      file.name.var2 = paste0(dir.exe,"/Recession_h_pooling/Config_b1_VAR.txt")
      cat(ncurves,     file =file.name.var2, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var2, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b1_',    file =file.name.var2, append = TRUE,sep="")
        cat(i,         file =file.name.var2, append = TRUE,sep="")  
        cat('"',       file =file.name.var2, append = TRUE,sep="\n")  
        cat(prior[9],  file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[8],  file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[6],  file =file.name.var2, append = TRUE,sep=",")
        cat(",",       file =file.name.var2, append = TRUE, sep=",")
        cat(prior[7] , file =file.name.var2, append = TRUE, sep="\n")
      }
    }
    #-------------------------------------------------------------------------
    if (prior[15] == "var"){
      file.name.var3 = paste0(dir.exe,"/Recession_h_pooling/Config_n1_VAR.txt")
      cat(ncurves,      file =file.name.var3, sep="\n")  #Number of paramaters
      cat(4,            file =file.name.var3, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"n1_',     file =file.name.var3, append = TRUE, sep="")
        cat(i,          file =file.name.var3, append = TRUE, sep="")  
        cat('"',        file =file.name.var3, append = TRUE, sep="\n")  
        cat(prior[14],  file =file.name.var3, append = TRUE, sep="\n")
        cat(prior[13],  file =file.name.var3, append = TRUE, sep="\n")
        cat(prior[11],  file =file.name.var3, append = TRUE, sep=",")
        cat(",",        file =file.name.var3, append = TRUE, sep=",")
        cat(prior[12],  file =file.name.var3, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[20] == "var"){
      file.name.var4 = paste0(dir.exe,"/Recession_h_pooling/Config_a2_VAR.txt")
      cat(ncurves,      file =file.name.var4, sep="\n")  #Number of paramaters
      cat(4,            file =file.name.var4, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a2_',     file =file.name.var4, append = TRUE, sep="")
        cat(i,          file =file.name.var4, append = TRUE, sep="")  
        cat('"',        file =file.name.var4, append = TRUE, sep="\n")  
        cat(prior[19],  file =file.name.var4, append = TRUE, sep="\n")
        cat(prior[18],  file =file.name.var4, append = TRUE, sep="\n")
        cat(prior[16],  file =file.name.var4, append = TRUE, sep=",")
        cat(",",        file =file.name.var4, append = TRUE, sep=",")
        cat(prior[17] , file =file.name.var4, append = TRUE, sep="\n")
      }
    }
    #-----------------------------------------------------------------
    # file extra:
    file.xtra = paste0(dir.exe,"/Recession_h_pooling/Config_xtra.txt")
    cat("expexp", file = file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1,        file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 

    
    
  } else if ((model == "hyperb")|(model == "hyperb_bis")) { # ======================================================> 1 double exp
    cat(4,          file =file.name, append = TRUE, sep="\n")    #Number of paramaters  5
    cat('"a1"',     file =file.name, append = TRUE, sep="\n")
    cat(prior[4],   file =file.name, append = TRUE, sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"',  file =file.name, append = TRUE, sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[3], file =file.name, append = TRUE,sep="\n")
      cat(prior[1], file =file.name, append = TRUE,sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[2], file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------
    cat('"b1"',     file =file.name, append = TRUE, sep="\n")
    cat(prior[9],   file =file.name, append = TRUE, sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"',  file =file.name, append = TRUE, sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE, sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------
    cat('"n1"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[14],     file =file.name, append = TRUE, sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"',     file =file.name, append = TRUE, sep="\n")
      cat('"Config_n1_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE, sep="\n")
      cat(prior[11],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #---------------------------------------------------------
    cat('"a2"',        file =file.name, append = TRUE, sep="\n")
    cat(prior[19],     file =file.name, append = TRUE, sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"', file =file.name, append = TRUE, sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE, sep="\n")
    } else {
      cat(prior[18],   file =file.name, append = TRUE, sep="\n")
      cat(prior[16],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------------------------
    if (prior[5] == "var"){
      file.name.var1 = paste0(dir.exe,"/Recession_h_pooling/Config_a1_VAR.txt")
      cat(ncurves,     file =file.name.var1, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var1, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a1_',    file =file.name.var1, append = TRUE, sep="")
        cat(i,         file =file.name.var1, append = TRUE, sep="")  
        cat('"',       file =file.name.var1, append = TRUE, sep="\n")  
        cat(prior[4],  file =file.name.var1, append = TRUE, sep="\n")
        cat(prior[3],  file =file.name.var1, append = TRUE, sep="\n")
        cat(prior[1],  file =file.name.var1, append = TRUE, sep=",")
        cat(",",       file =file.name.var1, append = TRUE, sep=",")
        cat(prior[2] , file =file.name.var1, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[10] == "var"){
      file.name.var2 = paste0(dir.exe,"/Recession_h_pooling/Config_b1_VAR.txt")
      cat(ncurves, file =file.name.var2, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var2, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b1_', file =file.name.var2, append = TRUE,sep="")
        cat(i, file =file.name.var2, append = TRUE,sep="")  
        cat('"', file =file.name.var2, append = TRUE,sep="\n")  
        cat(prior[9] , file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[8], file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[6], file =file.name.var2, append = TRUE,sep=",")
        cat(",",file =file.name.var2, append = TRUE, sep=",")
        cat(prior[7] , file =file.name.var2, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[15] == "var"){
      file.name.var3 = paste0(dir.exe,"/Recession_h_pooling/Config_n1_VAR.txt")
      cat(ncurves, file =file.name.var3, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var3, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"n1_', file =file.name.var3, append = TRUE,sep="")
        cat(i, file =file.name.var3, append = TRUE,sep="")  
        cat('"', file =file.name.var3, append = TRUE,sep="\n")  
        cat(prior[14] , file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[13], file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[11], file =file.name.var3, append = TRUE,sep=",")
        cat(",",file =file.name.var3, append = TRUE, sep=",")
        cat(prior[12] , file =file.name.var3, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[20] == "var"){
      file.name.var4 = paste0(dir.exe,"/Recession_h_pooling/Config_a2_VAR.txt")
      cat(ncurves, file =file.name.var4, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var4, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a2_', file =file.name.var4, append = TRUE,sep="")
        cat(i, file =file.name.var4, append = TRUE,sep="")  
        cat('"', file =file.name.var4, append = TRUE,sep="\n")  
        cat(prior[19] , file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[18], file =file.name.var4, append = TRUE,sep="\n")
        cat(prior[16], file =file.name.var4, append = TRUE,sep=",")
        cat(",",file =file.name.var4, append = TRUE, sep=",")
        cat(prior[17] , file =file.name.var4, append = TRUE, sep="\n")
      }
    }
    #-----------------------------------------------------------------
    # file extra:
    file.xtra = paste0(dir.exe,"/Recession_h_pooling/Config_xtra.txt")
    cat("hyperb", file =file.xtra,  sep="\n")             # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1, file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
    
    
    
    
    
    
  } else if ((model == "Coutagne")|(model == "Coutagne_bis")) { # =====================================================> Coutagne
    cat(4, file =file.name, append = TRUE,sep="\n")       #Number of paramaters  5
    cat('"a1"',   file =file.name, append = TRUE,sep="\n")
    cat(prior[4], file =file.name, append = TRUE,sep="\n")
    if (prior[5] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[3],  file =file.name, append = TRUE,sep="\n")
      cat(prior[1],  file =file.name, append = TRUE,sep=",")
      cat(",",       file =file.name, append = TRUE, sep=",")
      cat(prior[2] , file =file.name, append = TRUE, sep="\n")
    }
    #------------------------------------------------------
    cat('"b1"',     file =file.name, append = TRUE,sep="\n")
    cat(prior[9],   file =file.name, append = TRUE,sep="\n")
    if (prior[10] == "var"){
      cat('"VAR"',  file =file.name, append = TRUE,sep="\n")
      cat('"Config_b1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[8], file =file.name, append = TRUE,sep="\n")
      cat(prior[6], file =file.name, append = TRUE, sep=",")
      cat(",",      file =file.name, append = TRUE, sep=",")
      cat(prior[7], file =file.name, append = TRUE, sep="\n")
    }
    #-----------------------------------------------------
    cat('"n1"',        file =file.name, append = TRUE,sep="\n")
    cat(prior[14],     file =file.name, append = TRUE,sep="\n")
    if (prior[15] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_n1_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[13],   file =file.name, append = TRUE,sep="\n")
      cat(prior[11],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[12],   file =file.name, append = TRUE, sep="\n")
    }
    #--------------------------------------------------------
    cat('"a2"',        file =file.name, append = TRUE,sep="\n")
    cat(prior[19],     file =file.name, append = TRUE,sep="\n")
    if (prior[20] == "var"){
      cat('"VAR"', file =file.name, append = TRUE,sep="\n")
      cat('"Config_a2_VAR.txt"', file =file.name, append = TRUE,sep="\n")
    } else {
      cat(prior[18],   file =file.name, append = TRUE,sep="\n")
      cat(prior[16],   file =file.name, append = TRUE, sep=",")
      cat(",",         file =file.name, append = TRUE, sep=",")
      cat(prior[17],   file =file.name, append = TRUE, sep="\n")
    }
    #-------------------------------------------------------------------------
    if (prior[5] == "var"){
      file.name.var1 = paste0(dir.exe,"/Recession_h_pooling/Config_a1_VAR.txt")
      cat(ncurves,    file =file.name.var1, sep="\n")  #Number of paramaters
      cat(4,          file =file.name.var1, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a1_',   file =file.name.var1, append = TRUE,sep="")
        cat(i,        file =file.name.var1, append = TRUE,sep="")  
        cat('"',      file =file.name.var1, append = TRUE,sep="\n")  
        cat(prior[4], file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[3], file =file.name.var1, append = TRUE,sep="\n")
        cat(prior[1], file =file.name.var1, append = TRUE,sep=",")
        cat(",",      file =file.name.var1, append = TRUE, sep=",")
        cat(prior[2], file =file.name.var1, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[10] == "var"){
      file.name.var2 = paste0(dir.exe,"/Recession_h_pooling/Config_b1_VAR.txt")
      cat(ncurves, file =file.name.var2, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var2, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"b1_', file =file.name.var2, append = TRUE,sep="")
        cat(i, file =file.name.var2, append = TRUE,sep="")  
        cat('"', file =file.name.var2, append = TRUE,sep="\n")  
        cat(prior[9] , file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[8], file =file.name.var2, append = TRUE,sep="\n")
        cat(prior[6], file =file.name.var2, append = TRUE,sep=",")
        cat(",",file =file.name.var2, append = TRUE, sep=",")
        cat(prior[7] , file =file.name.var2, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[15] == "var"){
      file.name.var3 = paste0(dir.exe,"/Recession_h_pooling/Config_n1_VAR.txt")
      cat(ncurves, file =file.name.var3, sep="\n")  #Number of paramaters
      cat(4, file =file.name.var3, append = TRUE,sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"n1_', file =file.name.var3, append = TRUE,sep="")
        cat(i, file =file.name.var3, append = TRUE,sep="")  
        cat('"', file =file.name.var3, append = TRUE,sep="\n")  
        cat(prior[14] , file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[13], file =file.name.var3, append = TRUE,sep="\n")
        cat(prior[11], file =file.name.var3, append = TRUE,sep=",")
        cat(",",file =file.name.var3, append = TRUE, sep=",")
        cat(prior[12] , file =file.name.var3, append = TRUE, sep="\n")
      }
    }
    #--------------------------------------------------------------------------
    if (prior[20] == "var"){
      file.name.var4 = paste0(dir.exe,"/Recession_h_pooling/Config_a2_VAR.txt")
      cat(ncurves,     file =file.name.var4, sep="\n")  #Number of paramaters
      cat(4,           file =file.name.var4, append = TRUE, sep="\n") #column with period index    
      for (i in 1:ncurves) {
        cat('"a2_',    file =file.name.var4, append = TRUE, sep="")
        cat(i,         file =file.name.var4, append = TRUE, sep="")  
        cat('"',       file =file.name.var4, append = TRUE, sep="\n")  
        cat(prior[19], file =file.name.var4, append = TRUE, sep="\n")
        cat(prior[18], file =file.name.var4, append = TRUE, sep="\n")
        cat(prior[16], file =file.name.var4, append = TRUE, sep=",")
        cat(",",       file =file.name.var4, append = TRUE, sep=",")
        cat(prior[17], file =file.name.var4, append = TRUE, sep="\n")
      }
    }
    #----------------------------------------------------------------
    # file extra:
    file.xtra = paste0(dir.exe,"/Recession_h_pooling/Config_xtra.txt")
    cat("Coutagne", file =file.xtra,  sep="\n")        # recession model: "exp", "expexp", "Boussinesq", "hyperb"
    cat(1,          file = file.xtra, append = TRUE, sep="\n")  # number of exp terms 
  }
  
  
  
  
  
  ###################################################################################################################
  file.name2 = paste0(dir.exe,"/Recession_h_pooling/Config_Data.txt")
  cat("'Recession_h_pooling\\Curves_Data.txt'", file =file.name2, sep="\n")# path to data file
  cat(1, file =file.name2,sep="\n",append = TRUE)      	# number of header lines
  cat(nobs, file =file.name2,sep="\n",append = TRUE)    # Nobs, number of rows in data file (excluding header lines)
  cat(4, file =file.name2,sep="\n",append = TRUE)   		# number of columns in the data file
  cat(1, file =file.name2,sep="\n",append = TRUE)      	# columns for X (observed inputs) in data file - comma-separated if several (order: t, h, T, T_smooth)
  cat(0, file =file.name2,sep="\n",append = TRUE)       # 8,9,10,11,12 !!! columns for Xu (random uncertainty in X, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="\n",append = TRUE)    	  # columns for Xb (systematic uncertainty in X, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="\n",append = TRUE)       # columns for Xb_indx (index of systematic errors in X - use 0 for a no-error assumption)
  cat(2, file =file.name2,sep="\n",append = TRUE)       # 16,18,20,21 columns for Y (observed outputs) in data file - comma-separated if several (order: Q, g0, cos(pi*g0), sin(pi*g0))
  cat(3, file =file.name2,sep="\n",append = TRUE)       # 17,19,22,23 columns for Yu (uncertainty in Y, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="\n",append = TRUE)     	# columns for Yb (systematic uncertainty in Y, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0, file =file.name2,sep="",append = TRUE)     		# columns for Yb_indx (index of systematic errors in Y - use 0 for a no-error assumption)
  #---------------------------------------------------------
  file.mcmc = paste0(dir.exe,"/Recession_h_pooling/Config_MCMC.txt")
  cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
  cat(nmcmc,    file =file.mcmc, append = TRUE,sep="\n")    #Nadapt
  cat(ncycles,  file =file.mcmc, append = TRUE,sep="\n")    #Ncycles
  cat(0.1,      file =file.mcmc, append = TRUE,sep="\n")    #minMoveRate
  cat(0.5,      file =file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
  cat(jump.neg, file =file.mcmc, append = TRUE,sep="\n")    #DownMult
  cat(jump.pos, file =file.mcmc, append = TRUE,sep="\n")    #UpMult
  cat(0,        file =file.mcmc, append = TRUE,sep="\n")    #mode for init jump distr
  cat("****",   file =file.mcmc, append = TRUE,sep="\n") 
  cat(0.1,      file =file.mcmc, append = TRUE,sep="\n")    #MultFact
  cat(0.1,      file =file.mcmc, append = TRUE,sep=",")     #RC MultiFact
  cat(0.1,      file =file.mcmc, append = TRUE,sep=",")     
  cat(0.1,      file =file.mcmc, append = TRUE,sep="\n")
  cat(0.1,      file =file.mcmc, append = TRUE,sep=",")     #Remnant MultiFact
  cat(0.1,      file =file.mcmc, append = TRUE,sep="\n")
  #------------------------------------------------------------------------
  file.Pred1 = paste0(dir.exe,"/Recession_h_pooling/Config_Pred_Master.txt")
  cat('3', file =file.Pred1,sep="\n")
  cat("'Config_Pred_Maxpost.txt'", file = file.Pred1, append = TRUE,sep="\n")
  cat("'Config_Pred_ParamU.txt'", file =file.Pred1, append = TRUE,sep="\n")
  cat("'Config_Pred_TotalU.txt'", file =file.Pred1, append = TRUE,sep="\n")
  #------------------------------------------------------------------------- 
  file.Pred3 = paste0(dir.exe,"/Recession_h_pooling/Config_Pred_Maxpost.txt")
  cat('"Recession_h\\tgrid.txt"', file =file.Pred3, sep="\n")
  cat(ngrid, file =file.Pred3,sep="\n", append = TRUE)
  cat("1", file = file.Pred3, append = TRUE,sep="\n")                                   #!!! n of spaghetti
  cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                             #!!! Propagate parametric uncertainty?
  cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1", file = file.Pred3, append = TRUE,sep="\n")                                  #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'ht_Maxpost.spag'", file = file.Pred3, append = TRUE,sep="\n")                   #!!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
  cat("'ht_Maxpost.env'", file = file.Pred3, append = TRUE,sep="\n")                    #!!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred3, append = TRUE,sep="\n")                            #!!! Do state prediction? (size nState)
  #------------------------------------------------------------------------- 
  file.Pred4 = paste0(dir.exe,"/Recession_h_pooling/Config_Pred_ParamU.txt")
  cat("'Recession_h\\tgrid.txt'", file =file.Pred4, sep="\n")
  cat(ngrid, file =file.Pred4,sep="\n", append = TRUE)
  cat("1", file = file.Pred4, append = TRUE,sep="\n")                                   #!!! n of spaghetti
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
  cat(".false.", file = file.Pred4, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1", file = file.Pred4, append = TRUE,sep="\n")                                  #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'ht_ParamU.spag'", file = file.Pred4, append = TRUE,sep="\n")                    #!!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
  cat("'ht_ParamU.env'", file = file.Pred4, append = TRUE,sep="\n")                     #!!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred4, append = TRUE,sep="\n")                            #!!! Do state prediction? (size nState)
  #------------------------------------------------------------------------- 
  file.Pred5 = paste0(dir.exe,"/Recession_h_pooling/Config_Pred_TotalU.txt")
  cat("'Recession_h\\tgrid.txt'", file =file.Pred5,sep="\n")
  cat(ngrid, file =file.Pred5, append = TRUE, sep="\n")
  cat("1", file = file.Pred5, append = TRUE,sep="\n")                                  #!!! n of spaghetti
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Propagate parametric uncertainty?
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1", file = file.Pred5, append = TRUE,sep="\n")                                 #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'ht_TotalU.spag'", file = file.Pred5, append = TRUE,sep="\n")                   #!!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Post-processing: create envelops? (size nY)
  cat("'ht_TotalU.env'", file = file.Pred5, append = TRUE,sep="\n")                    #!!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred5, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  #---------------------------------------------------------------------------- 
  file.remnant = paste0(dir.exe,"/Recession_h_pooling/Config_RemnantSigma.txt")
  if (gamma.model == "constant"){
    cat("'Constant'", file = file.remnant, sep="\n")                          #! Function f used in sdev=f(Qrc) 
    cat(1, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
    cat("gamma", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(prior.gamma[4], file = file.remnant, append = TRUE, sep="\n")       #! Initial Guess
    cat(prior.gamma[3], file = file.remnant, append = TRUE, sep="\n")       #! Prior distribution
    cat(prior.gamma[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(prior.gamma[2],file =file.remnant, append = TRUE, sep="\n")
    
  } else if (gamma.model == "linear"){
    cat("'Linear'", file = file.remnant, sep="\n")                          #! Function f used in sdev=f(Qrc) 
    cat(2, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
    cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(prior.gamma[4], file = file.remnant, append = TRUE, sep="\n")       #! Initial Guess
    cat(prior.gamma[3], file = file.remnant, append = TRUE, sep="\n")       #! Prior distribution
    cat(prior.gamma[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(prior.gamma[2],file =file.remnant, append = TRUE, sep="\n")
    #----------------------------------------------------------------
    cat("gamma2", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(prior.gamma[8], file = file.remnant, append = TRUE, sep="\n")       #! Initial Guess
    cat(prior.gamma[7], file = file.remnant, append = TRUE, sep="\n")       #! Initial Guess
    cat(prior.gamma[5],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(prior.gamma[6],file =file.remnant, append = TRUE, sep="\n")
    # cat(gamma[2], file =file.remnant, append = TRUE,sep="\n")
    # cat('"VAR"', file =file.remnant, append = TRUE,sep="\n")
    # cat('"Config_g2_VAR.txt"', file =file.remnant, append = TRUE,sep="\n")
  }
  
  #--------------------------------------------------------------------------------  RUNNING OPTIONS
  file.Pred6 = paste0(dir.exe,"/Recession_h_pooling/Config_RunOptions.txt")
  cat(".true.", file = file.Pred6, sep="\n")                 # Do MCMC?
  cat(".true.", file = file.Pred6, append = TRUE, sep="\n")  # Do MCMC summary?
  cat(".true.", file = file.Pred6, append = TRUE,sep="\n")   # Do Residual diagnostics?
  cat(".false.", file = file.Pred6, append = TRUE,sep="\n")  # Do Predictions?
  #--------------------------------------------------------------------------------  COOKING CONFIG
  file.cooking = paste0(dir.exe,"/Recession_h_pooling/Config_Cooking.txt")
  cat("'Results_MCMC_Cooked.txt'" , file =file.cooking ,sep="\n")  #Result file
  cat(burn, file =file.cooking, append = TRUE, sep="\n")           #Burn factor
  cat(slim, file =file.cooking, append = TRUE, sep="\n")           #Nslim
  #--------------------------------------------------------------------------------  RESIDUALS CONFIG
  file.residuals = paste0(dir.exe,"/Recession_h_pooling//Config_Residuals.txt")
  cat("'Results_Residuals.txt'" , file =file.residuals ,sep="\n")    #Result file
  #--------------------------------------------------------------------------------  SUMMARY CONFIG
  file.summary = paste0(dir.exe,"/Recession_h_pooling//Config_Summary.txt")
  cat("'Results_Summary.txt'" , file =file.summary ,sep="\n")    #Result file
}


































#############################################################################################################
read.results.regression.rec = function(dir.recess,
                                       dir.segm.recessions, 
                                       dir.BaM.rec.pool,
                                       rec.model, 
                                       BayesianOption,
                                       which.recession, 
                                       time.rec,
                                       is.b1.var,
                                       time.limits,
                                       asymptote.limits,
                                       stage.scale.shift ,
                                       chi) {
  ############################################################################################################## 
  # Author: Matteo Darienzo, Inrae
  # Objective: this function read the results of the multi stage-recession estimation,
  #            In particular it provides statistics of temporal behavior of the recession model parameters 
  #            (e.g. the asymptotica stage)
  
  ############################
  if (BayesianOption ==1) {
  ############################
    summary.rec =NULL
    h.infinity=NULL; recess.mcmc = NULL; recess.env=NULL; recess.data=NULL
    for (r in which.recession) {
      if (long.recessions == TRUE) {
        recess.mcmc[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                  r,"_long/BaM/Results_MCMC_cooked.txt"),header=TRUE)
        summary.rec[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C", r,
                                                  "_long/BaM/Results_Summary.txt"),header=TRUE)
        h.infinity[[r]] = c(quantile(recess.mcmc[[r]][,7], p = c(0.025, 0.5, 0.975)), 
                            mean =summary.rec[[r]]$a4[5],  
                            stdev =  summary.rec[[r]]$a4[11], 
                            MAP=  summary.rec[[r]]$a4[16], 
                            t =rec.selection[[2]][r])
        recess.data[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                  r,"_long/BaM/Curves_Data.txt"),header=TRUE)
        recess.data[[r]] = cbind(recess.data[[r]], t= rep(r,length(recess.data[[r]][1]) ))
        recess.env[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                 r,"_long/BaM/ht_TotalU.env"),header=TRUE)  
      } else {
        recess.mcmc[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                  r,"/BaM/Results_MCMC_cooked.txt"),header=TRUE)
        summary.rec[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                  r,"/BaM/Results_Summary.txt"),header=TRUE)
        h.infinity[[r]] = c(quantile(recess.mcmc[[r]][,7], p = c(0.025, 0.5, 0.975)), 
                            mean =summary.rec[[r]]$a4[5],  
                            stdev =  summary.rec[[r]]$a4[11], 
                            MAP=  summary.rec[[r]]$a4[16], 
                            t =rec.selection[[2]][r])
        recess.data[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                  r,"/BaM/Curves_Data.txt"),header=TRUE)
        recess.data[[r]] = cbind(recess.data[[r]], t= rep(r,length(recess.data[[r]][1]) ))
        recess.env[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                 r,"/BaM/ht_TotalU.env"),header=TRUE)
      }
    }
    df.h.infinity = data.frame(t(sapply(h.infinity, function(x) x[1:max(lengths(h.infinity))])))
    write.table(df.h.infinity, file =paste0(dir.segm.recessions,"/df.h.infinity.txt"),
                sep="\t", row.names=FALSE, col.names = TRUE)
    return(list(recess.mcmc=recess.mcmc, 
                summary.rec=summary.rec,
                h.infinity=h.infinity,
                recess.data=recess.data,
                recess.env=recess.env,
                df.h.infinity=df.h.infinity))
    
    

    
  #########################################################################################
  } else if (BayesianOption==2) {
  #########################################################################################
    # directories where to save the results of var parameters for the segmentation:
    dir.rec.pool.test               = dir.recess
    dir.create(paste0(dir.segm.recessions,"/",rec.model))
    dir.segm.recessions.model       = paste0(dir.segm.recessions,"/", rec.model)
    dir.create(paste0(dir.segm.recessions.model,"/chi_",chi))
    dir.segm.recessions.model.param = paste0(dir.segm.recessions.model,"/chi_",chi)
    #copy results of regression to the folder for the segmentation:
    list.files.pool <- c( paste0(dir.rec.pool.test,"/Results_MCMC_Cooked.txt"),
                          paste0(dir.rec.pool.test,"/Results_Residuals.txt"),
                          paste0(dir.rec.pool.test,"/Results_Summary.txt"))
    for (i in 1:length(list.files.pool)) {
             file.copy(list.files.pool[i], dir.segm.recessions.model.param ,overwrite = TRUE)
    }
    # read results BaM files :
    summary.rec   = read.table(file = paste0(dir.segm.recessions.model.param, "/Results_Summary.txt"), header=TRUE)
    mcmc.rec      = read.table(file = paste0(dir.segm.recessions.model.param, "/Results_MCMC_cooked.txt"), header=TRUE)
    residuals.rec = read.table(file = paste0(dir.segm.recessions.model.param, "/Results_Residuals.txt"), header=TRUE)

    # Info:
    #######
    # For each model there is a different parameterization (a1,b1, a2,b2, a3, ...). 
    # The last parameter is the asymptotic stage.
    # and (k) indicate wheter the parameter is recession-specific or stationary (in common between all recessions)
    
    if (rec.model =="3expWithAsympt"){
      #################################################################################################     
      # read single parameters results :   a1(k) , b1, a2, b2, a3, b3, a4(k) 
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean    = mean(a1.mcmc[,1]), 
                               stdev   = std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
      b1.summary = summary.rec[,tail(which.recession,1)+1]
      b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                               stdev =std(b1.mcmc), maxpost = b1.summary[16])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      #a2:
      a2.mcmc = mcmc.rec[,tail(which.recession,1) +2]
      a2.summary = summary.rec[,tail(which.recession,1) +2]
      a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc), 
                               stdev =std(a2.mcmc), maxpost = a2.summary[16])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      # write.table(a2.df, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
      #             sep="\t", row.names=FALSE, col.names = TRUE)
      
      #b2:
      b2.mcmc = mcmc.rec[,tail(which.recession,1)+3]
      b2.summary = summary.rec[,tail(which.recession,1) +3]
      b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b2.mcmc), 
                               stdev =std(b2.mcmc), maxpost = b2.summary[16])))
      names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #a3:
      a3.mcmc = mcmc.rec[,tail(which.recession,1)+4]
      a3.summary = summary.rec[,tail(which.recession,1)+4]
      a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a3.mcmc), 
                               stdev =std(a3.mcmc), maxpost = a3.summary[16])))
      names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      # write.table(a3.df, file =paste0(dir.segm.recessions.model.param,"/df.a3.txt"),
      #             sep="\t", row.names=FALSE, col.names = TRUE)
      #b3:
      b3.mcmc = mcmc.rec[,tail(which.recession,1)+5]
      b3.summary = summary.rec[,tail(which.recession,1)+5]
      b3.df = data.frame(  t(c(quantile(b3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b3.mcmc), 
                               stdev =std(b3.mcmc), maxpost = b3.summary[16])))
      names(b3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #a4: asymptote !!!
      a4.mcmc = mcmc.rec[, tail(which.recession,1) + 5 + which.recession]
      a4.summary = summary.rec[, tail(which.recession,1) + 5 + which.recession]
      a4.df = data.frame(  t(c(quantile(a4.mcmc[,1], p = c(0.025, 0.5, 0.975)), 
                               mean=mean(a4.mcmc[,1]), 
                               stdev =std(a4.mcmc[,1]), 
                               maxpost = a4.summary[16,1])))
      names(a4.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a4.df = rbind(a4.df, t(c(quantile(a4.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a4.mcmc[,i]), 
                                 stdev =std(a4.mcmc[,i]), maxpost = a4.summary[16,i]) ))
      }    
      a4.df[, c(1,2,3,4,6)] = a4.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a4.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a4.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      # plots:
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
             theme_bw(base_size = 20) + 
               geom_point()+
               geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
               xlab("Time [day]") + ylab("a1") +
               theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
                      ,panel.grid.major  = element_blank()
                      ,panel.grid.minor  = element_blank())
      # p.a2 = ggplot(data = a2.df.time, aes(x=t, y=mean))+ 
      #   theme_bw() + 
      #   geom_point()+
      #   geom_errorbar(aes(x=a2.df.time$t, ymin=a2.df.time$`2.5%`, ymax=a2.df.time$`97.5%`))
      # 
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1, 
                            p.asympt, 
                            labels = c('alpha_1', 'beta'), 
                            label_size = 25,
                            label_fontface = "bold",
                            ncol = 1, nrow=2)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()
      
      
      } else if (rec.model =="3expWithAsympt_bis"){
      #################################################################################################     
        # read single parameters results :   a1(k) , b1, a2(k), b2, a3, b3, a4(k) 
        # a1:
        a1.mcmc = mcmc.rec[,which.recession]
        a1.summary = summary.rec[,which.recession]
        a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                                 mean=mean(a1.mcmc[,1]), 
                                 stdev =std(a1.mcmc[,1]), 
                                 maxpost = a1.summary[16,1])))
        names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        for (i in which.recession[-1]){
          a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                   stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
        }
        a1.df.time = cbind(a1.df ,   t =  time.rec)
        write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                    sep="\t", row.names=FALSE, col.names = TRUE)
        #b1:
        b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
        b1.summary = summary.rec[,tail(which.recession,1)+1]
        b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                                 stdev =std(b1.mcmc), maxpost = b1.summary[16])))
        names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        
        #a2:
        a2.mcmc = mcmc.rec[,tail(which.recession,1) +1 + which.recession]
        a2.summary = summary.rec[,tail(which.recession,1) +1 + which.recession]
        a2.df = data.frame(  t(c(quantile(a2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                                 mean  = mean(a2.mcmc[,1]), 
                                 stdev = std(a2.mcmc[,1]), 
                                 maxpost = a2.summary[16,1])))
        names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        for (i in which.recession[-1]){
          a2.df = rbind(a2.df, t(c(quantile(a2.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                                   mean    = mean(a2.mcmc[,i]), 
                                   stdev   = std(a2.mcmc[,i]), 
                                   maxpost = a2.summary[16,i]) ))
        }
        a2.df.time = cbind(a2.df ,   t =  time.rec)
        write.table(a2.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
                     sep="\t", row.names=FALSE, col.names = TRUE)
        
        #b2:
        b2.mcmc = mcmc.rec[, 2*tail(which.recession,1)+2]
        b2.summary = summary.rec[, 2*tail(which.recession,1)+2]
        b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b2.mcmc), 
                                 stdev =std(b2.mcmc), maxpost = b2.summary[16])))
        names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        
        
        #a3:
        a3.mcmc = mcmc.rec[, 2*tail(which.recession,1)+3]
        a3.summary = summary.rec[, 2*tail(which.recession,1)+2]
        a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a3.mcmc), 
                                 stdev =std(a3.mcmc), maxpost = a3.summary[16])))
        names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        # write.table(a3.df, file =paste0(dir.segm.recessions.model.param,"/df.a3.txt"),
        #             sep="\t", row.names=FALSE, col.names = TRUE)
        
        #b3:
        b3.mcmc = mcmc.rec[, 2*tail(which.recession,1)+4]
        b3.summary = summary.rec[, 2*tail(which.recession,1)+4]
        b3.df = data.frame(  t(c(quantile(b3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b3.mcmc), 
                                 stdev =std(b3.mcmc), maxpost = b3.summary[16])))
        names(b3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        
        #a4: asymptote !!!
        a4.mcmc = mcmc.rec[, 2*tail(which.recession,1)+4 + which.recession]
        a4.summary = summary.rec[, 2*tail(which.recession,1)+4 + which.recession]
        a4.df = data.frame(  t(c(quantile(a4.mcmc[,1], p = c(0.025, 0.5, 0.975)), 
                                 mean=mean(a4.mcmc[,1]), 
                                 stdev =std(a4.mcmc[,1]), 
                                 maxpost = a4.summary[16,1])))
        names(a4.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        for (i in which.recession[-1]){
          a4.df = rbind(a4.df, t(c(quantile(a4.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a4.mcmc[,i]), 
                                   stdev =std(a4.mcmc[,i]), maxpost = a4.summary[16,i]) ))
        }    
        a4.df[, c(1,2,3,4,6)] = a4.df[,c(1,2,3,4,6)]  -  stage.scale.shift
        h.infinity = cbind(a4.df ,   t =  time.rec)
        write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a4.txt"),
                    sep="\t", row.names=FALSE, col.names = TRUE)
        
        p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
          theme_bw(base_size = 20) + 
          geom_point()+
          geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
          xlab("Time [day]") + ylab("a1") +
          theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
                 ,panel.grid.major  = element_blank()
                 ,panel.grid.minor  = element_blank())
        
        p.a2 = ggplot(data = a2.df.time, aes(x=t, y=mean))+
          theme_bw(base_size = 20) +
          geom_point()+
          geom_errorbar(aes(x=a2.df.time$t, ymin=a2.df.time$`2.5%`, ymax=a2.df.time$`97.5%`))+
          xlab("Time [day]") + ylab("a2") +
          theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
                 ,panel.grid.major  = element_blank()
                 ,panel.grid.minor  = element_blank())

        p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
          theme_bw(base_size = 20) + 
          geom_point()+
          geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
          xlab("Time [day]") + ylab("Asymptotic stage") +
          theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
                 ,panel.grid.major  = element_blank()
                 ,panel.grid.minor  = element_blank())
        
        bind.plt =  plot_grid(p.a1, 
                              p.a2,
                              p.asympt, 
                              labels = c('alpha_1', 'alpha_2', 'beta'), 
                              label_size = 25,
                              label_fontface = "bold",
                              ncol = 1, nrow=3)  
        pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
        print(bind.plt)
        dev.off()
        
      

      } else if (rec.model =="2expWithAsympt_rel"){
        ################################################################################################# 
        # read single parameters results : a1(k) , b1, a2, b2, a3(k) 
        # a1:
        a1.mcmc = mcmc.rec[,which.recession]
        a1.summary = summary.rec[,which.recession]
        a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                                 mean=mean(a1.mcmc[,1]), 
                                 stdev =std(a1.mcmc[,1]), 
                                 maxpost = a1.summary[16,1])))
        names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        for (i in which.recession[-1]){
          a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                   stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
        }
        a1.df.time = cbind(a1.df ,   t =  time.rec)
        write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                    sep="\t", row.names=FALSE, col.names = TRUE)
        #b1:
        b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
        b1.summary = summary.rec[,tail(which.recession,1)+1]
        b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                                 stdev =std(b1.mcmc), maxpost = b1.summary[16])))
        names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        
        #a2:
        a2.mcmc    = mcmc.rec[,tail(which.recession,1)+2]
        a2.summary = summary.rec[,tail(which.recession,1)+2]
        a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), 
                                 mean=mean(a2.mcmc), 
                                 stdev =std(a2.mcmc),
                                 maxpost = a2.summary[16])))
        names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        
        #b2:
        b2.mcmc = mcmc.rec[, tail(which.recession,1) + 3]
        b2.summary = summary.rec[, tail(which.recession,1) + 3]
        b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b2.mcmc), 
                                 stdev =std(b2.mcmc), maxpost = b2.summary[16])))
        names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        
        #a3:asymptote !!!
        a3.mcmc = mcmc.rec[, tail(which.recession,1) + 3 + which.recession]
        a3.summary = summary.rec[, tail(which.recession,1) + 3 + which.recession]
        a3.df = data.frame(  t(c(quantile(a3.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                                 mean=mean(a3.mcmc[,1]), 
                                 stdev =std(a3.mcmc[,1]), 
                                 maxpost = a3.summary[16,1])))
        names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
        for (i in which.recession[-1]){
          a3.df = rbind(a3.df, t(c(quantile(a3.mcmc[,i], p = c(0.025, 0.5, 0.975)),
                                   mean=mean(a3.mcmc[,i]), 
                                   stdev =std(a3.mcmc[,i]), 
                                   maxpost = a3.summary[16,i]) ))
        }
        a3.df[, c(1,2,3,4,6)] = a3.df[,c(1,2,3,4,6)]  -  stage.scale.shift
        h.infinity = cbind(a3.df ,   t =  time.rec)
        write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a3.txt"),
                    sep="\t", row.names=FALSE, col.names = TRUE)
        
        p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
          theme_bw(base_size = 20) + 
          geom_point()+
          geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
          xlab("Time [day]") + ylab("a1") +
          theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
                 ,panel.grid.major  = element_blank()
                 ,panel.grid.minor  = element_blank())
        
        p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
          theme_bw(base_size = 20) + 
          geom_point()+
          geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
          xlab("Time [day]") + ylab("Asymptotic stage") +
          theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
                 ,panel.grid.major  = element_blank()
                 ,panel.grid.minor  = element_blank())
        
        bind.plt =  plot_grid(p.a1,
                              p.asympt, 
                              labels = c('alpha_1', 'beta'), 
                              label_size = 25,
                              label_fontface = "bold",
                              ncol = 1, nrow=2)  
        pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
        print(bind.plt)
        dev.off()
        
        
        
        
      
    } else if (rec.model =="2expWithAsympt_bis"){
      ################################################################################################# 
      # read single parameters results : a1(k) , b1, a2(k), b2, a3(k) 
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
      b1.summary = summary.rec[,tail(which.recession,1)+1]
      b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                               stdev =std(b1.mcmc), maxpost = b1.summary[16])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #a2:
      a2.mcmc = mcmc.rec[,tail(which.recession,1) +1 + which.recession]
      a2.summary = summary.rec[,tail(which.recession,1) +1 + which.recession]
      a2.df = data.frame(  t(c(quantile(a2.mcmc[,1], 
                                        p = c(0.025, 0.5, 0.975)),
                               mean=mean(a2.mcmc[,1]), 
                               stdev =std(a2.mcmc[,1]), 
                               maxpost = a2.summary[16,1])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a2.df = rbind(a2.df, t(c(quantile(a2.mcmc[,i], p = c(0.025, 0.5, 0.975)),
                                 mean=mean(a2.mcmc[,i]), 
                                 stdev =std(a2.mcmc[,i]), 
                                 maxpost = a2.summary[16,i]) ))
      }
      a2.df.time = cbind(a2.df ,   t =  time.rec)
      write.table(a2.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #b2:
      b2.mcmc = mcmc.rec[,2*tail(which.recession,1)+2]
      b2.summary = summary.rec[,2*tail(which.recession,1)+2]
      b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b2.mcmc), 
                               stdev =std(b2.mcmc), maxpost = b2.summary[16])))
      names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #a3:asymptote !!!
      a3.mcmc = mcmc.rec[, 2*tail(which.recession,1)+2 + which.recession]
      a3.summary = summary.rec[, 2*tail(which.recession,1)+2 + which.recession]
      a3.df = data.frame(  t(c(quantile(a3.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a3.mcmc[,1]), 
                               stdev =std(a3.mcmc[,1]), 
                               maxpost = a3.summary[16,1])))
      names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a3.df = rbind(a3.df, t(c(quantile(a3.mcmc[,i], p = c(0.025, 0.5, 0.975)),
                                 mean=mean(a3.mcmc[,i]), 
                                 stdev =std(a3.mcmc[,i]), 
                                 maxpost = a3.summary[16,i]) ))
      }
      a3.df[, c(1,2,3,4,6)] = a3.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a3.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a3.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a1") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())

      p.a2 = ggplot(data = a2.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a2.df.time$t, ymin=a2.df.time$`2.5%`, ymax=a2.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a2") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1,
                            p.a2,
                            p.asympt, 
                            labels = c('alpha_1', "alpha_2", 'beta'), 
                            label_size = 25,
                            label_fontface = "bold",
                            ncol = 1, nrow=3)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()
      
      
      
      
      
      
    } else if (rec.model =="2expWithAsympt"){
      ################################################################################################# 
      # read single parameters results : a1(k) , b1, a2, b2, a3(k) 
      # a1:
      a1.mcmc    = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df      = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
      b1.summary = summary.rec[,tail(which.recession,1)+1]
      b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                               stdev =std(b1.mcmc), maxpost = b1.summary[16])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #a2:
      a2.mcmc    = mcmc.rec[,tail(which.recession,1)+2]
      a2.summary = summary.rec[,tail(which.recession,1)+2]
      a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), 
                               mean=mean(a2.mcmc), 
                               stdev =std(a2.mcmc),
                               maxpost = a2.summary[16])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #b2:
      b2.mcmc = mcmc.rec[, tail(which.recession,1) + 3]
      b2.summary = summary.rec[, tail(which.recession,1) + 3]
      b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b2.mcmc), 
                               stdev =std(b2.mcmc), maxpost = b2.summary[16])))
      names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #a3:asymptote !!!
      a3.mcmc = mcmc.rec[, tail(which.recession,1) + 3 + which.recession]
      a3.summary = summary.rec[, tail(which.recession,1) + 3 + which.recession]
      a3.df = data.frame(  t(c(quantile(a3.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a3.mcmc[,1]), 
                               stdev =std(a3.mcmc[,1]), 
                               maxpost = a3.summary[16,1])))
      names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a3.df = rbind(a3.df, t(c(quantile(a3.mcmc[,i], p = c(0.025, 0.5, 0.975)),
                                 mean=mean(a3.mcmc[,i]), 
                                 stdev =std(a3.mcmc[,i]), 
                                 maxpost = a3.summary[16,i]) ))
      }
      a3.df[, c(1,2,3,4,6)] = a3.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a3.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a3.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a1")+
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1,
                            p.asympt, 
                            labels = c("alpha_1", "beta"), 
                            label_size = 25,
                            label_fontface = "bold",
                            ncol = 1, nrow=2)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()
      
      
      
      
      
    } else if (rec.model =="1expWithAsympt"){
      ################################################################################################# 
      # read single parameters results :  a1(k) , b1,  a2(k) 
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
      b1.summary = summary.rec[,tail(which.recession,1)+1]
      b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                               stdev =std(b1.mcmc), maxpost = b1.summary[16])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #a3:asymptote !!!
      a2.mcmc = mcmc.rec[, tail(which.recession,1) + 1 + which.recession]
      a2.summary = summary.rec[, tail(which.recession,1) + 1 + which.recession]
      a2.df = data.frame(  t(c(quantile(a2.mcmc[,1], p = c(0.025, 0.5, 0.975)), 
                               mean=mean(a2.mcmc[,1]), 
                               stdev =std(a2.mcmc[,1]), 
                               maxpost = a2.summary[16,1])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a2.df = rbind(a2.df, t(c(quantile(a2.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                                 mean= mean(a2.mcmc[,i]), 
                                 stdev = std(a2.mcmc[,i]), 
                                 maxpost = a2.summary[16,i]) ))
      }
      a2.df[, c(1,2,3,4,6)] = a2.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a2.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a1") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1, 
                            p.asympt, 
                            labels = c('alpha_1', 'beta'), 
                            label_size = 25,
                            label_fontface = "plain",
                            ncol = 1, nrow=2)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()  
      
      
      
      
      
      
    } else if (rec.model =="expexp"){
      #################################################################################################
      # read single parameters results : a1(k) , b1, n1,  a2(k) 
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
      b1.summary = summary.rec[,tail(which.recession,1)+1]
      b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                               stdev =std(b1.mcmc), maxpost = b1.summary[16])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      write.table(b1.df, file =paste0(dir.segm.recessions.model.param,"/df.b1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #n1:
      n1.mcmc = mcmc.rec[,tail(which.recession,1) +2]
      n1.summary = summary.rec[,tail(which.recession,1) +2]
      n1.df = data.frame(  t(c(quantile(n1.mcmc, p = c(0.025, 0.5, 0.975)),
                               mean=mean(n1.mcmc), 
                               stdev =std(n1.mcmc), 
                               maxpost = n1.summary[16])))
      names(n1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      write.table(n1.df, file =paste0(dir.segm.recessions.model.param,"/df.n1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #a2:asymptote !!!
      a2.mcmc = mcmc.rec[, tail(which.recession,1) + 2 + which.recession]
      a2.summary = summary.rec[, tail(which.recession,1) + 2 + which.recession]
      a2.df = data.frame(  t(c(quantile(a2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a2.mcmc[,1]), 
                               stdev =std(a2.mcmc[,1]), 
                               maxpost = a2.summary[16,1])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a2.df = rbind(a2.df, t(c(quantile(a2.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc[,i]), 
                                 stdev =std(a2.mcmc[,i]), maxpost = a2.summary[16,i]) ))
      }
      a2.df[, c(1,2,3,4,6)] = a2.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a2.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a1") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1, 
                            p.asympt, 
                            labels = c('alpha_1', 'beta'), 
                            label_size = 25,
                            label_fontface = "plain",
                            ncol = 1, nrow=2)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()
      
      
      
      
      
    } else if ((rec.model =="expexp_bis") | (rec.model =="hyperb_bis")| (rec.model =="Coutagne_bis")){
      #################################################################################################
      # read single parameters results :   a1(k) , b1(k), n1,  a2(k) 
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1) + which.recession]
      b1.summary = summary.rec[,tail(which.recession,1) + which.recession]
      b1.df = data.frame(  t(c(quantile(b1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean    = mean(b1.mcmc[,1]), 
                               stdev   = std(b1.mcmc[,1]), 
                               maxpost = b1.summary[16,1])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        b1.df = rbind(b1.df, t(c(quantile(b1.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                                 mean  = mean(b1.mcmc[,i]), 
                                 stdev = std(b1.mcmc[,i]), 
                                 maxpost = b1.summary[16,i]) ))
      }
      b1.df.time = cbind(b1.df ,   t =  time.rec)
      write.table(b1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.b1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      #n1:
      n1.mcmc = mcmc.rec[, 2*tail(which.recession,1) +1]
      n1.summary = summary.rec[, 2*tail(which.recession,1) +1]
      n1.df = data.frame(  t(c(quantile(n1.mcmc, p = c(0.025, 0.5, 0.975)),
                               mean=mean(n1.mcmc), 
                               stdev =std(n1.mcmc), 
                               maxpost = n1.summary[16])))
      names(n1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      write.table(n1.df, file =paste0(dir.segm.recessions.model.param,"/df.n1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #a2:asymptote !!!
      a2.mcmc = mcmc.rec[, 2*tail(which.recession,1) + 1 + which.recession]
      a2.summary = summary.rec[, 2*tail(which.recession,1) + 1 + which.recession]
      a2.df = data.frame(  t(c(quantile(a2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a2.mcmc[,1]), 
                               stdev =std(a2.mcmc[,1]), 
                               maxpost = a2.summary[16,1])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a2.df = rbind(a2.df, t(c(quantile(a2.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc[,i]), 
                                 stdev =std(a2.mcmc[,i]), maxpost = a2.summary[16,i]) ))
      }
      a2.df[, c(1,2,3,4,6)] = a2.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a2.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a1") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.b1 = ggplot(data = b1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=b1.df.time$t, ymin=b1.df.time$`2.5%`, ymax=b1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("b1") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1,
                            p.b1,
                            p.asympt, 
                            labels = c('alpha_1', 'lambda1', 'beta'), 
                            label_size = 25,
                            label_fontface = "bold",
                            ncol = 1, nrow=3)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()
      
    
      
      
    } else if (rec.model =="hyperb"){
      ################################################################################################# 
      # read single parameters results : a1(k) , b1, n1,  a2(k) 
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
      b1.summary = summary.rec[,tail(which.recession,1)+1]
      b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                               stdev =std(b1.mcmc), maxpost = b1.summary[16])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      write.table(b1.df, file =paste0(dir.segm.recessions.model.param,"/df.b1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #n1:
      n1.mcmc = mcmc.rec[,tail(which.recession,1) +2]
      n1.summary = summary.rec[,tail(which.recession,1) +2]
      n1.df = data.frame(  t(c(quantile(n1.mcmc, p = c(0.025, 0.5, 0.975)),
                               mean=mean(n1.mcmc), 
                               stdev =std(n1.mcmc), 
                               maxpost = n1.summary[16])))
      names(n1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      write.table(n1.df, file =paste0(dir.segm.recessions.model.param,"/df.n1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #a2:asymptote !!!
      a2.mcmc = mcmc.rec[, tail(which.recession,1) + 2 + which.recession]
      a2.summary = summary.rec[, tail(which.recession,1) + 2 + which.recession]
      a2.df = data.frame(  t(c(quantile(a2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a2.mcmc[,1]), 
                               stdev =std(a2.mcmc[,1]), 
                               maxpost = a2.summary[16,1])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a2.df = rbind(a2.df, t(c(quantile(a2.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc[,i]), 
                                 stdev =std(a2.mcmc[,i]), maxpost = a2.summary[16,i]) ))
      }
      a2.df[, c(1,2,3,4,6)] = a2.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a2.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a1") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1, 
                            p.asympt, 
                            labels = c('alpha_1', 'beta'), 
                            label_size = 25,
                            label_fontface = "bold",
                            ncol = 1, nrow=2)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()
      
      
      
      
    } else if (rec.model =="Coutagne"){
      ################################################################################################# 
      # read single parameters results : a1(k) , b1, n1 , a2(k)
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      a1.df.time = cbind(a1.df ,   t =  time.rec)
      write.table(a1.df.time, file =paste0(dir.segm.recessions.model.param,"/df.a1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+1]
      b1.summary = summary.rec[,tail(which.recession,1)+1]
      b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                               stdev =std(b1.mcmc), maxpost = b1.summary[16])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      write.table(b1.df, file =paste0(dir.segm.recessions.model.param,"/df.b1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #n1:
      n1.mcmc = mcmc.rec[,tail(which.recession,1) +2]
      n1.summary = summary.rec[,tail(which.recession,1) +2]
      n1.df = data.frame(  t(c(quantile(n1.mcmc, p = c(0.025, 0.5, 0.975)),
                               mean=mean(n1.mcmc), 
                               stdev =std(n1.mcmc), 
                               maxpost = n1.summary[16])))
      names(n1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      write.table(n1.df, file =paste0(dir.segm.recessions.model.param,"/df.n1.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      
      #a2:asymptote !!!
      a2.mcmc = mcmc.rec[, tail(which.recession,1) + 2 + which.recession]
      a2.summary = summary.rec[, tail(which.recession,1) + 2 + which.recession]
      a2.df = data.frame(  t(c(quantile(a2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a2.mcmc[,1]), 
                               stdev =std(a2.mcmc[,1]), 
                               maxpost = a2.summary[16,1])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession[-1]){
        a2.df = rbind(a2.df, t(c(quantile(a2.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc[,i]), 
                                 stdev =std(a2.mcmc[,i]), maxpost = a2.summary[16,i]) ))
      }
      a2.df[, c(1,2,3,4,6)] = a2.df[,c(1,2,3,4,6)]  -  stage.scale.shift
      h.infinity = cbind(a2.df ,   t =  time.rec)
      write.table(h.infinity, file =paste0(dir.segm.recessions.model.param,"/df.a2.txt"),
                  sep="\t", row.names=FALSE, col.names = TRUE)
      
      p.a1 = ggplot(data = a1.df.time, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=a1.df.time$t, ymin=a1.df.time$`2.5%`, ymax=a1.df.time$`97.5%`))+
        xlab("Time [day]") + ylab("a1") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      p.asympt = ggplot(data = h.infinity, aes(x=t, y=mean))+ 
        theme_bw(base_size = 20) + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`)) +
        xlab("Time [day]") + ylab("Asymptotic stage") +
        theme( plot.margin       = unit(c(1.5, 0.5, 1.5, 1),"cm")
               ,panel.grid.major  = element_blank()
               ,panel.grid.minor  = element_blank())
      
      bind.plt =  plot_grid(p.a1, 
                            p.asympt, 
                            labels = c('alpha_1', 'beta'), 
                            label_size = 25,
                            label_fontface = "bold",
                            ncol = 1, nrow=2)  
      pdf(paste0(dir.rec.pool.test,"/Figure_parameters_timeseries.pdf"), 12, 18 ,useDingbats=F)
      print(bind.plt)
      dev.off()
    }
    ############################################################################################
    return(list(mcmc.rec      = mcmc.rec,
                summary.rec   = summary.rec,
                h.infinity    = h.infinity,
                df.h.infinity = h.infinity,
                pl            = p.asympt))
    ############################################################################################
  }
}































#############################################################################################################
bt.from.gaugings <- function(nperiods,
                             df.limni, 
                             dir.gaugings.results,   
                             t.shift.for.b, 
                             h_G, t_G, color_G,                                       
                             times.uncert         ) {                                              
#############################################################################################################
# Author: Matteo Darienzo  
# Objective: this function provides information on the river bed estimated by gaugings for each stable period
  
  #Reading and plotting b1 and b2:
  #*******************************
  # Boxplots of b1 and b2:
  SPD.summary = read.table(file=paste0(dir.gaugings.results, "/Results_Summary.txt"))
  SPD.mcmc.cooked = read.table(file=paste0(dir.gaugings.results, "/Results_MCMC_Cooked.txt"), header=TRUE)
  ###################################################
  # b1:
  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked[,1:nperiods]; 
  names(SPD.mcmc.cooked.b1) = seq(1,nperiods,1)
  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked.b1 %>% gather(period, b1, na.rm = FALSE, convert = FALSE)
  
  bt1.mcmc = SPD.mcmc.cooked[,1:nperiods]
  bt1.summary = SPD.summary[,1:nperiods]
  bt1.df = data.frame(  t(c(quantile(bt1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                            mean=    mean(bt1.mcmc[,1]), 
                            stdev =   std(bt1.mcmc[,1]), 
                            maxpost = bt1.summary[16,1])))
  names(bt1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
  for (i in 2:nperiods){
    bt1.df = rbind(bt1.df, t(c(quantile(bt1.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                               mean = mean(bt1.mcmc[,i]), 
                               stdev =std(bt1.mcmc[,i]), 
                               maxpost = bt1.summary[16,i]) ))
  }
  write.table(bt1.df, paste0(dir.gaugings.results,"/bt1_df.txt"), 
              sep ="\t", row.names=FALSE)
  ##################
  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked[,(nperiods+3 ):(2*nperiods+2)]; 
  names(SPD.mcmc.cooked.b2) = seq(1,nperiods,1)
  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked.b2 %>% gather(period, b2, na.rm = FALSE, convert = FALSE)
  
  bt2.mcmc = SPD.mcmc.cooked[,(nperiods+3 ):(2*nperiods+2)]; 
  bt2.summary = SPD.summary[,(nperiods+3 ):(2*nperiods+2)]; 
  bt2.df = data.frame(  t(c(quantile(bt2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                            mean =    mean(bt2.mcmc[,1]), 
                            stdev =   std(bt2.mcmc[,1]), 
                            maxpost = bt2.summary[16,1])))
  names(bt2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
  for (i in 2:nperiods){
    bt2.df = rbind(bt2.df, t(c(quantile(bt2.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                               mean = mean(bt2.mcmc[,i]), 
                               stdev =std(bt2.mcmc[,i]), 
                               maxpost = bt2.summary[16,i]) ))
  }
  write.table(bt2.df, paste0(dir.gaugings.results,"/bt2_df.txt"), 
              sep ="\t", row.names=FALSE)

  
  # Read shift times:
  shifts.all = t.shift.for.b$t.adj
  if (times.uncert ==TRUE) {   # /!\ this has to be changed !!!!!!!!!!!!!!
    t.shifts.before = c(0,t.shift.for.b$t.adj)
    t.shifts.plus = c(t.shift.for.b$t.adj, tail(df.limni$t_limni,1))
    bt2.df =cbind(bt2.df, 
                  t.shifts.before = t.shifts.before , 
                  t.shifts.plus=t.shifts.plus)
    
    bt1.df =cbind(bt1.df, 
                  t.shifts.before = t.shifts.before , 
                  t.shifts.plus=    t.shifts.plus)
  } else {
    t.shifts = shifts[-length(shifts)]
    t.shifts.before = c(0,t.shifts)
    t.shifts.plus = sort(c(t.shifts, tail(df.limni$t_limni,1)))
  }
  
  # return the parameter b statistics:
  return(list(bt1.df        = bt1.df, 
              bt2.df        = bt2.df,
              t.shift.for.b = t.shift.for.b))
}














































#############################################################################################################
read.results.regression.rec.long = function(dir.recess, dir.segm.recessions,  
                                            rec.model, 
                                            rec.selection,
                                            BayesianOption,
                                            which.recession,
                                            time.rec) {
############################################################################################################## 
  ##########################
  if (BayesianOption ==1) {
  ##########################
    summary.rec =NULL
    h.infinity=NULL; recess.mcmc = NULL; recess.env=NULL; recess.data=NULL;
    a1.df = NULL;   b1.df = NULL;   a1.new.prior =NULL;    b1.new.prior;
    a2.df = NULL;   b2.df = NULL;   a2.new.prior = NULL ;  b2.new.prior =NULL ;
    a3.df = NULL;   b3.df = NULL;   a3.new.prior = NULL;   b3.new.prior = NULL;
    
    for (r in which.recession) {
      recess.mcmc[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                r,"_long/BaM/Results_MCMC_cooked.txt"),header=TRUE)
      summary.rec[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                r,"_long/BaM/Results_Summary.txt"),header=TRUE)
      a1.df[[r]]       = c(quantile(recess.mcmc[[r]][,1], p = c(0.025, 0.5, 0.975)), 
                           mean =summary.rec[[r]]$a1[5],  
                           stdev =  summary.rec[[r]]$a1[11], 
                           MAP=  summary.rec[[r]]$a1[16], 
                           t =rec.selection[[2]][r])
      b1.df[[r]]       = c(quantile(recess.mcmc[[r]][,2], p = c(0.025, 0.5, 0.975)), 
                           mean =summary.rec[[r]]$b1[5],  
                           stdev =  summary.rec[[r]]$b1[11], 
                           MAP=  summary.rec[[r]]$b1[16], 
                           t =rec.selection[[2]][r])
      a2.df[[r]]       = c(quantile(recess.mcmc[[r]][,3], p = c(0.025, 0.5, 0.975)), 
                           mean =summary.rec[[r]]$a2[5],  
                           stdev =  summary.rec[[r]]$a2[11], 
                           MAP=  summary.rec[[r]]$a2[16], 
                           t =rec.selection[[2]][r])
      b2.df[[r]]       = c(quantile(recess.mcmc[[r]][,4], p = c(0.025, 0.5, 0.975)), 
                           mean =summary.rec[[r]]$b2[5],  
                           stdev =  summary.rec[[r]]$b2[11], 
                           MAP=  summary.rec[[r]]$b2[16], 
                           t =rec.selection[[2]][r])
      a3.df[[r]]       = c(quantile(recess.mcmc[[r]][,5], p = c(0.025, 0.5, 0.975)), 
                           mean =summary.rec[[r]]$a3[5],  
                           stdev =  summary.rec[[r]]$a3[11], 
                           MAP=  summary.rec[[r]]$a3[16], 
                           t =rec.selection[[2]][r])
      b3.df[[r]]       = c(quantile(recess.mcmc[[r]][,6], p = c(0.025, 0.5, 0.975)), 
                           mean =summary.rec[[r]]$b3[5],  
                           stdev =  summary.rec[[r]]$b3[11], 
                           MAP=  summary.rec[[r]]$b3[16], 
                           t =rec.selection[[2]][r])
      #*********************  
      h.infinity[[r]] = c(quantile(recess.mcmc[[r]][,7], p = c(0.025, 0.5, 0.975)), 
                          mean =summary.rec[[r]]$a4[5],  
                          stdev =  summary.rec[[r]]$a4[11], 
                          MAP=  summary.rec[[r]]$a4[16], 
                          t =rec.selection[[2]][r])
      recess.data[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                                r,"_long/BaM/Curves_Data.txt"),header=TRUE)
      recess.data[[r]] = cbind(recess.data[[r]], t= rep(r,length(recess.data[[r]][1]) ))
      recess.env[[r]] = read.table(file=paste0(dir.recess,"/curves_regression/C",
                                               r,"_long/BaM/ht_TotalU.env"),header=TRUE)
    }
    
    # Unlist into a dataframe:
    df.h.infinity = data.frame(t(sapply(h.infinity, function(x) x[1:max(lengths(h.infinity))])))
    write.table(df.h.infinity, file =paste0(dir.segm.recessions,"/df.h.infinity_long.txt"),
                sep="\t", row.names=FALSE, col.names = TRUE)
    #
    a1.new.prior = c(mean(data.frame(t(sapply(a1.df, function(x) x[1:max(lengths(a1.df))])))$mean),
                     # std(data.frame(t(sapply(b2.df, function(x) x[1:max(lengths(b2.df))])))$mean) +
                     #     sqrt(sum(data.frame(t(sapply(b2.df, function(x) (x[1:max(lengths(b2.df))])^2)))$stdev))) 
                     sqrt(sum(data.frame(t(sapply(a1.df, function(x) (x[1:max(lengths(a1.df))])^2)))$stdev))) 
    
    b1.new.prior = c(mean(data.frame(t(sapply(b1.df, function(x) x[1:max(lengths(b1.df))])))$mean),
                     # std(data.frame(t(sapply(b2.df, function(x) x[1:max(lengths(b2.df))])))$mean) +
                     #     sqrt(sum(data.frame(t(sapply(b2.df, function(x) (x[1:max(lengths(b2.df))])^2)))$stdev))) 
                     sqrt(sum(data.frame(t(sapply(b1.df, function(x) (x[1:max(lengths(b1.df))])^2)))$stdev))) 
    
    a2.new.prior = c(mean(data.frame(t(sapply(a2.df, function(x) x[1:max(lengths(a2.df))])))$mean),
                     # std(data.frame(t(sapply(b2.df, function(x) x[1:max(lengths(b2.df))])))$mean) +
                     #     sqrt(sum(data.frame(t(sapply(b2.df, function(x) (x[1:max(lengths(b2.df))])^2)))$stdev))) 
                     sqrt(sum(data.frame(t(sapply(a2.df, function(x) (x[1:max(lengths(a2.df))])^2)))$stdev))) 
    
    b2.new.prior = c(mean(data.frame(t(sapply(b2.df, function(x) x[1:max(lengths(b2.df))])))$mean),
                     # std(data.frame(t(sapply(b2.df, function(x) x[1:max(lengths(b2.df))])))$mean) +
                     #     sqrt(sum(data.frame(t(sapply(b2.df, function(x) (x[1:max(lengths(b2.df))])^2)))$stdev))) 
                     sqrt(sum(data.frame(t(sapply(b2.df, function(x) (x[1:max(lengths(b2.df))])^2)))$stdev))) 
    
    a3.new.prior = c(mean(data.frame(t(sapply(a3.df, function(x) x[1:max(lengths(a3.df))])))$mean),
                     # std(data.frame(t(sapply(b2.df, function(x) x[1:max(lengths(b2.df))])))$mean) +
                     #     sqrt(sum(data.frame(t(sapply(b2.df, function(x) (x[1:max(lengths(b2.df))])^2)))$stdev))) 
                     sqrt(sum(data.frame(t(sapply(a3.df, function(x) (x[1:max(lengths(a3.df))])^2)))$stdev))) 
    b3.new.prior = c(mean(data.frame(t(sapply(b3.df, function(x) x[1:max(lengths(b3.df))])))$mean),
                     # std(data.frame(t(sapply(b2.df, function(x) x[1:max(lengths(b2.df))])))$mean) +
                     #     sqrt(sum(data.frame(t(sapply(b2.df, function(x) (x[1:max(lengths(b2.df))])^2)))$stdev))) 
                     sqrt(sum(data.frame(t(sapply(b3.df, function(x) (x[1:max(lengths(b3.df))])^2)))$stdev))) 
    return(list(recess.mcmc   = recess.mcmc, 
                summary.rec   = summary.rec,
                h.infinity    = h.infinity,
                recess.data   = recess.data,
                recess.env    = recess.env,
                df.h.infinity = df.h.infinity,
                a1.df         = a1.df,
                b1.df         = b1.df, 
                a2.df         = a2.df,
                b2.df         = b2.df,
                a3.df         = a3.df,
                b3.df         = b3.df,
                #####################
                a1.new.prior = a1.new.prior,
                b1.new.prior = b1.new.prior,
                a2.new.prior = a2.new.prior,
                b2.new.prior = b2.new.prior,
                a3.new.prior = a3.new.prior,
                b3.new.prior = b3.new.prior))
     
    
    ###############################
  } else if (BayesianOption==2) {
    ###############################
    dir.regression = paste0(dir.recess, "/curves_regression")
    dir.rec.pool = paste0(dir.regression,"/Pooling")
    dir.rec.pool.test = paste0(dir.rec.pool,"/test",3)
    # read results BaM files :
    summary.rec = read.table(file=paste0(dir.rec.pool.test,"/Results_Summary.txt"), header=TRUE)
    mcmc.rec = read.table(file=paste0(dir.rec.pool.test,"/Results_MCMC_cooked.txt"), header=TRUE)
    residuals.rec = read.table(file=paste0(dir.rec.pool.test,"/Results_Residuals.txt"), header=TRUE)
    # read single parameters results :
    if (rec.model == "3expWithAsympt") {
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+which.recession]
      b1.summary = summary.rec[,tail(which.recession,1)+which.recession]
      b1.df = data.frame(  t(c(quantile(b1.mcmc[,1], p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc[,1]), 
                               stdev =std(b1.mcmc[,1]), maxpost = b1.summary[16,1])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession){
        b1.df = rbind(b1.df, t(c(quantile(b1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc[,i]), 
                                 stdev =std(b1.mcmc[,i]), maxpost = b1.summary[16,i]) ))
      }
      #a2:
      a2.mcmc = mcmc.rec[,tail(which.recession,1)*2 +1]
      a2.summary = summary.rec[,tail(which.recession,1)*2 +1]
      a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc), 
                               stdev =std(a2.mcmc), maxpost = a2.summary[16])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #b2:
      b2.mcmc = mcmc.rec[,tail(which.recession,1)*2 +2]
      b2.summary = summary.rec[,tail(which.recession,1)*2 +2]
      b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b2.mcmc), 
                               stdev =std(b2.mcmc), maxpost = b2.summary[16])))
      names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      
      #a3:
      a3.mcmc = mcmc.rec[,tail(which.recession,1)*2+3]
      a3.summary = summary.rec[,tail(which.recession,1)*2+3]
      a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a3.mcmc), 
                               stdev =std(a3.mcmc), maxpost = a3.summary[16])))
      names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      
      #b3:
      b3.mcmc = mcmc.rec[,tail(which.recession,1)*2+4]
      b3.summary = summary.rec[,tail(which.recession,1)*2+4]
      b3.df = data.frame(  t(c(quantile(b3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b3.mcmc), 
                               stdev =std(b3.mcmc), maxpost = b3.summary[16])))
      names(b3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      
      #a4: asymptote !!!
      a4.mcmc = mcmc.rec[, tail(which.recession,1)*2 + 4+ which.recession]
      a4.summary = summary.rec[, tail(which.recession,1)*2 + 4 + which.recession]
      a4.df = data.frame(  t(c(quantile(a4.mcmc[,1], p = c(0.025, 0.5, 0.975)), mean=mean(a4.mcmc[,1]), 
                               stdev =std(a4.mcmc[,1]), maxpost = a4.summary[16,1])))
      names(a4.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in 2:tail(which.recession,1)){
        a4.df = rbind(a4.df, t(c(quantile(a4.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a4.mcmc[,i]), 
                                 stdev =std(a4.mcmc[,i]), maxpost = a4.summary[16,i]) ))
      }
      
      h.infinity = cbind(a4.df,  t =   time.rec[which.recession])
      
      pl = ggplot(data=h.infinity, aes(x=t, y=mean))+ theme_bw() + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`))+
        scale_y_continuous(limits = c(-150, 50))
      
    }
    
    return(list(mcmc.rec=mcmc.rec,
                summary.rec=summary.rec,
                h.infinity=h.infinity,
                df.h.infinity=h.infinity,
                pl=pl
                ))
  }
}





























































#############################################################################################################
read.results.regression.rec.multiple = function(dir.recess, dir.segm.recessions,  rec.model,
                                                which.group, data.to.cancel, time.rec) {
############################################################################################################## 
    dir.regression = paste0(dir.recess, "/curves_regression")
    dir.rec.pool = paste0(dir.regression,"/Pooling")
    dir.rec.pool.test = paste0(dir.rec.pool,"/test",3)
    nn = length(which.group) #number of datasets to merge
    summary.rec=NULL
    mcmc.rec=NULL
    residuals.rec=NULL
    h.infinity= NULL
    
    for (j in 1:nn){
    which.recession = which.group[[j]]
    # read results BaM files :
    summary.rec[[j]] = read.table(file=paste0(dir.rec.pool.test,"/Results_Summary.txt"), header=TRUE)
    mcmc.rec[[j]] = read.table(file=paste0(dir.rec.pool.test,"/Results_MCMC_cooked.txt"), header=TRUE)
    residuals.rec[[j]] = read.table(file=paste0(dir.rec.pool.test,"/Results_Residuals.txt"), header=TRUE)
    # read single parameters results :
    if (rec.model == "3expWithAsympt") {
      # a1:
      a1.mcmc = mcmc.rec[,which.recession]
      a1.summary = summary.rec[,which.recession]
      a1.df = data.frame(  t(c(quantile(a1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                               mean=mean(a1.mcmc[,1]), 
                               stdev =std(a1.mcmc[,1]), 
                               maxpost = a1.summary[16,1])))
      names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession){
        a1.df = rbind(a1.df, t(c(quantile(a1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a1.mcmc[,i]), 
                                 stdev =std(a1.mcmc[,i]), maxpost = a1.summary[16,i]) ))
      }
      #b1:
      b1.mcmc = mcmc.rec[,tail(which.recession,1)+which.recession]
      b1.summary = summary.rec[,tail(which.recession,1)+which.recession]
      b1.df = data.frame(  t(c(quantile(b1.mcmc[,1], p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc[,1]), 
                               stdev =std(b1.mcmc[,1]), maxpost = b1.summary[16,1])))
      names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in which.recession){
        b1.df = rbind(b1.df, t(c(quantile(b1.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc[,i]), 
                                 stdev =std(b1.mcmc[,i]), maxpost = b1.summary[16,i]) ))
      }
      #a2:
      a2.mcmc = mcmc.rec[,tail(which.recession,1)*2 +1]
      a2.summary = summary.rec[,tail(which.recession,1)*2 +1]
      a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc), 
                               stdev =std(a2.mcmc), maxpost = a2.summary[16])))
      names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      #b2:
      b2.mcmc = mcmc.rec[,tail(which.recession,1)*2 +2]
      b2.summary = summary.rec[,tail(which.recession,1)*2 +2]
      b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b2.mcmc), 
                               stdev =std(b2.mcmc), maxpost = b2.summary[16])))
      names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      
      #a3:
      a3.mcmc = mcmc.rec[,tail(which.recession,1)*2+3]
      a3.summary = summary.rec[,tail(which.recession,1)*2+3]
      a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a3.mcmc), 
                               stdev =std(a3.mcmc), maxpost = a3.summary[16])))
      names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      
      #b3:
      b3.mcmc = mcmc.rec[,tail(which.recession,1)*2+4]
      b3.summary = summary.rec[,tail(which.recession,1)*2+4]
      b3.df = data.frame(  t(c(quantile(b3.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b3.mcmc), 
                               stdev =std(b3.mcmc), maxpost = b3.summary[16])))
      names(b3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      
      
      #a4: asymptote !!!
      a4.mcmc = mcmc.rec[, tail(which.recession,1)*2 + 4+ which.recession]
      a4.summary = summary.rec[, tail(which.recession,1)*2 + 4 + which.recession]
      a4.df = data.frame(  t(c(quantile(a4.mcmc[,1], p = c(0.025, 0.5, 0.975)), mean=mean(a4.mcmc[,1]), 
                               stdev =std(a4.mcmc[,1]), maxpost = a4.summary[16,1])))
      names(a4.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
      for (i in 2:tail(which.recession,1)){
        a4.df = rbind(a4.df, t(c(quantile(a4.mcmc[,i], p = c(0.025, 0.5, 0.975)), mean=mean(a4.mcmc[,i]), 
                                 stdev =std(a4.mcmc[,i]), maxpost = a4.summary[16,i]) ))
      }
      h.infinity = rbind(h.infinity, a4.df)
      h.infinity = cbind(a4.df,  t =  time.rec[which.recession])
    } 
    h.infinity = h.infinity[-data.to.cancel, ]
    pl = ggplot(data=h.infinity, aes(x=t, y=mean))+ theme_bw() + 
        geom_point()+
        geom_errorbar(aes(x=h.infinity$t, ymin=h.infinity$`2.5%`, ymax=h.infinity$`97.5%`))+
        scale_y_continuous(limits = c(-150, 50))
    ggsave(pl, filename =paste0(dir.rec.pool,"/time_series_asymptote.png"), 
                device = "png", width = 16, height =13, dpi = 400, units = "in")
    
    return(list(mcmc.rec=mcmc.rec,
                summary.rec=summary.rec,
                h.infinity=h.infinity,
                df.h.infinity=h.infinity,
                pl=pl))
  }
}






























  
  
  











#########################################################################################################
realtime.rec   <-  function(            t_limni, h_limni, 
                                        rec.model, 
                                        data.recess.realtime,
                                        which.recession.realtime,
                                        dir.case_study, dir_code, dir_exe, 
                                        station.name, data.period,  
                                        prior.rec.realtime, 
                                        prior.gamma.realtime,
                                        ncycl, ncycles.max, nmcmc, nslim, jump.pos, jump.neg, nburn,  
                                        pred, 
                                        ncurves.h, curves.h, index.h, Qpeak.h, tpeak.h, Qf.h, 
                                        limits.y, limits.x,
                                        long.recessions,
                                        stage.scale.shift,
                                        starting.time) {
  ##########################################################################################################
  start_time <- Sys.time()
  message("REAL TIME Recession estimation - using BaM !!!"); flush.console()
  # function for the regression OF RECESSIONS (separately):
  dir.create(paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression"))
  dir.create(paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression/Real_Time"))
  dir.regression <- paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression/Real_Time")
  colfunc = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))   
  #Initialisation of files and plots:
  output_file_h = paste0(dir.regression,"/Param_rec_h.csv")
  cat(c("a1","b1","a2","b2","a3","b3","a4"), file = output_file_h, sep=";", append = FALSE)
  cat("", file = output_file_h, sep="\n", append = TRUE)
  reg.plot <-  ggplot() + 
               theme_bw(base_size = 15)+
               scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
                         expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
               scale_y_continuous(name = "Stage h (cm)", limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
               theme(plot.title = element_text(hjust = 0.5),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank()) 
                     #axis.line = element_line(colour = "black"))            
              ggsave(reg.plot, filename =paste(dir.regression,"/regression_realtime.png", sep=""), 
                      bg = "transparent", width = 12, height =7, dpi = 200)

  #initialisation:
  asymptote.h = 0;    curve_good.h = 0; hpeakgood.h = 0; t.real.good.h =0; index.good.h = 0; 
  asym.h.maxpost = 0; asym.h.stdev= 0;  asym.h.Q10= 0;   asym.h.Q90 = 0; asym.h.mean = 0;
  asym.h.Q2 = 0;      asym.h.Q95 =0;    Resultss =NULL;  quantiles.rec = NULL;
  theta1.maxpost = 0; theta1.stdev= 0; theta1.Q10= 0; theta1.Q90 = 0; theta1.mean = 0;
  theta2.maxpost = 0; theta2.stdev= 0; theta2.Q10= 0; theta2.Q90 = 0; theta2.mean = 0;
  theta3.maxpost = 0; theta3.stdev= 0; theta3.Q10= 0; theta3.Q90 = 0; theta3.mean = 0;
  theta4.maxpost = 0; theta4.stdev= 0; theta4.Q10= 0; theta4.Q90 = 0; theta4.mean = 0;
  theta5.maxpost = 0; theta5.stdev= 0; theta5.Q10= 0; theta5.Q90 = 0; theta5.mean = 0;
  theta6.maxpost = 0; theta6.stdev= 0; theta6.Q10= 0; theta6.Q90 = 0; theta6.mean = 0;
  results.regression = 0; asymptote.df <- data.frame(NULL); asym.df.temp =NULL

  #Exponential regression of the real time recession:
    d.h = NULL
    starting.time  = starting.time #days
    tot.length.rec = length(data.recess.realtime[[1]][[which.recession.realtime]][,1])
    final.time = tail(data.recess.realtime[[1]][[which.recession.realtime]][,1], 1)
    print(paste0("tot length of the recess=", tot.length.rec, " final time =",final.time ))
    for (icurve in (1:( tot.length.rec - starting.time))){      
      ind = icurve + starting.time
      d.h[[icurve]] = data.recess.realtime[[1]][[which.recession.realtime]][1:ind,]
      setwd(dir_exe)
      print(icurve)
      dir.create(paste0(dir.regression,"/C",icurve))
      dir.create(paste0(dir.regression,"/C",icurve,"/BaM"))
      dir.rec = paste0(dir.regression,"/C",icurve,"/BaM")
      
      dir.BaM.rec = paste0(dir_code,"/BaM_exe/Recession_h")  
      hpeakgood.h[icurve] = Qpeak.h[icurve]; 
      t.real.good.h[icurve] = tpeak.h[icurve];
      #tpeakgood.h[curve_good.h] = tpeak.h[i];
      index.good.h[icurve] = index.h[icurve];
      #nobs.h <- index.h[i] - 1 - index.h[i-1]
      nobs.h =length(d.h[[icurve]]$t)
      curve_data.h = "Recession_h/Curves_Data.txt"
      write.table(d.h[[icurve]], file=curve_data.h, append = FALSE, sep = "\t", eol = "\n", 
                  na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh"))
      conv=FALSE
      iterat=0
      ncycles = ncycl
      while( (conv == FALSE)& (ncycles <= ncycles.max)){
        iterat=iterat+1
        if (iterat>1){  print(c("increasing Ncycles of mcmc =",ncycles)) }
        BaM_config1.h(tgrid = seq(limits.x[1], limits.x[2], 0.5),
                      nobs  = nobs.h, 
                      ncycles = ncycles, nmcmc =nmcmc, nslim=nslim, jump.pos=jump.pos,  jump.neg=jump.neg,  nburn=nburn,
                      model = rec.model, 
                      pred = pred, 
                      asympt.before = curves.h$Qcurve[index.h[icurve]-1], 
                      remnant ="Linear",
                      prior=prior.rec.realtime)
        
        #setwd(Dir_exe)
        #BAM:
        if (rec.model == "2expWithAsympt") {
          system2("BaM_2exp_pool2.exe",stdout =NULL, stderr =NULL)
        } else if (rec.model == "1expWithAsymptNorm") {
          system2("BaM_normalised_1exp_OK.exe", stdout =FALSE)
        } else if (rec.model == "2expWithoutAsymptNorm") { 
          system2("BaM_normalised_2exp_without_asym.exe", stdout =FALSE)
        } else if (rec.model == "3expWithAsympt") {
          system2("BaM_3exp.exe",stdout =NULL, stderr =NULL)
        }
        #save:
        list.of.files.rec <- c(paste0(dir.BaM.rec,"/Results_MCMC_Cooked.txt"),
                               paste0(dir.BaM.rec,"/Results_Residuals.txt"),
                               paste0(dir.BaM.rec,"/Results_Summary.txt"),
                               paste0(dir.BaM.rec,"/Config_Model.txt"),
                               paste0(dir.BaM.rec,"/Curves_Data.txt"),
                               paste0(dir.BaM.rec,"/tgrid.txt"),
                               paste0(dir.BaM.rec,"/ht_Maxpost.spag"),
                               paste0(dir.BaM.rec,"/ht_ParamU.spag"),
                               paste0(dir.BaM.rec,"/ht_ParamU.env"),
                               paste0(dir.BaM.rec,"/ht_TotalU.env"),
                               paste0(dir.BaM.rec,"/ht_TotalU.spag"))
        for (i in 1:length(list.of.files.rec)) {
          file.copy(list.of.files.rec[i], dir.rec ,overwrite = TRUE)
        }
        plot.mcmc(workspace= dir.rec, iter= icurve)  #plot density and traceplots of MCMC
        conv = Convergence.test(dir.seg = dir.rec ,  npar=7 , dir.plot = dir.rec)
        ncycles=ncycles+1000
        setwd(dir_exe)
      }
      if (conv==FALSE) {
        print(paste0("recess ",icurve," Not converged !!! "))
      }
      
      #reading results of BaM:
      tgrid.file <- paste0(dir_exe,"/Recession_h/tgrid.txt")
      # env.tot.h <- EnvelopLayer(file=paste(dir_exe,"/Recession_h/ht_TotalU.env", sep=""), 
      #             color=colfunc(80)[curve_good.h], alpha = 0.3)
      env = read.table(paste0(dir.BaM.rec,"/ht_TotalU.env"), header =TRUE)
      env.param = read.table(paste0(dir.BaM.rec,"/ht_ParamU.env"), header =TRUE)
      RecCurve.h = read.table("Recession_h/Curves_Data.txt", header =TRUE)
      Results.h = read.table("Recession_h/Results_Summary.txt", header =TRUE)
      #results.regression[[curve_good.h]] <- Results.h
      mcmc.cooked = read.table(paste0(dir_exe,"/Recession_h/Results_MCMC_Cooked.txt"), header =TRUE)
      if (rec.model == "2expWithAsympt") {
        for (param in 1:7) {
          quantiles.rec[[param]] <- quantile(x = mcmc.cooked[,param], p = c(0.025, 0.1, 0.5, 0.9, 0.975)) 
        }
        quant = data.frame(a1 = quantiles.rec[[1]],
                           b1 = quantiles.rec[[2]],
                           a2 = quantiles.rec[[3]],
                           b2 = quantiles.rec[[4]],
                           a3 = quantiles.rec[[5]],
                           Y1_gamma1 = quantiles.rec[[6]], 
                           Y1_gamma2 = quantiles.rec[[7]])
      } else if (rec.model == "3expWithAsympt") { 
        for (param in 1:9) {
          quantiles.rec[[param]] <- quantile(x = mcmc.cooked[,param], p = c(0.025,0.1, 0.5, 0.9, 0.975)) 
        }
        quant = data.frame(a1 = quantiles.rec[[1]], 
                           b1 = quantiles.rec[[2]],
                           a2 = quantiles.rec[[3]],
                           b2 = quantiles.rec[[4]],
                           a3 = quantiles.rec[[5]],
                           b3 = quantiles.rec[[6]],
                           a4 = quantiles.rec[[7]],
                           Y1_gamma1 = quantiles.rec[[8]],
                           Y1_gamma2 = quantiles.rec[[9]])
      }
      #print(quantiles.rec)
      Resultss[[icurve]] = rbind(Results.h, quant)
      write.table(Resultss[[icurve]], file=paste0(dir.rec, "/results.txt"), append = FALSE, sep = "\t", eol = "\n", 
                  na = "NA", dec = ".", row.names = FALSE)
      h_bam = 0
      ttt = seq(limits.x[1],limits.x[2], 0.5)
      if (rec.model == "2expWithAsympt") {
        for (ii in 1:length(ttt)) {
          h_bam[ii] = Resultss[[icurve]][16,1]*exp(- Resultss[[icurve]][16,2]*ttt[ii]) + 
            Resultss[[icurve]][16,3]*exp(- Resultss[[icurve]][16,4]*ttt[ii]) + 
            Resultss[[icurve]][16,5]
        }
        boxplot.df <- data.frame(x=1,
                                 y0= min(mcmc.cooked[,5]), 
                                 y2= quantile(mcmc.cooked[,5],0.025), 
                                 y50 = median(mcmc.cooked[,5]), 
                                 y97 =quantile(mcmc.cooked[,5],0.975),
                                 y100 = max(mcmc.cooked[,5])) 
      } else if (rec.model == "1expWithAsymptNorm") {
        for (ii in 1:length(ttt)) {
          h_bam[ii] = (1- Maxpost.h[[1]])*exp(- Maxpost.h[[2]]*ttt[ii]) + Maxpost.h[[1]]
        }
      } else if (rec.model == "3expWithAsympt") { 
        for (ii in 1:length(ttt)) {   # MAXPOST:
          h_bam[ii] = Resultss[[icurve]][16,1]*exp(- Resultss[[icurve]][16,2]*ttt[ii])+ 
            Resultss[[icurve]][16,3]*exp(- Resultss[[icurve]][16,4]*ttt[ii])+ 
            Resultss[[icurve]][16,5]*exp(- Resultss[[icurve]][16,6]*ttt[ii])+ 
            Resultss[[icurve]][16,7]
        }
        # boxplot.df <- data.frame(x=1, y0= min(mcmc.cooked[,7]), 
        #                          y2= quantile(mcmc.cooked[,7],0.025),
        #                          y50 = median(mcmc.cooked[,7]), 
        #                          y97 =quantile(mcmc.cooked[,7],0.975),
        #                          y100 = max(mcmc.cooked[,7])) 
      }
      #***************************************************************************************************
      if (pred ==TRUE) {
        temp <- data.frame(x = RecCurve.h[,1] , 
                           y = RecCurve.h[,2] - stage.scale.shift, 
                           z = RecCurve.h[,3])
        temp2 <- data.frame(xx = ttt , 
                            yy = h_bam   - stage.scale.shift, 
                            zz = env[,2] - stage.scale.shift,
                            kk = env[,3] - stage.scale.shift)
        temp3 <- data.frame(xx = ttt , 
                            yy = h_bam         - stage.scale.shift, 
                            zz = env.param[,2] - stage.scale.shift, 
                            kk = env.param[,3] - stage.scale.shift)
        
        write.table(temp, file=paste0(dir.rec, "/temp.txt"), append = FALSE, sep = "\t", eol = "\n", 
                    na = "NA", dec = ".", row.names = FALSE, col.names=c("x", "y", "z"))
        write.table(temp2, file=paste0(dir.rec, "/temp2.txt"), append = FALSE, sep = "\t", eol = "\n", 
                    na = "NA", dec = ".", row.names = FALSE, col.names=c("xx", "yy", "zz", "kk"))
        write.table(temp3, file=paste0(dir.rec, "/temp3.txt"), append = FALSE, sep = "\t", eol = "\n", 
                    na = "NA", dec = ".", row.names = FALSE, col.names=c("xx", "yy", "zz", "kk"))
        #****************************************************************************************
        
        reg.plot =  reg.plot  +
          # geom_point(data = temp , aes(x = x, y= y), color=colfunc(80)[curve_good.h] , size= 0.5)+
          #geom_line(data = temp2, aes(x=xx, y=yy), color= "gray", size= 0.1) + #colfunc(Ncurves)[icurve], size = 0.1)+
          geom_point(data = temp , aes(x = x, y= y), color= "black", size= 2)+  #colfunc(Ncurves)[icurve], size= 2)+
          geom_ribbon(data = temp2, aes(x = xx, ymin = zz , ymax = kk), 
                      fill = "gray40", alpha=0.1) #colfunc(Ncurves)[icurve], alpha = 0.1)
        # geom_ribbon(data = temp3, aes(x = xx, ymin = zz , ymax = kk), 
        #             fill = colfunc(Ncurves)[icurve], alpha = 0.3)
        
        ggsave(reg.plot, filename =paste(dir.regression,"/regression_realtime.png", sep=""),
               bg = "transparent", width = 12, height =7, dpi = 200)
      }
      #**************************************************************************************
      asym.df.temp[[icurve]] = data.frame(  t = data.recess.realtime[[2]][icurve], 
                                            h = Resultss[[icurve]][16,7] - stage.scale.shift, 
                                            q2 = Resultss[[icurve]][17,7] - stage.scale.shift, 
                                            q97= Resultss[[icurve]][21,7] - stage.scale.shift) 
      write.table(asym.df.temp[[icurve]], 
                  file=paste0(dir.rec, "/df.asymptote.txt"), append = FALSE, sep = "\t", eol = "\n", 
                  na = "NA", dec = ".", row.names = FALSE, col.names=c("t", "h", "q2", "q97"))            
      asymptote.df = rbind(asymptote.df, asym.df.temp[[icurve]])
      # asympt.plot =  asympt.plot +
      #                geom_point(data = asym.df.temp[[icurve]], aes(x = t, y = h), colour = "black",size = 2) +
      #                #geom_boxplot(boxplot.df, aes(x = x, ymin=y0, lower =y2, middle= y50, 
      #                              #upper=y97, ymax= y100), stat = "identity") +
      #                             #geom_line(data = asym.df.temp[[i]], aes(x =t, y =h), colour = "gray",size = 1) +
      #                geom_errorbar(data = asym.df.temp[[icurve]], aes(x= t, ymin = q2, ymax =q97), width=50, size = 0.3)
      #                ggsave(asympt.plot, filename =paste0(dir.regression,"/asymptote.png"),
      #                       bg = "transparent", width = 12, height =6, dpi = 200)
      #***************************************************************************************
      setwd(dir_code)
      
    }

    # save:
    if (rec.model == "2expWithAsympt") {
      Data.segm.rec <- data.frame("t"=t.real.good.h, "Y"=unlist(lapply(Resultss, "[[", 16,5), use.names=FALSE), 
                                  "uY"= unlist(lapply(Resultss, "[[",11,5), use.names=FALSE))
    } else if (rec.model == "3expWithAsympt") {
      Data.segm.rec <- data.frame("t"=t.real.good.h, "Y"=unlist(lapply(Resultss, "[[", 16,7), use.names=FALSE), 
                                  "uY"= unlist(lapply(Resultss, "[[",11,7), use.names=FALSE))
    }
    end_time <- Sys.time()
    print(c("computat. time for regression of recessions=", end_time-start_time ))
    return(list(asymptote.df, Data.segm.rec))
}


































#########################################################################################################
recession.regression.pooling <- function(t_limni, h_limni, rec.model, dir.case_study, dir_code, dir_exe,
                                         station.name, data.period, 
                                         ncycles.pool,
                                         nmcmc.pool,
                                         jumps,
                                         slim,
                                         burn, minMoveRate, maxMoveRate,
                                         tgood, pred, ncurves.h,
                                         curves.h,  index.h ,   hpeak.h , tpeak.h , hf.h , tf.h,
                                         hmin,  hmax,  tmax,
                                         limits.y, limits.x, length.avg,
                                         theta.b.prior, sd.theta.b.prior,
                                         theta.a.prior,
                                         sd.theta.a.prior,
                                         gamma,
                                         stdev.gamma) {
  ##########################################################################################################
  start_time = Sys.time()
  # Regression OF RECESSIONS:
  dir.create(paste(dir.case_study,"/Results/segmentation_recessions/curves_regression", sep=""))
  dir.regression <- paste(dir.case_study,"/Results/segmentation_recessions/curves_regression", sep="")
  colfunc = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))   
  
  output_file_h = paste(dir.regression,"/Param_rec_h.csv", sep="")
  if (rec.model == "2expWithAsympt") {
          cat(c("a1","b1","a2","b2","a3"), file = output_file_h, sep=";", append = FALSE)
          cat("", file = output_file_h, sep="\n", append = TRUE)
          reg.plot <- ggplot() + 
            scale_x_continuous(name = "time (day)", limits =c(0, tmax), expand = c(0,0), breaks=seq(0,tmax,10))+
            scale_y_continuous(name = "Stage h (cm)", limits = c(hmin, hmax), expand = c(0,0)) +
            theme(plot.title = element_text(hjust = 0.5),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            theme_light(base_size = 20)
          ggsave(reg.plot, filename =paste(dir.regression,"/regression.png", sep=""), 
                 bg = "transparent", width = 8, height =4, dpi = 400)
  } else if(rec.model == "3expWithAsympt") {
          cat(c("a1","b1","a2","b2","a3","b3","a4"), file = output_file_h, sep=";", append = FALSE)
          cat("", file = output_file_h, sep="\n", append = TRUE)
          #first plot:
          reg.plot <- ggplot() + 
            scale_x_continuous(name = "time (day)", limits =c(0,60), expand = c(0,0), breaks=seq(0,60,5))+
            scale_y_continuous(name = "Stage h (cm)", limits = c(-100,200), expand = c(0,0)) +
            theme(plot.title = element_text(hjust = 0.5),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
          ggsave(reg.plot, filename =paste(dir.regression,"/regression.png", sep=""), 
                 bg = "transparent", width = 8, height =4, dpi = 400)
  } else if (rec.model == "1expWithAsymptNorm") {
          cat(c("a1","b1"), file = output_file_h, sep=";", append = FALSE)
          cat("", file = output_file_h, sep="\n", append = TRUE)
          reg.plot <- ggplot() + 
            scale_x_continuous(name = "time (day)", limits =c(0,40), expand = c(0,0), breaks=seq(0,40,5))+
            scale_y_continuous(name = "Stage h (cm)", limits = c(0,1), expand = c(0,0)) +
            theme(text=element_text(size=16,  family="Serif"))+
            theme_light(base_size = 10)
  } else if (rec.model == "2expWithAsymptNorm") { 
          cat(c("a1","b1","a2","b2"), file = output_file_h, sep=";", append = FALSE)
          cat("", file = output_file_h, sep="\n", append = TRUE)
          reg.plot <- ggplot() + 
            scale_x_continuous(name = "time (day)", limits =c(0,40), expand = c(0,0), breaks=seq(0,40,5))+
            scale_y_continuous(name = "Stage h (cm)", limits = c(0,1), expand = c(0,0)) +
            theme(plot.title = element_text(hjust = 0.5),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            theme_light(base_size = 10)
  }
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
  results.regression = 0;
  asymptote.h = 0; curve_good.h = 0; hpeakgood.h = 0; t.real.good.h =0; index.good.h = 0;
  asymptote.h[0] = 1 #first asymptote starting point
  d.h.pooling <- NULL
  
  
  
  # POOLING Recession h, all curves together
  #########################################################################################
  #directories:
  dir.create(paste(dir.regression,"/Pooling", sep=""))
  dir.rec.pool <- paste(dir.regression,"/Pooling", sep="")
  dir.create(paste(dir.rec.pool,"/test",2, sep=""))
  dir.rec.pool.test <- paste(dir.rec.pool,"/test",2, sep="")
  dir.BaM.rec.pool <- paste(dir_code,"/BaM_exe/Recession_h_pooling",sep="")
  #***********************************************************************
  asympt.plot =  ggplot()+
    theme_light(base_size = 15)+   
    scale_x_continuous(name="time [day]",expand = c(0,0),limits = c(0,tail(t_limni,1))) +
    scale_y_continuous(name="Asymptote h [cm]",
                       limits = c(-100,50),expand = c(0,2)) +
    xlab("Time [day]")+ ylab("Asymptote [cm]")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  ggsave(asympt.plot, filename =paste0(dir.rec.pool.test,"/asymptote.png"),
         bg = "transparent", width = 8, height =4, dpi = 200)
  
  
  #create the data file for BaM (all curves appended in same file):
  #*******************************************************************
  for (i in 2:ncurves.h) {
    d.h.pooling.new = data.frame(t=curves.h$tcurve[(index.h[i-1]+1):(index.h[i]-1)],
                                 h=curves.h$Qcurve[(index.h[i-1]+1):(index.h[i]-1)],
                                 uh = 0/100*curves.h$Qcurve[(index.h[i-1]+1):(index.h[i]-1)])
    if ((tf.h[i] >= tgood) &(length(d.h.pooling.new$t) >= 10) ) {
      setwd(dir_exe)
      curve_good.h = curve_good.h + 1
      # print(curve_good.h)
      dir.rec <- paste(dir.regression,"/C",curve_good.h,"/BaM", sep="")
      hpeakgood.h[curve_good.h] = hpeak.h[i];
      t.real.good.h[curve_good.h] = tpeak.h[i];
      index.good.h[curve_good.h] = index.h[i];
      d.h.pooling.new = data.frame(t =curves.h$tcurve[(index.h[i-1]+1):(index.h[i]-1)],
                                   h = curves.h$Qcurve[(index.h[i-1]+1):(index.h[i]-1)],
                                   uh = 0/100*curves.h$Qcurve[(index.h[i-1]+1):(index.h[i]-1)],
                                   Period = curve_good.h)
      
      #average the recession stage:
      d.h.pooling.avg =moving.average.h(stage.h=d.h.pooling.new$h, time.h=d.h.pooling.new$t , Period.h=curve_good.h, length.avg)
      d.h.pooling <- rbind(d.h.pooling, d.h.pooling.avg)
    }
  }
  nobs.pooling <- length(d.h.pooling$t)
  ncurves.pooling <- curve_good.h
  print(c("ncurves=",ncurves.pooling))
  curve_data = "Recession_h_pooling/Curves_Data.txt"
  write.table(d.h.pooling, file=curve_data, append = FALSE, sep = "\t", eol = "\n",
              na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh", "Period"))
  
  BaM_config.pooling(
    dir_code = dir_code , model = rec.model, nobs=nobs.pooling, 
    ncycles = ncycles.pool, nmcmc =  nmcmc.pool, pred = pred, ncurves = ncurves.pooling,
    jumps,
    slim,
    burn, minMoveRate, maxMoveRate,
    theta.b.prior, sd.theta.b.prior,
    theta.a.prior,
    sd.theta.a.prior,
    gamma,
    stdev.gamma)
  message("****************************************************************"); flush.console()
  message("Applying BaM Regression by pooling method !!!  Please, wait ... "); flush.console()
  message("****************************************************************"); flush.console()
  if (rec.model == "2expWithAsympt") {
    system2("BaM_2exp_pool.exe",stdout =NULL, stderr =NULL)
  } else if (rec.model == "3expWithAsympt")  {
    system2("BaM_3exp.exe",stdout =NULL, stderr =NULL)
  }
  
  
  #**********************
  #save and copy results:
  #**********************
  list.files.pool <- c( paste(dir.BaM.rec.pool,"/Results_MCMC_Cooked.txt", sep=""),
                        paste(dir.BaM.rec.pool,"/Results_Residuals.txt", sep=""),
                        paste(dir.BaM.rec.pool,"/Results_Summary.txt", sep=""),
                        paste(dir.BaM.rec.pool,"/Config_Model.txt", sep=""),
                        paste(dir.BaM.rec.pool,"/Curves_Data.txt", sep=""))
  for (i in 1:length(list.files.pool)) {
    file.copy(list.files.pool[i], dir.rec.pool.test ,overwrite = TRUE)
  }
  setwd(dir.rec.pool.test)
  Results_summary.pool <- read.table(paste(dir.rec.pool.test,"/Results_Summary.txt",sep=""), header =TRUE)
  Results_mcmc.pool <- read.table(paste(dir.rec.pool.test,"/Results_MCMC_Cooked.txt",sep=""), header =TRUE)
  Results_residuals.pool <- read.table(paste(dir.rec.pool.test,"/Results_Residuals.txt",sep=""), header =TRUE)
  
  
  
  #*****************
  # Multi plot MCMC:
  #*****************
  nparam = ncol(Results_mcmc.pool)
  multiplot.mcmc =MCMCplot(MCMCfile='Results_MCMC_Cooked.txt', doPar=T,doLogPost=T,doDPar=F, type="trace",
                           xlab='',ylab='', ncol=5, prior=NULL, burn=0, slim=1)
  ggsave(multiplot.mcmc, filename =paste(dir.rec.pool.test,"/mcmc_multiplot.png", sep=""),
         bg = "transparent", width = 15, height =20, dpi = 300)
  
  #*****************************************
  # Convergence Gelman and autocorrelations:
  #*****************************************
  # conv <- Convergence.test(dir.seg = dir.rec.pool.test , npar=nparam, dir.plot =dir.rec.pool.test )
  
  
  #***********
  # Residuals:
  #***********
  residuals.plot <- ggplot(Results_residuals.pool)+
    geom_point(aes(x=Results_residuals.pool$X1_obs,y=Results_residuals.pool$Y1_sim), color = "red", size = 3)+
    geom_point(aes(x=Results_residuals.pool$X1_obs,y=Results_residuals.pool$Y1_obs), color = "black")+
    theme_bw()
  ggsave(residuals.plot, filename =paste(dir.rec.pool.test,"/residuals.png", sep=""),
         bg = "transparent", width = 8, height =4, dpi = 200)
  
  ############################################### save:
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
    col.asympt = seq((3*ncurves.pooling+3  + 1), (3*ncurves.pooling+3 +  ncurves.pooling), 1)
    asymptote.h = Results_summary.pool[16,col.asympt]
    u.asymptote.h = Results_summary.pool[11,col.asympt]
    asymptote.df = data.frame(
      t = t.real.good.h, 
      h2 = quantile(Results_mcmc.pool[,col.asympt], 0.025),
      h97 = quantile(Results_mcmc.pool[,col.asympt], 0.975),
      hMAP = asymptote.h
    )
    Data.segm.rec <- data.frame("t"=t.real.good.h, 
                                "Y"=asymptote.h, 
                                "uY"= u.asymptote.h
    )
    #*****************************************************************************
  }
  ggplot(data=asymptote.df)+geom_point(aes(x=time, y=hMAP)) + 
    geom_errorbar(aes(x=time, ymin =h2, ymax =h97), width=.02, size = 0.2)
  end_time <- Sys.time()
  print(c("computat. time for regression of recessions=", end_time-start_time ))
  return(list(asymptote.df, Data.segm.rec))
}




















#####################################################################
moving.average.h = function(stage.h, time.h , Period.h, length.avg) {
#####################################################################  
  #df.averaged.h = movavg(stage.h, 10, type="t")
  n = 0; averaged.h =NULL; averaged.t =NULL; 
  for (obs in 1:length(stage.h)) {
    if ((obs %% length.avg) == 0) {
      n=n+1
      averaged.h[n] = mean(stage.h[(obs-length.avg):obs])
      averaged.t[n] = mean(time.h[(obs-length.avg):obs])
    }
  }
  df.averaged.h = data.frame(
    t =averaged.t,
    h = averaged.h,
    uh =rep(0,n),
    Period = Period.h
  )
  return(df.averaged.h)
}






































##########################################################################################################################
# SEGMENTATION Of RECESSIONS (IN PARTICULAR OF THE ASYMPTOTIC STAGE PARAMETER):
recession.segmentation <- function(dir.exe,
                                   dir.segment.rec,
                                   file.options.general,
                                   file.options.recess,
                                   stage.record,
                                   gaugings,
                                   initial.time,
                                   colors.period) {
##########################################################################################################################
  # Segmentation of recessions dynamic parameters:
  # read inputs and options for computation:
  source(file.options.general)
  source(file.options.recess)
  BayesianOption         = 2                                   # Type of Bayesian model (1 = simple, 2 = pooling). only 2 is available!
  plot.b.from.gaugings   = FALSE                               # plot the river bed estimation obtained by using gaugings
  x.name                 = "Time [days]"                       # x-axis label for plot of time series
  y.name                 = "Asymptotic stage [m]"              # y-axis label for plot of time series of asymptotic stage
  save.all.results       = TRUE                                # TRUE = save all segmentation computations
  plot.gamma.uncertainty = TRUE                                # plot structural uncertainty on the segmentation.
  segm_all_parameters    = FALSE
  stage.scale.shift      = 1000                               # is in [cm]: parameter used to shift stage values and avoid negative values: h = h + stage.scale.shift
  plot.recession.uncert  = TRUE                           # plot the total uncertainty ribbons of the recession curves.
  plot.dates             = FALSE                             # plot tha shift times in date format

  # to remove in future:
  limits.Y.alpha         = c(0, 1000, 200)                     # Limits for the alpha (initial stage) plot.
  limits.Y.lambda        = c(0, 40, 10)                        # Limits for the lambda (recession rate) plot.
  stage.limits           = c(-150, 50, 50)                     # limits for recession stage is in "cm" !!!
  limits.y               = stage.limits
  limits.x.recess        = c(0, 150, 30)                       # [c(min, max, step)] limits for the recession period in days. by default = c(0, 150, 30).
  limits.x               = limits.x.recess
  asymptote.limits       = c(-150, 50, 50)                # limits for beta (asymptotic stages) in "cm"  !!!
  

  # directories :
  dir.segment.rec.test1    = paste0(dir.segment.rec,"/",name.folder.results.recession)
  dir.config.segmentation  = paste0(dir.exe, "/Segmentation")
  dir.estim.recessions     = paste0(dir.segment.rec.test1,"/2_curves_estimation")
  dir.segm.recessions      = paste0(dir.segment.rec.test1,"/3_curves_segmentation")
  dir.create(paste0(dir.segment.rec,"/", name.folder.results.recession))
  dir.create( paste0(dir.segment.rec.test1,"/3_curves_segmentation"))
  
  seg.rec = NULL; read.res.rec = NULL; results.segment = list()

  
  # # Read results of gaugings if any:
  # # Read river bed b (or b1 and b2) estimated from gaugings:
  # read.bt.gaug = bt.from.gaugings(nperiods              = nperiod.from.segm.gaugings.1,
  #                                 stage.record          = stage.record,
  #                                 dir.gaugings.results  = dir.SPD.segm.gaug.results,
  #                                 t.shift.for.b         = shift.times.gaugings.1,
  #                                 gaugings.period       = gaugings,
  #                                 times.uncert          = TRUE)
  
  
  
  #########################################
  for (iter.model in 1:length(rec.model)) {
  #########################################
    message("##############################################################################"); flush.console()    
    message(paste0("SEGMENTATION OF PARAMETERS OF STAGE-RECESSION MODEL '", rec.model[iter.model], "'"));  flush.console()
    message("##############################################################################"); flush.console()
    message("Method used for segmenting the time series: Darienzo et al., 2021"); flush.console()
    message("
*****************************************************************
A few information:
*****************************************************************
- Please, notice that it may take some time (e.g. 1 hour). 
- It is a single-pass multi change point detection based on
  Bayesian approach:
- The segmentation of the time series is performed for increasing 
  number of segments.
- then the choice of the optimal segmentation is based on a
  criterion (e.g. min BIC).
- it is possible to adjust the detected shift times (if any).
- the segmentation accounts for data uncertainties (combining both 
  time series data uncertainties and structural uncertainty). 
*****************************************************************
")
    
  if (rec.model[[iter.model]] == "1expWithAsympt" ){
    param.var.model = c("a1", "a2")

  } else if (rec.model[[iter.model]] =="2expWithAsympt" ){
    param.var.model = c("a1", "a3")

  } else if (rec.model[[iter.model]] =="2expWithAsympt_rel" ){
    param.var.model = c("a1", "a3")

  } else if (rec.model[[iter.model]] =="3expWithAsympt" ){
    param.var.model = c("a1", "a4")

  } else if (rec.model[[iter.model]] =="3expWithAsympt_bis" ){
    param.var.model = c("a1", "a2", "a4")

  } else if (rec.model[[iter.model]] =="expexp" ){
    param.var.model = c("a1", "a2")

  } else if (rec.model[[iter.model]] =="hyperb" ){
    param.var.model = c("a1", "a2")

  } else if ((rec.model[[iter.model]] =="hyperb_bis" )|
             (rec.model[[iter.model]] =="expexp_bis" )|
             (rec.model[[iter.model]] =="Coutagne_bis" )) {
    param.var.model = c("a1", "b1", "a2")

  } else if (rec.model[[iter.model]] =="2expWithAsympt_bis" ){
    param.var.model = c("a1", "a2", "a3")
  }

  dir.segm.recessions.model        = paste0(dir.segm.recessions,"/", rec.model[[iter.model]])
  dir.estim.recessions.model       = paste0(dir.estim.recessions,"/Pooling/model_", rec.model[[iter.model]])
  dir.segm.recessions.model.param  = paste0(dir.segm.recessions.model ,"/chi_",chi)
  dir.estim.recessions.model.param = paste0(dir.estim.recessions.model,"/chi_",chi)
  dir.segm.recessions.param        = NULL
  results.segment[[iter.model]]    = list()
  dir.create(paste0(dir.segm.recessions,"/", rec.model[[iter.model]]))
  dir.create(paste0(dir.segm.recessions.model,"/chi_",chi))
  
  
  
  
  
  
  #######################################
  if (seg.plot.results.only == TRUE) {
  #######################################
    if((file.exists(paste0(dir.estim.recessions.model.param,"/Results_Summary.txt")) ==FALSE)|
       (file.exists(paste0(dir.estim.recessions.model.param,"/Results_MCMC_cooked.txt")) == FALSE)|
       (file.exists(paste0(dir.estim.recessions.model.param,"/Results_Residuals.txt")) == FALSE)|
       (file.exists(paste0(dir.estim.recessions.model.param,"/Curves_Data.txt")) == FALSE)){
      
      plot.results.new= as.logical(dlgInput(c("You have selected the option to read results only",
                                              "(seg.plot.results.only= T) in the Options_Recession file", 
                                              "But there are no results to read and process.", " ",
                                              "--> If you want to perform also the recession estimation please press 'F'",
                                              "--> Otherwise press 'T'."), 
                                            Sys.info()[" "])$res)
      if (plot.results.new == FALSE){
          seg.plot.results.only = FALSE
      } else {
          err1 = print("******* ERROR: there are no results to process, please put seg.plot.results.only = FALSE")
          return(err1)
      }
    }
  }
  
  
  
  #####################################
  if (seg.plot.results.only == FALSE) {
  ####################################
  if (segm_all_parameters == TRUE){
      parm = seq(1,length(param.var.model),1)
  } else {
      parm = length(param.var.model) # only the asymptot parameter
  }
  #####################
  for (param in parm) {
  #####################
       message("  "); flush.console()
       message("*****************************************************************"); flush.console()
       message(paste0("Segmentation of parameter ", param.var.model[param])); flush.console()
       message("*****************************************************************"); flush.console()
       message("Options selected:"); flush.console()
       message(paste0("- Prior mean for the segments means ~ ",prior.mu.rec.segment[1], "(", prior.mu.rec.segment[2],",", prior.mu.rec.segment[3],")"));  flush.console()
       message(paste0("- Prior for structural uncertainty of the segments ~ ",gamma.prior.rec[1], "(", gamma.prior.rec[2], ",",  gamma.prior.rec[3],")"));  flush.console()
       message(paste0("- Maximum number of segments: ", nSmax.rec));  flush.console()
       message(paste0("- Criterion for the choice of optimal segmentation: ", criterion.rec));  flush.console()
       message(paste0("- Number of cycles of mcmc during segmentation: ", Ncycle.rec.segment));  flush.console()
       message(paste0("- Number of mcmc per cycle during segmentation: ", Nmcmc.rec.segment));  flush.console()
       message(paste0("- Fraction of mcmc that is burned: ", Nburn.rec.segment));  flush.console()
       message(paste0("- Minimum time between two shifts: ", tmin.rec));  flush.console()
       if (shift.time.adjustment.type.rec == 1){
         message("- Shift time adjustment: always the Maximum A Posterior estimate.")
       } else if  (shift.time.adjustment.type.rec == 2){
         message("- Shift time adjustment: always the largest flood peak in the 95% CI.") 
       } else {
         message("- Shift time adjustment: manual selection at each iteration") 
       }
       if (plot.b.from.gaugings == TRUE) {
         message("- You have chosen to perform a comparison with gaugings segmentation.");  flush.console()
       } 
       message("Start segmentation ...")
 
       dir.create(paste0(dir.segm.recessions.model.param,"/", param.var.model[param]))
       dir.segm.recessions.param[[param]] = paste0(dir.segm.recessions.model.param,"/",param.var.model[param])
       file.dir.with.data                 = paste0(dir.segm.recessions.model.param,"/df.",param.var.model[param],".txt")
       
       #initialisation of recessions segmentation :
       seg.period.rec       = 1; 
       end.end              = FALSE
       i_init.rec           = i_final.rec  = BIC.rec = AIC.rec = npar.rec  = maxpost.rec  = 0
       i_init.rec[1]        = 1
       times.of.shift.rec   =  mean.of.segments.rec = t.q10.rec = t.q90.rec = final.period.rec = NULL
       
         
       # read data to segment from file:
       if (file.exists(file.dir.with.data) == FALSE) {
           errr = print(paste0("********* ERROR: cannot find file '", file.dir.with.data, 
                               "'!! It seems that you have not performed the recessions estimation yet!",
                               " Please check."))
           return(list(err = errr))
       }
       
       #########################
       if (BayesianOption ==1) { # old version of the method (without pooling approach)
       #########################
             data.segm.rec  = read.table(file = paste0(dir.segm.recessions,"/df.h.infinity.txt"),header=TRUE)
             tss_tot_ns.rec = c(1,length(data.segm.rec$t));
             #write data file for BaM :
             data.segm.rec.BaM  = data.frame(t  = data.segm.rec$t,  h  = data.segm.rec$mean,  uh = data.segm.rec$stdev)
             write.table(data.segm.rec.BaM, paste0(dir.config.segmentation,"/Segm_data.txt"),sep="\t",row.names=FALSE)
       ######## 
       } else { # Bayesian pooling approach:
       ########
             data.segm.rec  = read.table(file = file.dir.with.data, header=TRUE)
             tss_tot_ns.rec = c(1,length(data.segm.rec$t));
             data.segm.rec.BaM  = data.frame(t  = data.segm.rec$t, h  = data.segm.rec$maxpost, uh = data.segm.rec$stdev)
             if (param.var.model[param] == tail(param.var.model,1)){
                      data.segm.rec.BaM  = data.frame(t  = data.segm.rec$t, 
                                                      h  = data.segm.rec$maxpost, # - stage.scale.shift, 
                                                      uh = data.segm.rec$stdev)
             }   
             write.table(data.segm.rec.BaM, paste0(dir.config.segmentation,"/Segm_data.txt"),sep="\t",row.names=FALSE)
       }
       ################################################################################################################
       # initialisation of criteria:  
       #                   AIC = Aikake Information Criterion
       #                   BIC = Bayesian Information Criterion or SIC (Schwarz)
       #                   mBIC = modified Bayesianb Information Criterion
       #                   HQC = Hannan-Quinn information Criterion
       #                   DIC = Deviance Informatiojn Criterion
       AIC = AICc= BIC= mBIC = HQC = DIC=0;
       npar = maxpost = varLogpost= loglikelihood = maxLikelihood = varLogLikelihood=0;
       N = length(data.segm.rec.BaM$h) #Number of observations
       end.seg = FALSE
       i = 1
       # start Bayesian segmentation (Darienzo et al. 2021): 
       ########################################
       while ((i <= nSmax.rec) & ((N-i) >=2)) {
      #########################################
              nS = i
              print(c("n.segm = ",nS))
              npar[i] = 2*nS
              setwd(dir.exe)
              dir.create(paste0(dir.segm.recessions.param[[param]],"/nS",nS))
              dir.nS = paste0(dir.segm.recessions.param[[param]],"/nS",nS)
              # Launch BaM exe segmentation:
              Segmentation_app.rec(dir.exe          = dir.exe,
                                   nobs             = N, 
                                   nS               = nS, 
                                   tmin             = tmin.rec,
                                   tfirst           = data.segm.rec.BaM$t[1],
                                   tfin             = tail(data.segm.rec.BaM$t,1),
                                   ncycles          = Ncycle.rec.segment,
                                   nmcmc            = Nmcmc.rec.segment,
                                   Nslim            = Nslim.rec.segment,
                                   Nburn            = Nburn.rec.segment,
                                   tP               = data.segm.rec.BaM$t,
                                   resid.uncertaint = TRUE,
                                   prior.mu.segm    = prior.mu.rec.segment,
                                   gamma_prior      = gamma.prior.rec)

             list.of.files.segment <- c(paste0(dir.config.segmentation,"/Config_model.txt"),
                                        paste0(dir.config.segmentation,"/Segm_data.txt"),
                                        paste0(dir.config.segmentation,"/Results_MCMC_cooked.txt"))
             for (ll in 1:length(list.of.files.segment)) {
                   file.copy(list.of.files.segment[ll], dir.nS, overwrite = TRUE)
             }
             mcmc.segm    = read.table(file=paste0(dir.config.segmentation,"/Results_MCMC_cooked.txt"),header=TRUE)
             resid.segm   = read.table(file=paste0(dir.config.segmentation,"/Results_Residuals.txt"),header=TRUE)
             summary.segm = read.table(file=paste0(dir.config.segmentation,"/Results_Summary.txt"),header=TRUE)
             
             write.table(mcmc.segm,    file=paste0(dir.nS,"/Results_MCMC_cooked_","nS",nS,".csv"), sep=";")
             write.table(resid.segm,   file=paste0(dir.nS,"/Results_Residuals_","nS",nS,".csv"), sep=";")
             write.table(summary.segm, file=paste0(dir.nS,"/Results_Summary_","nS",nS,".csv"), sep=";")
             
             if (save.all.results == TRUE) {
                 plot.mcmc.segment(workspace = dir.nS, seg.iter=1, nS=nS)
                 converg = Convergence.test.segment(dir.seg  = dir.nS,  npar = npar[nS], dir.plot = dir.nS)
                 if (converg==FALSE) {
                    print(paste0("Segmentation ",nS," does not converged !!! Need to increse Ncycles ..."))
                 }
             }
             logpost           = mcmc.segm[,(2*nS+1)]
             loglikelihood     = singlelikelihoods = single.prior.mu = single.prior.tau = 0;
             len.mcmc          = length(mcmc.segm[,1])
             priors.mu         = matrix(NA, nrow = length(mcmc.segm[,1]), ncol = nS)
             if (prior.mu.rec.segment[1] == "Gaussian"){
               for (j in 1:nS) {
                 priors.mu[,j] = dnorm(mcmc.segm[,j],  mean = as.numeric(prior.mu.rec.segment[2]), sd = as.numeric(prior.mu.rec.segment[3]), log = TRUE)
               }
             } else if (prior.mu.rec.segment[1] == "Uniform") {
               for (j in 1:nS) {
                 priors.mu[,j] = dunif(mcmc.segm[,j], min = as.numeric(prior.mu.rec.segment[2]), max = as.numeric(prior.mu.rec.segment[3]), log = TRUE)
               }  
             }
             if (nS > 1) {
               priors.tau = matrix(0, nrow = len.mcmc, ncol = nS-1)
             } else {
               priors.tau = matrix(0, nrow = len.mcmc, ncol = 1)
             }
             priors.gamma = 0
             priors.gamma = dunif(mcmc.segm[, 2*nS], min = as.numeric(gamma.prior.rec[2]), max = as.numeric(gamma.prior.rec[3]), log = TRUE)
             #
             logprior = 0
             for (ll in 1:len.mcmc){
               logprior[ll] = sum(priors.mu[ll,]) + sum(priors.tau[ll,])  + priors.gamma[ll]   
             }
             loglikelihood        = logpost - logprior
             df.mcmc.LL           = data.frame(loglikelihood = loglikelihood,  logprior = logprior, logposterior  = logpost)
             maxpost[nS]          = max(logpost)
             maxLikelihood[nS]    = max(loglikelihood)
             Likelihood.maxpost   = loglikelihood[which.max(logpost)]
             varLogpost[nS]       = var(logpost)
             varLogLikelihood[nS] = var(loglikelihood)
             MeanDev              = -2*mean(loglikelihood)
             write.table(df.mcmc.LL, file=paste0(dir.nS,"/likelihood_","nS",nS,".csv"), sep=";")
             #---------------------------------------------------------------------------------------------
             #AIC[i] = 2*npar[nS]- 2*Likelihood.maxpost[nS]          # Aikake,1974
             AIC[i] = 2*npar[nS]- 2*maxLikelihood[nS]                # Aikake,1974
             #AICc[i] = AIC[i] +  2*(npar[nS])*(npar[nS]+1)/(N-npar[nS]-1) #Hurvich and Tsai,1989
             #BIC[i] = (npar[nS])*log(N)-2*Likelihood.maxpost[nS]    # Schwarz,1978
             BIC[i] = log(N)*(npar[nS]) -2*maxLikelihood[nS]         # Schwarz,1978
             HQC[i] = 2*(npar[nS])*log(log(N)) - 2*maxLikelihood[nS] # Hannan-Quinn
             #DIC[i] = 2*MeanDev - Dev.of.mean                       # Spiegelhalter et al.,2002
             #DIC[i] = MeanDev + 2*varLogpost[nS]                    # Gelman 2004 "Bayesian data analysis"
             DIC[i] = MeanDev + 2*varLogLikelihood[nS]               # Gelman 2004 "Bayesian data analysis"
             #---------------------------------------------------------------------------------------------
             i = i+1
  }
  
       
  ######################################
  # Analysis of results: CHOICE OF MODEL
  ######################################
  # BIC / AIC / AICc / HQC / DIC:
  #min values  of the criteria:
  BICmin = which.min(BIC);
  AICmin = which.min(AIC);
  HQCmin = which.min(HQC);
  DICmin = which.min(DIC);
  #AICcmin = which.min(AICc, na.rm=TRUE);
  # dataframe with all criteria (for the plot):
  criteria.df =  data.frame(BIC = BIC, AIC = AIC, HQC = HQC, DIC = DIC, x = seq(1, nS, 1), BICmin = BICmin, AICmin=AICmin, HQCmin = HQCmin, DICmin = DICmin)
  BICplot  = model.selection.plot(criteria.df, dir.segm.recessions.param[[param]],  seg.iter =1)
  if (criterion.rec == "AIC") {
    nS.ok = AICmin
  } else if (criterion.rec == "AICc") {
    nS.ok = AICcmin
  } else if (criterion.rec == "BIC") {
    nS.ok = BICmin
  } else if (criterion.rec == "DIC") {
    nS.ok = DICmin
  } else if (criterion.rec == "HQC") {
    nS.ok = HQCmin
  }
  print(paste0("=========> Optimal number of segments (considering the minimum ", criterion.rec,") = ", nS.ok))
  dir.nS.ok        = paste0(dir.segm.recessions.param[[param]],"/nS", nS.ok)
  
  
  
  
  
  
  # READ RESULTS OF THE OPTIMAL SEGMENTATION:
  ###########################################
  Residuals        = read.table(file= paste0(dir.nS.ok,"/Results_Residuals_nS",nS.ok,".csv"),sep=";",header=TRUE)
  mu.s             = as.numeric(Residuals[,5])
  Results.seg      = read.table(file= paste0(dir.nS.ok,"/Results_Summary_nS",nS.ok,".csv"),sep=";",header=TRUE)
  mcmc.segment     = read.table(file= paste0(dir.nS.ok,"/Results_MCMC_cooked_nS",nS.ok,".csv"),sep=";",header=TRUE)
  #initialisation of results:
  Q2.ts = NULL; Q97.ts = NULL; Q2.mu = NULL; Q97.mu = NULL;
  if ( nS.ok > 1) {
    # change point times "tau":
    ts.res         <- as.numeric(c(Results.seg[16,(nS.ok+1):(npar[nS.ok]-1)]))  #the maxpost of all mcmc
    ts.mean        <- as.numeric(c(Results.seg[5,(nS.ok+1):(npar[nS.ok]-1)]))  #the maxpost of all mcmc
    ts.median      <- as.numeric(c(Results.seg[6,(nS.ok+1):(npar[nS.ok]-1)]))  #the maxpost of all mcmc
    ts.res.plus    <- as.numeric(c(Results.seg[16,(nS.ok+1):(npar[nS.ok]-1)], tail(data.segm.rec$t,1)))
    ts.res.before  <- as.numeric(c(data.segm.rec$t[1],Results.seg[16,((nS.ok+1):(npar[nS.ok]-1))]))
    stdev.ts       <- as.numeric(c(Results.seg[11,(nS.ok+1):(npar[nS.ok]-1)]))
    for (i in 1:(nS.ok-1)) {
         Q2.ts[i]  <- c(quantile(mcmc.segment[,(nS.ok+i)], p = c(0.025)))
         Q97.ts[i] <- c(quantile(mcmc.segment[,(nS.ok+i)], p = c(0.975)))
    }
    Q10.ts         <- as.numeric(c(Results.seg[7,(nS.ok+1):(npar[nS.ok]-1)]))
    Q90.ts         <- as.numeric(c(Results.seg[10,(nS.ok+1):(npar[nS.ok]-1)]))
    # segments mean "mu":
    mu.res         <- as.numeric(c(Results.seg[16,1:nS.ok]))
    mu.mean        <- as.numeric(c(Results.seg[5,1:nS.ok]))  #the maxpost of all mcmc
    mu.median      <- as.numeric(c(Results.seg[6,1:nS.ok]))  #the maxpost of all mcmc
    for (j in 1:nS.ok) {
         Q2.mu[j]  <- c(quantile(mcmc.segment[,j], p = c(0.025)))
         Q97.mu[j] <- c(quantile(mcmc.segment[,j], p = c(0.975)))
    }
    Q10.mu.res     <- as.numeric(Results.seg[7,1:(nS.ok)])
    Q90.mu.res     <- as.numeric(Results.seg[10,1:(nS.ok)])
    stdev.mu       <- as.numeric(Results.seg[11,1:(nS.ok)])
    #structural error parameter "gamma":
    gamma.MAP      <- as.numeric(Results.seg[16,npar[nS.ok]])
    gamma.stdev    <- as.numeric(Results.seg[11,npar[nS.ok]])
    gamma.mean     <- as.numeric(Results.seg[5,npar[nS.ok]])
    Q2.gamma       <- quantile(mcmc.segment[,npar[nS.ok]], p = c(0.025))
    Q97.gamma      <- quantile(mcmc.segment[,npar[nS.ok]], p = c(0.975))
    pdf.ts.rec     =  mcmc.segment[, (nS.ok+1):(2*nS.ok-1)]
    
    
  } else { # if no change points ==> only one segment !!
    ts.res         <- NULL; stdev.ts <- NULL; ts.res.before <- NULL;  ts.res.plus <- NULL; ts.mean <- NULL;
    ts.median      <- NULL; stdev.ts <- NULL; Q2.ts <- NULL; Q97.ts <- NULL; Q10.ts <- NULL; Q90.ts <-NULL;
    mu.res         <- as.numeric(Results.seg[16,1])
    mu.mean        <- as.numeric(Results.seg[5,1])  #the maxpost of all mcmc         
    mu.median      <- as.numeric(Results.seg[6,1])  #the maxpost of all mcmc
    Q2.mu          <- quantile(mcmc.segment[,1], p = c(0.025))
    Q97.mu         <- quantile(mcmc.segment[,1], p = c(0.975))
    Q10.mu.res     <- as.numeric(Results.seg[7,1])
    Q90.mu.res     <- as.numeric(Results.seg[10,1])
    stdev.mu       <- as.numeric(Results.seg[11,1])
    gamma.MAP      <- as.numeric(Results.seg[16,2])
    gamma.stdev    <- as.numeric(Results.seg[11,2])
    gamma.mean     <- as.numeric(Results.seg[5,2])
    Q2.gamma       <- quantile(mcmc.segment[,2], p = c(0.025))
    Q97.gamma      <- quantile(mcmc.segment[,2], p = c(0.975))
    pdf.ts.rec     = NULL
  }
  
  #saving to a data.frame statistics of segmentation parameters tau, mu and gamma:
  tau.results.df = data.frame(# change point times "tau":
                              tau.MAP= ts.res, tau.q2=Q2.ts, tau.q10 = Q10.ts, tau.q90 = Q90.ts, tau.q97 = Q97.ts, 
                              tau.stdev = stdev.ts , tau.mean = ts.mean,  tau.median = ts.median)
  mu.results.df= data.frame(# Segment mean "mu":
                            mu.MAP= mu.res, mu.q2= Q2.mu, mu.q10 = Q10.mu.res, mu.q90 = Q90.mu.res, mu.q97 = Q97.mu,
                            mu.stdev = stdev.mu, mu.mean = mu.mean,  mu.median = mu.median)
  gamma.results.df= data.frame(# Segment mean "mu":
                                gamma.MAP = gamma.MAP, gamma.q2= Q2.gamma, gamma.q97 = Q97.gamma,
                                gamma.stdev = gamma.stdev, gamma.mean = gamma.mean)

  times.of.shift.MAP <- NULL; mean.of.segments <- NULL; t.q10 <- NULL; t.q90 <- NULL;  t.q2 <- NULL; t.q97 <- NULL;
  ts.all.real <- NULL; ts.all.MAP <- NULL; ts.all.q2 = NULL; ts.all.q10 = NULL; ts.all.q90 = NULL; ts.all.q97 = NULL;
  ts.morpho.real = NULL; ts.morpho.MAP = NULL; ts.morpho.q2 = NULL; ts.morpho.q97 = NULL;
 
   # saving to vectors:
  times.of.shift.MAP <- c(times.of.shift.MAP, ts.res)
  mean.of.segments <- c(mean.of.segments,mu.res)
  t.q2  = c(t.q2, Q2.ts)
  t.q10 = c(t.q10, Q10.ts)
  t.q90 = c(t.q90, Q90.ts)
  t.q97 = c(t.q97, Q97.ts)
  times.of.shift.MAP = sort(times.of.shift.MAP)
  mean.of.segments   = sort(mean.of.segments)
  t.q2  = sort(t.q2)
  t.q10 = sort(t.q10)
  t.q90 = sort(t.q90)
  t.q97 = sort(t.q97)
  t.q2  = t.q2[c(TRUE,  !t.q2[-length(t.q2)] == t.q2[-1])]
  t.q10 = t.q10[c(TRUE, !t.q10[-length(t.q10)] == t.q10[-1])]
  t.q90 = t.q90[c(TRUE, !t.q90[-length(t.q90)] == t.q90[-1])]
  t.q97 = t.q97[c(TRUE, !t.q97[-length(t.q97)] == t.q97[-1])]
  times.of.shift.MAP <- times.of.shift.MAP[c(TRUE, !times.of.shift.MAP[-length(times.of.shift.MAP)] == times.of.shift.MAP[-1])]
  mean.of.segments   <-  mean.of.segments[c(TRUE, !mean.of.segments[-length( mean.of.segments)] ==  mean.of.segments[-1])]
  
  
  
  
  
  ###############################################################################################################
  #Interval of the shift time:     
  interval = list(); ts.real = hflood = tflood= hflood2 = tflood2=NULL;
  if (nS.ok > 1) { # if at least one shift has been detected:
    for (i in 1:length(tau.results.df$tau.MAP)) {
      interval[[i]] = which((stage.record$t_limni >= min(tau.results.df$tau.q2[i],  tau.results.df$tau.MAP[i])) &
                            (stage.record$t_limni <= max(tau.results.df$tau.MAP[i], tau.results.df$tau.q97[i])))
      #maximum stage value (if known):
      if (!is.null(interval[[i]])) {
        hflood[i] = max(h_limni[interval[[i]]])
        tflood[i] = stage.record$t_limni[which.max(stage.record$h_limni[interval[[i]]])  +  interval[[i]][1]]
      } else {
        tflood[i] = tau.MAP[i]
      }
    }
    tau.results.df  = cbind(tau.results.df, tflood = tflood) 
    # gaugings <- data.frame("h"      = h_Gaug,
    #                        "Q"      = Q_BaRatin,
    #                        "uQ"     = uQ_Gaug,
    #                        "Period" = 1 , 
    #                        "t"      = t_Gaug,
    #                        "t.true" = t_gaug.true)
    if (!is.null(gaugings)){
        CdT.tot         = data.frame(t = gaugings$t,  Q = gaugings$Q,  h = gaugings$h)      
    } else {
        CdT.tot = NULL
    }
    initial.tsplot  = initial.ts.plot.rec(CdT.P           = CdT.tot, 
                                          df.limni        = stage.record, 
                                          tshift          = tau.results.df,
                                          limni.labels    = limni.labels, 
                                          grid_limni.ylim = grid_limni.ylim, 
                                          dir.seg.gaug    = dir.segm.recessions.param[[param]],
                                          seg.iter        = 1, 
                                          mcmc.segment    = mcmc.segment, 
                                          nS              = nS.ok)
    
    # pupup plot of stage record with proposed segmentation
    X11(); print(initial.tsplot)
    # Options:  1) if there is a flood, assign the shift time to the max peak !
    #           2) if the shift is due to other causes (vegetation, works, ice ...)
    #              then assign the shift time to the maxpost
    #           3) if the shift time is known then fix it.
    user.choice.ts=NULL
    i =0
    while (i < length(tau.results.df$tau.MAP)) {
      i = i +1
      
      # 3 Options:
      ############
      if (shift.time.adjustment.type.rec == 1) {   # Option 1: always chose the MAP of the shift time !!!
          ts.real[i] = tau.results.df$tau.MAP[i] 
        
          
      } else if (shift.time.adjustment.type.rec == 2) {  # Option 2: select largest flood peak.
      ############
        if (any(ts.real==tflood[i]) | any(ts.all.real==tflood[i])){
          print("Same flood selected twice...")
          hflood2[i] = hflood[i]
          delay.flood = 0.01 # days  ####### to change this ??????!!!!!
          print(paste0("...delay the second flood of: ", delay.flood))
          tflood2[i] = tflood[i] + delay.flood 
          # hflood2[i] = sort(h_limni[interval[[i]]], TRUE)[2]
          # tflood2[i] = t_limni[which_nth_highest_vaalue(x=h_limni[interval[[i]]], n=2)[1]
          #             + interval[[i]][1]]
          ts.real[i] = tau.results.df$tau.MAP[i]
        } else {
          hflood2[i] = hflood[i]
          tflood2[i] = tflood[i]
          ts.real[i] = tflood2[i]
        }
        ts.morpho.real = c(ts.morpho.real, ts.real[i])
        ts.morpho.MAP  = c(ts.morpho.MAP,  ts.res[i])
        ts.morpho.q2   = c(ts.morpho.q2,   Q2.ts[i])
        ts.morpho.q97  = c(ts.morpho.q97,  Q97.ts[i])
        
        
      } else {    # Option 3: select manually the shift time.
      ##############
        # pop-up window for user choice of the TRUE shift time option:
        user.choice.ts[i] <- dlgInput(paste0("Which Adjusted shift time do you chose for ts",i," ? \n",
                                             "1 = MAP shift time =  ", ts.res[i], " days     interval=[ ", Q2.ts[i], " ; ", Q97.ts[i], " ]  \n",
                                             "2 = stage max (morphogenic flood at t_flood =", tflood[i], "   ;   with stage h_flood = ", hflood[i], " ) \n",
                                             "3 = other time, e.g. earthquake, cyclon, ..."), Sys.info()[" "])$res
        #choices for the adjustment:
        if (user.choice.ts[i] == 1) {  # MAP
             ts.real[i]        = ts.res[i]
             
        } else if (user.choice.ts[i] == 2) {  # flood event
             ts.real[i]        = tflood[i]
             ts.morpho.real    = c(ts.morpho.real, ts.real[i])
             ts.morpho.MAP     = c(ts.morpho.MAP,  ts.res[i])
             ts.morpho.q2      = c(ts.morpho.q2,   Q2.ts[i])
             ts.morpho.q97     = c(ts.morpho.q97,  Q97.ts[i])
             
        } else if (user.choice.ts[i] == 3) {  # other events, insert time manually.
             ts.real[i]        = dlgInput("insert the date (in days) ", Sys.info()[" "])$res
             ts.real[i]        = as.numeric(ts.real[i])
             
             # define morphological shifts (if known):
             morpho.ask        = dlgInput("is it related to sediment transport dynamics ? [Y/N]", Sys.info()[" "])$res
             if (morpho.ask== "Y") {
             #----------------------
                ts.morpho.real = c(ts.morpho.real, ts.real[i])
                ts.morpho.MAP  = c(ts.morpho.MAP,  ts.res[i])
                ts.morpho.q2   = c(ts.morpho.q2,   Q2.ts[i])
                ts.morpho.q97  = c(ts.morpho.q97,  Q97.ts[i])
             }
        }    
      }
    }
    dev.off() 
    
    
    # Update shift times df:
    ts.real             = as.numeric(ts.real)
    ts.res.plus         = c( as.numeric(ts.real), tail(data.segm.rec$t,1))
    ts.res.before       = c(data.segm.rec$t[1],  ts.real)
    df.shift.times      = data.frame(Q2.ts, Q97.ts, ts.res,  ts.real)
    df.shift.times.plus = data.frame(ts.res.before, ts.res.plus, Q2.mu, Q10.mu.res, Q90.mu.res, Q97.mu,  mu.res)
    df.shift.times2     = data.frame(tMAP	= ts.res, treal =	ts.real, t2 =	Q2.ts, t10 = Q10.ts, t90 = Q90.ts, t97 = Q97.ts)
    
  ##########
  } else { # no shifts detected:  
  ##########
    
    df.shift.times         = NULL
    df.shift.times$ts.real = NULL
    df.shift.times.plus    = NULL
    ts.res                 = NULL
    stdev.ts               = NULL 
    ts.res.before          = NULL
    ts.res.plus            = NULL
    ts.mean                = NULL 
    ts.median              = NULL
    stdev.ts               = NULL
    Q2.ts                  = NULL
    Q97.ts                 = NULL
    Q10.ts                 = NULL
    Q90.ts                 = NULL
  }
  
  #######################################################################
  
  # Gaugings:
  message("Determining periods of stability of gaugings ..."); flush.console()
  if (!is.null(gaugings)){
      color= colors.period
      c_Gaug = 0; P_Gaug = 0;
      if (!is.null(df.shift.times$ts.real)) {
          for (i in 1:length(gaugings$t)) {
            if(gaugings$t[i] <= df.shift.times$ts.real[1]) {
                #points(x=hP[i], y=QP[i], log ="y", col = color[1],pch=1,lwd=4)
                c_Gaug[i] = color[1]
                P_Gaug[i] = 1
            }
          }
          for (j in 2:nS.ok) {
            for (i in 1:length(gaugings$t)) {
              if ((gaugings$t[i] <= tail(gaugings$t,1)) & 
                 (gaugings$t[i] >  df.shift.times$ts.real[j-1])) {
                  #points(x=hP[i], y=QP[i], log ="y", col = color[j],pch=1,lwd=4)
                  c_Gaug[i] = color[j]
                  P_Gaug[i] = j
              }
            }
          }
      } else {
          for (i in 1:length(gaugings$t)) {
            #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
            c_Gaug[i] = color[1]
            P_Gaug[i] = 1
          }
      }
      df.RC <- data.frame(gaugings$h, gaugings$Q, gaugings$uQ,  P_Gaug, gaugings$t, c_Gaug)
      names(df.RC) = c("h","Q", "uQ", "Period", "t","color")
  } else {
      message("No gaugings provided !!!"); flush.console()
      df.RC = NULL
  }
  
  
  ############################################################################
  # plot results of segmentation:
  message("Saving and plotting results ..."); flush.console()
  if (param.var.model[param] == tail(param.var.model,1)){
    limits.Y = asymptote.limits #stage.limits
    Y.name   = y.name
  } else {
    limits.Y = "automatic"
    Y.name   = param.var.model[param] 
  }
  
  
  limni.time.limits = c(stage.record$t_limni[1],  tail(stage.record$t_limni,1))
  # plotting results of segmentation:
  plot1 = plot.segmentation.results.rec(dir.results          = dir.segm.recessions.param[[param]],
                                        df                   = data.segm.rec.BaM,
                                        df.shifts            = tau.results.df,
                                        ts.res.before        = ts.res.before,
                                        ts.res               = ts.res,
                                        ts.real              = ts.real,
                                        ts.res.plus          = ts.res.plus,
                                        df.mu                = mu.results.df,
                                        mu.res               = mu.res,
                                        mcmc.segment         = mcmc.segment,
                                        nS                   = nS.ok,
                                        x.name               = x.name ,
                                        limits.X             = limni.time.limits,
                                        y.name               = Y.name,
                                        limits.Y             = limits.Y,
                                        error.type           = 3,
                                        plot.axis.x          = TRUE,
                                        plot.axis.y          = TRUE,
                                        plot.title           = param.var.model[param],
                                        plot.mu.uncertainty  = TRUE,
                                        plot.tau.uncertainty = FALSE,
                                        plot.shift.dates     = TRUE,
                                        plot.pdf             = TRUE,
                                        points.size          = 1,
                                        plot.file.name       = paste0("Segment_", param.var.model[param]))
  # save into files:
  write.table(df.RC, paste0(dir.segm.recessions.param[[param]],"/data_with_periods.txt"), sep ="\t", row.names=FALSE)
  write.table(mcmc.segment, paste0(dir.segm.recessions.param[[param]],"/mcmc_segmentation.txt"),
              sep ="\t", row.names=FALSE)

  write.table(df.shift.times, paste0(dir.segm.recessions.param[[param]],"/df.shift.times.txt"),
              sep ="\t", row.names=FALSE)
  write.table(df.shift.times2, paste0(dir.segm.recessions.param[[param]],"/shift_times.txt"),
              sep ="\t", row.names=FALSE)
  write.table(pdf.ts.rec, paste0(dir.segm.recessions.param[[param]],"/pdf_ts.txt"),
              sep ="\t", row.names=FALSE)
  
  write.table(tau.results.df, paste0(dir.segm.recessions.param[[param]],"/tau.results.df.txt"),
              sep ="\t", row.names=FALSE)
  write.table(mu.results.df, paste0(dir.segm.recessions.param[[param]],"/mu.results.df.txt"),
              sep ="\t", row.names=FALSE)
  write.table(gamma.results.df, paste0(dir.segm.recessions.param[[param]],"/gamma.results.df.txt"),
              sep ="\t", row.names=FALSE)

  write.table(df.shift.times.plus, paste0(dir.segm.recessions.param[[param]],"/df.shift.times.plus.txt"),
              sep ="\t", row.names=FALSE)
  write.table(ts.res.before, paste0(dir.segm.recessions.param[[param]],"/ts.res.before.txt"),
              sep ="\t", row.names=FALSE)
  write.table(ts.res, paste0(dir.segm.recessions.param[[param]],"/ts.res.txt"),
              sep ="\t", row.names=FALSE)
  write.table(ts.res.plus, paste0(dir.segm.recessions.param[[param]],"/ts.res.plus.txt"),
              sep ="\t", row.names=FALSE)
  write.table(Q2.ts, paste0(dir.segm.recessions.param[[param]],"/Q2.ts.txt"),
              sep ="\t", row.names=FALSE)
  write.table(Q97.ts, paste0(dir.segm.recessions.param[[param]],"/Q97.ts.txt"),
              sep ="\t", row.names=FALSE)
  write.table(Q2.mu, paste0(dir.segm.recessions.param[[param]],"/Q2.mu.txt"),
              sep ="\t", row.names=FALSE)
  write.table(mu.res, paste0(dir.segm.recessions.param[[param]],"/mu.res.txt"),
              sep ="\t", row.names=FALSE)
  write.table(Q97.mu, paste0(dir.segm.recessions.param[[param]],"/Q97.mu.txt"),
              sep ="\t", row.names=FALSE)
  write.table(data.segm.rec, paste0(dir.segm.recessions.param[[param]],"/Data.segm.rec.txt"),
              sep ="\t", row.names=FALSE)
  # final results for the rating shift times detection method:
  data.annotate.recess = data.frame(q2=Q2.ts,  MAP=ts.res, q97=Q97.ts, t.adj=ts.res)
  # list of results for all parameters and models:
  results.segment[[iter.model]][[param]] = list(model.and.par.name   = c("model"= rec.model[[iter.model]], "par" = param.var.model[param]),
                                                nS.ok                = nS.ok, 
                                                tau.results.df       = tau.results.df, 
                                                mu.results.df        = mu.results.df,
                                                ts.res.before        = ts.res.before, 
                                                ts.res               = ts.res, 
                                                ts.real              = ts.real,
                                                ts.res.plus          = ts.res.plus, 
                                                Q2.ts                = Q2.ts, 
                                                Q97.ts               = Q97.ts,
                                                Q2.mu                = Q2.mu, 
                                                mu.res               = mu.res, 
                                                Q97.mu               = Q97.mu,
                                                data.segm.rec.BaM    = data.segm.rec.BaM, 
                                                mcmc.segment         = mcmc.segment,
                                                df.shift.times       = df.shift.times, 
                                                df.shift.times.plus  = df.shift.times.plus,
                                                df.RC                = df.RC, 
                                                df.gamma             = gamma.results.df,
                                                data.annotate.recess = data.annotate.recess)
    }
  }
  
  
  ############################################
  ############################################
  # plot figures for results of only one model:
  plot.segmentation.recession.one.model( dir.rec.pool.test       = dir.estim.recessions.model.param,  
                                         dir.rec.segm.test       = dir.segm.recessions.model.param,
                                         rec.mod                 = rec.model[[iter.model]],
                                         stage.limits            = stage.limits,
                                         stage.scale.shift       = stage.scale.shift,
                                         asymptote.limits        = asymptote.limits,
                                         limits.Y.lambda         = limits.Y.lambda,
                                         limits.Y.alpha          = limits.Y.alpha,
                                         limits.x.recess         = limits.x.recess,
                                         x.name                  = x.name,
                                         plot.b.from.gaugings    = plot.b.from.gaugings,
                                         plot.recession.uncert   = plot.recession.uncert,
                                         plot.gamma.uncertainty  = plot.gamma.uncertainty,
                                         plot.dates              = plot.dates,
                                         date_origin             = date_origin,
                                         initial.time            = initial.time,
                                         gaugings                = gaugings,
                                         limni.labels            = limni.labels,
                                         grid_limni.ylim         = c(stage.limits[1]/100, grid_limni.ylim[2], grid_limni.ylim[3]),
                                         df.limni                = stage.record)

  }
  print("****************")
  print("   All done!    ")
  print("****************")
#####################################################################
  return(results.segment)
####################################################################
}   







































####################################################################################################################
read.results.segment.recess = function(dir.segm.recessions.par, 
                                       officialShiftsTime, 
                                       Gaugings, 
                                       plot.dates) {
####################################################################################################################
# This function has the objective to read the results of the segmentation of the time series of recession parameters  
  mcmc.seg.rec         = read.table(file=paste0(dir.segm.recessions.par, "/mcmc_segmentation.txt"), header=TRUE)
  mu.results.df.rec    = read.table(file=paste0(dir.segm.recessions.par, "/mu.results.df.txt"), header=TRUE)
  gamma.results.df.rec = read.table(file=paste0(dir.segm.recessions.par, "/gamma.results.df.txt"), header=TRUE)
  gamma_segm_recess    = c(mean(mcmc.seg.rec$Y1_gamma1), std(mcmc.seg.rec$Y1_gamma1))
  Q2.mu.rec            = read.table(file=paste0(dir.segm.recessions.par, "/Q2.mu.txt"), header=TRUE)
  mu.res.rec           = read.table(file=paste0(dir.segm.recessions.par, "/mu.res.txt"), header=TRUE)
  Q97.mu.rec           = read.table(file=paste0(dir.segm.recessions.par, "/Q97.mu.txt"), header=TRUE)
  Data.segm.rec        = read.table(file=paste0(dir.segm.recessions.par, "/Data.segm.rec.txt"), header=TRUE)
  
  if (length(mu.results.df.rec$mu.MAP) > 1){ # at least one shift detected !!!!!
    tau.results.df.rec      = read.table(file=paste0(dir.segm.recessions.par, "/tau.results.df.txt"), header=TRUE)
    df.shift.times.rec      = read.table(file=paste0(dir.segm.recessions.par, "/df.shift.times.txt"), header=TRUE)
    shift.times.recessions  = read.table(file=paste0(dir.segm.recessions.par,"/df.shift.times.txt"), header = TRUE)
    df.shift.times.plus.rec = read.table(file=paste0(dir.segm.recessions.par, "/df.shift.times.plus.txt"), header=TRUE)
    ts.res.before.rec       = read.table(file=paste0(dir.segm.recessions.par, "/ts.res.before.txt"), header=TRUE)
    ts.res.rec              = read.table(file=paste0(dir.segm.recessions.par, "/ts.res.txt"), header=TRUE)  
    ts.res.plus.rec         = read.table(file=paste0(dir.segm.recessions.par, "/ts.res.plus.txt"), header=TRUE)
    Q2.ts.rec               = read.table(file=paste0(dir.segm.recessions.par, "/Q2.ts.txt"), header=TRUE)
    Q10.ts.rec              = tau.results.df.rec$tau.q10
    Q90.ts.rec              =  tau.results.df.rec$tau.q90
    Q97.ts.rec              = read.table(file=paste0(dir.segm.recessions.par, "/Q97.ts.txt"), header=TRUE)
    nS.ok.rec               = nrow(df.shift.times.rec)+1
    pdf.ts.rec              = mcmc.seg.rec[, (nS.ok.rec+1):(2*nS.ok.rec-1)]
    data.annotate.recess    = data.frame(t      = shift.times.recessions$ts.res,
                                         tstart = shift.times.recessions$Q2.ts,
                                         tend   = shift.times.recessions$Q97.ts,
                                         t.adj  = shift.times.recessions$ts.real)
    if (plot.dates == TRUE){
      dates.recess = 0
      for (dat in 1:length(data.annotate.recesst$t.adj)) {
        dates.recess[dat]  =  data.annotate.recess$t.adj[dat]  +  Gaugings$X.3[1]
      }
      RealDates.recess = as.Date(dates.recess,  origin = "1899-12-30")
      RealDates.recess.new = strptime(as.character(RealDates.recess), "%Y-%m-%d" )
      RealDates.recess.newnew = format( RealDates.recess.new, "%d/%m/%Y")
    } else {
      RealDates.recess.newnew = NULL
    }
  } else {
    
    nS.ok.rec  = 1
    pdf.ts.rec = NULL
    data.annotate.recess = NULL
    data.annotate.recess.adjust = NULL
    RealDates.recess.newnew = NULL
    tau.results.df.rec = NULL
    ts.res.before.rec = NULL 
    ts.res.rec = NULL
    ts.res.plus.rec = NULL 
    Q2.ts.rec = NULL
    Q10.ts.rec = NULL
    Q90.ts.rec = NULL
    Q97.ts.rec = NULL
  }
  
  ############## Official shift dates:
  if (is.null(officialShiftsTime)==FALSE) {
    data.annotate.off <- data.frame(xeffect = officialShiftsTime,
                                    xpotent = officialShiftsTime)
  } else {
    data.annotate.off =NULL
  }
  
  if (!is.null(Gaugings)){
      gaugings.df.recess = read.table(file =paste0(dir.segm.recessions.par, "/data_with_periods.txt"), header = TRUE)
  } else {
      gaugings.df.recess = NULL
  }
  

    
  return(list(
    nS.ok                       = nS.ok.rec, 
    mcmc.segment                = mcmc.seg.rec,
    gaugings.df                 = gaugings.df.recess,
    data.annotate.recess        = data.annotate.recess,
    pdf.ts.rec                  = pdf.ts.rec,
    tau.results.df              = tau.results.df.rec, 
    mu.results.df               = mu.results.df.rec,
    gamma.results.df            = gamma.results.df.rec,
    ts.res.before               = unlist(ts.res.before.rec), 
    ts.res                      = unlist(ts.res.rec), 
    ts.res.plus                 = unlist(ts.res.plus.rec), 
    Q2.ts                       = unlist(Q2.ts.rec), 
    Q10.ts                      = unlist(Q10.ts.rec), 
    Q90.ts                      = unlist(Q90.ts.rec), 
    Q97.ts                      = unlist(Q97.ts.rec), 
    Q2.mu                       = unlist(Q2.mu.rec), 
    mu.res                      = unlist(mu.res.rec), 
    Q97.mu                      = unlist(Q97.mu.rec),
    data.annotate.off           = data.annotate.off,
    dates                       = RealDates.recess.newnew,
    gamma                       = gamma_segm_recess))
}































###########################################################################################
Segmentation_app.rec <- function(dir.exe, nobs,nS, tmin, tfirst, tfin,
                                 ncycles , nmcmc, Nslim, Nburn, tP, resid.uncertaint, 
                                 prior.mu.segm, gamma_prior) {
###########################################################################################
  Segmentation_config.rec(dir.exe, nobs, nS, tmin, tfirst, tfin, 
                          ncycles, nmcmc, Nslim, Nburn, tP, resid.uncertaint, 
                          prior.mu.segm, gamma_prior)
  message("Applying Segmentation - BaM !!!  Wait ... "); flush.console()
  # system2(paste(dir_code,"/BaM_exe/BaM_Segmentation.exe",sep=""),stdout =NULL, stderr = NULL);
  system2(paste0(dir.exe,"/BaM_Segmentation2.exe"), stdout =NULL, stderr = NULL); 
}
















####################################################################################################
Segmentation_config.rec <- function(dir.exe, nobs, nS, tmin, tfirst, tfin,
                                    ncycles, nmcmc, Nslim, Nburn, tP, resid.uncertaint, 
                                    prior.mu.segm, gamma_prior) {
####################################################################################################
  npar = nS + nS - 1  # Number of model parameters
  N = length(tP)      # number of observations to segment
  # Config model file for BaM:
  file.name1 = paste0(dir.exe,"/Segmentation/Config_Model.txt")
  cat('"Segmentation2"', file =file.name1,sep="\n")
  cat(1,    file =file.name1, append = TRUE,sep="\n")
  cat(1,    file =file.name1, append = TRUE,sep="\n")
  cat(npar, file =file.name1, append = TRUE,sep="\n")
  for (i in 1:nS) {
    cat(paste('"k',i,'"',sep=""), file =file.name1, append = TRUE,sep="\n")
    cat(prior.mu.segm[2],         file =file.name1, append = TRUE,sep="\n")
    cat(prior.mu.segm[1],         file =file.name1, append = TRUE,sep="\n")            #! Prior distribution
    cat(prior.mu.segm[2],         file =file.name1, append = TRUE,sep=",")
    cat(",",                      file =file.name1, append = TRUE,sep=",")
    cat(prior.mu.segm[3],         file =file.name1, append = TRUE,sep="\n")
    # cat("'FlatPrior'", file =file.name1, append = TRUE,sep="\n")
    # cat(" ", file =file.name1, append = TRUE,sep="\n")
  }
  #************************************************************************
  if (nS>1) {
    tP_half=0
    for (ii in 1:length(tP)) {
      tP_half[ii] = (tP[ii+1]+tP[ii])/2
    }
    if ((length(tP)-nS) <= 2) {
      tstart=0
      for (ee in 1: (nS-1)) {
        tstart[ee] = tP_half[ee]
      }
    } else {
      step = trunc((length(tP)-1)/(nS-1), digits = 0)
      init = trunc((length(tP)-1)/2)  
      tstart=0
      aaa = 0
      for (ee in 1: (nS-1)) {
        if (ee %% 2 == 0){
          tstart[ee] = tP_half[init + step*aaa]
        } else {
          tstart[ee] = tP_half[init - step*aaa]
          aaa=aaa+1
        }
      }
      tstart =sort(tstart)
      #tstart = start.values.time(nS, tfirst, tfin)
    }
    for (i in 1:(nS-1)) {
      cat(paste('"tau',i,'"',sep=""), file =file.name1, append = TRUE,sep="\n")
      cat(tstart[i],      file = file.name1, append = TRUE,sep="\n")
      cat("'FlatPrior+'", file = file.name1, append = TRUE,sep="\n")
      cat(" ",            file = file.name1, append = TRUE,sep="\n")
      # cat('"Gaussian"', file = file.name1, append = TRUE, sep="\n")            #! Prior distribution
      # cat((tfirst+(tfin-tfirst)/nS*i), file =file.name1, append = TRUE, sep=",")
      # cat(",",file =file.name1, append = TRUE, sep=",")
      # cat((tfin-tfirst)/4, file =file.name1, append = TRUE, sep="\n")
    }
  }
  #####################################################
  file.name2 = paste0(dir.exe,"/Segmentation/Config_Data.txt")
  cat("'Segmentation\\Segm_data.txt'", file =file.name2,sep="\n")
  cat(1,    file = file.name2, append = TRUE,sep="\n")
  cat(nobs, file = file.name2, append = TRUE,sep="\n") 
  cat(3,    file = file.name2, append = TRUE,sep="\n")
  cat(1,    file = file.name2, append = TRUE,sep="\n")
  cat(0,    file = file.name2, append = TRUE,sep="\n")
  cat(0,    file = file.name2, append = TRUE,sep="\n")
  cat(0,    file = file.name2, append = TRUE,sep="\n")
  cat(2,    file = file.name2, append = TRUE,sep="\n")
  cat(3,    file = file.name2, append = TRUE,sep="\n")
  cat(0,    file = file.name2, append = TRUE,sep="\n")
  cat(0,    file = file.name2, append = TRUE,sep="\n")
  ####################################################
  file_BaM = paste0(dir.exe,"/Config_BaM.txt")
  #creation of Config_BaM.txt
  cat('"Segmentation/"',           file = file_BaM , sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"',   file = file_BaM , sep="\n", append = TRUE)   
  cat('"Config_Model.txt"',        file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_xtra.txt"',         file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Data.txt"',         file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file_BaM , sep="\n", append = TRUE)                                      
  cat('"Config_MCMC.txt"',         file = file_BaM , sep="\n", append = TRUE)                                           
  cat('"Config_Cooking.txt"',      file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Summary.txt"',      file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"',    file = file_BaM , sep="\n", append = TRUE)
  #cat('"Config_Pred_Master.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('" "',                       file = file_BaM , sep="\n", append = TRUE)
  ###########################################################################
  # file extra:
  # file.xtra = paste(dir_code,"/BaM_exe/Segmentation/Config_xtra.txt",sep="")
  # cat(nS, file =file.xtra,sep="\n")
  # cat(tmin, file = file.xtra, append = TRUE,sep="\n")
  file.xtra = paste0(dir.exe,"/Segmentation/Config_xtra.txt")
  cat(nS,  file = file.xtra,sep="\n")
  cat(tmin, file = file.xtra, append = TRUE, sep="\n") # tmin
  cat(1,   file = file.xtra, append = TRUE, sep="\n")   # Nmin
  cat(1,   file = file.xtra, append = TRUE, sep="\n")   # option 
  #####################################################################  RUN OPTIONS
  file.run = paste0(dir.exe,"/Segmentation/Config_RunOptions.txt")
  cat(".true.",  file =file.run,sep="\n")                     #Do MCMC?
  cat(".true.",  file =file.run, append = TRUE, sep="\n")     #Do MCMC summary?
  cat(".true.",  file =file.run, append = TRUE, sep="\n")     #Do Residual diagnostics?
  cat(".false.", file =file.run, append = TRUE, sep="\n")    #Do Predictions?
  ######################################################################   RESIDUALS CONFIG
  file.residuals = paste0(dir.exe,"/Segmentation/Config_Residuals.txt")
  cat("Results_Residuals.txt" , file =file.residuals ,sep="\n")     #Result file
  ######################################################################   SUMMARY CONFIG
  file.summary = paste0(dir.exe,"/Segmentation/Config_Summary.txt")
  cat("Results_Summary.txt" , file =file.summary ,sep="\n")          #Result file
  ######################################################################   COOKING CONFIG
  file.cooking = paste0(dir.exe,"/Segmentation/Config_Cooking.txt")
  cat("Results_MCMC_Cooked.txt" , file =file.cooking ,sep="\n")      #Result file
  cat(Nburn, file =file.cooking, append = TRUE, sep="\n")            #Burn factor
  cat(Nslim, file =file.cooking, append = TRUE, sep="\n")            #Nslim
  ###################################################################################   REMNANT ERROR CONFIG
  file.remnant = paste0(dir.exe,"/Segmentation/Config_RemnantSigma.txt")
  cat("'Constant'",   file = file.remnant, sep="\n")                    #! Function f used in sdev = f(Qrc)
  cat(1,              file = file.remnant, append = TRUE, sep="\n")     #! Number of parameters gamma for f
  cat("gamma1",       file = file.remnant, append = TRUE, sep="\n")     #! Parameter Name
  cat(gamma_prior[4], file = file.remnant, append = TRUE, sep="\n")     #! Initial Guess
  cat(gamma_prior[1], file = file.remnant, append = TRUE, sep="\n")     #! Prior distribution
  cat(gamma_prior[2], file = file.remnant, append = TRUE, sep=",")
  cat(",",            file = file.remnant, append = TRUE, sep=",")
  cat(gamma_prior[3], file = file.remnant, append = TRUE, sep="\n")
  
  # cat('"LogNormal"', file = file.remnant, append = TRUE, sep="\n")   #! Prior distribution
  # cat(-10,file =file.remnant, append = TRUE, sep=",")
  # cat(",",file =file.remnant, append = TRUE, sep=",")
  # cat(1,file =file.remnant, append = TRUE, sep="\n")
  
  ##########################################################################   MCMC
  file.mcmc = paste0(dir.exe,"/Segmentation/Config_MCMC.txt")
  cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
  cat(nmcmc,   file = file.mcmc, append = TRUE,sep="\n")    # Nadapt
  cat(ncycles, file = file.mcmc, append = TRUE,sep="\n")  # Ncycles
  cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")       # minMoveRate
  cat(0.5,     file = file.mcmc, append = TRUE,sep="\n")       # maxMoveRate
  cat(0.9,     file = file.mcmc, append = TRUE,sep="\n")       # DownMult
  cat(1.1,     file = file.mcmc, append = TRUE,sep="\n")       # UpMult
  cat(0,       file = file.mcmc, append = TRUE,sep="\n")       # mode for init jump distr
  cat("****",  file = file.mcmc, append = TRUE,sep="\n")
  cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")       # MultFact
  cat(0.1,     file = file.mcmc, append = TRUE,sep=",")       # RC MultiFact
  cat(0.1,     file = file.mcmc, append = TRUE,sep=",")    
  cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")
  cat(0.1,     file = file.mcmc, append = TRUE,sep=",")       # Remnant MultiFact
  cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")
}






















  
  
  
  
###################################################################################################
performance.recess.analysis = function() {
###################################################################################################  
  #plot DIC and other performance criteria to compare all models and all chi:
  plot.performance.model.comparison(model.names     =   c("1expWithAsympt", 
                                                          "2expWithAsympt", 
                                                          "2expWithAsympt_bis",
                                                          "3expWithAsympt", 
                                                          "3expWithAsympt_bis", 
                                                          #"2expWithAsympt_rel",
                                                          "expexp", 
                                                          "expexp_bis",
                                                          "hyperb", 
                                                          "hyperb_bis", 
                                                          "Coutagne", 
                                                          "Coutagne_bis"),
                                    
                                    model.titles     = c("M1", "M2", "M3", "M4", "M5", "M6", 
                                                         "M7", "M8", "M9", "M10", "M11"),
                                    chi.test        = c(10, 30, 50),  #for different values of chi
                                    dir.case_study  = dir.case_study, 
                                    priors          = list(prior.param.rec[[1]],     #priors for the each model
                                                           prior.param.rec[[2]],
                                                           prior.param.rec[[7]],
                                                           prior.param.rec[[3]],
                                                           prior.param.rec[[8]],
                                                           #prior.param.rec[[12]],
                                                           prior.param.rec[[4]],
                                                           prior.param.rec[[9]],
                                                           prior.param.rec[[5]],
                                                           prior.param.rec[[10]],
                                                           prior.param.rec[[6]],
                                                           prior.param.rec[[11]]
                                    ), 
                                    which.recession     = which.recession,
                                    bt.from.gaugingsss  = read.bt,  #river bed estimation
                                    pdf.ts.gaugings     = pdf.ts.results.1,
                                    data.annotate.off   = data.annotate.off,
                                    time.limits         = limni.time.limits,
                                    grid_limni.ylim     = c(-2 , 4.5 , 1))
  
  
  # final results for the rating shift times detection method:
  data.annotate.recess <- data.frame(q2    = read.res.rec$Q2.ts, 
                                     MAP   = read.res.rec$data.annotate.recess$t ,
                                     q97   = read.res.rec$Q97.ts,
                                     t.adj = read.res.rec$data.annotate.recess.adjust$t.adj)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



### REAL TIME analysis :
###########################################################################################################
real.time.recession = function(){
###########################################################################################################
dddf = bind_rows( rec.selection[[1]], .id = "column_label")
rt.rec = ggplot()+
         geom_point(aes(x=dddf$t,y=dddf$h),   size=1)
#         ggplot()+
#         geom_point(aes(x=rec.selection[[1]][[15]]$t,
#                        y=rec.selection[[1]][[15]]$h), 
#                    size=1)+
#         geom_point(aes(x=rec.selection[[1]][[40]]$t,
#                        y=rec.selection[[1]][[40]]$h), 
#                    size=1, color="green")+
#         geom_point(aes(x=rec.selection[[1]][[50]]$t,
#                        y=rec.selection[[1]][[50]]$h), 
#                    size=1, color="red")+
#         scale_y_continuous(limits=c(900, 1200))+
#         scale_x_continuous(limits=c(0,150))
# 
recess.realtime = realtime.recess.2(       rec.model.rt             = rec.model[2],
                                           data.recess.realtime     = rec.selection,
                                           which.recession.realtime = which.recession.realtime,
                                           dir.case_stud   = dir.case_study,
                                           dir_code        = dir_code, 
                                           dir_exe         = dir.exe,
                                           station.name    = station.name,
                                           data.period     = data.period,
                                           prior.rec.realtime     = prior.param.rec[[2]],
                                           prior.gamma.realtime   = c(0, 100, "'Uniform'", 0.1,   # gamma1
                                                                      0, 100, "'Uniform'", 1),  # gamma2
                                           ncycl               = 1000, 
                                           ncycles.max         = 2000,# Ncycles.mcmc.rec,
                                           nmcmc               = 100,
                                           nslim               = 50,
                                           jump.pos            = 1.2,
                                           jump.neg            = 0.8,
                                           nburn               = 0.5, 
                                           pred                = TRUE,
                                           chi.rt              = c(50),
                                           delta.t.min.rt      = 0.3,
                                           delta.t.max.rt      = 100,
                                           uh.rec.rt           = 0.5,
                                           limits.y            = c(-100, 350, 50), #stage.limits,
                                           limits.x            = limits.x.recess,  # recession time is in "days"
                                           long.recessions     = FALSE,
                                           stage.scale.shift   = stage.scale.shift,
                                           starting.time       = 3, # min number of days
                                           realtime.period     = c(5860, 6500),
                                           frames.time         = c(6, 15, 23))


recess.realtime.animation = realtime.recess.3(rec.model.rt             = rec.model,
                                              data.recess.realtime     = rec.selection,
                                              which.recession.realtime = which.recession.realtime,
                                              dir.case_stud            = dir.case_study,
                                              dir_code                 = dir_code, 
                                              dir_exe                  = dir.exe,
                                              station.name             = station.name,
                                              data.period              = data.period,
                                              prior.rec.realtime       = prior.param.rec,
                                              prior.gamma.realtime     = c(0, 100, "'Uniform'", 0.1,   # gamma1
                                                                           0, 100, "'Uniform'", 1),  # gamma2
                                              ncycl                    = 1000, 
                                              ncycles.max              = 2000,   # Ncycles.mcmc.rec,
                                              nmcmc                    = 100,
                                              nslim                    = 50,
                                              jump.pos            = 1.2,
                                              jump.neg            = 0.8,
                                              nburn               = 0.5, 
                                              pred                = TRUE,
                                              chi.rt              = c(50),
                                              delta.t.min.rt      = 0.3,
                                              delta.t.max.rt      = 100,
                                              uh.rec.rt           = 0.5,
                                              limits.y            = c(-100, 200, 50), #stage.limits,
                                              limits.x            = c(0, 100, 20),  # recession time is in "days"
                                              long.recessions     = FALSE,
                                              stage.scale.shift   = stage.scale.shift,
                                              starting.time       = 3, # min number of days
                                              realtime.period     = c(5860, 6500))
 return()
}








  








######################################################################################################################
realtime.recession <- function(data.rec.rt,
                               df.asymptotes.rt.old,
                               t_limni.rt.rec,
                               rec.model.rt,
                               param.var.model.rt,
                               Nburn.rec.rt, 
                               Nslim.rec, 
                               jump.pos.rec, jump.neg.rec,
                               Nmin.rec,
                               dir_code, 
                               dir_exe, 
                               dir.real.time,
                               dir.rec.retrosp.rt,
                               station.name,
                               data.period,
                               ncycl.rec, # Ncycles.mcmc.rec,
                               nmcmc.rec,
                               tgood,    #tgood ,
                               pred,
                               prior.param.rec.rt,
                               prior.gamma.rec.rt,
                               limits.y, limits.x,
                               time.step,
                               stage.scale.shift,
                               data4BaRatin.rt) {   # recession time in "days", stage in "cm"
  #####################################################################################################################
  start_time <- Sys.time()
  #Initialisation of files and plots:
  output_file_h = paste0(dir.real.time,"/Param_rec_h.csv")
  colfunc = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple"))   
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
  #**********************************************************************************
  
  
  #BaM application (Benjamin Renard, Irstea):
  setwd(dir_exe)
  dir.BaM.rec            = paste0(dir_code,"/BaM_exe/Recession_h")
  curve_data.h           = "Recession_h/Curves_Data.txt"
  data.rec.rt.modified   = data.rec.rt
  data.rec.rt.modified$h = data.rec.rt.modified$h + stage.scale.shift
  data.rec.rt.modified   = data.rec.rt.modified[-c(1:Nburn.rec.rt),]
  nobs.h                 = length(data.rec.rt.modified$t)
  write.table(data.rec.rt.modified, 
              file= curve_data.h, append = FALSE, sep = "\t", eol = "\n", 
              na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh"))
  other.grid = seq(round(tail(data.rec.rt.modified$t,1),0), 
                   #round(tail(data.rec.rt.modified$t,1),0) + 100, 2)
                   120, 2)
  tgrid =c(data.rec.rt.modified$t, 
           other.grid)
  
  # Recover the prior from past recessions:
  curves_retrospect.rt = read.table(paste0(dir.rec.retrosp.rt,"/Pooling/test_",
                                           rec.model.rt, "/chi_",
                                           chi.rt,"/Curves_Data.txt"), header =TRUE)
  Results_retrospect.rt = read.table(paste0(dir.rec.retrosp.rt,"/Pooling/test_",
                                            rec.model.rt, "/chi_",
                                            chi.rt,"/Results_Summary.txt"), header =TRUE)
  Nrecess = tail(curves_retrospect.rt$Period,1)
  
  if ((rec.model.rt== "2expWithAsympt")){
    #######################################################
    priors.rt = c(prior.param.rec.rt[1:4], "static",
                  Results_retrospect.rt[5, Nrecess+1], Results_retrospect.rt[11, Nrecess+1], "Gaussian", Results_retrospect.rt[5, Nrecess+1], "static" ,
                  Results_retrospect.rt[5, Nrecess+2], Results_retrospect.rt[11, Nrecess+2], "Gaussian", Results_retrospect.rt[5,Nrecess+2], "static",
                  Results_retrospect.rt[5, Nrecess+3], Results_retrospect.rt[11, Nrecess+3], "Gaussian", Results_retrospect.rt[5,Nrecess+3], "static",
                  prior.param.rec.rt[21:24], "static")
    priors.gamma.rt = c(Results_retrospect.rt[5, 2*Nrecess+4], Results_retrospect.rt[11,2*Nrecess+4], "Gaussian", Results_retrospect.rt[5,2*Nrecess+4],
                        Results_retrospect.rt[5, 2*Nrecess+5], Results_retrospect.rt[11,2*Nrecess+5], "Gaussian", Results_retrospect.rt[5,2*Nrecess+5])
    
  } else if ((rec.model.rt == "2expWithAsympt_bis")){
    ######################################################
    priors.rt = c(prior.param.rec.rt[1:4], "static",
                  Results_retrospect.rt[5, Nrecess+1], Results_retrospect.rt[11, Nrecess+1], "Gaussian", Results_retrospect.rt[5, Nrecess+1], "static" ,
                  prior.param.rec.rt[11:14], "static",
                  Results_retrospect.rt[5,2*Nrecess+2], Results_retrospect.rt[11,2*Nrecess+2], "Gaussian", Results_retrospect.rt[5,2*Nrecess+2], "static",
                  prior.param.rec.rt[21:24], "static")
    priors.gamma.rt = c(Results_retrospect.rt[5,3*Nrecess + 3], Results_retrospect.rt[11,3*Nrecess+3], "Gaussian", Results_retrospect.rt[5,3*Nrecess+3],
                        Results_retrospect.rt[5,3*Nrecess + 4], Results_retrospect.rt[11,3*Nrecess+4], "Gaussian", Results_retrospect.rt[5,3*Nrecess+4])
    
  } else if ((rec.model.rt == "3expWithAsympt_bis")){
    ######################################################
    priors.rt = c(prior.param.rec.rt[1:4], "static",
                  Results_retrospect.rt[5, Nrecess+1], Results_retrospect.rt[11, Nrecess+1], "Gaussian", Results_retrospect.rt[5, Nrecess+1], "static" ,
                  prior.param.rec.rt[11:14], "static",
                  Results_retrospect.rt[5,2*Nrecess+2], Results_retrospect.rt[11,2*Nrecess+2], "Gaussian", Results_retrospect.rt[5,2*Nrecess+2], "static",
                  Results_retrospect.rt[5,2*Nrecess+3], Results_retrospect.rt[11,2*Nrecess+3], "Gaussian", Results_retrospect.rt[5,2*Nrecess+3], "static",
                  Results_retrospect.rt[5,2*Nrecess+4], Results_retrospect.rt[11,2*Nrecess+4], "Gaussian", Results_retrospect.rt[5,2*Nrecess+4], "static",
                  prior.param.rec.rt[31:34], "static"
    )
    priors.gamma.rt = c(Results_retrospect.rt[5,3*Nrecess + 5], Results_retrospect.rt[11,3*Nrecess+5], "Gaussian", Results_retrospect.rt[5,3*Nrecess+5],
                        Results_retrospect.rt[5,3*Nrecess + 6], Results_retrospect.rt[11,3*Nrecess+6], "Gaussian", Results_retrospect.rt[5,3*Nrecess+6])
    
  } else if ( (rec.model.rt == "3expWithAsympt")){
    
  }
  
  #######################################################
  # Bam config for recession analysis:
  BaM_config1.h(  tgrid         = seq(limits.x[1], 
                                      limits.x[2], 0.5),
                  nobs          = nobs.h,
                  ncycles       = ncycl.rec,
                  nmcmc         = nmcmc.rec,
                  nslim         = Nslim.rec,
                  jump.pos      = jump.pos.rec, 
                  jump.neg      = jump.neg.rec,  
                  nburn         = 0.5,
                  model         = rec.model.rt,
                  pred          = pred,
                  asympt.before = NULL, #curves.h$Qcurve[index.h[icurve]-1],
                  remnant       = "Linear",
                  prior         =  priors.rt ,
                  prior.gamma   =  priors.gamma.rt)
  # prior    =     prior.rec.realtime[[mod]] ,
  #prior.gamma =  prior.gamma.realtime)
  #Launch BAM exe:
  system2("BaM_recession_multi_model_final.exe", stdout =NULL, stderr =NULL)
  #save results files:
  list.of.files.rec <- c(paste0(dir.BaM.rec,"/Results_MCMC_Cooked.txt"),
                         paste0(dir.BaM.rec,"/Results_Residuals.txt"),
                         paste0(dir.BaM.rec,"/Results_Summary.txt"),
                         paste0(dir.BaM.rec,"/Config_Model.txt"),
                         paste0(dir.BaM.rec,"/Curves_Data.txt"),
                         paste0(dir.BaM.rec,"/tgrid.txt"),
                         paste0(dir.BaM.rec,"/ht_Maxpost.spag"),
                         paste0(dir.BaM.rec,"/ht_ParamU.spag"),
                         paste0(dir.BaM.rec,"/ht_ParamU.env"),
                         paste0(dir.BaM.rec,"/ht_TotalU.env"),
                         paste0(dir.BaM.rec,"/ht_TotalU.spag"))
  for (i in 1:length(list.of.files.rec)) {
    file.copy(list.of.files.rec[i], dir.real.time ,overwrite = TRUE)
  }
  # plot.mcmc(workspace = dir.real.time, 
  #           iter= 2)  #plot density and traceplots of MCMC
  # conv = Convergence.test(dir.seg = dir.real.time ,  npar=7 , dir.plot = dir.real.time)
  setwd(dir.exe)
  # if (conv==FALSE) {
  #   print(paste0("recess ",icurve," Not converged !!! "))
  # }
  Results_residuals.rt = read.table("Recession_h/Results_Residuals.txt", header =TRUE)
  residuals.plot = ggplot(Results_residuals.rt) +
    geom_point(aes(x = Results_residuals.rt$X1_obs,
                   y = Results_residuals.rt$Y1_sim), color = "red", size = 2)+
    geom_point(aes(x = Results_residuals.rt$X1_obs,
                   y = Results_residuals.rt$Y1_obs), color = "black", size=1)+
    theme_bw()
  ggsave(residuals.plot, filename =paste0(dir.real.time,"/residuals.png"),
         bg = "transparent", width = 8, height =4, dpi = 200)
  
  #reading results of BaM:
  tgrid.file <- paste0(dir.exe,"/Recession_h/tgrid.txt")
  maxpost = read.table(paste0(dir.BaM.rec,"/ht_Maxpost.spag"), header =FALSE)
  env     = read.table(paste0(dir.BaM.rec,"/ht_TotalU.env"), header =TRUE)
  env.param = read.table(paste0(dir.BaM.rec,"/ht_ParamU.env"), header =TRUE)
  RecCurve.h = read.table(paste0(dir.BaM.rec,"/Curves_Data.txt"), header =TRUE)
  summary.rec = read.table(paste0(dir.BaM.rec,"/Results_Summary.txt"), header =TRUE)
  mcmc.rec = read.table(paste0(dir.BaM.rec,"/Results_MCMC_Cooked.txt"), header =TRUE)
  
  if ((rec.model.rt == "2expWithAsympt")| (rec.model.rt == "2expWithAsympt_bis")){
    # a1:
    a1.mcmc = mcmc.rec[,1]
    a1.summary = summary.rec[,1]
    a1.df = data.frame(  t(c(quantile(a1.mcmc, p = c(0.025, 0.5, 0.975)),
                             mean=mean(a1.mcmc),
                             stdev =std(a1.mcmc),
                             maxpost = a1.summary[16])))
    names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    #b1:
    b1.mcmc = mcmc.rec[,2]
    b1.summary = summary.rec[,2]
    b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)),
                             mean=mean(b1.mcmc),
                             stdev =std(b1.mcmc),
                             maxpost = b1.summary[16])))
    names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    #a2:
    a2.mcmc = mcmc.rec[,3]
    a2.summary = summary.rec[,3]
    a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)),
                             mean=mean(a2.mcmc),
                             stdev =std(a2.mcmc), maxpost = a2.summary[16])))
    names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    
    #b2:
    b2.mcmc = mcmc.rec[,4]
    b2.summary = summary.rec[,4]
    b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)),
                             mean=mean(b2.mcmc),
                             stdev =std(b2.mcmc), maxpost = b2.summary[16])))
    names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    
    
    #a3:
    a3.mcmc = mcmc.rec[,5]
    a3.summary = summary.rec[,5]
    a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)),
                             mean=mean(a3.mcmc),
                             stdev =std(a3.mcmc), maxpost = a3.summary[16])))
    names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    a3.df[, c(1,2,3,4,6)] = a3.df[,c(1,2,3,4,6)]  -  stage.scale.shift
    h.infinity = cbind(a3.df ,   t =  t_limni.rt.rec)
    
  } else if ((rec.model.rt =="3expWithAsympt")| (rec.model.rt =="3expWithAsympt_bis")){
    ###############################################################################################     
    # read single parameters results :   a1(k) , b1, a2, b2, a3, b3, a4(k) 
    # a1:
    a1.mcmc = mcmc.rec[,1]
    a1.summary = summary.rec[,1]
    a1.df = data.frame(  t(c(quantile(a1.mcmc, p = c(0.025, 0.5, 0.975)), 
                             mean=mean(a1.mcmc), 
                             stdev =std(a1.mcmc), 
                             maxpost = a1.summary[16])))
    names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    #b1:
    b1.mcmc = mcmc.rec[,2]
    b1.summary = summary.rec[,2]
    b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), 
                             mean=mean(b1.mcmc), 
                             stdev =std(b1.mcmc), 
                             maxpost = b1.summary[16])))
    names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    #a2:
    a2.mcmc = mcmc.rec[,3]
    a2.summary = summary.rec[,3]
    a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), 
                             mean=mean(a2.mcmc), 
                             stdev =std(a2.mcmc),
                             maxpost = a2.summary[16])))
    names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    
    #b2:
    b2.mcmc = mcmc.rec[,4]
    b2.summary = summary.rec[,4]
    b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)),
                             mean=mean(b2.mcmc), 
                             stdev =std(b2.mcmc),
                             maxpost = b2.summary[16])))
    names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    
    
    #a3:
    a3.mcmc = mcmc.rec[,5]
    a3.summary = summary.rec[,5]
    a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)), 
                             mean=mean(a3.mcmc), 
                             stdev =std(a3.mcmc),
                             maxpost = a3.summary[16])))
    names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    
    #b3:
    b3.mcmc = mcmc.rec[,6]
    b3.summary = summary.rec[,6]
    b3.df = data.frame(  t(c(quantile(b3.mcmc, p = c(0.025, 0.5, 0.975)), 
                             mean=mean(b3.mcmc), 
                             stdev =std(b3.mcmc), 
                             maxpost = b3.summary[16])))
    names(b3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    
    #a4: asymptote !!!
    a4.mcmc = mcmc.rec[,7]
    a4.summary = summary.rec[,7]
    a4.df = data.frame(  t(c(quantile(a4.mcmc, p = c(0.025, 0.5, 0.975)), 
                             mean=mean(a4.mcmc), 
                             stdev =std(a4.mcmc), 
                             maxpost = a4.summary[16])))
    names(a4.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
    a4.df[, c(1,2,3,4,6)] = a4.df[,c(1,2,3,4,6)]  -  stage.scale.shift
    h.infinity = cbind(a4.df ,   t =  t_limni.rt.rec)
    
    
  }
  #############################################################################################
  #plot regression:
  ttttt = seq(limits.x[1],limits.x[2], 0.5)
  temp <- data.frame(x = data.rec.rt.modified$time   -  data.rec.rt$time[1],
                     y = data.rec.rt.modified$h      - stage.scale.shift,
                     z = data.rec.rt.modified$uh)
  temp2 <- data.frame(xx = ttttt ,
                      yy = maxpost - stage.scale.shift,
                      zz = env[,2] - stage.scale.shift,
                      kk = env[,3] - stage.scale.shift)
  write.table(temp, file=paste0(dir.real.time, "/temp.txt"), append = FALSE, sep = "\t", eol = "\n",
              na = "NA", dec = ".", row.names = FALSE, col.names=c("x", "y", "z"))
  write.table(temp2, file=paste0(dir.real.time, "/temp2.txt"), append = FALSE, sep = "\t", eol = "\n",
              na = "NA", dec = ".", row.names = FALSE, col.names=c("xx", "yy", "zz", "kk"))
  
  reg.rt.plot     =  ggplot() +
    theme_bw(base_size = 15)+
    scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
                       expand = c(0,0),     breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
    scale_y_continuous(name = "Stage h (cm)", limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin=unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
    geom_point(data = temp , aes(x = x, y= y), color= "black", size= 2)+  #colfunc(Ncurves)[icurve], size= 2)+
    geom_ribbon(data = temp2, aes(x = xx,    ymin = zz , ymax = kk),
                fill = "green", alpha=0.1)  # colfunc(Ncurves)[icurve], alpha = 0.1)
  ggsave(reg.rt.plot, filename =paste0(dir.real.time,"/regression_realtime.png"),
         bg = "transparent", width = 12, height =7, dpi = 200)
  
  
  
  
  #################################################################################################################
  #   asympt.plot =  geom_plot() +
  #                  geom_point(data = h.infinity, aes(x = t, y = maxpost), colour = "black",size = 3) +
  #                  #geom_boxplot(boxplot.df, aes(x = x, ymin=y0, lower =y2, middle= y50,
  #                                #upper=y97, ymax= y100), stat = "identity") +
  #                  #geom_line(data = asym.df.temp[[i]], aes(x =t, y =h), colour = "gray",size = 1) +
  #                  geom_errorbar(data = h.infinity, aes(x= t, ymin = `2.5%`, ymax =`97.5%`), width=5, size = 0.2)
  #                  ggsave(asympt.plot[[chi.i]][[mod]], filename =paste0(dir.regression.chi.mod,"/asymptote.png"),
  #                         bg = "transparent", width = 12, height =6, dpi = 200)
  #################################################################################################################
  # segmentation of the asymptotes (only one change point searched):
  seg.rec = NULL;
  read.res.rec=NULL;
  dir.create(paste0(dir.real.time, "/segmentation"))                    
  dir.segm.recessions = paste0(dir.real.time, "/segmentation")
  dir.segmentation <- paste0(dir_code,"/BaM_exe/Segmentation")
  param.var.model.rt = c( "a3")
  dir.segm.recessions.model = paste0(dir.segm.recessions,"/",rec.model.rt)
  dir.create(dir.segm.recessions.model)
  dir.segm.recessions.model.param = paste0(dir.segm.recessions.model,"/chi_",chi.rt)
  dir.create(dir.segm.recessions.model.param)
  dir.segm.recessions.param = paste0(dir.segm.recessions.model.param,"/",param.var.model.rt)
  dir.create(dir.segm.recessions.param)
  # dir.create(paste0(dir.c, "/segmentation"))
  # dir.c.segm = paste0(dir.c, "/segmentation")
  hasympt.new   = h.infinity
  hasympt.old   = df.asymptotes.rt.old
  #data.segm.rec = read.table(file=paste0(dir.segm.recessions.param,"/Data.segm.rec.txt"),    header=TRUE)
  # data.segm.rec1 = data.frame(  X2.5.   = hasympt.new$`2.5%`,
  #                               X50.    = hasympt.new$`50%`,
  #                               X97.5.  = hasympt.new$`97.5%`,
  #                               mean    = hasympt.new$mean,
  #                               stdev   = hasympt.new$stdev,
  #                               maxpost = hasympt.new$maxpost,
  #                               t       = hasympt.new$t)
  #hasympt$t =  h.infinity$t  + df.limni$t_limni[realtime.period[1]]
  data.segm.rec2  = rbind( hasympt.old, 
                           hasympt.new) #data.segm.rec1)
  write.table(data.segm.rec2, paste0(dir.segm.recessions.param,"/Segm_data.txt"),sep="\t",row.names=FALSE)
  
  seg.rec = recession.segmentation(
    dir.segmentation     = dir.segmentation,
    dir_code             = dir_code,
    dir.case_study       = dir.case_study,
    dir_exe              = dir.exe,
    dir.segm.recessions  = dir.segm.recessions.param,
    param.to.segment     = param.var.model.rt,
    all.param.to.segment = param.var.model.rt,
    file.dir.with.data   = paste0(dir.segm.recessions.param,"/Segm_data.txt"),
    Ncycle.segment       = 1500,
    Nmcmc.segment        = 100,
    Nslim.segment        = 50,
    Nburn.segment        = 0.5,
    df.limni             = df.limni.rt,
    nSmax                = 2,  # look for 1 shift only !!!!
    tmin                 = 1,
    criterion            = "DIC",
    limits.X             = c(limni.time.limits[1], limni.time.limits[2], 1000),
    limits.Y             = c(asymptote.limits[1], asymptote.limits[2], 50),
    x.name               = "Time [days]",
    y.name               = "Asymptotic stage [m]",
    obs.uncertainty.y    = TRUE,
    prior.mu.segm        = c("Uniform", asymptote.limits),
    gamma_prior          = c("Uniform", 0, 50, 0.1) ,
    save.all.results     = TRUE,
    floodpeakalways      = FALSE,
    MAPalways            = TRUE,
    FinalColors          = colo,
    gaugings             = data4BaRatin.rt,
    BayesianOption       = BayesianOption,
    stage.scale.shift    = stage.scale.shift
  )
  
  read.res.rec = read.results.segment.recess(dir.segm.recessions = dir.segm.recessions.param ,
                                             officialShiftsTime  = officialShiftsTime,
                                             Gaugings            = Gaugings,
                                             plot.dates          = FALSE)
  
  
  
  shift.true = TRUE
  if ((read.res.rec$nS.ok >1)& (shift.true == FALSE)){
    ####################################################
    
    print(paste0("shift detected after ", h.infinity$t , " days from flood peak"))
    shift.time.rt = h.infinity$t
    asympt.plot      =  asympt.plot +
      geom_vline(aes(xintercept = shift.time.rt), color="red", linetype = "dashed", size = 1)
    ggsave(asympt.plot, filename =paste0(dir.regression.chi.mod,"/asymptote.png"),
           bg = "transparent", width = 12, height =6, dpi = 200)
  } else {
    print(paste0("Recession is stable !"))
  }
  setwd(dir_code)
  
  
  
  #   # save the 3 frames that you want:
  #   if ((ttt == frames.time[1])|(ttt == frames.time[2])|(ttt == frames.time[3])) {
  #     ##############################################################################
  #     fram = fram+1
  #     reg.rt.plot[[fram]] =  ggplot() + 
  #       theme_bw(base_size = 25)+
  #       scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
  #                          expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
  #       scale_y_continuous(name = "Stage h (cm)", limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
  #       theme(plot.title = element_text(hjust = 0.5),
  #             panel.grid.major = element_blank(), 
  #             panel.grid.minor = element_blank(),
  #             panel.background = element_blank(),
  #             plot.margin=unit(c(0.5, 1, 0.5, 0.5),"cm")) +
  #       coord_cartesian(clip="off")+
  #       # geom_ribbon(data= PltData.rec,
  #       #             aes(x=t,
  #       #                 ymin=inf,
  #       #                 ymax= sup,
  #       #                 group=Period), color="green", alpha=0.1) +
  #       # geom_path(data= PltData.rec,
  #       #             aes(x=t, 
  #       #                 y = maxpost,
  #       #                 group=Period), color="green")+       
  #       geom_point(data=data.rec.obs.1,
  #                  aes(x=time, 
  #                      y = h- stage.scale.shift,
  #                      group=Period), color="gray90", size= 2)+
  #       geom_vline(xintercept = h.infinity$t, color="blue", linetype ="8f", size=0.3)+
  #       geom_point(data = temp , aes(x = x, y= y), color= "blue", size= 4) +
  #       geom_ribbon(data = temp2, aes(x = xx,    ymin = zz , ymax = kk), 
  #                   fill = "blue", alpha=0.1) +
  #       geom_line(data = temp2, aes(x = xx,  y= V1), color = "blue")+
  #       annotate("text", 
  #                x =  70,
  #                y = 200, 
  #                label = paste0("after  ", round(h.infinity$t,1)," days"), 
  #                color = "red",  fontface= 1,
  #                size=15)+
  #       annotate("text", 
  #                x =  40,
  #                y = 50, 
  #                label = "past \n recessions", 
  #                color = "gray80",  fontface= 1,
  #                size=8)+
  #       annotate("text", 
  #                x =  30,
  #                y = -70, 
  #                label = "actual \n recession", 
  #                color = "blue",  fontface= 1,
  #                size=8)
  #     
  ###############################################################################################       
  # asympt.rt.plot <- ggplot()+
  #                   theme_bw(base_size = 25)+
  #                   scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
  #                                      expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
  #                                      scale_y_continuous(name = TeX("Asymptotic level  $\\; \\beta$ (cm)"), 
  #                                      limits = c(limits.y[1], 100), expand = c(0,0)) +
  #                   theme(plot.title = element_text(hjust = 0.5),
  #                         panel.grid.major = element_blank(), 
  #                         panel.grid.minor = element_blank(),
  #                         panel.background = element_blank(),
  #                         plot.margin=unit(c(0.5, 1, 0.5, 0.5),"cm")) +    
  #                   geom_errorbar(data = h.infinity, aes(x= t, ymin = `2.5%`, ymax =`97.5%`),
  #                                 width=4, size = 1, color="blue")+     
  #                   geom_point(data = h.infinity, aes(x = t, y = maxpost),
  #                              colour = "blue",
  #                              size = 4) 
  # 
  # pdf(paste0(dir.regression.chi.mod,"/Figure_realtime_regression.pdf"), 12, 7,  useDingbats=F)
  # print(reg.plot[[chi.i]][[mod]])
  # dev.off()
  
  
  
  #       geom_vline(xintercept = h.infinity$t, color="blue", linetype ="8f", size=0.3)+
  #       annotate(geom = "rect",
  #                xmin = 0, xmax= 150,
  #                ymin = read.res.rec.old$mu.results.df$mu.q2, 
  #                ymax = read.res.rec.old$mu.results.df$mu.q97,
  #                fill= "gray60", alpha= 0.1)+
  #       geom_segment(aes(x= 0, xend=150, 
  #                        y= read.res.rec.old$mu.results.df$mu.mean, 
  #                        yend= read.res.rec.old$mu.results.df$mu.mean), 
  #                    color="gray80")+
  #       annotate("text", 
  #                x =  70,
  #                y = read.res.rec.old$mu.results.df$mu.mean + 12, 
  #                label =TeX("$\\beta$ of previous recessions (mean)"), 
  #                color = "gray80",  fontface= 1,
  #                size=8)
  #     if (read.res.rec$nS.ok >1) {
  #       asympt.rt.plot[[fram]] = asympt.rt.plot[[fram]]+
  #         annotate("text", 
  #                  x =  70,
  #                  y = -70, 
  #                  label ="Shift  \n detected !!!", 
  #                  color = "red",  fontface= 2,
  #                  size=7)
  ###########################################################################
  return(list(#Data.segm.rec           = Data.segm.rec, 
    temp                    = temp,
    temp2                   = temp2,
    Results.segmentation    = read.res.rec,
    df.asymptotes.rt.new    = hasympt.new))
  ###########################################################################
}





























#########################################################################################################
realtime.recess.2   <-  function(            t_limni, 
                                             h_limni, 
                                             rec.model.rt, 
                                             data.recess.realtime,
                                             which.recession.realtime,
                                             dir.case_study, dir_code, dir_exe, 
                                             station.name, data.period,  
                                             prior.rec.realtime, 
                                             prior.gamma.realtime,
                                             ncycl, ncycles.max, nmcmc, nslim, jump.pos, jump.neg, nburn,  
                                             pred, 
                                             chi.rt,
                                             delta.t.min.rt,
                                             delta.t.max.rt,
                                             uh.rec.rt,
                                             limits.y, limits.x,
                                             long.recessions,
                                             stage.scale.shift,
                                             starting.time,
                                             realtime.period,
                                             frames.time) {
  ##########################################################################################################
  start_time <- Sys.time()
  print("Uploading results of retrospective analysis on past recessions ...")
  
  # function for the regression OF RECESSIONS (separately):
  dir.create(paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression"))
  dir.create(paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression/Real_Time"))
  dir.regression <- paste0(dir.case_study,"/Results/segmentation_recessions/curves_regression/Real_Time")
  
  
  #selection of the period to study:
  ###################################
  t.limni.rt = df.limni$t_limni[realtime.period[1]:realtime.period[2]]
  h.limni.rt = df.limni$h_limni[realtime.period[1]:realtime.period[2]]
  min_loc.h.rt <- localmin(t.limni.rt, h.limni.rt)  #stage in centimeters !!!
  recess.h.rt <- rec(min_loc.h.rt$tmin, 
                     min_loc.h.rt$Qmin,  
                     chi=50, 
                     delta.t.max = 100)    #chi.rt/100)
  curves.h.rt.1 <- extract_curve(recess.h.rt$trec, recess.h.rt$Qrec)
  curves.h.rt.1$tcurve = curves.h.rt.1$tcurve  - curves.h.rt.1$tcurve[1]
  eliminate.index =0
  for (ccc in 2:length( curves.h.rt.1$tcurve)) {
    if ((curves.h.rt.1$tcurve[ccc] - curves.h.rt.1$tcurve[ccc-1]) < delta.t.min.rt) {
      eliminate.index = c(eliminate.index, ccc)
    }
  }
  curves.h.rt = curves.h.rt.1[-eliminate.index,]
  starting.time.index = which( curves.h.rt$tcurve >starting.time)[1] 
  reg.plot =NULL;
  asympt.plot=NULL;
  plot.bind =NULL;
  plot.bind.chi=NULL;
  
  # plot the points to test the dataset:
  plot(x=curves.h.rt$tcurve, y=curves.h.rt$Qcurve)
  
  
  
  
  # do a specific computation for the results of each model and of each chi value:  
  ########################################################################################################
  for (chi.i in 1:length(chi.rt)){
  ########################################################################################################
    print(paste0("CHI :  ", chi.rt[chi.i], "  #########################################################################################"))
    reg.plot[[chi.i]]= list()
    asympt.plot[[chi.i]] = list()
    plot.bind[[chi.i]] = list()
    dir.create(paste0(dir.regression,"/chi_",chi.rt[chi.i]))
    dir.regression.chi <- paste0(dir.regression,"/chi_",chi.rt[chi.i])
    
    # mu.results.df.rec = NULL; gamma.results.df.rec = NULL;  Q2.mu.rec  = NULL; mu.res.rec  = NULL; 
    # Q97.mu.rec  = NULL;  Data.segm.rec = NULL;  mcmc.seg.rec = NULL;
    # gamma_segm_recess = NULL;  tau.results.df.rec = NULL; df.shift.times.rec = NULL;
    # df.shift.times.plus.rec = NULL; ts.res.before.rec = NULL; ts.res.rec = NULL; ts.res.plus.rec = NULL; Q2.ts.rec = NULL;
    # Q97.ts.rec = NULL; nS.ok.rec = NULL; pdf.ts.rec = NULL;shift.times.recessions = NULL;
    # data.annotate.recess = NULL; data.annotate.recess.adjust = NULL;
    
    
    ########################################################################################################
    for (mod in 1:length(rec.model.rt)) {
    ########################################################################################################
      # read the priors from retrospective analysis:
      dir.rec.retrosp = paste0(dir.case_study,
                               "/Results/segmentation_recessions/curves_regression/Pooling/test_",
                               rec.model.rt[mod], "/chi_",
                               chi.rt[chi.i]
      )
      dir.segm.retrosp= paste0(dir.case_study,
                               "/Results/segmentation_recessions/segmentation/",
                               rec.model.rt[mod], "/chi_",
                               chi.rt[chi.i]
      )
      curves_retrospect = read.table(paste0(dir.rec.retrosp,
                                            "/Curves_Data.txt"), header =TRUE)
      Results_retrospect = read.table(paste0(dir.rec.retrosp,
                                             "/Results_Summary.txt"), header =TRUE)
      Nrecess = tail(curves_retrospect$Period,1)
      model.title =0
      
      # # recess model with tot uncertainty:
      # #*******************************************
      # Recession.model = function(theta, model, t){ 
      #   #*******************************************
      #   h=0*t
      #   if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      #     h = theta[1]*exp(-theta[2]*t) +
      #       theta[3]*exp(-theta[4]*t) + 
      #       theta[5]*exp(-theta[6]*t) +
      #       theta[7] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[8],theta[9]))
      #   } else if (model =="2expWithAsympt"){
      #     h = theta[1]*exp(-theta[2]*t) +
      #       theta[3]*exp(-theta[4]*t) + 
      #       theta[5] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[6],theta[7]))
      #   } else if (model =="2expWithAsympt_bis"){
      #     h = theta[1]*exp(-theta[2]*t) +
      #       theta[3]*exp(-theta[4]*t) + 
      #       theta[5] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[6],theta[7]))
      #   } else if (model =="2expWithAsympt_rel"){
      #     h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
      #       theta[5] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[6],theta[7]))
      #   } else if (model =="1expWithAsympt"){
      #     h = theta[1]*exp(-theta[2]*t) +
      #       theta[3] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[4],theta[5]))
      #   } else if ((model =="expexp")|(model =="expexp_bis")){
      #     h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[5],theta[6]))
      #   } else if ((model =="hyperb")|(model =="hyperb_bis")){
      #     h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[5],theta[6]))
      #   } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      #     h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
      #     res.h = sapply(h, 
      #                    function(h, theta){h + rnorm(1, mean=0, sd=theta[1] + theta[2]*h)}, 
      #                    theta=c(theta[5],theta[6]))
      #   } 
      #   return(res.h)
      # }
      # 
      # # recess model for maxpost:
      # #************************************************
      # Recess.Maxpost.model = function(theta, model, t){
      #   #************************************************
      #   h=0*t
      #   if ((model =="3expWithAsympt")|(model =="3expWithAsympt_bis")){
      #     h = theta[1]*exp(-theta[2]*t) +
      #       theta[3]*exp(-theta[4]*t) + 
      #       theta[5]*exp(-theta[6]*t) +
      #       theta[7] 
      #   }  else if (model =="2expWithAsympt"){
      #     h = theta[1]*exp(-theta[2]*t) +
      #       theta[3]*exp(-theta[4]*t) + 
      #       theta[5] 
      #   }  else if (model =="2expWithAsympt_bis"){
      #     h = theta[1]*exp(-theta[2]*t) +
      #       theta[3]*exp(-theta[4]*t) + 
      #       theta[5] 
      #   }  else if (model =="2expWithAsympt_rel"){
      #     h = theta[1]*(exp(-theta[2]*t) +  theta[3]*exp(-theta[4]*t)) + 
      #       theta[5] 
      #   } else if ((model =="1expWithAsympt")|(model =="1expWithAsympt_bis")){
      #     h = theta[1]*exp(-theta[2]*t) + theta[3] 
      #   } else if ((model =="expexp")|(model =="expexp_bis")){
      #     h = theta[1]*exp(-theta[2]*t^theta[3]) + theta[4] 
      #   } else if ((model =="hyperb")|(model =="hyperb_bis")){
      #     h = theta[1]/((1 + theta[2]*t)^theta[3]) + theta[4] 
      #   } else if ((model =="Coutagne")|(model =="Coutagne_bis")){
      #     h = theta[1]*(1 + (theta[3]-1)*theta[2]*t)^(theta[3]*(1 - theta[3])) + theta[4] 
      #   }
      #   return(h)
      # }
      # 
      
      #initialisation of the lists of plot objects:
      # title.model=NULL
      # dir.rec.segm.test.param=NULL;
      # tau.results.df.rec =NULL; mu.results.df.rec=NULL; gamma.results.df.rec=NULL; df.shift.times.rec=NULL
      # df.shift.times.plus.rec=NULL; ts.res.before.rec=NULL; ts.res.rec=NULL; ts.res.plus.rec=NULL; Q2.ts.rec=NULL
      # Q97.ts.rec=NULL; Q2.mu.rec=NULL; mu.res.rec=NULL; Q97.mu.rec=NULL; Data.segm.rec=NULL; nS.ok.rec=NULL;
      # mcmc.seg.rec=NULL; pdf.ts.rec=NULL; gamma_segm_recess=NULL; gaugings.df.recess=NULL;
      # shift.times.recessions=NULL; data.annotate.recess=NULL; data.annotate.recess.adjust=NULL; 
      # X=NULL; X1=NULL; data.tmp= NULL;data.tmp.2= NULL; Ysim=NULL; time.adjust.before =NULL; time.adjust.plus=NULL; parameters =NULL; parameters.names =NULL 
      # read data mcmc:  
      data.MCMC.cooked=as.matrix(read.table(paste0(dir.rec.retrosp,"/Results_MCMC_Cooked.txt"), header=TRUE,dec=".", sep=""))
      data.MCMC.MaxPost = as.numeric(read.table(paste0(dir.rec.retrosp,"/Results_Summary.txt"), row.names=1,dec=".",sep="", skip = 16))
      summary.rec   = read.table(file=paste0(dir.rec.retrosp,"/Results_Summary.txt"), header=TRUE)
      mcmc.rec      = read.table(file=paste0(dir.rec.retrosp,"/Results_MCMC_cooked.txt"), header=TRUE)
      residuals.rec = read.table(file=paste0(dir.rec.retrosp,"/Results_Residuals.txt"), header=TRUE)
      curves.data.rec = read.table(file=paste0(dir.rec.retrosp,"/Curves_Data.txt"), header=TRUE)
      nsample = length(data.MCMC.cooked[,1])
      tgrid = seq(limits.x.recess[1],  limits.x.recess[2],  0.5) 
      Ncurves.pool = tail(curves.data.rec$Period,1)
      
      #Initialisation:
      if (rec.model.rt[mod] =="3expWithAsympt"){
        #######################################
        b1.var =FALSE
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a4", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] -2
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])

        MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=10) # 9 param + # of period
        MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=10)      # 9 param + # of period
        for(i in 1:Ncurves.pool){
          # col.num = c( i, Ncurves.pool+1,                                  #a1, b1,
          #              Ncurves.pool +1+i, Ncurves.pool*2 +1+1,            #a2, b2
          #              Ncurves.pool*2 +1+1+ i, Ncurves.pool*3+1+1+1,       #a3, b3,
          #              Ncurves.pool*3 +1+1+1+ i,                            #a4 (asymptotic level parameter)
          #              Ncurves.pool*4 + 1 +1+1+1, Ncurves.pool*4 +1+1+1+2)  #gamma1, gamma2
          col.num = c( i, Ncurves.pool+1,                                  #a1(var), b1,
                       Ncurves.pool +2, Ncurves.pool + 3,                  #a2, b2
                       Ncurves.pool +4, Ncurves.pool + 5,                  #a3, b3,
                       Ncurves.pool +5 + i,                                #a4(var) (asymptotic level parameter)
                       Ncurves.pool*2 + 5+1, Ncurves.pool*2 + 5+2)         #gamma1, gamma2
          MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
          MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
        }

      } else if (rec.model.rt[mod] =="3expWithAsympt_bis"){
        #####################################################
        b1.var =FALSE
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a4", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 2
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])

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

      }  else if (rec.model.rt[mod] =="2expWithAsympt_bis"){
        ####################################################
        b1.var =FALSE
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a3", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 2
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])
        MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=8) # 7 param + # of period
        MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=8)      # 7 param + # of period
        for(i in 1:Ncurves.pool){
          col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                       Ncurves.pool +2 + i, Ncurves.pool*2 +1+1,   #a2, b2
                       Ncurves.pool*2 +1+1+i,                       #a3,
                       Ncurves.pool*3 + 2 + 1, Ncurves.pool*3 +2+2)  #gamma1, gamma2
          MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
          MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
        }

      }  else if ((rec.model.rt[mod] =="2expWithAsympt")|(rec.model.rt[mod] =="2expWithAsympt_rel")){
        #################################################################################################
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a3", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 2
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])
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


      } else if (rec.model.rt[mod] =="1expWithAsympt"){
        ################################################
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a2", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 2
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])
        MCMC.save    = matrix(NA, nrow=Ncurves.pool*nsample, ncol=6) # 5 param + # of period
        MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=6)      # 5 param + # of period
        for(i in 1:Ncurves.pool){
          col.num = c( i, Ncurves.pool+1,                      #a1, b1,
                       Ncurves.pool +1+ i,                      #a2
                       Ncurves.pool*2 +1+1, Ncurves.pool*2 +1 + 2)  #gamma1, gamma2
          MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
          MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
        }


      } else if (rec.model.rt[mod] =="expexp"){
        #######################################
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a2", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 1
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])
        MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
        MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period
        for(i in 1:Ncurves.pool){
          col.num = c( i, Ncurves.pool + 1,                      #a1, b1,
                       Ncurves.pool + 2,                      #n1
                       Ncurves.pool + 2 + i,                  #a2 (asymptotic level parameter)
                       Ncurves.pool*2 + 2 +1, Ncurves.pool*2 +2 +2)  #gamma1, gamma2
          MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
          MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
        }

      } else if ((rec.model.rt[mod] =="expexp_bis")|(rec.model.rt[mod] =="hyperb_bis")|(rec.model.rt[mod] =="Coutagne_bis")){
        #####################################################################################################################
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a2", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))

        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 1
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])

        MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
        MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period
        for(i in 1:Ncurves.pool){
          col.num = c( i, Ncurves.pool + i,                   #a1, b1,
                       Ncurves.pool*2 + 1,                      #n1
                       Ncurves.pool*2 + 1 + i,                  #a2 (asymptotic level parameter)
                       Ncurves.pool*3 + 2, Ncurves.pool*3 +3)  #gamma1, gamma2
          MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
          MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
        }

      } else if (rec.model.rt[mod] =="hyperb"){
        #######################################
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a2", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 2
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])
        MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
        MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period
        for(i in 1:Ncurves.pool){
          col.num = c( i, Ncurves.pool + 1,                      #a1, b1,
                       Ncurves.pool + 2,                      #n1
                       Ncurves.pool + 2 + i,                  #a2 (asymptotic level parameter)
                       Ncurves.pool*2 + 2 + 1, Ncurves.pool*2 +2 +2)  #gamma1, gamma2
          MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
          MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
        }



      } else if (rec.model.rt[mod] =="Coutagne"){
        ##########################################
        Data.segm.rec = read.table(file=paste0(dir.segm.retrosp,"/a2", "/Data.segm.rec.txt"), header=TRUE)
        mcmc.new     = rnorm(10000, mean= mean(Data.segm.rec$mean), sd=mean(Data.segm.rec$stdev))
        nRec = which(Data.segm.rec$t > df.limni$t_limni[realtime.period[1]])[1] - 2
        asymp.old.min = min(Data.segm.rec$X2.5.[1:nRec])
        asymp.old.max = max(Data.segm.rec$X97.5.[1:nRec])
        asymp.old.mean = mean(Data.segm.rec$mean[1:nRec])
        MCMC.save = matrix(NA, nrow=Ncurves.pool*nsample, ncol=7) # 6 param + # of period
        MaxPost.save = matrix(NA, nrow=Ncurves.pool, ncol=7)      # 6 param + # of period
        for(i in 1:Ncurves.pool){
          col.num = c( i, Ncurves.pool + 1,                   #a1, b1,
                       Ncurves.pool + 2,                      #n1
                       Ncurves.pool + 2 + i,                  #a2 (asymptotic level parameter)
                       Ncurves.pool*2 + 2+1, Ncurves.pool*2 +2 +2)  #gamma1, gamma2
          MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
          MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num],i)
        }
      }
      # ########################################################################################################
      # # Apply recession model (total uncertainty and maxpost):
      # Rec.Post    = apply(MCMC.save,    MARGIN=1, Recession.model,      model=rec.model.rt,  t=tgrid) #add structural error:
      # Rec.MaxPost = apply(MaxPost.save, MARGIN=1, Recess.Maxpost.model, model=rec.model.rt,  t=tgrid) # Maximum posterior 
      # # Quantiles: 
      # message("Plotting all Recession curves !!!  Wait ... "); flush.console()
      # List.Rec.quants = list(NULL)
      # for(i in 1:nRec){
      #   data.tmp = apply(Rec.Post[,(nsample*(i-1)+1):(nsample*i)], 
      #                    MARGIN=1, quantile, probs=c(0.025, 0.975),  na.rm=TRUE)
      #   List.Rec.quants[[i]] = data.frame(cbind(tgrid, 
      #                                           t(data.tmp - stage.scale.shift),
      #                                           Rec.MaxPost[,i] - stage.scale.shift))  #!!!!!!!!!!!!!! 
      #   colnames(List.Rec.quants[[i]]) = c("t", "inf", "sup", "maxpost")
      # }
      # #prepare plot:
      # inter.per=seq(1,nRec,1) 
      # data.plot.Rec = data.frame(do.call("rbind", 
      #                                    List.Rec.quants[inter.per]), 
      #                            Period=rep(inter.per,each=length(tgrid)))
      # inter.null=which(data.plot.Rec$maxpost==0)
      data.rec.obs = read.table(paste0(dir.rec.retrosp,"/Curves_Data.txt"),
                                header=TRUE,dec=".",sep="") # Gaugings loading
      # pos.num =function(x.int){
      #   inter=which(x.int==data.rec.obs$period);
      #   return(inter)}
      # #inter.rec.obs =unlist(sapply(inter.per, pos.num), recursive = TRUE, use.names = TRUE)
      # data.Rec = data.plot.Rec #[-inter.null,]
      # data.Rec$Period = factor(data.Rec$Period)
      # #data.rec.obs$Period = as.factor(data.rec.obs$Period)
      # #write.table( data.Rec, paste0(dir.rec.pool.test,"/Rec_SPD_env.txt"), sep ="\t", row.names=FALSE)
      # 
      data.Rec = read.table(paste0(dir.rec.retrosp,"/Rec_SPD_env.txt"), header =TRUE)
      PltData.rec.1 <- data.Rec #data.Rec[data.Rec$inf > -200,] 
      PltData.rec = PltData.rec.1[1:(nRec*length(tgrid)),]
      data.rec.obs.1 = data.rec.obs[which(data.rec.obs$Period< nRec),]
      dir.create(paste0(dir.regression.chi,"/", rec.model.rt[mod]))
      dir.regression.chi.mod <- paste0(dir.regression.chi,"/", rec.model.rt[mod])
      
      # plot:
      ###########################################################################
      reg.plot[[chi.i]][[mod]] <-  ggplot() + 
        theme_bw(base_size = 15)+
        scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
                           expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
        scale_y_continuous(name = "Stage h (cm)", limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.margin=unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
        geom_ribbon(data= PltData.rec,
                  aes(x=t,
                      ymin=inf,
                      ymax= sup,
                      group=Period),
                  fill="green", alpha =0.2,size=0.1)
      # geom_smooth(data= PltData.rec, 
      #             aes(x=t,
      #                 y=maxpost,
      #                 ymax=sup,
      #                 ymin=inf),
      #             fill="green", color="green",
      #             size=0.1, stat='identity', alpha=0.1)
      # geom_linerange(data= PltData.rec,
      #                aes(x=time, 
      #                    ymax= (h + 2*uh) - stage.scale.shift, 
      #                    ymin=(h-2*uh) - stage.scale.shift), color = "green",
      #                data=data.rec.obs, size=0.1)
      # geom_point(aes(x=time, 
      #                y=h-stage.scale.shift), color ="green", 
      #                data=data.rec.obs, shape=16, size=1)+
      
      # title.model[[mod]] <- ggdraw() + 
      #   draw_label( model.title[mod],
      #               x = 0, hjust = 0,
      #               size = 15) +
      #   theme(plot.margin = margin(0, 0, 0, 7))
      #   reg.pool.plot2[[mod]] =  plot_grid(title.model[[mod]],
      #                                    reg.pool.plot[[mod]] ,  
      #                                    ncol = 1,
      #                                    # rel_heights values control vertical title margins
      #                                    rel_heights = c(0.2, 1))
      
      #axis.line = element_line(colour = "black"))            
      ggsave(reg.plot[[chi.i]][[mod]], filename =paste(dir.regression.chi.mod,"/regression_realtime.png", sep=""), 
             bg = "transparent", width = 12, height =7, dpi = 200)
      
      asympt.plot[[chi.i]][[mod]] <- ggplot()+
        theme_bw(base_size = 15)+
        scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
                           expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
        scale_y_continuous(name = TeX("Asymptotic level  $\\; \\beta$ (cm)"), limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.margin=unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
        annotate(geom = "rect",
                 xmin = 0, xmax= 150,
                ymin = asymp.old.min, ymax= asymp.old.max,
                 fill= "green", alpha= 0.2)+
        geom_segment(aes(x= 0, xend=150, y= asymp.old.mean, yend=asymp.old.mean), color="green")
      
      
      ggsave(asympt.plot[[chi.i]][[mod]], filename =paste0(dir.regression.chi.mod,"/asymptote.png"),
             bg = "transparent", width = 12, height =6, dpi = 200)
      
      if (rec.model.rt[mod] == "3expWithAsympt"){
      #############################################
        priors.rt = c(prior.rec.realtime[[mod]][1:4], "static",
                      Results_retrospect[5, Nrecess+1], Results_retrospect[11, Nrecess+1], "Gaussian", Results_retrospect[5, Nrecess+1], "static" ,
                      
                      Results_retrospect[5,Nrecess+2], Results_retrospect[11,Nrecess+2], "Gaussian", Results_retrospect[5,Nrecess+2], "static",
                      Results_retrospect[5,Nrecess+3], Results_retrospect[11,Nrecess+3], "Gaussian", Results_retrospect[5,Nrecess+3], "static",
                      
                      Results_retrospect[5,Nrecess+4], Results_retrospect[11,Nrecess+4], "Gaussian", Results_retrospect[5,Nrecess+4], "static",
                      Results_retrospect[5,Nrecess+5], Results_retrospect[11,Nrecess+5], "Gaussian", Results_retrospect[5,Nrecess+5], "static",
                      
                      prior.rec.realtime[[mod]][31:34], "static"
        )
        priors.gamma.rt = c(Results_retrospect[5,2*Nrecess+6], Results_retrospect[11,2*Nrecess+6], "Gaussian", Results_retrospect[5,2*Nrecess+6],
                            Results_retrospect[5,2*Nrecess+7], Results_retrospect[11,2*Nrecess+7], "Gaussian", Results_retrospect[5,2*Nrecess+7])
        
      } else if ((rec.model.rt[mod] == "2expWithAsympt")|(rec.model.rt[mod] == "2expWithAsympt_rel")){
      #################################################################################################
        priors.rt = c(prior.rec.realtime[[mod]][1:4], "static",
                      Results_retrospect[5, Nrecess+1], Results_retrospect[11, Nrecess+1], "Gaussian", Results_retrospect[5, Nrecess+1], "static" ,
                      
                      Results_retrospect[5,Nrecess+2], Results_retrospect[11,Nrecess+2], "Gaussian", Results_retrospect[5,Nrecess+2], "static",
                      Results_retrospect[5,Nrecess+3], Results_retrospect[11,Nrecess+3], "Gaussian", Results_retrospect[5,Nrecess+3], "static",
                      
                      prior.rec.realtime[[mod]][21:24], "static"
        )
        priors.gamma.rt = c(Results_retrospect[5,2*Nrecess+4], Results_retrospect[11,2*Nrecess+4], "Gaussian", Results_retrospect[5,2*Nrecess+4],
                            Results_retrospect[5,2*Nrecess+5], Results_retrospect[11,2*Nrecess+5], "Gaussian", Results_retrospect[5,2*Nrecess+5])
        
      } else if ((rec.model.rt[mod] == "1expWithAsympt")){
      ###################################################
        priors.rt = c(prior.rec.realtime[[mod]][1:4], "static",
                      Results_retrospect[5, Nrecess+1], Results_retrospect[11, Nrecess+1], "Gaussian", Results_retrospect[5, Nrecess+1], "static" ,
                      prior.rec.realtime[[mod]][11:14], "static"
        )
        priors.gamma.rt = c(Results_retrospect[5,2*Nrecess+2], Results_retrospect[11,2*Nrecess+2], "Gaussian", Results_retrospect[5,2*Nrecess+2],
                            Results_retrospect[5,2*Nrecess+3], Results_retrospect[11,2*Nrecess+3], "Gaussian", Results_retrospect[5,2*Nrecess+3])
        
      } else if ((rec.model.rt[mod] == "2expWithAsympt_bis")){
        ######################################################
        priors.rt = c(prior.rec.realtime[[mod]][1:4], "static",
                      Results_retrospect[5, Nrecess+1], Results_retrospect[11, Nrecess+1], "Gaussian", Results_retrospect[5, Nrecess+1], "static" ,
                      
                      prior.rec.realtime[[mod]][11:14], "static",
                      Results_retrospect[5,2*Nrecess+2], Results_retrospect[11,2*Nrecess+2], "Gaussian", Results_retrospect[5,2*Nrecess+2], "static",
                      
                      prior.rec.realtime[[mod]][21:24], "static"
        )
        priors.gamma.rt = c(Results_retrospect[5,3*Nrecess + 3], Results_retrospect[11,3*Nrecess+3], "Gaussian", Results_retrospect[5,3*Nrecess+3],
                            Results_retrospect[5,3*Nrecess + 4], Results_retrospect[11,3*Nrecess+4], "Gaussian", Results_retrospect[5,3*Nrecess+4])
      }
      
      
      #########################
      data.segm.rec = read.table(file=paste0(dir.segm.recessions.param,"/Data.segm.rec.txt"),    header=TRUE)
      data.segm.rec2 = data.segm.rec[1:nRec, ]
      dir.create(paste0(dir.regression.chi.mod,"/segment_old_recessions"))
      dir.segment.old.recess <- paste0(dir.regression.chi.mod,"/segment_old_recessions")
      write.table(data.segm.rec2, paste0(dir.segment.old.recess,"/Segm_data.txt"),sep="\t",row.names=FALSE)
      seg.rec = recession.segmentation(
        dir.segmentation   = dir.segmentation, 
        dir_code            = dir_code,
        dir.case_study      = dir.case_study, 
        dir_exe             = dir.exe, 
        dir.segm.recessions = dir.segment.old.recess,
        param.to.segment    = param.var.model,
        all.param.to.segment = param.var.model,
        file.dir.with.data  = paste0(dir.segment.old.recess,"/Segm_data.txt"), 
        Ncycle.segment      = 1500,
        Nmcmc.segment       = 100,
        Nslim.segment       = 50,
        Nburn.segment       = 0.5,
        df.limni            = df.limni,
        nSmax               = 1,       # 5, change this for each case study !!!!!!! 
        tmin                = 1, 
        criterion           = "DIC", 
        limits.X            = c(limni.time.limits[1], limni.time.limits[2], 1000), 
        limits.Y            = c(asymptote.limits[1], asymptote.limits[2], 50),
        x.name              = "Time [days]", 
        y.name              = "Asymptotic stage [m]",
        obs.uncertainty.y   = TRUE,
        prior.mu.segm       = c("Uniform", asymptote.limits),
        gamma_prior         = c("Uniform", 0, 50, 0.1) ,
        save.all.results    = TRUE,
        floodpeakalways     = FALSE,
        MAPalways           = TRUE,
        FinalColors         = colo,
        gaugings            = data4BaRatin,
        BayesianOption      = BayesianOption,
        stage.scale.shift   = stage.scale.shift)
      read.res.rec.old  = read.results.segment.recess(dir.segm.recessions = dir.segment.old.recess , 
                                                      officialShiftsTime  = officialShiftsTime, 
                                                      Gaugings            = Gaugings,
                                                      plot.dates          = FALSE)
      setwd(dir_code)
      d.h = NULL
      # rt.period = seq(realtime.period[1], realtime.period[2], 10)
      message("REAL TIME Recession estimation - using BaM ....!!!"); flush.console()
      shift.true =FALSE
      reg.rt.plot =NULL
      asympt.rt.plot = NULL
      combined.plot.rt=NULL
      fram=0
      # REAL TIME LOOP:  
      ################################################################################################
      for (ttt in  frames.time){ #    #starting.time.index:length(curves.h.rt$tcurve)){
      ################################################################################################
        icurve = ttt - (starting.time.index -1)
        d.h[[icurve]] = data.frame(t = curves.h.rt$tcurve[1:ttt],
                                   h = curves.h.rt$Qcurve[1:ttt]*100 +
                                       stage.scale.shift,
                                   uh = uh.rec.rt)
        setwd(dir.exe)
        print(icurve)
        dir.create(paste0(dir.regression.chi.mod,"/C",icurve))
        dir.c = paste0(dir.regression.chi.mod,"/C",icurve)
        dir.create(paste0(dir.regression.chi.mod,"/C",icurve,"/BaM"))
        dir.rec = paste0(dir.regression.chi.mod,"/C",icurve,"/BaM")
        dir.BaM.rec = paste0(dir_code,"/BaM_exe/Recession_h")  
        nobs.h =length(d.h[[icurve]]$t)
        curve_data.h = "Recession_h/Curves_Data.txt"
        write.table(d.h[[icurve]], 
                    file=curve_data.h, append = FALSE, sep = "\t", eol = "\n", 
                    na = "NA", dec = ".", row.names = FALSE, col.names=c("time", "h", "uh"))
        BaM_config1.h(  tgrid    = seq(limits.x[1], limits.x[2], 0.5),
                        nobs     = nobs.h, 
                        ncycles  = ncycl, 
                        nmcmc    = nmcmc, 
                        nslim    = nslim, 
                        jump.pos = jump.pos,  jump.neg=jump.neg,  nburn=nburn,
                        model    = rec.model.rt[mod], 
                        pred     = pred, 
                        asympt.before = curves.h$Qcurve[index.h[icurve]-1], 
                        remnant  ="Linear",
                        prior    =     priors.rt ,
                        prior.gamma =  priors.gamma.rt )
        # prior    =     prior.rec.realtime[[mod]] ,
        #prior.gamma =  prior.gamma.realtime)
        #Launch BAM exe:
        system2("BaM_recession_multi_model_final.exe", stdout =NULL, stderr =NULL)
        #save results files:
        list.of.files.rec <- c(paste0(dir.BaM.rec,"/Results_MCMC_Cooked.txt"),
                               paste0(dir.BaM.rec,"/Results_Residuals.txt"),
                               paste0(dir.BaM.rec,"/Results_Summary.txt"),
                               paste0(dir.BaM.rec,"/Config_Model.txt"),
                               paste0(dir.BaM.rec,"/Curves_Data.txt"),
                               paste0(dir.BaM.rec,"/tgrid.txt"),
                               paste0(dir.BaM.rec,"/ht_Maxpost.spag"),
                               paste0(dir.BaM.rec,"/ht_ParamU.spag"),
                               paste0(dir.BaM.rec,"/ht_ParamU.env"),
                               paste0(dir.BaM.rec,"/ht_TotalU.env"),
                               paste0(dir.BaM.rec,"/ht_TotalU.spag"))
        for (i in 1:length(list.of.files.rec)) {
          file.copy(list.of.files.rec[i], dir.rec ,overwrite = TRUE)
        }
        plot.mcmc(workspace= dir.rec, iter= icurve)  #plot density and traceplots of MCMC
        #conv = Convergence.test(dir.seg = dir.rec ,  npar=7 , dir.plot = dir.rec)
        setwd(dir.exe)
        # if (conv==FALSE) {
        #   print(paste0("recess ",icurve," Not converged !!! "))
        # }
        Results_residuals.rt = read.table("Recession_h/Results_Residuals.txt", header =TRUE)
        residuals.plot = ggplot(Results_residuals.rt)+
          geom_point(aes(x = Results_residuals.rt$X1_obs, 
                         y = Results_residuals.rt$Y1_sim), color = "red", size = 2)+
          geom_point(aes(x = Results_residuals.rt$X1_obs, 
                         y = Results_residuals.rt$Y1_obs), color = "black", size=1)+
          theme_bw()
        ggsave(residuals.plot, filename =paste(dir.rec,"/residuals.png", sep=""),
               bg = "transparent", width = 8, height =4, dpi = 200)
        
        #reading results of BaM:
        tgrid.file <- paste0(dir.exe,"/Recession_h/tgrid.txt")
        maxpost = read.table(paste0(dir.BaM.rec,"/ht_Maxpost.spag"), header =FALSE)
        env     = read.table(paste0(dir.BaM.rec,"/ht_TotalU.env"), header =TRUE)
        env.param = read.table(paste0(dir.BaM.rec,"/ht_ParamU.env"), header =TRUE)
        RecCurve.h = read.table(paste0(dir.BaM.rec,"/Curves_Data.txt"), header =TRUE)
        summary.rec = read.table(paste0(dir.BaM.rec,"/Results_Summary.txt"), header =TRUE)
        mcmc.rec = read.table(paste0(dir.BaM.rec,"/Results_MCMC_Cooked.txt"), header =TRUE)
        
        
        
        
        if ((rec.model.rt[mod] =="3expWithAsympt")| (rec.model.rt[mod] =="3expWithAsympt_bis")){
          ###############################################################################################     
          # read single parameters results :   a1(k) , b1, a2, b2, a3, b3, a4(k) 
          # a1:
          a1.mcmc = mcmc.rec[,1]
          a1.summary = summary.rec[,1]
          a1.df = data.frame(  t(c(quantile(a1.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a1.mcmc), 
                                   stdev =std(a1.mcmc), 
                                   maxpost = a1.summary[16])))
          names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #b1:
          b1.mcmc = mcmc.rec[,2]
          b1.summary = summary.rec[,2]
          b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(b1.mcmc), 
                                   stdev =std(b1.mcmc), 
                                   maxpost = b1.summary[16])))
          names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #a2:
          a2.mcmc = mcmc.rec[,3]
          a2.summary = summary.rec[,3]
          a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a2.mcmc), 
                                   stdev =std(a2.mcmc),
                                   maxpost = a2.summary[16])))
          names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          
          #b2:
          b2.mcmc = mcmc.rec[,4]
          b2.summary = summary.rec[,4]
          b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)),
                                   mean=mean(b2.mcmc), 
                                   stdev =std(b2.mcmc),
                                   maxpost = b2.summary[16])))
          names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          
          
          #a3:
          a3.mcmc = mcmc.rec[,5]
          a3.summary = summary.rec[,5]
          a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a3.mcmc), 
                                   stdev =std(a3.mcmc),
                                   maxpost = a3.summary[16])))
          names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          
          #b3:
          b3.mcmc = mcmc.rec[,6]
          b3.summary = summary.rec[,6]
          b3.df = data.frame(  t(c(quantile(b3.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(b3.mcmc), 
                                   stdev =std(b3.mcmc), 
                                   maxpost = b3.summary[16])))
          names(b3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          
          #a4: asymptote !!!
          a4.mcmc = mcmc.rec[,7]
          a4.summary = summary.rec[,7]
          a4.df = data.frame(  t(c(quantile(a4.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a4.mcmc), 
                                   stdev =std(a4.mcmc), 
                                   maxpost = a4.summary[16])))
          names(a4.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          a4.df[, c(1,2,3,4,6)] = a4.df[,c(1,2,3,4,6)]  -  stage.scale.shift
          h.infinity = cbind(a4.df ,   t =  curves.h.rt$tcurve[ttt])
          
          
          
          
          
        } else if ((rec.model.rt[mod] =="2expWithAsympt")|(rec.model.rt[mod] =="2expWithAsympt_bis")|(rec.model.rt[mod] =="2expWithAsympt_rel")){
          ################################################################################################# 
          # a1:
          a1.mcmc = mcmc.rec[,1]
          a1.summary = summary.rec[,1]
          a1.df = data.frame(  t(c(quantile(a1.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a1.mcmc), 
                                   stdev =std(a1.mcmc), 
                                   maxpost = a1.summary[16])))
          names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #b1:
          b1.mcmc = mcmc.rec[,2]
          b1.summary = summary.rec[,2]
          b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(b1.mcmc), 
                                   stdev =std(b1.mcmc), 
                                   maxpost = b1.summary[16])))
          names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #a2:
          a2.mcmc = mcmc.rec[,3]
          a2.summary = summary.rec[,3]
          a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a2.mcmc), 
                                   stdev =std(a2.mcmc), maxpost = a2.summary[16])))
          names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          
          #b2:
          b2.mcmc = mcmc.rec[,4]
          b2.summary = summary.rec[,4]
          b2.df = data.frame(  t(c(quantile(b2.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(b2.mcmc), 
                                   stdev =std(b2.mcmc), maxpost = b2.summary[16])))
          names(b2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          
          
          #a3:
          a3.mcmc = mcmc.rec[,5]
          a3.summary = summary.rec[,5]
          a3.df = data.frame(  t(c(quantile(a3.mcmc, p = c(0.025, 0.5, 0.975)),
                                   mean=mean(a3.mcmc), 
                                   stdev =std(a3.mcmc), maxpost = a3.summary[16])))
          names(a3.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          a3.df[, c(1,2,3,4,6)] = a3.df[,c(1,2,3,4,6)]  -  stage.scale.shift
          h.infinity = cbind(a3.df ,   t =  curves.h.rt$tcurve[ttt])
          
          
          
        } else if (rec.model.rt[mod] =="1expWithAsympt"){
          ################################################################################################# 
          # read single parameters results :  a1 , b1,  a2(k) 
          # a1:
          a1.mcmc = mcmc.rec[,1]
          a1.summary = summary.rec[,1]
          a1.df = data.frame(  t(c(quantile(a1.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a1.mcmc), 
                                   stdev =std(a1.mcmc), maxpost = a1.summary[16])))
          names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #b1:
          b1.mcmc = mcmc.rec[,2]
          b1.summary = summary.rec[,2]
          b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                                   stdev =std(b1.mcmc), maxpost = b1.summary[16])))
          names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #a2:
          a2.mcmc = mcmc.rec[,3]
          a2.summary = summary.rec[,3]
          a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc), 
                                   stdev =std(a2.mcmc), maxpost = a2.summary[16])))
          names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          a2.df[, c(1,2,3,4,6)] = a2.df[,c(1,2,3,4,6)]  -  stage.scale.shift
          h.infinity = cbind(a2.df ,   t =  curves.h.rt$tcurve[ttt])
          
          
          
        } else if ((rec.model.rt[mod] == "expexp")|(rec.model.rt[mod] == "expexp_bis")|(rec.model.rt[mod] == "hyperb")| 
                   (rec.model.rt[mod] == "Coutagne_bis")|(rec.model.rt[mod] == "hyperb_bis")|(rec.model.rt[mod] =="Coutagne")) {
          #################################################################################################
          # read single parameters results : a1, b1, n1,  a2
          # a1:
          a1.mcmc = mcmc.rec[,1]
          a1.summary = summary.rec[,1]
          a1.df = data.frame(  t(c(quantile(a1.mcmc, p = c(0.025, 0.5, 0.975)), 
                                   mean=mean(a1.mcmc), 
                                   stdev =std(a1.mcmc), maxpost = a1.summary[16])))
          names(a1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #b1:
          b1.mcmc = mcmc.rec[,2]
          b1.summary = summary.rec[,2]
          b1.df = data.frame(  t(c(quantile(b1.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(b1.mcmc), 
                                   stdev =std(b1.mcmc), maxpost = b1.summary[16])))
          names(b1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          
          #n1:
          n1.mcmc = mcmc.rec[, 3]
          n1.summary = summary.rec[, 3]
          n1.df = data.frame(  t(c(quantile(n1.mcmc, p = c(0.025, 0.5, 0.975)),
                                   mean=mean(n1.mcmc), 
                                   stdev =std(n1.mcmc), 
                                   maxpost = n1.summary[16])))
          names(n1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          #a2:
          a2.mcmc = mcmc.rec[,4]
          a2.summary = summary.rec[,4]
          a2.df = data.frame(  t(c(quantile(a2.mcmc, p = c(0.025, 0.5, 0.975)), mean=mean(a2.mcmc), 
                                   stdev =std(a2.mcmc), maxpost = a2.summary[16])))
          names(a2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
          a2.df[, c(1,2,3,4,6)] = a2.df[,c(1,2,3,4,6)]  -  stage.scale.shift
          h.infinity = cbind(a2.df ,   t =  curves.h.rt$tcurve[ttt])
          
        }    
        #############################################################################################
        ttttt = seq(limits.x[1],limits.x[2], 0.5)
        temp <- data.frame(x = RecCurve.h$time - RecCurve.h$time[1], 
                           y = RecCurve.h$h - stage.scale.shift, 
                           z = RecCurve.h$uh)
        temp2 <- data.frame(xx = ttttt , 
                            yy = maxpost - stage.scale.shift, 
                            zz = env[,2] - stage.scale.shift,
                            kk = env[,3] - stage.scale.shift)
        # temp3 <- data.frame(xx = ttttt , 
        #                     yy = h_bam, 
        #                     zz = env.param[,2] - stage.scale.shift, 
        #                     kk = env.param[,3] - stage.scale.shift)
        
        write.table(temp, file=paste0(dir.rec, "/temp.txt"), append = FALSE, sep = "\t", eol = "\n", 
                    na = "NA", dec = ".", row.names = FALSE, col.names=c("x", "y", "z"))
        write.table(temp2, file=paste0(dir.rec, "/temp2.txt"), append = FALSE, sep = "\t", eol = "\n", 
                    na = "NA", dec = ".", row.names = FALSE, col.names=c("xx", "yy", "zz", "kk"))
        # write.table(temp3, file=paste0(dir.rec, "/temp3.txt"), append = FALSE, sep = "\t", eol = "\n", 
        #             na = "NA", dec = ".", row.names = FALSE, col.names=c("xx", "yy", "zz", "kk"))
        
        ############################################################################################
        reg.plot[[chi.i]][[mod]] =  reg.plot[[chi.i]][[mod]]  +
          # scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
          #                    expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
          # geom_point(data = temp , aes(x = x, y= y), color=colfunc(80)[curve_good.h] , size= 0.5)+
          # geom_line(data = temp2, aes(x=xx, y=yy), color= "gray", size= 0.1) + #colfunc(Ncurves)[icurve], size = 0.1)+
          geom_point(data = temp , aes(x = x, y= y), color= "black", size= 2)+  #colfunc(Ncurves)[icurve], size= 2)+
          geom_ribbon(data = temp2, aes(x = xx,    ymin = zz , ymax = kk), 
                      fill = "blue", alpha=0.1)  # colfunc(Ncurves)[icurve], alpha = 0.1)
        # geom_ribbon(data = temp3, aes(x = xx, ymin = zz , ymax = kk), 
        #             fill = colfunc(Ncurves)[icurve], alpha = 0.3)
        ggsave(reg.plot[[chi.i]][[mod]], filename =paste0(dir.regression.chi.mod,"/regression_realtime.png"),
               bg = "transparent", width = 12, height =7, dpi = 200)
        ##############################################################################################
        asympt.plot[[chi.i]][[mod]] =  asympt.plot[[chi.i]][[mod]] +
                                       geom_point(data = h.infinity, aes(x = t, y = maxpost), colour = "black",size = 3) +
          #geom_boxplot(boxplot.df, aes(x = x, ymin=y0, lower =y2, middle= y50,
          #upper=y97, ymax= y100), stat = "identity") +
          #geom_line(data = asym.df.temp[[i]], aes(x =t, y =h), colour = "gray",size = 1) +
          geom_errorbar(data = h.infinity, aes(x= t, ymin = `2.5%`, ymax =`97.5%`), width=5, size = 0.2)
        ggsave(asympt.plot[[chi.i]][[mod]], filename =paste0(dir.regression.chi.mod,"/asymptote.png"),
               bg = "transparent", width = 12, height =6, dpi = 200)
        ###############################################################################################
        # segmentation of the asymptotes (only one change point searched):
        seg.rec = NULL; read.res.rec=NULL;
        dir.segm.recessions = paste0(dir.case_study,"/Results/segmentation_recessions/segmentation")
        dir.segmentation <- paste0(dir_code,"/BaM_exe/Segmentation")
        
        if (rec.model.rt[[mod]] =="1expWithAsympt" ){
          param.var.model = c( "a2") 
          
        } else if (rec.model.rt[[mod]] =="2expWithAsympt" ){
          param.var.model = c( "a3") 
          
        } else if (rec.model.rt[[mod]] =="2expWithAsympt_rel" ){
          param.var.model = c( "a3") 
          
        } else if (rec.model.rt[[mod]] =="3expWithAsympt" ){
          param.var.model = c("a4")
          
        } else if (rec.model.rt[[mod]] =="3expWithAsympt_bis" ){
          param.var.model = c("a4")
          
        } else if (rec.model.rt[[mod]]=="expexp" ){
          param.var.model = c( "a2")
          
        } else if (rec.model.rt[[mod]] =="hyperb" ){
          param.var.model = c("a2")
          
        } else if ((rec.model.rt[[mod]] =="hyperb_bis" )|
                   (rec.model.rt[[mod]] =="expexp_bis" )|
                   (rec.model.rt[[mod]] =="Coutagne_bis" )) { 
          param.var.model = c("a2")
          
        } else if (rec.model.rt[[mod]] =="2expWithAsympt_bis" ){
          param.var.model = c("a3")
        }
        dir.segm.recessions.model = paste0(dir.segm.recessions,"/",rec.model.rt[[mod]])
        dir.segm.recessions.model.param = paste0(dir.segm.recessions.model,"/chi_",chi.rt[chi.i])
        dir.segm.recessions.param = paste0(dir.segm.recessions.model.param,"/",param.var.model)
        dir.create(paste0(dir.c, "/segmentation"))
        dir.c.segm = paste0(dir.c, "/segmentation")
      
        hasympt = h.infinity
        hasympt$t =  h.infinity$t + df.limni$t_limni[realtime.period[1]]
        data.segm.rec3  = rbind(data.segm.rec2,   c(X2.5. = hasympt$`2.5%`,
                                                    X50. = hasympt$`50%`,
                                                    X97.5. = hasympt$`97.5%`,
                                                    mean = hasympt$mean,
                                                    stdev = hasympt$stdev,
                                                    maxpost = hasympt$maxpost,
                                                    t = hasympt$t))
        write.table(data.segm.rec3, paste0(dir.c.segm,"/Segm_data.txt"),sep="\t",row.names=FALSE)
        seg.rec = recession.segmentation(
          dir.segmentation   = dir.segmentation, 
          dir_code            = dir_code,
          dir.case_study      = dir.case_study, 
          dir_exe             = dir.exe, 
          dir.segm.recessions = dir.c.segm,
          param.to.segment    = param.var.model,
          all.param.to.segment = param.var.model,
          file.dir.with.data  = paste0(dir.c.segm,"/Segm_data.txt"), 
          Ncycle.segment      = 1500,
          Nmcmc.segment       = 100,
          Nslim.segment       = 50,
          Nburn.segment       = 0.5,
          df.limni            = df.limni,
          nSmax               = 2,       # 5, change this for each case study !!!!!!! 
          tmin                = 1, 
          criterion           = "DIC", 
          limits.X            = c(limni.time.limits[1], limni.time.limits[2], 1000), 
          limits.Y            = c(asymptote.limits[1], asymptote.limits[2], 50),
          x.name              = "Time [days]", 
          y.name              = "Asymptotic stage [m]",
          obs.uncertainty.y   = TRUE,
          prior.mu.segm       = c("Uniform", asymptote.limits),
          gamma_prior         = c("Uniform", 0, 50, 0.1) ,
          save.all.results    = TRUE,
          floodpeakalways     = FALSE,
          MAPalways           = TRUE,
          FinalColors         = colo,
          gaugings            = data4BaRatin,
          BayesianOption      = BayesianOption,
          stage.scale.shift   = stage.scale.shift)
        read.res.rec = read.results.segment.recess(dir.segm.recessions = dir.c.segm , 
                                                    officialShiftsTime  = officialShiftsTime, 
                                                    Gaugings            = Gaugings,
                                                    plot.dates          = FALSE)
        
        if ((read.res.rec$nS.ok >1)& (shift.true == FALSE)){
                             shift.true = TRUE
                             print(paste0("shift detected after ", h.infinity$t , " days from flood peak")) 
                             shift.time.rt = h.infinity$t 
                             asympt.plot[[chi.i]][[mod]] =  asympt.plot[[chi.i]][[mod]] +
                                                            geom_vline(aes(xintercept = shift.time.rt), 
                                                                       color="red",
                                                                       linetype = "dashed", size = 1)
                             ggsave(asympt.plot[[chi.i]][[mod]], filename =paste0(dir.regression.chi.mod,"/asymptote.png"),
                             bg = "transparent", width = 12, height =6, dpi = 200)
        }
        setwd(dir_code)
        
        
        
        # save the 3 frames that you want:
        if ((ttt == frames.time[1])|(ttt == frames.time[2])|(ttt == frames.time[3])) {
        ##############################################################################
            fram = fram+1
            reg.rt.plot[[fram]] =  ggplot() + 
                                   theme_bw(base_size = 25)+
                                   scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
                                                      expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
                                     scale_y_continuous(name = "Stage h (cm)", limits = c(limits.y[1], limits.y[2]), expand = c(0,0)) +
                                     theme(plot.title = element_text(hjust = 0.5),
                                           panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank(),
                                           panel.background = element_blank(),
                                           plot.margin=unit(c(0.5, 1, 0.5, 0.5),"cm")) +
                                    coord_cartesian(clip="off")+
                                     # geom_ribbon(data= PltData.rec,
                                     #             aes(x=t,
                                     #                 ymin=inf,
                                     #                 ymax= sup,
                                     #                 group=Period), color="green", alpha=0.1) +
                                     # geom_path(data= PltData.rec,
                                     #             aes(x=t, 
                                     #                 y = maxpost,
                                     #                 group=Period), color="green")+       
                                     geom_point(data=data.rec.obs.1,
                                                aes(x=time, 
                                                    y = h- stage.scale.shift,
                                                    group=Period), color="gray90", size= 2)+
                                     geom_vline(xintercept = h.infinity$t, color="blue", linetype ="8f", size=0.3)+
                                     geom_point(data = temp , aes(x = x, y= y), color= "blue", size= 4) +
                                     geom_ribbon(data = temp2, aes(x = xx,    ymin = zz , ymax = kk), 
                                                 fill = "blue", alpha=0.1) +
                                     geom_line(data = temp2, aes(x = xx,  y= V1), color = "blue")+
                                     annotate("text", 
                                              x =  70,
                                              y = 200, 
                                              label = paste0("after  ", round(h.infinity$t,1)," days"), 
                                              color = "red",  fontface= 1,
                                              size=15)+
                                     annotate("text", 
                                              x =  40,
                                              y = 50, 
                                              label = "past \n recessions", 
                                              color = "gray80",  fontface= 1,
                                              size=8)+
                                     annotate("text", 
                                                x =  30,
                                                y = -70, 
                                                label = "actual \n recession", 
                                                color = "blue",  fontface= 1,
                                                size=8)

            ##############################################################################################       
            asympt.rt.plot[[fram]] <- ggplot()+
                                           theme_bw(base_size = 25)+
                                             scale_x_continuous(name = "time (day)", limits =c(limits.x[1], limits.x[2]),
                                                                expand = c(0,0), breaks=seq(limits.x[1], limits.x[2], limits.x[3]))+
                                             scale_y_continuous(name = TeX("Asymptotic level  $\\; \\beta$ (cm)"), 
                                                                limits = c(limits.y[1], 100), expand = c(0,0)) +
                                             theme(plot.title = element_text(hjust = 0.5),
                                                   panel.grid.major = element_blank(), 
                                                   panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(),
                                                   plot.margin=unit(c(0.5, 1, 0.5, 0.5),"cm")) +
                                             geom_vline(xintercept = h.infinity$t, color="blue", 
                                                        linetype ="8f", size=0.3)+
                                             annotate(geom = "rect",
                                                      xmin = 0, xmax= 150,
                                                      ymin = read.res.rec.old$mu.results.df$mu.q2, 
                                                      ymax = read.res.rec.old$mu.results.df$mu.q97,
                                                      fill= "gray60", alpha= 0.1)+
                                             geom_segment(aes(x= 0, xend=150, 
                                                              y= read.res.rec.old$mu.results.df$mu.mean, 
                                                              yend= read.res.rec.old$mu.results.df$mu.mean), 
                                                          color="gray80")+
                                             geom_point(data = h.infinity, aes(x = t, y = maxpost),
                                                        colour = "blue",
                                                        size = 4) +
                                             geom_errorbar(data = h.infinity, aes(x= t, ymin = `2.5%`, ymax =`97.5%`), 
                                                           width=4, size = 1, color="blue")+
              annotate("text", 
                       x =  70,
                       y = read.res.rec.old$mu.results.df$mu.mean + 12, 
                       label =TeX("$\\beta$ of previous recessions (mean)"), 
                       color = "gray80",  fontface= 1,
                       size=8)
            if (read.res.rec$nS.ok >1) {
              asympt.rt.plot[[fram]] = asympt.rt.plot[[fram]]+
              annotate("text", 
                       x =  70,
                       y = -70, 
                       label ="Shift  \n detected !!!", 
                       color = "red",  fontface= 2,
                       size=7)
            }
            combined.plot.rt[[fram]] = plot_grid( reg.rt.plot[[fram]] , 
                                                 asympt.rt.plot[[fram]],
                                                 ncol = 1, nrow =2,
                                                 rel_heights = c(1,0.5)) 
           
        }        
        
      } # end of loop in real time
      
      
      
      ggsave(reg.plot[[chi.i]][[mod]], filename =paste0(dir.regression.chi.mod,"/regression_realtime.png"),
             bg = "transparent", width = 12, height =7, dpi = 200)
      plot.bind[[chi.i]][[mod]] = plot_grid(reg.plot[[chi.i]][[mod]] , 
                                            asympt.plot[[chi.i]][[mod]],
                                            ncol = 1, nrow =2,
                                            rel_heights = c(1,0.5)) 
      pdf(paste0(dir.regression.chi.mod,"/Figure_realtime_regression.pdf"), 12, 7,  useDingbats=F)
      print(reg.plot[[chi.i]][[mod]])
      dev.off()
      pdf(paste0(dir.regression.chi.mod,"/Figure_realtime_asymptote.pdf"), 12, 7 , useDingbats=F)
      print(asympt.plot[[chi.i]][[mod]])
      dev.off()
      pdf(paste0(dir.regression.chi.mod,"/Figure_realtime_recession.pdf"), 12, 12 , useDingbats=F)
      print(plot.bind[[chi.i]][[mod]])
      dev.off()
      
      plot.figure.paper.rt = plot_grid( plot.bind[[chi.i]][[1]]  , 
                                        plot.bind[[chi.i]][[2]] ,
                                        ncol = 1, nrow =2,
                                        rel_heights = c(1,1))   
      ggsave(plot.figure.paper.rt, 
             filename =paste0(dir.regression,"/figure_paper_realtime.png"),
             bg = "transparent", width = 12, height =7, dpi = 200)
    
      
      
    } # end of loop in models!
    plot.bind.chi[[chi.i]] = plot_grid( plot.bind[[chi.i]][[1]]  , 
                                        plot.bind[[chi.i]][[2]] ,
                                        ncol = 1, nrow =2,
                                        rel_heights = c(1,1))   
  } # end of loop in chi!
  ###############################################################################################    
  plot.bind.chi.tot = plot_grid( combined.plot.rt[[1]]  , 
                                 combined.plot.rt[[2]] ,
                                 combined.plot.rt[[3]], 
                                 ncol = 3, nrow =1,
                                 rel_widths = c(1,1,1)
                                 )  
  pdf(paste0(dir.regression,"/Figure_realtime_final.pdf"), 23, 16 , useDingbats=F)
  print(plot.bind.chi.tot)
  dev.off()
  # # save:
  # if (rec.model == "2expWithAsympt") {
  #   Data.segm.rec <- data.frame("t"=t.real.good.h, "Y"=unlist(lapply(Resultss, "[[", 16,5), use.names=FALSE), 
  #                               "uY"= unlist(lapply(Resultss, "[[",11,5), use.names=FALSE))
  # } else if (rec.model == "3expWithAsympt") {
  #   Data.segm.rec <- data.frame("t"=t.real.good.h, "Y"=unlist(lapply(Resultss, "[[", 16,7), use.names=FALSE), 
  #                               "uY"= unlist(lapply(Resultss, "[[",11,7), use.names=FALSE))
  # }
  end_time <- Sys.time()
  print(c("computat. time for regression of recessions=", end_time-start_time ))
  return(list(asymptote.df))
}











