##########################################################################################################
read_all_results_shift_detection = function(dir.segment.g, dir.segment.rec, dir.segment.gaug.rec, 
                                            stage.record,        
                                            gaugings,       
                                            official.dates, 
                                            colors.period,
                                            file.options.general,
                                            file.options.segment,
                                            file.options.SPD,
                                            file.options.recessions){
##########################################################################################################
    dir.create(dir.segment.gaug.rec)
    source(file.options.general)
    source(file.options.recessions)
    source(file.options.segment)
    
  
    # Gaugings:
    ############
    user.choice.gaugings  = dlgInput("Do you have results from Gaugings segmentation ? [Y/N]",  Sys.info()[" "])$res
    if (user.choice.gaugings == "Y"){
    list_all_tests_gaugings     = list.dirs(path = dir.segment.g, recursive=FALSE)
    if (length(list_all_tests_gaugings) > 0) {
       if (length(list_all_tests_gaugings) >1) {
          user.choice.folder.gaugings  = list_all_tests_gaugings[ as.numeric(
                                      dlgInput(c("Insert the folder with the results from gaugings segmentation: \n",
                                      paste0(seq(1, length(list_all_tests_gaugings), 1), " = ",
                                             list_all_tests_gaugings)),  Sys.info()[" "])$res)]
          dir.segment.gaug.1           = user.choice.folder.gaugings
       } else {
          dir.segment.gaug.1           = list_all_tests_gaugings
       }
       list_selected_results_gaug   = list.files(path = dir.segment.gaug.1, recursive=FALSE)
       
       if (length(list_selected_results_gaug) > 0) {
          if (any(list_selected_results_gaug == "shift_times.txt")) {
             g.s.results.1                = read.table(paste0(dir.segment.gaug.1,"/data_with_periods.txt"), sep ="\t", header =TRUE) 
             pdf.ts.results.1             = read.table(paste0(dir.segment.gaug.1,"/pdf_ts.txt"), sep ="\t", header =TRUE) 
             nperiod.from.segm.gaugings.1 = tail(g.s.results.1[[4]],1)
             shift.times.gaugings.1       = read.table(file =paste0(dir.segment.gaug.1,"/shift_times.txt"), header = TRUE)
             data.annotate.gaug.1 <- data.frame( q2   = shift.times.gaugings.1$t2,
                                                 q10  = shift.times.gaugings.1$t10,
                                                 MAP  = shift.times.gaugings.1$tMAP, 
                                                 q90  = shift.times.gaugings.1$t90,
                                                 q97  = shift.times.gaugings.1$t97)
             data.annotate.gaug.1.sort <- data.annotate.gaug.1[  with( data.annotate.gaug.1, order(MAP)),]
             data.annotate.gaug.1.sort = cbind(data.annotate.gaug.1.sort,  t.adj = shift.times.gaugings.1$treal)
             pdf.ts.results.1.sort = pdf.ts.results.1[,  with( data.annotate.gaug.1, order(MAP)) ]
             
          } else {
             g.s.results.1                = NULL
             pdf.ts.results.1             = NULL
             nperiod.from.segm.gaugings.1 = NULL
             shift.times.gaugings.1       = NULL
             data.annotate.gaug.1.sort    = NULL
             pdf.ts.results.1.sort        = NULL
          }
       } else {
          g.s.results.1                = NULL
          pdf.ts.results.1             = NULL
          nperiod.from.segm.gaugings.1 = NULL
          shift.times.gaugings.1       = NULL
          data.annotate.gaug.1.sort    = NULL
          pdf.ts.results.1.sort        = NULL
       }
    } else {
       print('****** Gaugings segmentation has NOT been performed yet!  Please check.')
       g.s.results.1                = NULL
       pdf.ts.results.1             = NULL
       nperiod.from.segm.gaugings.1 = NULL
       shift.times.gaugings.1       = NULL
       data.annotate.gaug.1.sort    = NULL
       pdf.ts.results.1.sort        = NULL
    }
    
    } else {
      print('You choose to NOT accounting for gaugings segmentation results.')
      g.s.results.1                = NULL
      pdf.ts.results.1             = NULL
      nperiod.from.segm.gaugings.1 = NULL
      shift.times.gaugings.1       = NULL
      data.annotate.gaug.1.sort    = NULL
      pdf.ts.results.1.sort        = NULL
      
    }
    
    
    
    
    
    
    # Recessions:
    #############
    user.choice.recessions  = dlgInput("Do you have results from stage-recession analysis ? [Y/N]",  Sys.info()[" "])$res
    if (user.choice.recessions == "Y"){
    user.choice.folder.recess =1
    folder = dir.segment.rec
    while (user.choice.folder.recess >0){
      list_folders = list.dirs(path = folder, recursive=FALSE)
      if (length(list_folders) > 0) {
        user.choice.folder.recess  =dlgInput(c("Looking for file recession segmentation results", 
                                               "(e.g., mcmc_segmentation.txt):",
                                               "*********************************************",
                                            "move to folder: [index]", " ",
                                            paste0(0, "   = ", folder, " (STOP, stay in this folder)"),
                                            " ",
                                            paste0(seq(1, length(list_folders), 1), "   = ",
                                                   list_folders)),  Sys.info()[" "])$res
        if (user.choice.folder.recess != "0"){
          folder = list_folders[as.numeric(user.choice.folder.recess)]
        }
      } else {
        user.choice.folder.recess =0 
        print("*****No folders in recession analysis. Please check!")
        break
      }
    } 
      
      print(paste0("Opening folder: ", folder))
      list_selected_results_rec   = list.files(path = folder, recursive=FALSE)
      print(list_selected_results_rec)
      if (length(list_selected_results_rec)>0) {
         dir.segment.rec = folder
         read.res.rec = read.results.segment.recess(dir.segm.recessions = dir.segment.rec, 
                                                    officialShiftsTime  = officialShiftsTime, 
                                                    Gaugings            = Gaugings,
                                                    plot.dates          = FALSE)
         data.annotate.recess <- data.frame(q2    = read.res.rec$Q2.ts, 
                                            q10   = read.res.rec$Q10.ts, 
                                            MAP   = read.res.rec$data.annotate.recess$t,
                                            q90   = read.res.rec$Q90.ts, 
                                            q97   = read.res.rec$Q97.ts,
                                            t.adj = read.res.rec$data.annotate.recess$t) # not the t.adj !
      } else {
        read.res.rec = NULL
        data.annotate.recess = NULL
      }
    } else {
      print('You choose to NOT accounting for stage-recession analysis results.')
      read.res.rec = NULL
      data.annotate.recess = NULL
    }
      
    
    
      
      
    
      
    # OFFICIAL SHIFT TIMES:
    ######################
    if (!is.null(officialShiftsTime)) {
      data.annotate.off =  data.frame(xeffect = officialShiftsTime,
                                      xpotent = officialShiftsTime)
    } else {
      data.annotate.off =NULL
    }
    
    

    

    
    # Combining results of recession and gaugings:
    if (!is.null(data.annotate.gaug.1.sort)){
      if (!is.null(data.annotate.recess)){
       data.combined.gaug.recess = rbind(cbind(data.annotate.gaug.1.sort), 
                                           cbind(data.annotate.recess))
       minlength = min(length(pdf.ts.results.1.sort$V1), length(read.res.rec$pdf.ts.rec$tau1))
       pdf.ts.combined = cbind(pdf.ts.results.1.sort[1:minlength,],  
                               read.res.rec$pdf.ts.rec[1:minlength,])
       data.combined.gaug.recess.sort = data.combined.gaug.recess[
                                        with( data.combined.gaug.recess, order(t.adj)),
                                          ]
       pdf.ts.combined.sort = pdf.ts.combined[ ,with( data.combined.gaug.recess, order(t.adj))]
      } else {
        data.combined.gaug.recess.sort = data.annotate.gaug.1.sort
        pdf.ts.combined.sort           = pdf.ts.results.1.sort
      }
    } else {
       if (!is.null(data.annotate.recess)){
          data.combined.gaug.recess.sort = data.annotate.recess
          pdf.ts.combined.sort                = read.res.rec$pdf.ts.rec
       } else {
          data.combined.gaug.recess.sort = NULL
          pdf.ts.combined.sort                = NULL
       }  
    }

    
    

    if (!is.null(pdf.ts.combined.sort)) {
      if (!is.null(dev.list())){  #clear all devices if any
        while (!is.null(dev.list())) {
          dev.off()
        }
      }   
       #plot all results:
      X11() ;   
      plot.shifts =  plot.time.shifts.step1(           dir                  = dir.segment.gaug.rec ,
                                                       gaug.1               = gaugings, 
                                                       rec.1                = NULL, 
                                                       data.annotate.off    = data.annotate.off, 
                                                       data.annotate.gaug.1 = data.annotate.gaug.1.sort ,
                                                       data.annotate.rec.1  = data.annotate.recess ,
                                                       data.annotate.step1  = data.combined.gaug.recess.sort,
                                                       color                = colors.period, 
                                                       df.limni             = stage.record, 
                                                       limni.labels         = limni.labels , 
                                                       grid_limni.ylim      = grid_limni.ylim,
                                                       pdf.ts.1             = pdf.ts.results.1.sort, 
                                                       pdf.ts.2             = read.res.rec$pdf.ts.rec,
                                                       pdf.ts.3             = pdf.ts.combined.sort,
                                                       save_file_name       = "/STEP1_time_series.pdf")
       print(plot.shifts)
       
                    user.choice.step1 <- dlgInput(paste0("do you want to discard any of these shift times ? [Y / N] \n" ), 
                                           Sys.info()[" "])$res
                    if (user.choice.step1 == "Y") {  
                     # eliminate shift times ?
                     user.remove.ts.step1 =1; user.eliminate.ts.step1=c()
                     while (user.remove.ts.step1 != 0){                          
                        user.remove.ts.step1 = as.numeric(dlgInput(c("Select the shift time/s to remove (enter 0 when finished):",
                                                                     "0   = STOP", " ",
                                                                     paste0(seq(1, length(data.combined.gaug.recess.sort$t.adj), 1), "   = ",  
                                                                     c(data.combined.gaug.recess.sort$t.adj))),  Sys.info()[" "])$res)
                        user.eliminate.ts.step1 = c(user.eliminate.ts.step1, user.remove.ts.step1)
                     }
                     user.eliminate.ts.step1 = user.eliminate.ts.step1[-c(length(user.eliminate.ts.step1))]
                     data.annotate.step1     = data.combined.gaug.recess.sort[ - user.eliminate.ts.step1, ]
                     message("You have removed these shift times: ")
                     print(data.combined.gaug.recess.sort[user.eliminate.ts.step1, ])
                     pdf.ts.step1            = pdf.ts.combined.sort[ , - user.eliminate.ts.step1]
                     plot.shifts =  plot.time.shifts.step1(dir                  = dir.segment.gaug.rec ,
                                                           gaug.1               = gaugings, 
                                                           rec.1                = NULL, 
                                                           data.annotate.off    = data.annotate.off, 
                                                           data.annotate.gaug.1 = data.annotate.gaug.1.sort ,
                                                           data.annotate.rec.1  = data.annotate.recess ,
                                                           data.annotate.step1  = data.annotate.step1,
                                                           color                = colors.period, 
                                                           df.limni             = stage.record, 
                                                           limni.labels         = limni.labels , 
                                                           grid_limni.ylim      = grid_limni.ylim,
                                                           pdf.ts.1             = pdf.ts.results.1.sort, 
                                                           pdf.ts.2             = read.res.rec$pdf.ts.rec,
                                                           pdf.ts.3             = pdf.ts.step1,
                                                           save_file_name       = "/STEP1_time_series_adjusted.pdf")
                                             
                                             
                  } else {
                            data.annotate.step1 = data.combined.gaug.recess.sort
                            pdf.ts.step1        = pdf.ts.combined.sort                     
                  }  
                  # Sorting the final dataframe of shift times:  
                  if (!is.null(data.annotate.step1)){              
                    if (length(data.annotate.step1$t.adj) >1 ) {               
                      data.annotate.step1.sort = data.annotate.step1[with( data.annotate.step1, order(t.adj)), ]
                      pdf.ts.step1.sort        = pdf.ts.step1[ ,with( data.annotate.step1, order(t.adj))]
                    } else {
                      data.annotate.step1.sort  = data.annotate.step1
                      pdf.ts.step1.sort         = pdf.ts.step1
                    }
                  } else {
                    data.annotate.step1.sort  = NULL
                    pdf.ts.step1.sort         = NULL
                  }
    } else {
        data.annotate.step1.sort  = NULL
        pdf.ts.step1.sort         = NULL
    }
      
  
   if (!is.null(dev.list())){
     while (!is.null(dev.list())) {
       dev.off()
     }
   }   
   
   print("****************")
   print("   All done!    ")
   print("****************")
   return(list( df_gaug      = data.annotate.gaug.1.sort,
                pdf_gaug     = pdf.ts.results.1.sort,
                df_rec       = data.annotate.recess,
                pdf_rec      = read.res.rec$pdf.ts.rec,
                df_gaug_rec  = data.annotate.step1.sort,
                pdf_gaug_rec = pdf.ts.step1.sort,
                df_offic     = data.annotate.off))
}
































##########################################################################################################
read_all_results_shift_detection2 = function(dir.segment.g, dir.segment.rec, dir.segment.gaug.rec, dir.Sedim.transp,
                                             stage.record,        
                                             gaugings,       
                                             official.dates, 
                                             colors.period,
                                             file.options.general,
                                             file.options.segment,
                                             file.options.SPD,
                                             file.options.recessions,
                                             file.options.ST){
##########################################################################################################
  dir.create(dir.segment.gaug.rec)
  source(file.options.general)
  source(file.options.recessions)
  source(file.options.segment)
  source(file.options.ST)
  
  
  # Gaugings:
  ############
  list_all_tests_gaugings     = list.dirs(path = dir.segment.g, recursive=FALSE)
  if (length(list_all_tests_gaugings) > 0) {
    if (length(list_all_tests_gaugings) >1) {
      user.choice.folder.gaugings  = list_all_tests_gaugings[ as.numeric(
        dlgInput(c("Insert the folder with the results from gaugings segmentation: \n",
                   paste0(seq(1, length(list_all_tests_gaugings), 1), " = ",
                          list_all_tests_gaugings)),  Sys.info()[" "])$res)]
      dir.segment.gaug.1           = user.choice.folder.gaugings
    } else {
      dir.segment.gaug.1           = list_all_tests_gaugings
    }
    list_selected_results_gaug   = list.files(path = dir.segment.gaug.1, recursive=FALSE)
    
    if (length(list_selected_results_gaug) > 0) {
      if (any(list_selected_results_gaug == "shift_times.txt")) {
        g.s.results.1                = read.table(paste0(dir.segment.gaug.1,"/data_with_periods.txt"), sep ="\t", header =TRUE) 
        pdf.ts.results.1             = read.table(paste0(dir.segment.gaug.1,"/pdf_ts.txt"), sep ="\t", header =TRUE) 
        nperiod.from.segm.gaugings.1 = tail(g.s.results.1[[4]],1)
        shift.times.gaugings.1       = read.table(file =paste0(dir.segment.gaug.1,"/shift_times.txt"), header = TRUE)
        data.annotate.gaug.1 <- data.frame( q2   = shift.times.gaugings.1$t2,
                                            q10  = shift.times.gaugings.1$t10,
                                            MAP  = shift.times.gaugings.1$tMAP, 
                                            q90  = shift.times.gaugings.1$t90,
                                            q97  = shift.times.gaugings.1$t97)
        data.annotate.gaug.1.sort <- data.annotate.gaug.1[  with( data.annotate.gaug.1, order(MAP)),]
        data.annotate.gaug.1.sort = cbind(data.annotate.gaug.1.sort,  t.adj = shift.times.gaugings.1$treal)
        pdf.ts.results.1.sort = pdf.ts.results.1[,  with( data.annotate.gaug.1, order(MAP)) ]
        
      } else {
        g.s.results.1                = NULL
        pdf.ts.results.1             = NULL
        nperiod.from.segm.gaugings.1 = NULL
        shift.times.gaugings.1       = NULL
        data.annotate.gaug.1.sort    = NULL
        pdf.ts.results.1.sort        = NULL
      }
    } else {
      g.s.results.1                = NULL
      pdf.ts.results.1             = NULL
      nperiod.from.segm.gaugings.1 = NULL
      shift.times.gaugings.1       = NULL
      data.annotate.gaug.1.sort    = NULL
      pdf.ts.results.1.sort        = NULL
    }
  } else {
    print('****** Gaugings segmentation has not been performed yet!  Please check.')
    g.s.results.1                = NULL
    pdf.ts.results.1             = NULL
    nperiod.from.segm.gaugings.1 = NULL
    shift.times.gaugings.1       = NULL
    data.annotate.gaug.1.sort    = NULL
    pdf.ts.results.1.sort        = NULL
  }
  
  
  
  
  
  
  # Recessions:
  #############
  user.choice.folder.recess =1
  folder = dir.segment.rec
  while (user.choice.folder.recess >0){
    list_folders = list.dirs(path = folder, recursive=FALSE)
    if (length(list_folders) > 0) {
      user.choice.folder.recess  =dlgInput(c("Looking for file recession segmentation results", 
                                             "(e.g., mcmc_segmentation.txt):",
                                             "*********************************************",
                                             "move to folder: [index]", " ",
                                             paste0(0, "   = ", folder, " (STOP, stay in this folder)"),
                                             " ",
                                             paste0(seq(1, length(list_folders), 1), "   = ",
                                                    list_folders)),  Sys.info()[" "])$res
      if (user.choice.folder.recess != "0"){
        folder = list_folders[as.numeric(user.choice.folder.recess)]
      }
    } else {
      user.choice.folder.recess =0 
      print("*****No folders in recession analysis. Please check!")
      break
    }
  } 
  
  print(paste0("Opening folder: ", folder))
  list_selected_results_rec   = list.files(path = folder, recursive=FALSE)
  print(list_selected_results_rec)
  if (length(list_selected_results_rec)>0) {
    dir.segment.rec = folder
    read.res.rec = read.results.segment.recess(dir.segm.recessions = dir.segment.rec, 
                                               officialShiftsTime  = officialShiftsTime, 
                                               Gaugings            = Gaugings,
                                               plot.dates          = FALSE)
    data.annotate.recess <- data.frame(q2    = read.res.rec$Q2.ts, 
                                       q10   = read.res.rec$Q10.ts, 
                                       MAP   = read.res.rec$data.annotate.recess$t,
                                       q90   = read.res.rec$Q90.ts, 
                                       q97   = read.res.rec$Q97.ts,
                                       t.adj = read.res.rec$data.annotate.recess$t) # not the t.adj !
  } else {
    read.res.rec = NULL
    data.annotate.recess = NULL
  }
  
  
  

  
  
  
  # OFFICIAL SHIF TIMES:
  ######################
  if (!is.null(officialShiftsTime)) {
    data.annotate.off =  data.frame(xeffect = officialShiftsTime,
                                    xpotent = officialShiftsTime)
  } else {
    data.annotate.off =NULL
  }
  
  
  
  
  
  # Combining results of recession and gaugings:
  if (!is.null(data.annotate.gaug.1.sort)){
    if (!is.null(data.annotate.recess)){
      data.combined.gaug.recess = rbind(cbind(data.annotate.gaug.1.sort), 
                                        cbind(data.annotate.recess))
      minlength = min(length(pdf.ts.results.1.sort$V1), length(read.res.rec$pdf.ts.rec$tau1))
      pdf.ts.combined = cbind(pdf.ts.results.1.sort[1:minlength,],  
                              read.res.rec$pdf.ts.rec[1:minlength,])
      data.combined.gaug.recess.sort = data.combined.gaug.recess[
        with( data.combined.gaug.recess, order(t.adj)),
        ]
      pdf.ts.combined.sort = pdf.ts.combined[ ,with( data.combined.gaug.recess, order(t.adj))]
    } else {
      data.combined.gaug.recess.sort = data.annotate.gaug.1.sort
      pdf.ts.combined.sort           = pdf.ts.results.1.sort
    }
  } else {
    if (!is.null(data.annotate.recess)){
      data.combined.gaug.recess.sort = data.annotate.recess
      pdf.ts.combined.sort                = read.res.rec$pdf.ts.rec
    } else {
      data.combined.gaug.recess.sort = NULL
      pdf.ts.combined.sort                = NULL
    }  
  }
  
  
  
  
  
  
  # Sediment transport (potential shifts):
  #########################################
  user.choice.folder.ST =1
  folder = dir.Sedim.transp
  while (user.choice.folder.ST >0){
    list_folders = list.dirs(path = folder, recursive=FALSE)
    if (length(list_folders) > 0) {
      user.choice.folder.ST  =dlgInput(c("Looking for file with sediment transport results", 
                                             "(e.g., mcmc_segmentation.txt):",
                                             "*********************************************",
                                             "move to folder: [index]", " ",
                                             paste0(0, "   = ", folder, " (STOP, stay in this folder)"),
                                             " ",
                                             paste0(seq(1, length(list_folders), 1), "   = ",
                                                    list_folders)),  Sys.info()[" "])$res
      if (user.choice.folder.recess != "0"){
        folder = list_folders[as.numeric(user.choice.folder.recess)]
      }
    } else {
      user.choice.folder.recess =0 
      print("*****No folders in recession analysis. Please check!")
      break
    }
  } 
  
  print(paste0("Opening folder: ", folder))
  list_selected_results_rec   = list.files(path = folder, recursive=FALSE)
  print(list_selected_results_rec)
  if (length(list_selected_results_rec)>0) {
    dir.segment.rec = folder
    read.res.rec = read.results.segment.recess(dir.segm.recessions = dir.segment.rec, 
                                               officialShiftsTime  = officialShiftsTime, 
                                               Gaugings            = Gaugings,
                                               plot.dates          = FALSE)
    data.annotate.recess <- data.frame(q2    = read.res.rec$Q2.ts, 
                                       q10   = read.res.rec$Q10.ts, 
                                       MAP   = read.res.rec$data.annotate.recess$t,
                                       q90   = read.res.rec$Q90.ts, 
                                       q97   = read.res.rec$Q97.ts,
                                       t.adj = read.res.rec$data.annotate.recess$t) # not the t.adj !
  } else {
    read.res.rec = NULL
    data.annotate.recess = NULL
  }
  
  
  
  
  
  
  
  
  
  
  
  if (!is.null(pdf.ts.combined.sort)) {
    if (!is.null(dev.list())){  #clear all devices if any
      while (!is.null(dev.list())) {
        dev.off()
      }
    }   
    #plot all results:
    X11() ;   
    plot.shifts =  plot.time.shifts.step1(           dir                  = dir.segment.gaug.rec ,
                                                     gaug.1               = gaugings, 
                                                     rec.1                = NULL, 
                                                     data.annotate.off    = data.annotate.off, 
                                                     data.annotate.gaug.1 = data.annotate.gaug.1.sort ,
                                                     data.annotate.rec.1  = data.annotate.recess ,
                                                     data.annotate.step1  = data.combined.gaug.recess.sort,
                                                     color                = colors.period, 
                                                     df.limni             = stage.record, 
                                                     limni.labels         = limni.labels , 
                                                     grid_limni.ylim      = grid_limni.ylim,
                                                     pdf.ts.1             = pdf.ts.results.1.sort, 
                                                     pdf.ts.2             = read.res.rec$pdf.ts.rec,
                                                     pdf.ts.3             = pdf.ts.combined.sort,
                                                     save_file_name       = "/STEP1_time_series.pdf")
    print(plot.shifts)
    
    user.choice.step1 <- dlgInput(paste0("do you want to discard any of these shift times ? [Y / N] \n" ), 
                                  Sys.info()[" "])$res
    if (user.choice.step1 == "Y") {  
      # eliminate shift times ?
      user.remove.ts.step1 =1; user.eliminate.ts.step1=c()
      while (user.remove.ts.step1 != 0){                          
        user.remove.ts.step1 = as.numeric(dlgInput(c("Select the shift time/s to remove (enter 0 when finished):",
                                                     "0   = STOP", " ",
                                                     paste0(seq(1, length(data.combined.gaug.recess.sort$t.adj), 1), "   = ",  
                                                            c(data.combined.gaug.recess.sort$t.adj))),  Sys.info()[" "])$res)
        user.eliminate.ts.step1 = c(user.eliminate.ts.step1, user.remove.ts.step1)
      }
      user.eliminate.ts.step1 = user.eliminate.ts.step1[-c(length(user.eliminate.ts.step1))]
      data.annotate.step1     = data.combined.gaug.recess.sort[ - user.eliminate.ts.step1, ]
      message("You have removed these shift times: ")
      print(data.combined.gaug.recess.sort[user.eliminate.ts.step1, ])
      pdf.ts.step1            = pdf.ts.combined.sort[ , - user.eliminate.ts.step1]
      plot.shifts =  plot.time.shifts.step1(dir                  = dir.segment.gaug.rec ,
                                            gaug.1               = gaugings, 
                                            rec.1                = NULL, 
                                            data.annotate.off    = data.annotate.off, 
                                            data.annotate.gaug.1 = data.annotate.gaug.1.sort ,
                                            data.annotate.rec.1  = data.annotate.recess ,
                                            data.annotate.step1  = data.annotate.step1,
                                            color                = colors.period, 
                                            df.limni             = stage.record, 
                                            limni.labels         = limni.labels , 
                                            grid_limni.ylim      = grid_limni.ylim,
                                            pdf.ts.1             = pdf.ts.results.1.sort, 
                                            pdf.ts.2             = read.res.rec$pdf.ts.rec,
                                            pdf.ts.3             = pdf.ts.step1,
                                            save_file_name       = "/STEP1_time_series_adjusted.pdf")
      
      
    } else {
      data.annotate.step1 = data.combined.gaug.recess.sort
      pdf.ts.step1        = pdf.ts.combined.sort                     
    }  
    # Sorting the final dataframe of shift times:  
    if (!is.null(data.annotate.step1)){              
      if (length(data.annotate.step1$t.adj) >1 ) {               
        data.annotate.step1.sort = data.annotate.step1[with( data.annotate.step1, order(t.adj)), ]
        pdf.ts.step1.sort        = pdf.ts.step1[ ,with( data.annotate.step1, order(t.adj))]
      } else {
        data.annotate.step1.sort  = data.annotate.step1
        pdf.ts.step1.sort         = pdf.ts.step1
      }
    } else {
      data.annotate.step1.sort  = NULL
      pdf.ts.step1.sort         = NULL
    }
  } else {
    data.annotate.step1.sort  = NULL
    pdf.ts.step1.sort         = NULL
  }
  
  
  if (!is.null(dev.list())){
    while (!is.null(dev.list())) {
      dev.off()
    }
  }   
  
  print("****************")
  print("   All done!    ")
  print("****************")
  return(list( df_gaug      = data.annotate.gaug.1.sort,
               pdf_gaug     = pdf.ts.results.1.sort,
               df_rec       = data.annotate.recess,
               pdf_rec      = read.res.rec$pdf.ts.rec,
               df_gaug_rec  = data.annotate.step1.sort,
               pdf_gaug_rec = pdf.ts.step1.sort,
               df_offic     = data.annotate.off))
}
