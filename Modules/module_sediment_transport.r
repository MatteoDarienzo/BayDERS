#########################################
# MODULE FOR SEDIMENT TRANSPORT (BAYDERS)
#########################################
# Author: Matteo Darienzo, 
# Affiliation: Inrae (Lyon, France), Cima Foundation (Savona, Italy)
# Last update: 16/11/2021

# In this module you find several functions for the sediment transport analysis with the main objective of
# determining new potential morphological shifts in the rating curve.
# For more information please read PhD thesis manuscript of Matteo Darienzo, 2021.
#
# Main functions:
# 1) read.ref.morphogenic.events()     will define the reference past shift times that are assumed to be
#                                      origined by morphogenic floods.
# 2) sediment.transport.segmentation() will compute the sediment transport on the whole stage record
#
# 3) linear.estimation()               will perform a linear regression between the sediment transported 
#                                      volumes computed for each reference event vs. the estimated shift 
#                                      magnitude for each event (obtained from BaratinSPD analysis).
########################################











#########################################################
find_max = function(nmax, series,  interval, max_deltat){
#########################################################  
  
  int     = which(series$t_limni >= interval[1] & series$t_limni <= interval[2])
  vec     = series$h_limni[int[1]:tail(int,1)]
  t_vec   = series$t_limni[int[1]:tail(int,1)]
  indexes = 0
  N       = nmax 
  Nth     = sort(vec, partial = partial)[partial]
  indexes = which(vec >= Nth)
  h_indixes = vec[indexes]
  
  first_peak = which()
  
  indexes_remove = which.min(c(vec[2], vec[1]))
  t_vec_final    = t_vec[indexes_remove]
  h_vec_final    = vec[indexes_remove]
  
  
  while (length(indexes) < nmax){
     partial <- length(vec) - N + 1
     for (i in 3:length(indexes)){
        # limit of 2 days between peaks
        temp_index = c(i, i-1)
        if ((t_vec[indexes[i]] - tail(t_vec_final,1)) < max_deltat){  
           print("remove")
           indexes_remove = c(indexes_remove, temp_index[which.min(c(vec[indexes[i]], vec[indexes[i-1]]))])
        } else {
          t_vec_final = c(t_vec_final, t_vec[indexes[i]] )
          h_vec_final = c(h_vec_final, vec[indexes[i]] )
        }
     }
     indexes = indexes[-indexes_remove]
     #print(N)
     N = N + 1
  }
  return(list(hmax    = h_vec_final,
              indexes = indexes,
              t_hmax  = t_vec_final))
}
















###################################################
find_peaks <- function (x, m){
###################################################
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}













###########################################################################################################
read.refer.morphogenic.events = function(dir.exe, dir.case_study, dir.Sedim.transp, 
                                         file.options.general, file.options.segment, 
                                         file.options.recessions, file.options.ST,                                                            
                                         stage.record,  gaugings ,   official.dates, colors.period,
                                         res_gaug_recess){
###########################################################################################################
  source(file.options.general)
  source(file.options.recessions)
  source(file.options.segment)
  source(file.options.ST)
  dir.create(dir.Sedim.transp)
  dir.ST.test         = paste0(dir.Sedim.transp,"/",name.folder.results.ST)
  dir.create(dir.ST.test)
  dir.ST.test.step1   = paste0(dir.ST.test,"/1_reference_events")
  dir.create(dir.ST.test.step1)
  
  
  dir.create(paste0(dir.ST.test.step1,"/SPD"))
  dir.sed.transp.SPD = paste0(dir.ST.test.step1,"/SPD")
  dir.create(paste0(dir.exe, "/Linear"))
  dir.linear =paste0(dir.exe, "/Linear")
  
  
  # saving file with the selected options for the current test:
  file.options = paste0(dir.ST.test.step1,"/options.txt")
  cat(paste0("stage.limni.u.m.      =   " ,        u.m.limni,  "                 # unity of the stage record (df.limni$h_limni)"), file = file.options,  sep="\n")   
  cat(paste0("reference shift times =   ",    file_name_ref_shift_times,   "             # filename with the reference morphogenic shift times (*** FALSE if none***) "), file = file.options, append = TRUE, sep="\n")   
  
  user.choice.ST = dlgInput(paste0("REFERENCE SHIFT TIMES options: \n",
                                   "******************************\n",
                                   "1. Use your own .txt file with reference shif times.\n\n",
                                   "2. Use  shift times obtained from previous step (gaugings segmentation + recession analysis)"), 
                            Sys.info()[" "])$res

  
  #**************************
  if (user.choice.ST == "2"){
  #**************************
      # # Reading results of step 1 of the retro analysis:
      # file.copy(c(paste0(dir.retro.an,"/STEP1_SPD/bt1_df.txt"), 
      #             paste0(dir.retro.an,"/STEP1_SPD/bt2_df.txt"),
      #             paste0(dir.retro.an,"/STEP1_data_with_periods.txt"),
      #             paste0(dir.retro.an,"/STEP1_shift_times.txt")), 
      #           dir.sed.transp.SPD, overwrite = TRUE)

      
      ###################################################################################################
      # ask to user if she/he wants to delete some shift times because are possibly not related to floods
      # but to other causes:
      if (!is.null(dev.list())){
        while (!is.null(dev.list())) {
          dev.off()
        }
      }   
      X11(width = 20, height = 20)
      plot.shifts= plot.time.shifts.step1( dir                  = dir.ST.test.step1,
                                           gaug.1               = gaugings, 
                                           rec.1                = NULL, 
                                           data.annotate.off    = res_gaug_recess$df_offic, 
                                           data.annotate.gaug.1 = res_gaug_recess$df_gaug ,
                                           data.annotate.rec.1  = res_gaug_recess$df_rec ,
                                           data.annotate.step1  = res_gaug_recess$df_gaug_rec,
                                           colo                 = colors.period, 
                                           df.limni             = stage.record, 
                                           limni.labels         = limni.labels , 
                                           grid_limni.ylim      = grid_limni.ylim,
                                           pdf.ts.1             = res_gaug_recess$pdf_gaug, 
                                           pdf.ts.2             = res_gaug_recess$pdf_rec,
                                           pdf.ts.3             = res_gaug_recess$pdf_gaug_rec,
                                           save_file_name       = "/reference_shift_times.pdf")
      print(plot.shifts)
      user.choice.ST <- dlgInput(paste0("REFERENCE MORPHOGENIC EVENTS: \n",
                                        "*******************************\n",
                                        "Modify some shift times ? [Y / N] \n" ), 
                                        Sys.info()[" "])$res; dev.off();  
      if (user.choice.ST == "Y") {  
      #***************************
           user.eliminate.ST =   c()
           user.modify.ST    =   c()
           data.annotate.ST  =   res_gaug_recess$df_gaug_rec
           pdf.ts.ST         =   res_gaug_recess$pdf_gaug_rec
           if (!is.null(gaugings)){
               gaugings.P        =   gaugings; 
               names(gaugings.P) =   c("hP", "QP", "uQP", "Period", "tP", "t.trueP")
           } else {
               gaugings.P = NULL
           }
           stage.record.NA = na.omit(stage.record)
           
           for (ii in 1:length(res_gaug_recess$df_gaug_rec$t.adj)){ 
           #*******************************************************   
            # loop on all the selected events to find the flood peak:
             int     = which(stage.record.NA$t_limni >= res_gaug_recess$df_gaug_rec$q2[ii] &
                               stage.record.NA$t_limni <= res_gaug_recess$df_gaug_rec$q97[ii])
             # h_peaks = find_max(nmax       = 10, 
             #                    series     = stage.record, 
             #                    interval   = c(res_gaug_recess$df_gaug_rec$q2[ii],
             #                                   res_gaug_recess$df_gaug_rec$q97[ii]),
             #                    max_deltat = 20) # min days between consecutive peaks
             
            
             h_peaks = data.frame(hmax   = round(stage.record.NA$h_limni[int][find_peaks(x = stage.record.NA$h_limni[int], 
                                                                                         m=deltat_peaks)], digits =2),
                                  t_hmax = round(stage.record.NA$t_limni[int][find_peaks(x = stage.record.NA$h_limni[int],
                                                                                         m=deltat_peaks)], digits =2))
                      
             
             message(paste0("event ", ii, ":"))
             message("***************")
             message("A few peaks (timings) in the CI of the event: ")
             cat(paste0(seq(1, length(h_peaks$t_hmax),1) , ")  ", round(h_peaks$t_hmax,digits=2), " days   (h =", round(h_peaks$hmax, digits =3), " m)\n"))
             message("")
             
             X11(width = 20, height = 10); print(
               # call function from "module_gaugings_segmentation_plots.R"
               exemple_ts.plot(CdT.P        = gaugings.P, 
                               df.limni     = stage.record, 
                               Q10.ts       = res_gaug_recess$df_gaug_rec$q2, 
                               Q90.ts       = res_gaug_recess$df_gaug_rec$q97, 
                               ts.res       = res_gaug_recess$df_gaug_rec$MAP,
                               ts.real      = res_gaug_recess$df_gaug_rec$t.adj, 
                               station.name,
                               limni.labels, 
                               grid_limni.ylim,
                               dir.seg.gaug = dir.ST.test.step1,
                               seg.iter     = 1, 
                               t_Gaug       = gaugings$t, 
                               h_Gaug       = gaugings$h, 
                               pdf.ts       = res_gaug_recess$pdf_gaug_rec,
                               time_period  = c(res_gaug_recess$df_gaug_rec$q2[ii], 
                                                res_gaug_recess$df_gaug_rec$q97[ii]),
                               df_peaks     = h_peaks, 
                               save_name    = paste0("shift_event_", ii)))
             # remove or modify this shift time ?
             
             
             user.adjust.ts.ST =   dlgInput(c("Proposed shift time (days)= ", 
                                              res_gaug_recess$df_gaug_rec$t.adj[ii],
                                              "************************************* ",
                                              "-> ACCEPT ?   [Y]",
                                              "-> REMOVE?   [F]",
                                              "-> MODIFY MANUALLY? [insert the date]",
                                              "-> MODIFY TO ONE OF THE LIST BELOW (times of h peaks)? [insert index]",
                                              " ",
                                              paste0(seq(1, length(h_peaks$t_hmax), 1), "  = ", h_peaks$t_hmax, 
                                                     " days  (",  h_peaks$hmax, " m)")), 
                                            Sys.info()[" "])$res
             
             if (user.adjust.ts.ST == "Y"){
               print(paste0("time ", res_gaug_recess$df_gaug_rec$t.adj[ii], " accepted"))
               data.annotate.ST[ii,] = res_gaug_recess$df_gaug_rec[ii,]
             } else if (user.adjust.ts.ST == "F"){
               
               user.eliminate.ST = c(user.eliminate.ST, ii)
               print(paste0("time ", res_gaug_recess$df_gaug_rec$t.adj[ii], " removed"))
             } else {
               if (as.numeric(user.adjust.ts.ST) > 100){ 
                 data.annotate.ST[ii,]$t.adj = as.numeric(user.adjust.ts.ST)
                 print(paste0("time ", res_gaug_recess$df_gaug_rec$t.adj[ii], " adjusted to ", data.annotate.ST[ii,]$t.adj))
               } else {
                 #choose one of the peaks proposed in the list:
                 data.annotate.ST[ii,]$t.adj = h_peaks$t_hmax[as.numeric(user.adjust.ts.ST)]
                 print(paste0("time ", res_gaug_recess$df_gaug_rec$t.adj[ii], " adjusted to ", data.annotate.ST[ii,]$t.adj))
               }
             }
             dev.off()
           } 
           if (!is.null(user.eliminate.ST)){
             data.annotate.ST    = data.annotate.ST[-c(user.eliminate.ST),]
             pdf.ts.ST           = res_gaug_recess$pdf_gaug_rec[ , - user.eliminate.ST]
           }
           # plot new modified reference shift times:
           plot.shifts =  plot.time.shifts.step1(dir                  = dir.ST.test.step1 ,
                                                 gaug.1               = gaugings, 
                                                 rec.1                = NULL, 
                                                 data.annotate.off    = res_gaug_recess$df_offic, 
                                                 data.annotate.gaug.1 = res_gaug_recess$df_gaug,
                                                 data.annotate.rec.1  = res_gaug_recess$df_rec,
                                                 data.annotate.step1  = data.annotate.ST,
                                                 color                = colors.period, 
                                                 df.limni             = stage.record, 
                                                 limni.labels         = limni.labels , 
                                                 grid_limni.ylim      = grid_limni.ylim,
                                                 pdf.ts.1             = res_gaug_recess$pdf_gaug, 
                                                 pdf.ts.2             = res_gaug_recess$pdf_rec,
                                                 pdf.ts.3             = pdf.ts.ST,
                                                 save_file_name       = "/reference_shift_times_adjusted.pdf")
           
      } else {
        # no modifications:
        ###################
        data.annotate.ST = res_gaug_recess$df_gaug_rec
        pdf.ts.ST        = res_gaug_recess$pdf_gaug_rec                     
      }  
      # # save into files:
      # write.table(pdf.ts.ST, paste0(dir.sed.transp.SPD,"/STEP1_pdf_ts_for_ST.txt"),  sep ="\t", row.names=FALSE)
      # write.table(data.annotate.ST, paste0(dir.sed.transp.SPD,"/STEP1_shift_times_for_ST.txt"),  sep ="\t", row.names=FALSE)
      
      
      
      
  ############
  } else {   # if reference shift times are directed provided by user through a text file:
  ############
      print(c("You have chosen to consider as reference shift times for",
                   "the sediment transport analysis your .txt file ", file_name_ref_shift_times))
      user.file.ST = dlgInput(paste0("Is the file below ?  [Y/F] \n", 
                                       file_name_ref_shift_times),  Sys.info()[" "])$res
      if (user.file.ST == "Y"){
          ref.shift.times      = read.table(file =paste0(dir.case_study,"/shift_times.txt"), header = TRUE)
      } else {
          user.new.file.ST = dlgInput(paste0("Enter the .txt file name:"),  Sys.info()[" "])$res
          ref.shift.times      = read.table(file =paste0(dir.case_study,"/", user.new.file.ST, ".txt"), header = TRUE)
      }
        
      pdf.ts.file = matrix(NA, nrow = 10000, ncol=length(ref.shift.times$treal)); 
      names(ref.shift.times)[1] = "t.adj"
      for (i in 1:length(ref.shift.times$treal)){
        pdf.ts.file[,i] = runif(10000, min =ref.shift.times$t2, max=ref.shift.times$t97)
      }
      
      if (!is.null(dev.list())){
        while (!is.null(dev.list())) {
          dev.off()
        }
      }  
      X11(width = 20, height = 20); plot.shifts = plot.time.shifts.file( dir                 = dir.ST.test.step1 ,
                                                                        gaug.1               = gaugings, 
                                                                        data.annotate.off    = data.annotate.off, 
                                                                        data.annotate.file   = ref.shift.times,
                                                                        data.annotate.chosen = ref.shift.times,
                                                                        color                = colors.period, 
                                                                        df.limni             = stage.record, 
                                                                        limni.labels         = limni.labels , 
                                                                        grid_limni.ylim      = grid_limni.ylim,
                                                                        pdf.ts.file          = pdf.ts.file,
                                                                        pdf.ts.chosen        = pdf.ts.file,
                                                                        save_file_name       = "/reference_shift_times_fromfile.pdf"
      ); print(plot.shifts); user.choice.ST <- dlgInput(paste0("REFERENCE MORPHOGENIC EVENTS: \n",
                                                               "Do you want to remove or modify some shift times to define the reference morphological shifts ? [Y / N] \n" ), 
                                                        Sys.info()[" "])$res ;  dev.off();  if (user.choice.ST == "Y") {  
                                                        gaugings.P = gaugings; names(gaugings.P) = c("hP", "QP", "uQP", "Period", "tP", "t.trueP")
                                                        stage.record.NA = na.omit(stage.record)
                                                        
                                                        for (ii in 1:length(ref.shift.times$t.adj)){  
                                                          
                                                          # loop on all the selected events to find the flood peak:
                                                          int = which(stage.record.NA$t_limni >= ref.shift.times$q2[ii] &
                                                                      stage.record.NA$t_limni <= ref.shift.times$q97[ii])
                                                          h_peaks = data.frame(hmax   = round(stage.record.NA$h_limni[int][find_peaks(x = stage.record.NA$h_limni[int], m=deltat_peaks)], digits =2),
                                                                               t_hmax = round(stage.record.NA$t_limni[int][find_peaks(x = stage.record.NA$h_limni[int], m=deltat_peaks)], digits =2))
                                
                                                          message(paste0("event ", ii, ":"))
                                                          message("***************")
                                                          message("A few peaks (timings) in the CI of the event: ")
                                                          cat(paste0(seq(1, length(h_peaks$t_hmax),1) , ")  ", round(h_peaks$t_hmax,digits=2), " days   (h =", round(h_peaks$hmax, digits =3), " m)\n"))
                                                          message("")
                                                          
                                                          X11(width = 20, height = 10); print(exemple_ts.plot(CdT.P           = gaugings.P, 
                                                                                                              df.limni        = stage.record, 
                                                                                                              Q10.ts          = ref.shift.times$q2, 
                                                                                                              Q90.ts          = ref.shift.times$q97, 
                                                                                                              ts.res          = ref.shift.times$MAP,
                                                                                                              ts.real         = ref.shift.times$t.adj, 
                                                                                                              station.name    = station.name,
                                                                                                              limni.labels    = limni.labels, 
                                                                                                              grid_limni.ylim = grid_limni.ylim,
                                                                                                              dir.seg.gaug    = dir.ST.test.step1,
                                                                                                              seg.iter        = 1, 
                                                                                                              t_Gaug          = gaugings$t, 
                                                                                                              h_Gaug          = gaugings$h, 
                                                                                                              pdf.ts          = pdf.ts.file,
                                                                                                              time_period     = c(ref.shift.times$q2[ii], 
                                                                                                                                  ref.shift.times$q97[ii]),
                                                                                                              df_peaks        = h_peaks, 
                                                                                                              save_name       = paste0("shift_event_", ii)))
                                                          # remove or modify this shift time ?
                                                          user.eliminate.ST =   c() 
                                                          user.modify.ST    =   c()
                                                          data.annotate.ST  =   ref.shift.times
                                                          pdf.ts.ST         =   pdf.ts.file
                                                          user.adjust.ts.ST =   dlgInput(c("Proposed shift time (days)= ", 
                                                                                           ref.shift.times$t.adj[ii],
                                                                                           "************************************* ",
                                                                                           "-> Accept it ?   [Y]",
                                                                                           "-> Remove it ?   [F]",
                                                                                           "-> Modify it to your own date? [insert the date]",
                                                                                           "-> Modify it to one of the list below (times of h peaks)? [insert index]",
                                                                                           " ",
                                                                                           paste0(seq(1, length(h_peaks$t_hmax), 1), "  = ", h_peaks$t_hmax, 
                                                                                                  " days (",  h_peaks$hmax, " m)")), 
                                                                                         Sys.info()[" "])$res
                                                      
                                                          
                                                          if (user.adjust.ts.ST == "Y"){
                                                            print(paste0("time ", ref.shift.times$t.adj[ii], " accepted"))
                                                            data.annotate.ST[ii,] = ref.shift.times
                                                          } else if (user.adjust.ts.ST == "F"){
                                                            user.eliminate.ST = c(user.eliminate.ST, ii)
                                                            print(paste0("time ", ref.shift.times$t.adj[ii], " removed"))
                                                            
                                                          } else {
                                                            if (as.numeric(user.adjust.ts.ST) > 100){ 
                                                              data.annotate.ST[ii,]$t.adj = as.numeric(user.adjust.ts.ST)
                                                              print(paste0("time ", res_gaug_recess$df_gaug_rec$t.adj[ii], " adjusted to ", data.annotate.ST[ii,]$t.adj))
                                                            } else {
                                                              data.annotate.ST[ii,]$t.adj = h_peaks$t_hmax[as.numeric(user.adjust.ts.ST)]
                                                              print(paste0("time ", ref.shift.times$t.adj[ii], " adjusted to ", data.annotate.ST[ii,]$t.adj))
                                                            }
                                                          }  
                                                          dev.off()
                                                        }
                                                        
                                                        if (!is.null(user.eliminate.ST)){
                                                          data.annotate.ST    = data.annotate.ST[-c(user.eliminate.ST),]
                                                          pdf.ts.ST           = res_gaug_recess$pdf_gaug_rec[ , - user.eliminate.ST]
                                                        }
                                                        
                                                        plot.shifts =  plot.time.shifts.file( dir                  = dir.ST.test.step1 ,
                                                                                              gaug.1               = gaugings, 
                                                                                              data.annotate.off    = data.annotate.off, 
                                                                                              data.annotate.file   = ref.shift.times,
                                                                                              data.annotate.chosen = data.annotate.ST,
                                                                                              color                = colors.period, 
                                                                                              df.limni             = stage.record, 
                                                                                              limni.labels         = limni.labels , 
                                                                                              grid_limni.ylim      = grid_limni.ylim,
                                                                                              pdf.ts.file          = pdf.ts.file,
                                                                                              pdf.ts.chosen        = pdf.ts.ST,
                                                                                              save_file_name       = "/reference_shift_times_adjusted.pdf")
                                                        
                                                      } else {
                                                        data.annotate.ST = res_gaug_recess$df_gaug_rec
                                                        pdf.ts.ST        = res_gaug_recess$pdf_gaug_rec                     
                                                      }  
    
    
      # # save into files:
      # write.table(pdf.ts.ST, paste0(dir.sed.transp.SPD,"/STEP1_pdf_ts_for_ST.txt"),  sep ="\t", row.names=FALSE)
      # write.table(data.annotate.ST, paste0(dir.sed.transp.SPD,"/STEP1_shift_times_for_ST.txt"),  sep ="\t", row.names=FALSE)

      # data.annotate.ST = res_gaug_recess$df_gaug_rec
      # pdf.ts.ST        = res_gaug_recess$pdf_gaug_rec  
    
  }
  
  if (!is.null(dev.list())){
    while (!is.null(dev.list())) {
      dev.off()
    }
  }  
  
  print("****************")
  print("   All done!    ")
  print("****************")
  
  return(list(
               df_ST     = data.annotate.ST,
               pdf_ST    = pdf.ts.ST,
               df_offic  = data.annotate.off))
}   




















#############################################################################################
transp <- function(s,d50,h,hcr,beta,alpha, S0.sed.transp) {
#############################################################################################
  #Function that computes the sediment transport flux (per unit of time) ==> q = m3/s
  # riverbed elevation z = b2
  # (h-hcritic) or (y-ycritic) 
  # Application of formula of Meyer-Peter and Mueller (1948) :
  if (h > hcr) {
    q = alpha * ((s-1) * 9.81* d50^3)^0.5 * (S0.sed.transp/(d50*(s-1)))^beta * (h - hcr)^beta
  } else {
    q = 0
  }
  return(q)
}







###################################################################
hcritic = function(ycr, bt) {
###################################################################
  hcr = 0
  hcr = ycr + bt
  return(hcr)
}











###################################################################
ycritic = function(taucritic, d50, S0.sed.transp, s.sed.transp) {
###################################################################
  if (taucritic == "Brownlie")  {
    taucritic = 0.22 * d50^(-0.9) + 0.06* 10^(-7.7*d50^(-0.9))  
  } else if (taucritic == "Soulsby") {
    taucritic = 0.3 / (1 + 1.2* d50) + 0.055* (1 - exp(-0.02* d50))
  }
  ycritic = (s.sed.transp -1)*d50/S0.sed.transp * taucritic
  #####
  return(ycritic)
}           

















################################################################################################################################
sediment.transport.segmentation <- function(dir.sed.transp, 
                                            dir.sed.transp.SPD, 
                                            res.SPD,
                                            file.options.general,
                                            file.options.ST,
#                                            sequence.ST, 
                                            df.limni.ST, 
                                            officialShiftsTime) {
################################################################################################################################
   # Initialization:     
   message("################################################################")
   message("            Sediment transport analysis (Darienzo, 2021)        ")
   message("################################################################")
   
   df.limni.ST = na.omit(df.limni.ST)
   
   # message("For this application you are:
   # -  reading gaugings and periods from (segmentation results): \n  '../",  basename(dir.segm.results),"'
   # -  saving SPD results in: \n  '../", basename(dir.SPD.results),"'
   # -  using 'module_BaRatinSPD.r' and 'module_prior_propagation.r' 
   # 
   # *****************************************************************
   #   Configurating ...
   # *****************************************************************")
                  source(file.options.general)
                  source(file.options.ST)
                  ST.segm = NULL; 
                  setwd(dir.sed.transp)
                  setwd("./..")
                  dir.sed.transp.ST = getwd()
                  dir.sed.transp.d50 = paste0(dir.sed.transp.ST,"/2_transport_computation")
                  dir.create(dir.sed.transp.d50)
                  # dir.create(paste0(dir.sed.transp,"/d",d50))
                  # dir.sed.transp.d50 = paste0(dir.sed.transp,"/d",d50)
                  
                  
   # b evolution per period (results from Baratin SPD):
                  message("- Computing the bed evolution (parameter b) in time !!!  Wait ... "); flush.console()
                  # Dataframes with asymptote parameter (stage at time infinity):
                  bt1.df          = read.table(paste0(dir.sed.transp.SPD,"/bt", 1, "_df.txt"), header =TRUE)
                  if (control.main.channel >1){
                     bt2.df          = read.table(paste0(dir.sed.transp.SPD,"/bt", control.main.channel, "_df.txt"), header =TRUE)
                     # needs to be improved this to be generalized !!!!!!
                  } else {
                     bt2.df = bt1.df
                  } 
                  gaugings_SPD_ST = read.table(paste0(dir.sed.transp,"/data_with_periods.txt"), header =TRUE)
                  gaugings_SPD    = gaugings_SPD_ST
                  nperiods        = tail(gaugings_SPD_ST$Period,1)
                  #
                  ts.ST           = read.table(paste0(dir.sed.transp,"/shift_times.txt"), header =TRUE)
                  ts.before.ST    = sort(c(0, ts.ST$t.adj))
                  ts.plus.ST      = sort(c(ts.ST$t.adj, tail(df.limni.ST$t_limni,1)))
                  bt1.df.ST       = bt1.df  #[c(sequence.ST,tail(sequence.ST,1)+1),]
                  bt2.df.ST       = bt2.df  #[c(sequence.ST,tail(sequence.ST,1)+1),]
                  nperiods.ST     = length(bt2.df.ST$maxpost)
                  
                  
  # Defining hcr and ycr for each stable period:
                  bt            = NULL 
                  hcrit         = NULL
                  nsim          = 10000
                  S0.sed.transp = S0.st
                  s.sed.transp  = s.st
                  d50.st        = d50.st
                  phi.st        = d50.st/S0.st
                  
                  
                  
                  set.seed(bt2.df.ST$mean)  
                  for (pp in 1:nperiods.ST){
                        bt[[pp]]    = rnorm(nsim, mean=bt2.df.ST$mean[pp], sd=bt2.df.ST$stdev[pp])
                        hcrit[[pp]] = vector(mode='double',length=nsim); 
                  }
                  for(i in 1:nsim){
                      # Compute a for each control
                      for (pp in 1:nperiods.ST) {
                          hcrit[[pp]][i]= hcritic(ycr = ycritic(taucritic,
                                                                d50.st, 
                                                                S0.sed.transp, 
                                                                s.sed.transp), 
                                                  bt[[pp]][i])
                      }
                  }
                  sim   = data.frame(hcrit[[1]])
                  for (pp in 2:nperiods.ST) {
                      sim = cbind(sim,  hcrit[[pp]])
                  }
                  # check sizes match and if so, initialize:
                  n         = NCOL(sim)
                  hcr       = vector("list",n)
                  transform = matrix(NA,nsim,n)
                  # start computations for each margin:
                  for(i in 1:n){
                      # Marginal prior parameters
                      hcr[[i]] = switch("Gaussian",
                                        'Gaussian'  = {c(mean(sim[,i]),      sd(sim[,i]))},
                                        'LogNormal' = {c(mean(log(sim[,i])), sd(log(sim[,i])))},
                                        NA)
                      # Transform simulations into Gaussian space to estimate the correlation of the Gaussian copula
                      transform[,i] = switch("Gaussian",
                                             'Gaussian' = {qnorm(pnorm(sim[,i],mean=sim[[i]][1],sd=sim[[i]][2]))},
                                             'LogNormal'= {qnorm(plnorm(sim[,i],meanlog=sim[[i]][1],sdlog=sim[[i]][2]))},
                                             NA)
                  }

                  #corel=cor(transform)
                  hcr.df = data.frame(t(sapply(hcr, function(x) x[1:max(lengths(hcr))])))
                  names(hcr.df) = c("mean", "stdev")
                  
                  
                  # plot hcr.df per period:
                  hcr.plot = ggplot(data = hcr.df)+ 
                             geom_point(aes(x=seq(1,length(hcr.df$mean),1), y=mean), color="green")+
                             geom_line(aes(x=seq(1,length(hcr.df$mean),1), y=mean), color="green")+
                               geom_errorbar(aes(x=seq(1,length(hcr.df$mean),1), 
                                                 ymin=mean-2*stdev, 
                                                 ymax=mean+2*stdev), color="green")+
                               geom_point(aes(x=seq(1,length(bt2.df.ST$mean),1), y=bt2.df.ST$mean), color="red")+
                               geom_line(aes(x=seq(1,length(bt2.df.ST$mean),1), y=bt2.df.ST$mean), color="red")+
                               geom_errorbar(aes(x    = seq(1,length(bt2.df.ST$mean),1), 
                                                 ymin = bt2.df.ST$mean - 2*bt2.df.ST$stdev, 
                                                 ymax = bt2.df.ST$mean + 2*bt2.df.ST$stdev), 
                                                 color= "red")+
                               ylab("hcritical")+ xlab("Period")+
                               theme_bw()
                             pdf(paste0(dir.sed.transp.d50,"/hcr_per_period.pdf"), 14, 8 ,useDingbats=F)
                             print(hcr.plot)
                             dev.off()

                  
                    
 
  # Plot b(t) and hcrit per period
                  message("- Plotting the b and hcr per period against the stage record.")
                  bt_and_hcrit.plot(dir.sed.transp  = dir.sed.transp.d50, 
                                    ts              = ts.ST,
                                    ts.before       = ts.before.ST, 
                                    ts.plus         = ts.plus.ST, 
                                    bt2.df          = bt2.df.ST, 
                                    bt1.df          = bt1.df.ST,
                                    hcr.df          = hcr.df,
                                    df.limni        = df.limni.ST, 
                                    gaugings_SPD    = gaugings_SPD,
                                    ylimits         = ylimits, 
                                    d50             = d50.st)
  # plot a zoom on the first event as an example :
                  # bt_and_hcrit_zoom.plot(dir.sed.transp,
                  #                        ts.ST, 
                  #                        ts.before.ST,
                  #                        ts.plus.ST, 
                  #                        bt2.df.ST,
                  #                        hcr.df, 
                  #                        df.limni.ST, 
                  #                        gaugings_SPD,
                  #                        ylimits, 
                  #                        d50)
  # defining shifting periods:
                  # Define the index j of the shifting events
                  t.shifts.limni = NULL; t  = 1; j  = NULL;  i  = 1;  ev = 1
                  while (df.limni.ST$t_limni[t] < ts.ST$t.adj[length(ts.ST$t.adj)])  {
                  #*******************************************************************
                      if ((df.limni.ST$t_limni[t] <= ts.ST$t.adj[i]) &
                          (df.limni.ST$t_limni[t+1] >= ts.ST$t.adj[i])) {
                      #**************************************************
                          t.shifts.limni[i] = df.limni.ST$t_limni[t+1]
                          j[i] = t+1   # Define the indexes of the time vector for peak t.shift 
                                       # (it is the time data before the shift time)
                          i = i +1     # i is the number of inter-event stable periods
                      }
                      t = t+1          # t is the time index of the limni
                  }
                  
                  
                  
                  
                  
                  # for each reference shift time define the event period (left and right respect to the peak):  
                  ##############################################################################################
                  event.t = NULL; event.h =NULL; start.event.index=NULL; end.event.index =NULL;
                  message("- For each reference shift time define the whole event period")
                  message("  thus, the left (increasing phase) and the right (decreasing phase) respect to the peak.")
                  for (i in 1:length(ts.ST$t.adj)) { 
                  #*********************************
                      k =0 ; jj = j[i]; 
                      event.left.t = NULL; event.left.h = NULL
                      # from the peak go to left in the stage record until h < hcrit[i]:
                      while ((df.limni.ST$h_limni[jj] > hcr.df$mean[i])) {
                          #print(jj)
                          jj= jj-1
                          k = k + 1
                          event.left.t[k] = df.limni.ST$t_limni[jj]
                          event.left.h[k] = df.limni.ST$h_limni[jj]
                      }
                      start.event.index[i]=jj
                      #from the peak go to right in the stage record until h < hcrit[i+1]:
                      k =0; jj = j[i]; event.right.t = NULL; event.right.h = NULL;
                      while (df.limni.ST$h_limni[jj] >  hcr.df$mean[i+1]) {
                          jj = jj +1
                          k = k + 1
                          event.right.t[k] = df.limni.ST$t_limni[jj]
                          event.right.h[k] = df.limni.ST$h_limni[jj]
                      }
                      end.event.index[i] = jj
                      event.t[[ev]] = c(sort(event.left.t), event.right.t)
                      event.h[[ev]] = c(sort(event.left.h), event.right.h)
                      
                      if (is.null(event.right.t)==TRUE){
                        calibration.ON = FALSE
                        print(paste0("********* Event flood at = ", ts.ST$t.adj[i], "days has been missed !!! "))
                        message("- Please, check your input settings. You may need to change some parameters:")
                        message("  e.g.,")
                        message("  . decrease the characteristic sediment diameter d50")
                        message("  . or, increase the longitudinal slope of the bed S0")
                        ev = ev
                      } else {
                        ev = ev+1
                      }
                  }  
                  
                  
                  
                  
                    
                  
                  
    ############################            
    # if (calibration.ON == TRUE){
    # Evolution of b over  time  (Linear interpolation during the shifting events):
                  b = NULL; jj = 1;  # k=1
                  while (df.limni.ST$t_limni[jj] < tail(df.limni.ST$t_limni,1)) {
                  #**************************************************************
                      if (df.limni.ST$t_limni[jj] < event.t[[1]][1]) {
                          b[jj] = bt2.df.ST$mean[1]
                      } else if (df.limni.ST$t_limni[jj] > tail(event.t[[length(event.t)]],1)) {
                          b[jj] = bt2.df.ST$mean[length(bt2.df.ST$mean)]
                      }
                      if ((df.limni.ST$t_limni[jj] >= event.t[[1]][1]) & (df.limni.ST$t_limni[jj] <= tail(event.t[[1]],1))) {
                          b[jj] = bt2.df.ST$mean[1] + (df.limni.ST$t_limni[jj]-event.t[[1]][1]) / (tail(event.t[[1]],1) 
                                                   - event.t[[1]][1])* (bt2.df.ST$mean[2] - bt2.df.ST$mean[1]) 
                      } 
                      for (i in 2:length(event.t)) {
                          if ((df.limni.ST$t_limni[jj] >= event.t[[i]][1]) & (df.limni.ST$t_limni[jj] <= tail(event.t[[i]],1))) {
                              b[jj] = bt2.df.ST$mean[i] + (df.limni.ST$t_limni[jj]-event.t[[i]][1]) /
                                                       (tail(event.t[[i]],1)-event.t[[i]][1]) * (bt2.df.ST$mean[i+1] - bt2.df.ST$mean[i]) 
                              # event.b[[k]][] = 
                              # k = k+1
                              # j
                          } else if ((df.limni.ST$t_limni[jj] > tail(event.t[[i-1]],1)) &  (df.limni.ST$t_limni[jj] < event.t[[i]][1])) {
                              b[jj] =  bt2.df.ST$mean[i]
                          }
                      }
                      jj= jj+1
                  }
                  b[jj] = bt2.df.ST$mean[length(bt2.df.ST$mean)]    #last value !
                  
  #defining the new hcrit(t) over time:
                  hcr = NULL
                  message("- Compute the critical stage time series hcr(t)")
                  hcr = b + ycritic(taucritic, 
                                    d50.st, 
                                    S0.sed.transp,
                                    s.sed.transp)
                  
  # Plot the first morphogenic event with transient b(t) and transietn hcrit:
                  
                  bt_hcrit_interpol_zoom.plot(dir.sed.transp.d50, 
                                              b,
                                              hcr, 
                                              event.t, 
                                              event.h, 
                                              ts.ST, 
                                              ts.before.ST, 
                                              ts.plus.ST,
                                              bt2.df.ST, 
                                              hcr.df, 
                                              df.limni.ST, 
                                              gaugings_SPD, 
                                              hlimits.ST = ylimits,
                                              d50.st)
                  
  #plot the water depth y time series (y =h-b) :
                  message("- Compute the water depth y(t) time series: y(t) = h(t) - b(t) [m])")
                  y_limni     = df.limni.ST$h_limni - b
                  df.limni.ST = cbind(df.limni.ST)
                  ycr         = ycritic(taucritic, 
                                        d50.st,
                                        S0.sed.transp, 
                                        s.sed.transp)
                  
                  # estimate the b associated to each gauging to obtain y_Gaug:
                  b_Gaug =0; g=1 ; y_Gaug =0
                  message("- Estimate the parameter 'b' associated to each gauging to obtain y_Gaug (water depth).")
                  for (g in 1:length(t_Gaug)) {
                    for (ll in 1:(length(df.limni.ST$t_limni)-1)) {
                      if ((t_Gaug[g]>= df.limni.ST$t_limni[ll]) & (t_Gaug[g]<= df.limni.ST$t_limni[ll+1])) {
                        b_Gaug[g] = b[ll]
                        y_Gaug[g] = h_Gaug[g] - b_Gaug[g]
                      }
                    }
                  }
                  y_Gaug[which(y_Gaug <0)] =0
                  
                  # y_limni.plot(dir.sed.transp = dir.sed.transp.d50, 
                  #              gaugings_SPD, 
                  #              df.limni, 
                  #              y_limni,
                  #              y_Gaug, 
                  #              ycr, 
                  #              ylimits=c(0,5), 
                  #              d50) 
                  
                  
######################################################################################################
                  # for each shift time define the event period (left and right respect to the peak):  
                  event.t = NULL; event.h =NULL; start.event.index=NULL; end.event.index =NULL;
                  message("- For each shift time define the event period (left and right respect to the peak)")
                  message("  Please, wait ...")
                  pb <- txtProgressBar(min = 0,               # Minimum value of the progress bar
                                       max = length(ts.ST$t.adj), # Maximum value of the progress bar
                                       style = 3,             # Progress bar style (also available style = 1 and style = 2)
                                       width = 50,            # Progress bar width. Defaults to getOption("width")
                                       char = "=")            # Character used to create the bar
                  
                  
                  for (i in 1:length(ts.ST$t.adj)){   #length(sequence.ST)) { 
                      k =0 ; jj = j[i]; event.left.t = NULL; event.left.h = NULL
                      # from the peak go to left in the stage record until h < hcrit[i]:
                      while (df.limni.ST$h_limni[jj] > hcr[jj]) {
                         jj= jj-1
                         k = k + 1
                         event.left.t[k] = df.limni.ST$t_limni[jj]
                         event.left.h[k] = df.limni.ST$h_limni[jj]
                      }
                      start.event.index[i] = jj
                      #from the peak go to right in the stage record until h < hcrit[i+1]:
                      k =0; jj = j[i]; event.right.t = NULL; event.right.h = NULL;
                      while (df.limni.ST$h_limni[jj] >  hcr[jj]) {
                         jj = jj +1
                         k = k + 1
                         event.right.t[k] = df.limni.ST$t_limni[jj]
                         event.right.h[k] = df.limni.ST$h_limni[jj]
                      }
                      end.event.index[i] = jj
                      event.t[[i]] = c(sort(event.left.t), event.right.t)
                      event.h[[i]] = c(sort(event.left.h), event.right.h)
                      setTxtProgressBar(pb, i)
                  }  
                  close(pb)
                  
                  
                  
                  
  #Application of sediment transport model qs(t) ==> MPM model:
                  message("- Computing the bedload qs(t) at instant t !!!  Wait ... "); flush.console()
                  qs = NULL
                  pb <- txtProgressBar(min = 0,               # Minimum value of the progress bar
                                       max = length(df.limni.ST$t_limni), # Maximum value of the progress bar
                                       style = 3,             # Progress bar style (also available style = 1 and style = 2)
                                       width = 50,            # Progress bar width. Defaults to getOption("width")
                                       char = "=")            # Character used to create the bar
                  
                  for (i in 1:length(df.limni.ST$t_limni)) {
                      qs[i] = transp(s             = s.st, 
                                     d50           = d50.st, 
                                     h             = df.limni.ST$h_limni[i], 
                                     hcr           = hcr[i],
                                     beta          = beta.sed.transp, 
                                     alpha         = alpha.sed.transp, 
                                     S0.sed.transp = S0.st)
                      setTxtProgressBar(pb, i)
                  }  
                  close(pb)
                    
                  
                  
  # Cumulative sediment transport qsc(m) for each event m selected for the analysis:         
                  event.hcritic =NULL; b.event = NULL; y.event =NULL; t.event =NULL; h.event =NULL;
                  qs.event.tot =NULL; cum.qs= NULL;
                  message("- Computing the cumulative bedload 'qsc(k)' for each morphogenic event k!  Wait ... "); flush.console()
                  for (i in 1:length(ts.ST$t.adj)){   #length(sequence.ST)) { 
                      cumulat=0
                      b.event[[i]] = b[start.event.index[i] : (end.event.index[i]-1)]
                      h.event[[i]] = h_limni[start.event.index[i] : (end.event.index[i]-1)]
                      y.event[[i]] = y_limni[start.event.index[i] : (end.event.index[i]-1)]
                      t.event[[i]] = df.limni.ST$t_limni[start.event.index[i] : (end.event.index[i]-1)]
                      qs.event =NULL; deltat=0;
                      qs.event[1]=0
                      for (jjj in 2:length(y.event[[i]])) {
                          qs.effect = transp(s     = s.st, 
                                             d50   = d50.st, 
                                             h     = h.event[[i]][jjj], 
                                             hcr   = hcr[start.event.index[i]+jjj-1],
                                             beta  = beta.sed.transp, 
                                             alpha = alpha.sed.transp, 
                                             S0.sed.transp)   #[m3/s]
                          qs.event[jjj] = qs.effect
                          deltat = (t.event[[i]][jjj] - t.event[[i]][jjj-1])*86400   #to have the time step in seconds  ==> qs,cum [m3]
                          cumulat = cumulat + qs.effect*deltat*Bc.st
                      }
                      qs.event.tot[[i]] = qs.event 
                      cum.qs[i] = cumulat  #volume [m3]
                  }  

                   
  
# identify among all the known shift events the one that has the lowest qs cumulative:
                  qsc_crit   = min(cum.qs)
                  message("*********************************************")
                  print(paste0("V_crit (minimum sedim. volume observed) = ", round(qsc_crit, digits=3), " m^3"))
                  message("*********************************************")
                  event.crit = which.min(cum.qs)

                  
                  
                  
                  
                  
##########################################################################################################
# searching for all other events with y >= y_crit and computing qsc :
                  eve=0; cumqs=0; qs.cum=NULL; start.event.indexes =NULL; end.event.indexes =NULL;
                  qs.cum[1]=0; dtt =0; cumqs=0; dtt[1]=0; qscum.tot =0;  Vt = 0;
                  message("- searching for all other events with y >= y_crit and computing qsc.")
                  
                  for (i in  2:(length(df.limni.ST$t_limni))) {
                  #********************************************
                      dtt[i]       = (df.limni.ST$t_limni[i] - df.limni.ST$t_limni[i-1])*86400
                      qscum.tot[i] = qscum.tot[i-1] + qs[i]*dtt[i]*Bc.st
                      Vt[i]        = qscum.tot[i-1] + qs[i]*dtt[i]*Bc.st
                      Vt[i]        = Vt[i-1] + qs[i]*dtt[i]*Bc.st
                      
                      if ((qs[i]>0 & qs[i-1]==0)) {
                         eve       = eve+1
                         start.event.indexes[eve] =i
                         cumqs     = qs[i]*dtt[i]*Bc.st
                         
                      } else if (qs[i-1]> 0 & qs[i]> 0) {
                         eve       = eve
                         cumqs     = qs[i]*dtt[i]*Bc.st + cumqs
                         
                      } else if (((qs[i-1]>0 & qs[i]==0))) {
                         Vt[i]     =  0
                         eve       = eve
                         end.event.indexes[eve] = i
                         qs.cum[eve] = cumqs
                         
                      } else {
                         Vt[i] =  0
                         cumqs = 0
                      }
                  }
                  
# merge events with a distance fixed:
                  print(paste0("- Merging events very close (time distance < ", tmax_between_events, " days)"))
                  message("  If you prefer, change the parameter 'tmax_between_events' in the options file.")
                  # Initialization:
                  finalqs.cum        = NULL 
                  eve1               = 1 
                  kkk                = 2 
                  finalqs.cum[1]     = qs.cum[1];  
                  start.event.new    = NULL; 
                  end.event.new      = NULL;
                  start.event.new[1] = start.event.indexes[1]; 
                  end.event.new[1]   = end.event.indexes[1];
                  cum.qs.merged      = cum.qs
                  
                  while (kkk < length(start.event.indexes)) {
                  #******************************************
                    if ((abs(df.limni.ST$t_limni[end.event.indexes[kkk]] - df.limni.ST$t_limni[start.event.indexes[kkk+1]]) <= tmax_between_events)) { 
                      #max of ....  2 ....days between twp consecutive events !!!! 
                      # if ((any(abs(df.limni$t_limni[end.event.index] - df.limni$t_limni[start.event.indexes[kkk+1]]) <= tmax_between_events))) {
                      #      ind.ref.merge = tail(which(any(abs(df.limni$t_limni[end.event.index] - df.limni$t_limni[start.event.indexes[kkk+1]]) <=  tmax_between_events)), 1)
                      #      cum.qs.merged[ind.ref.merge] =  cum.qs[ind.ref.merge]  + qs.cum[kkk]
                      # }
                      finalqs.cum[eve1]     = finalqs.cum[eve1] + qs.cum[kkk]
                      start.event.new[eve1] = start.event.new[eve1]
                      end.event.new[eve1]   = end.event.new[eve1]
                    } else {
                      
                      eve1                  = eve1 + 1 
                      finalqs.cum[eve1]     = qs.cum[kkk]
                      start.event.new[eve1] = start.event.indexes[kkk]
                      end.event.new[eve1]   = end.event.indexes[kkk]
                    }
                    kkk=kkk+1
                  }

                  
                  ### reference events:
                  diff          = NULL
                  ind.ref.merge = NULL 
                  cum.qs.merged = 0
                  start.event.index.merge = start.event.index
                  end.event.index.merge   = end.event.index
                  
                  for (jkl in 1:length(end.event.index)){
                  #**************************************
                      diff[[jkl]]                  = abs(df.limni.ST$t_limni[end.event.index[jkl]] - df.limni.ST$t_limni[start.event.indexes])
                      ind.ref.merge[[jkl]]         = which(diff[[jkl]] < tmax_between_events)
                      cum.qs.merged[jkl]           = sum(qs.cum[ind.ref.merge[[jkl]]])
                      start.event.index.merge[jkl] = min(start.event.indexes[ind.ref.merge[[jkl]]])
                      end.event.index.merge[jkl]   = max(end.event.indexes[ind.ref.merge[[jkl]]])
                  }
                  start.event.indexes.merge        = start.event.indexes[-c(unlist(ind.ref.merge))]
                  end.event.indexes.merge          = end.event.indexes[-c(unlist(ind.ref.merge))]
                  qs.cum.merge                     = qs.cum[-c(unlist(ind.ref.merge))]
                  
                  
#Plot the sediment transport qs(t):  
                  #df.event.potent=NULL
                  # qs.plot(dir.sed.transp.d50 = dir.sed.transp.d50, 
                  #         df.limni           = df.limni,
                  #         y_limni            = y_limni, 
                  #         qs                 = qs, 
                  #         ycr                = ycr, 
                  #         hcr                = hcr,
                  #         ts                 = ts, 
                  #         t.event            = t.event, 
                  #         y.event            = y.event, 
                  #         h.event            = h.event,
                  #         start.event.new    = start.event.new,
                  #         end.event.new      = end.event.new, 
                  #         ylimits.st         = c(0,5), 
                  #         d50                = d50,
                  #         gaugings_SPD       = gaugings_SPD, 
                  #         y_Gaug             = y_Gaug,
                  #         qscum.tot          = qscum.tot)
                  
    
# testing different phi =d/S0 for chapter manuscript !!!!!!!!!!!!!!:
#*******************************************************************
                  # summarizing results:
                  phi1                = phi.st
                  ycr1                = ycr
                  hcr1                = hcr
                  qsc_crit1           = qsc_crit
                  qs1                 = qs
                  start.event.new1    = start.event.indexes
                  end.event.new1      = end.event.indexes
                  finalqs.cum1        = qs.cum     #finalqs.cum
                  t.event1            = t.event
                  y.event1            = y.event
                  h.event1            = h.event
                  start.new1          = start.event.indexes.merge
                  end.new1            = end.event.indexes.merge
                  qscum.tot1          = qscum.tot
                  cum.qs1             = cum.qs.merged
                  start.event.index1  = start.event.index.merge
                  end.event.index1    = end.event.index.merge
                  bt1.df.ST.merge     = bt1.df.ST
                  bt2.df.ST.merge     = bt1.df.ST
                  
              
                  
                  message("- Total number of potential shift times (with ST analysis):")
                  print(length(start.event.new1))
                  message("- Saving potential shift times to txt file:")
                  print(paste0(dir.sed.transp.d50,"/all_ts_potential.txt"))
                  potential.shift.times = data.frame(ts.pot.in  = start.event.new1,
                                                     ts.pot.fin = end.event.new1)
                  write.table(potential.shift.times, paste0(dir.sed.transp.d50,"/all_ts_potential.txt"), sep ="\t", row.names=FALSE)
                  
                  
                  
                  message("Plotting final results: detected shift times and Volumes V")
                  message("Please, wait (it might take a few minutes) ...")
                  

                  qsc.t.plot.one.choice(dir.sed.transp.d50  = dir.sed.transp.d50, 
                                        df.limni.ST         = df.limni.ST, 
                                        y_limni             = y_limni, 
                                        ylimits.st          = ylimits, 
                                        hlimits.st          = c(ylimits[1], ylimits[2]), 
                                        gaugings_SPD        = gaugings_SPD,
                                        t.event1            = t.event1, 
                                        y.event1            = y.event1, 
                                        h.event1            = h.event1,
                                        phi1                = phi1,
                                        hcr1                = hcr1,
                                        ycr1                = ycr1,
                                        start.event.index1  = start.event.index1,
                                        end.event.index1    = end.event.index1,
                                        start.event.new111  = start.event.new1, 
                                        end.event.new111    = end.event.new1,
                                        cum.qs1             = cum.qs1,           # cum.qs.merged
                                        finalqs.cum1        = finalqs.cum1,      # finalqs.cum,
                                        start.new111        = start.event.new1 , # start.new1, 
                                        end.new111          = end.event.new1,    # end.new1,
                                        qscum.tot1          = qscum.tot1,
                                        plot.event.index    = FALSE,
                                        V.limits            = V.limits,
                                        V.for.text          = V.limits[2]) 
                  
                  
# the Volumes V associated to each reference event:
##################################################
                  df.sed.transport = data.frame(index.tstart =  start.event.index1,
                                                inde.tend    =  end.event.index1,
                                                V            =  cum.qs1)
                  
                  # this part is specific to Meyras hydraulic configuration with 2 low controls:!!!!!!                
                  delta.b1.MAP =NULL; delta.b1.stdev = NULL; delta.b1.mean =NULL;
                  for (evento in 1:(length(bt1.df.ST.merge$maxpost)-1)) {
                        delta.b1.MAP[evento]   = (bt1.df.ST.merge$maxpost[evento+1] - bt1.df.ST.merge$maxpost[evento])
                        delta.b1.mean[evento]  = (bt1.df.ST.merge$mean[evento+1] - bt1.df.ST.merge$mean[evento])
                        delta.b1.stdev[evento] = (bt1.df.ST.merge$stdev[evento+1]^2 + bt1.df.ST.merge$stdev[evento]^2)^0.5
                  }
                  
                  delta.b2.MAP =NULL; delta.b2.stdev = NULL; delta.b2.mean =NULL;
                  for (evento in 1:(length(bt2.df.ST.merge$maxpost)-1)) {
                       delta.b2.MAP[evento]   = (bt2.df.ST.merge$maxpost[evento+1] - bt2.df.ST.merge$maxpost[evento])
                       delta.b2.mean[evento]  = (bt2.df.ST.merge$mean[evento+1]    - bt2.df.ST.merge$mean[evento])
                       delta.b2.stdev[evento] = (bt2.df.ST.merge$stdev[evento+1]^2 + bt2.df.ST.merge$stdev[evento]^2)^0.5
                  }
                  df.rel.deltab.V.TOT = data.frame(index.tstart  =  start.event.index1,
                                                   inde.tend     =  end.event.index1,
                                                   V             =  cum.qs1,
                                                   deltab1       = delta.b1.mean,
                                                   stdev_deltab1 = delta.b1.stdev,
                                                   deltab2       = delta.b2.mean,
                                                   stdev_deltab2 = delta.b2.stdev,
                                                   Period        = seq(1,length(delta.b1.mean),1))
                  
                  #*************************************************************************************
                  # df.rel.deltab.qscum = df.rel.deltab.V.TOT[which(data.annotate.step2.sort$index>0,),]
                  # df.rel.deltab.qscum = cbind(data.annotate.ST, df.rel.deltab.qscum)

                  db1_db2_V.plot = ggplot(data=df.rel.deltab.V.TOT) +
                                   geom_point(aes(y=deltab1, x= V), color ="red" , size=2)+
                                   geom_errorbar(aes(ymin = deltab1 - 2*stdev_deltab1,  ymax = deltab1 + 2*stdev_deltab1,  x = V),
                                                 color = "red",  size = 0.2, width=20)+
                                   geom_point(aes(y=deltab2, x = V),   color ="blue" , size=2)+
                                   geom_errorbar(aes(ymin = deltab2 - 2*stdev_deltab2,  ymax = deltab2 + 2*stdev_deltab2, x = V),
                                                 color = "blue",  size = 0.2, width=20)+
                                   annotate("text",  x=df.rel.deltab.V.TOT$V, y=df.rel.deltab.V.TOT$deltab1 + df.rel.deltab.V.TOT$stdev_deltab1+0.05,
                                            label= df.rel.deltab.V.TOT$Period, color = "black", size=2) +
                                   geom_point(data =data.frame(x=df.rel.deltab.V.TOT$V,  y=df.rel.deltab.V.TOT$deltab1 + df.rel.deltab.V.TOT$stdev_deltab1+0.05),
                                              aes(x=x, y=y),  size = 3, color="black", pch=21, fill="NA")+
                                   ylab("Delta b1, Delta b2 [m]") + xlab("Cumulative bedload, V [m3]")+
                                   theme_bw(base_size = 15)
                                   ggsave(db1_db2_V.plot, filename=paste0(dir.sed.transp.d50,"/db1_db2_V.png"),bg = "transparent", device="png", width=6, height=4, dpi=400)
                                   
                                   

                    
                    
                
                  print("****************")
                  print("   All done!    ")
                  print("****************")
                  
                  
                  
# Main Outputs:
###############
                  return(list(ycritic = ycr , 
                              hcritic = hcr,
                              Vt      = Vt,
                              db1_db2_V.plot  = db1_db2_V.plot,
                              df.rel.deltab.V.TOT  = df.rel.deltab.V.TOT,
                              potential.shift.time  = potential.shift.times
                             ))
}

                                
    






























#########################################################################################################
linear.estimation = function(dir_code, 
                             dir.sed.transp,
                             df.deltab.V,
                             file.options.general,
                             file.options.ST) {
#########################################################################################################
  # with the results of the SPD estimation , estimate a relation between Delta(b1,b2) and qs,cum
  # Use a linear regression with BaM.
  # Apply a BaM Linear Regression estimation (with data uncertainties).
  
  # Preparation of the dataframes:
  df.linear.regress.deltab1 = data.frame(X  = df.deltab.V$V, 
                                         uX = rep(0,length(df.deltab.V$V)), 
                                         y  = df.deltab.V$deltab1, 
                                         uY = df.deltab.V$stdev_deltab1)
  
  df.linear.regress.deltab2 = data.frame(X= df.deltab.V$V, 
                                         uX = rep(0,length(df.deltab.V$V)), 
                                         y = df.deltab.V$deltab2,
                                         uY = df.deltab.V$stdev_deltab2)
  #directories:     
  setwd(dir.sed.transp)
  setwd("./..")
  dir.sed.transp.ST = getwd()
  dir.sed.transp.linear = paste0(dir.sed.transp.ST,"/3_linear_relation")
  dir.create(dir.sed.transp.linear)
  dir.config.linear     = paste0(dir_code,"/BaM_exe/Linear")
  dir.qs.delta_b1       = paste0(dir.sed.transp.linear, "/delta_b1")
  dir.qs.delta_b2       = paste0(dir.sed.transp.linear, "/delta_b2")
  dir.results.qs.deltab = dir.sed.transp.linear
  dir.create(paste0(dir.sed.transp.linear, "/delta_b1"))
  dir.create(paste0(dir.sed.transp.linear, "/delta_b2"))
  
  
  # user options:                 
  source(file.options.general)
  source(file.options.ST)
  xgrid   = seq(V.limits[1], V.limits[2], V.limits[3]) #(V.limits[2] -V.limits[1])/20)
  Nmcmc   = Nmcmc.st
  Ncycles = Ncycles.st
  # delta b1:
  
  write.table(df.linear.regress.deltab1, file=paste0(dir.config.linear,"/data.txt"), sep="\t", row.names=FALSE )
  setwd(dir.exe)
  #launch BaM :
  linear_app(model.type         = model.type ,
             Nmcmc              = Nmcmc, 
             Ncycles            = Ncycles, 
             xgrid              = xgrid,
             nobs               = length(df.linear.regress.deltab1$X), 
             simMCMC            = TRUE, 
             prediction         = TRUE, 
             prediction.t       = FALSE, 
             nobs.t             = 1,
             prior.param.propor = prior.param.propor,
             prior.gamma.linear = prior.gamma.linear) 
  #read results :
  env.linear.regress.b1       = read.table(file = paste0(dir.config.linear,"/Lin_TotalU.env"), header=TRUE)
  env.linear.regress.b1       = cbind(env.linear.regress.b1, xgrid)
  env.param.linear.regress.b1 = read.table(file = paste0(dir.config.linear,"/Lin_ParamU.env"), header=TRUE)
  env.param.linear.regress.b1 = cbind(env.param.linear.regress.b1, xgrid)
  maxpost.linear.regress.b1   = read.table(file = paste0(dir.config.linear,"/Lin_Maxpost.spag"), header=FALSE)
  maxpost.linear.regress.b1   = cbind(maxpost.linear.regress.b1, xgrid)
  #write and save results of this linear regression Deltab1,b2 vs qscum :
  list.of.files <- c(
    paste0(dir.config.linear,"/Lin_maxpost.spag"),
    paste0(dir.config.linear,"/Lin_ParamU.spag"), 
    paste0(dir.config.linear,"/Lin_TotalU.spag"),
    paste0(dir.config.linear,"/Lin_ParamU.env"),  
    paste0(dir.config.linear,"/Lin_TotalU.env"),
    paste0(dir.config.linear,"/Results_MCMC_Cooked.txt"), 
    paste0(dir.config.linear,"/Results_Residuals.txt"),
    paste0(dir.config.linear,"/Results_Summary.txt"), 
    paste0(dir.config.linear,"/Config_Model.txt"),
    paste0(dir.config.linear,"/xgrid.txt"), 
    paste0(dir.config.linear,"/data.txt")
  )
  for (ll in 1:length(list.of.files)) {
    file.copy(list.of.files[ll], dir.qs.delta_b1, overwrite = TRUE) 
  }
  
  #converg = Convergence.test(dir.seg = dir.qs.delta_b1, npar = 2, dir.plot = dir.qs.delta_b1)
  # plot.mcmc.seg(workspace=dir.qs.delta_b1, seg.iter=1, nS=2)
  
  # plot results :
  lin.relat.deltab1.plot = ggplot(data=df.deltab.V) +
    geom_ribbon(data = env.linear.regress.b1, aes(x = xgrid, ymin =Q_q2.5, ymax = Q_q97.5), fill = "blue", alpha = 0.1)+
    geom_ribbon(data = env.param.linear.regress.b1, aes(x = xgrid, ymin =Q_q2.5, ymax = Q_q97.5), fill = "blue", alpha = 0.2)+
    geom_errorbar(aes(ymin= df.linear.regress.deltab1$y - 2*df.linear.regress.deltab1$uY, 
                      ymax= df.linear.regress.deltab1$y + 2*df.linear.regress.deltab1$uY, x   = df.linear.regress.deltab1$X ),
                  color = "black",  size = 0.3, width=2000)+
    geom_point(aes(y=deltab1, x= V), fill ="red" , size=3, pch = 21)+
    geom_line(data = maxpost.linear.regress.b1, aes(x=xgrid, y =V1), color="blue", size=1)+
    #geom_hline(yintercept = 0, color="black", size=0.3, linetype ="dashed")+
    geom_point(data =data.frame(x=df.deltab.V$V, y=df.deltab.V$deltab1 + df.deltab.V$stdev_deltab1+0.1), 
               aes(x=x, y=y),  size = 10, color="black", pch=21, fill="white")+
    annotate("text", x=df.deltab.V$V,  y=df.deltab.V$deltab1 + df.deltab.V$stdev_deltab1 + 0.1, label= df.deltab.V$Period, color = "black", size=7) +
    scale_x_continuous(name = TeX("$V \\; \\left[m^{3} \\right] $"), expand = c(0,0), limits=c(xgrid[1], tail(xgrid,1))) +
    scale_y_continuous(name = TeX("$\\Delta b_1 \\; \\left[m\\right]$"), expand = c(0,0)) +
    theme_bw()+ coord_cartesian(clip = 'off')+
    theme(plot.background    = element_rect(fill ="transparent", color = NA)
          ,panel.grid.major  = element_blank()
          ,panel.grid.minor  = element_blank()
          ,panel.background  = element_rect(fill ="transparent") 
          ,axis.ticks        = element_line(colour = "black")
          ,plot.margin       = unit(c(0.3,0.7,0.3,0.3),"cm")
          ,text              = element_text(size=14)
          ,axis.title        = element_text(size=20)
          ,legend.key        = element_rect(colour = "transparent", fill = "transparent")
          ,legend.background = element_rect(colour = "transparent", fill = "transparent")
          ,legend.position   ="none")
  
  ggsave(lin.relat.deltab1.plot, filename=paste0(dir.qs.delta_b1,"/Linear_qs_Deltab1.png"), 
         bg = "transparent", width = 7, height =7, dpi = 400)
  ##############################################################################################################
  # delta b2:
  write.table(df.linear.regress.deltab2, file=paste0(dir.config.linear,"/data.txt"), sep="\t", row.names=FALSE )
  setwd(dir.exe)
  #launch BaM :
  linear_app(model.type         = model.type, 
             Nmcmc              = Nmcmc, 
             Ncycles            = Ncycles, 
             xgrid              = xgrid,
             nobs               = length(df.linear.regress.deltab2$X), 
             simMCMC            = TRUE, 
             prediction         = TRUE, 
             prediction.t       = FALSE, 
             nobs.t             = 1,
             prior.param.propor = prior.param.propor,
             prior.gamma.linear = prior.gamma.linear) 
  #read results :
  env.linear.regress.b2       = read.table(file = paste0(dir.config.linear,"/Lin_TotalU.env"), header=TRUE)
  env.linear.regress.b2       = cbind(env.linear.regress.b2, xgrid)
  env.param.linear.regress.b2 = read.table(file = paste0(dir.config.linear,"/Lin_ParamU.env"), header=TRUE)
  env.param.linear.regress.b2 = cbind(env.param.linear.regress.b2, xgrid)
  maxpost.linear.regress.b2   = read.table(file = paste0(dir.config.linear,"/Lin_Maxpost.spag"), header=FALSE)
  maxpost.linear.regress.b2   = cbind(maxpost.linear.regress.b2, xgrid)
  #write and save results of this linear regression Deltab1,b2 vs qscum :
  list.of.files <- c(
    paste0(dir.config.linear,"/Lin_maxpost.spag"),
    paste0(dir.config.linear,"/Lin_ParamU.spag"), 
    paste0(dir.config.linear,"/Lin_TotalU.spag"),
    paste0(dir.config.linear,"/Lin_ParamU.env"),  
    paste0(dir.config.linear,"/Lin_TotalU.env"),
    paste0(dir.config.linear,"/Results_MCMC_Cooked.txt"), 
    paste0(dir.config.linear,"/Results_Residuals.txt"),
    paste0(dir.config.linear,"/Results_Summary.txt"), 
    paste0(dir.config.linear,"/Config_Model.txt"),
    paste0(dir.config.linear,"/xgrid.txt"),
    paste0(dir.config.linear,"/data.txt")
  )
  for (ll in 1:length(list.of.files)) {
    file.copy(list.of.files[ll], dir.qs.delta_b2, overwrite = TRUE)
  }
  
  # plot results :
  lin.relat.deltab2.plot = ggplot(data=df.deltab.V) +
    geom_ribbon(data = env.linear.regress.b2, aes(x = xgrid, ymin =Q_q2.5, ymax = Q_q97.5), fill = "blue", alpha = 0.1)+
    geom_ribbon(data = env.param.linear.regress.b2, aes(x = xgrid, ymin =Q_q2.5, ymax = Q_q97.5), fill = "blue", alpha = 0.2)+
    geom_errorbar(aes(ymin = deltab2 - 2*stdev_deltab2, ymax = deltab2 + 2*stdev_deltab2, x = V), color = "blue",  size = 0.3, width=2000)+
    geom_point(aes(y=deltab2,   x = V),  fill ="blue" , size=3, pch = 21)+
    geom_line(data = maxpost.linear.regress.b2, aes(x=xgrid, y =V1), color="blue", size=1)+
    geom_point(data =data.frame(x=df.deltab.V$V,  y=df.deltab.V$deltab2 + df.deltab.V$stdev_deltab2+0.1),
               aes(x=x, y=y),  size = 10, color="black", pch=21, fill="white")+
    annotate("text", x=df.deltab.V$V, y=df.deltab.V$deltab2 + df.deltab.V$stdev_deltab2 + 0.1,
                    label= df.deltab.V$Period, color = "black", size=7) +
    scale_x_continuous(name = TeX("$V \\; \\left[m^{3} \\right] $"), expand = c(0,0), limits=c(xgrid[1], tail(xgrid,1))) +
    scale_y_continuous(name = TeX("$\\Delta b_2 \\; \\left[m\\right]$"), expand = c(0,0)) +
    theme_bw()+ coord_cartesian(clip = 'off')+
    theme(plot.background    = element_rect(fill ="transparent", color = NA)
          ,panel.grid.major  = element_blank()
          ,panel.grid.minor  = element_blank()
          ,panel.background  = element_rect(fill ="transparent") 
          ,axis.ticks        = element_line(colour = "black")
          ,plot.margin       = unit(c(0.3,0.7,0.3,0.3),"cm")
          ,text              = element_text(size=14)
          ,axis.title        = element_text(size=20)
          ,legend.key        = element_rect(colour = "transparent", fill = "transparent")
          ,legend.background = element_rect(colour = "transparent", fill = "transparent")
          ,legend.position   ="none")
  ggsave(lin.relat.deltab2.plot, filename=paste0(dir.qs.delta_b2,"/Linear_qs_Deltab2.png"), bg = "transparent",
         width = 7, height =7, dpi = 400)
  
  
  
  #both figures:
  both.plot = plot_grid( lin.relat.deltab1.plot, lin.relat.deltab2.plot,
                         ncol = 2, nrow = 1, rel_widths = c(1,1), labels=c("a)", "b)"), label_size = 25)
  ggsave(both.plot, filename=paste0(dir.sed.transp.linear,"/Linear_V_Deltab.png"), bg = "transparent",
         width = 14, height =8, dpi = 400)
  
  
  print("****************")
  print("   All done!    ")
  print("****************")           
  
  return(list(df.linear.regress.deltab1 = df.linear.regress.deltab1, 
              df.linear.regress.deltab2 = df.linear.regress.deltab2))
}






















#########################################################################################
linear_app <- function(  model.type ,
                         Nmcmc, Ncycles,  
                         xgrid, nobs, simMCMC, 
                         prediction , 
                         prediction.t , 
                         nobs.t, 
                         prior.param.propor,
                         prior.gamma.linear) {               # Linear Model regression
  #########################################################################################
  BaM.linear.config(model.type , Nmcmc, Ncycles, xgrid, nobs, simMCMC, 
                    prediction, prediction.t, nobs.t, 
                    prior.param.propor, prior.gamma.linear)
  message("Applying Linear regression - BaM !!!  Wait ... "); flush.console()
  #system2(paste(dir_code,"/BaM_exe/BaM_Segmentation.exe",sep=""),stdout =NULL, stderr = NULL);
  system2(paste0(dir_code,"/BaM_exe/BaM_2exp_pool2.exe")) #, stdout =NULL, stderr = NULL); 
}
















# Linear Model regression
###############################################################################
BaM.linear.config <- function(model.type , 
                              Nmcmc,Ncycles, xgrid, nobs, 
                              simMCMC, prediction, prediction.t, nobs.t,
                              prior.param.propor, prior.gamma.linear) {  
###############################################################################
  write.table(xgrid, file =paste0(dir_code,"/BaM_exe/Linear/xgrid.txt"), 
              col.names = FALSE, row.names = FALSE)
  ngrid = length(xgrid)
  theta1 <- 1  ; St_theta1 <- 100
  
  
  #---------------------------------------------------------------------------
  #creation of Config_BaM.txt
  file.bam = paste0(dir_code,"/BaM_exe/Config_BaM.txt")
  cat('"Linear/"',                 file = file.bam , sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"',   file = file.bam , sep="\n", append = TRUE)    
  cat('"Config_Model.txt"',        file = file.bam , sep="\n", append = TRUE)
  cat('""',                        file = file.bam , sep="\n", append = TRUE)
  cat('"Config_Data.txt"',         file = file.bam , sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file.bam , sep="\n", append = TRUE)                                       
  cat('"Config_MCMC.txt"',         file = file.bam , sep="\n", append = TRUE)                                            
  cat('"Config_Cooking.txt"',      file = file.bam , sep="\n", append = TRUE)
  cat('"Config_Summary.txt"',      file = file.bam , sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"',    file = file.bam , sep="\n", append = TRUE)
  if( (prediction ==TRUE) | (prediction.t ==TRUE)) {
    cat('"Config_Pred_Master.txt"',file = file.bam , sep="\n", append = TRUE)
  } else {
    cat('""',                      file = file.bam , sep="\n", append = TRUE)
  }
  #-------------------------------------------------------------------------
  file.name = paste0(dir_code,"/BaM_exe/Linear/Config_Model.txt")
  cat('"Linear"',                  file = file.name,sep="\n")
  cat(1,                           file = file.name, append = TRUE, sep="\n")
  cat(1,                           file = file.name, append = TRUE, sep="\n")
  cat(2,                           file = file.name, append = TRUE, sep="\n") # n parameters
  
  if (model.type == "null"){
  ###########################
    cat('"a1"',                    file = file.name, append = TRUE, sep="\n")
    cat(prior.param.propor[2],     file = file.name, append = TRUE, sep="\n")
    cat(prior.param.propor[3],     file = file.name, append = TRUE, sep="\n")
    cat(prior.param.propor[4],     file = file.name, append = TRUE, sep=",")
    cat(",",                       file = file.name, append = TRUE, sep=",")
    cat(prior.param.propor[5],     file = file.name, append = TRUE, sep="\n")
    
    cat('"a2"',                    file = file.name, append = TRUE, sep="\n")
    cat(prior.param.propor[7],     file = file.name, append = TRUE, sep="\n")
    cat(prior.param.propor[8],     file = file.name, append = TRUE, sep="\n")
    cat(prior.param.propor[9],     file = file.name, append = TRUE, sep=",")
    cat(",",                       file = file.name, append = TRUE, sep=",")
    cat(prior.param.propor[10],    file = file.name, append = TRUE, sep="\n")
    
    
  } else if (model.type == "proportional"){
    #########################################
    cat('"a1"',                    file = file.name, append = TRUE,sep="\n")
    cat(prior.param.propor[2],     file = file.name, append = TRUE,sep="\n")
    cat(prior.param.propor[3],     file = file.name, append = TRUE,sep="\n")
    cat(prior.param.propor[4],     file = file.name, append = TRUE,sep=",")
    cat(",",                       file = file.name, append = TRUE,sep=",")
    cat(prior.param.propor[5],     file = file.name, append = TRUE,sep="\n")
    
    cat('"a2"',                    file = file.name, append = TRUE,sep="\n")
    cat(prior.param.propor[7],     file = file.name, append = TRUE,sep="\n")
    cat(prior.param.propor[8],     file = file.name, append = TRUE,sep="\n")
    cat(prior.param.propor[9],     file = file.name, append = TRUE,sep=",")
    cat(",",                       file = file.name, append = TRUE,sep=",")
    cat(prior.param.propor[10],    file = file.name, append = TRUE,sep="\n")
  }
  #-------------------------------------------------------------------------  
  file.name2 = paste0(dir_code,"/BaM_exe/Linear/Config_Data.txt")
  name = "'Linear/data.txt'"
  cat(name,                        file = file.name2, sep="\n") 	             # path to data file
  cat(1,                           file = file.name2, sep="\n",append = TRUE)  # number of header lines
  cat(nobs,                        file = file.name2, sep="\n",append = TRUE)  # Nobs, number of rows in data file (excluding header lines)
  cat(4,                           file = file.name2, sep="\n",append = TRUE)  # number of columns in the data file
  cat(1,                           file = file.name2, sep="\n",append = TRUE)  # columns for X (observed inputs) in data file - comma-separated if several (order: t, h, T, T_smooth)
  cat(0,                           file = file.name2, sep="\n",append = TRUE)  # 8,9,10,11,12 !!! columns for Xu (random uncertainty in X, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0,                           file = file.name2, sep="\n",append = TRUE)  # columns for Xb (systematic uncertainty in X, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0,                           file = file.name2, sep="\n",append = TRUE)  # columns for Xb_indx (index of systematic errors in X - use 0 for a no-error assumption)
  cat(3,                           file = file.name2, sep="\n",append = TRUE)  # 16,18,20,21 columns for Y (observed outputs) in data file - comma-separated if several (order: Q, g0, cos(pi*g0), sin(pi*g0))
  cat(4,                           file = file.name2, sep="\n",append = TRUE)  # 17,19,22,23 columns for Yu (uncertainty in Y, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0,                           file = file.name2, sep="\n",append = TRUE)  # columns for Yb (systematic uncertainty in Y, EXPRESSED AS A STANDARD DEVIATION - use 0 for a no-error assumption)
  cat(0,                           file = file.name2, sep="",  append = TRUE)  # columns for Yb_indx (index of systematic errors in Y - use 0 for a no-error assumption)
  #------------------------------------------------------------------------
  file.mcmc = paste0(dir_code,"/BaM_exe/Linear/Config_MCMC.txt")
  cat('"Results_MCMC.txt"',        file = file.mcmc, sep="\n")
  cat(Nmcmc,                       file = file.mcmc, append = TRUE,sep="\n")   #Nadapt
  cat(Ncycles,                     file = file.mcmc, append = TRUE,sep="\n")  #Ncycles
  cat(0.1,                         file = file.mcmc, append = TRUE,sep="\n")    #minMoveRate
  cat(0.5,                         file = file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
  cat(0.9,                         file = file.mcmc, append = TRUE,sep="\n")    #DownMult
  cat(1.1,                         file = file.mcmc, append = TRUE,sep="\n")    #UpMult
  cat(0,                           file = file.mcmc, append = TRUE,sep="\n")      #mode for init jump distr
  cat("****",                      file = file.mcmc, append = TRUE,sep="\n") 
  cat(0.1,                         file = file.mcmc, append = TRUE,sep="\n")    #MultFact
  cat(0.1,                         file = file.mcmc, append = TRUE, sep=",")     #RC MultiFact
  cat(0.1,                         file = file.mcmc, append = TRUE, sep=",")     
  cat(0.1,                         file = file.mcmc, append = TRUE,sep="\n")
  cat(0.1,                         file = file.mcmc, append = TRUE, sep=",")      #Remnant MultiFact
  cat(0.1,                         file = file.mcmc, append = TRUE,sep="\n")

  
  
  # PREDICTIONS:
  ##############
  if(prediction ==TRUE) {  
    file.Pred1 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_Master.txt")
    if(prediction.t ==TRUE) {      
      cat(6,                         file = file.Pred1,sep="\n")
    } else {
      cat(3,                         file = file.Pred1,sep="\n")
    }  
    cat("'Config_Pred_Maxpost.txt'", file = file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_ParamU.txt'",  file = file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_TotalU.txt'",  file = file.Pred1, append = TRUE,sep="\n")
    
    #------------------------------------------------------------------------ 
    # MAXPOST:
    file.Pred3 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_Maxpost.txt")
    cat("'Linear\\xgrid.txt'",       file = file.Pred3, sep="\n")
    cat(ngrid,                       file = file.Pred3, sep="\n", append = TRUE)
    cat("1",                         file = file.Pred3, append = TRUE,sep="\n")   # n of spaghetti
    cat(".false.",                   file = file.Pred3, append = TRUE,sep="\n")   # Propagate parametric uncertainty?
    cat(".false.",                   file = file.Pred3, append = TRUE,sep="\n")   # Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1",                        file = file.Pred3, append = TRUE,sep="\n")   # Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Lin_Maxpost.spag'",        file = file.Pred3, append = TRUE,sep="\n")   # Files containing spaghettis for each output variable (size nY)
    cat(".true.",                    file = file.Pred3, append = TRUE,sep="\n")   # Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.",                    file = file.Pred3, append = TRUE,sep="\n")   # Post-processing: create envelops? (size nY)
    cat("'Lin_Maxpost.env'",         file = file.Pred3, append = TRUE,sep="\n")   # Post-processing: name of envelop files (size nY)
    cat(".true.",                    file = file.Pred3, append = TRUE,sep="\n")   # Print progress in console during computations?
    cat(".false." ,                  file = file.Pred3, append = TRUE,sep="\n")   # Do state prediction? (size nState)
    #-------------------------------------------------------------------------
    # PARAMETRIC UNCERTAINTY:
    file.Pred4 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_ParamU.txt")
    cat("'Linear\\xgrid.txt'",       file = file.Pred4, sep="\n")
    cat(ngrid,                       file = file.Pred4, sep="\n", append = TRUE)
    cat("1",                         file = file.Pred4, append = TRUE,sep="\n")   # n of spaghetti
    cat(".true.",                    file = file.Pred4, append = TRUE,sep="\n")   # Propagate parametric uncertainty?
    cat(".false.",                   file = file.Pred4, append = TRUE,sep="\n")   # Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1",                        file = file.Pred4, append = TRUE,sep="\n")   # Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Lin_ParamU.spag'",         file = file.Pred4, append = TRUE,sep="\n")   # Files containing spaghettis for each output variable (size nY)
    cat(".true.",                    file = file.Pred4, append = TRUE,sep="\n")   # Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.",                    file = file.Pred4, append = TRUE,sep="\n")   # Post-processing: create envelops? (size nY)
    cat("'Lin_ParamU.env'",          file = file.Pred4, append = TRUE,sep="\n")   # Post-processing: name of envelop files (size nY)
    cat(".true.",                    file = file.Pred4, append = TRUE,sep="\n")   # Print progress in console during computations?
    cat(".false." ,                  file = file.Pred4, append = TRUE,sep="\n")   # Do state prediction? (size nState)
    #------------------------------------------------------------------------- 
    # TOTAL UNCERTAINTY
    file.Pred5 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_TotalU.txt")
    cat("'Linear\\xgrid.txt'",       file = file.Pred5, sep ="\n")
    cat(ngrid,                       file = file.Pred5, append = TRUE, sep="\n")
    cat("1",                         file = file.Pred5, append = TRUE, sep="\n")   # n of spaghetti
    cat(".true.",                    file = file.Pred5, append = TRUE, sep="\n")   # Propagate parametric uncertainty?
    cat(".true.",                    file = file.Pred5, append = TRUE, sep="\n")   # Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1",                        file = file.Pred5, append = TRUE, sep="\n")   # Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Lin_TotalU.spag'",         file = file.Pred5, append = TRUE, sep="\n")   # Files containing spaghettis for each output variable (size nY)
    cat(".true.",                    file = file.Pred5, append = TRUE, sep="\n")   # Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.",                    file = file.Pred5, append = TRUE, sep="\n")   # Post-processing: create envelops? (size nY)
    cat("'Lin_TotalU.env'",          file = file.Pred5, append = TRUE, sep="\n")   # Post-processing: name of envelop files (size nY)
    cat(".true.",                    file = file.Pred5, append = TRUE, sep="\n")   # Print progress in console during computations?
    cat(".false." ,                  file = file.Pred5, append = TRUE, sep="\n")   # Do state prediction? (size nState)
  }
  #---------------------------------------------------------------------- 
  if(prediction.t==TRUE) {
    file.Pred1 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_Master.txt")
    if(prediction==FALSE) {      
      cat(3,                            file =file.Pred1,sep="\n")
    }  
    cat("'Config_Pred_qs_Maxpost.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_qs_ParamU.txt'",  file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_qs_TotalU.txt'",  file =file.Pred1, append = TRUE,sep="\n")
    ###################################################################
    file.Pred6 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_qs_Maxpost.txt")
    cat("'Linear\\qs_t.txt'",            file = file.Pred6, sep="\n")
    cat(nobs.t,                          file = file.Pred6, append = TRUE, sep="\n")
    cat("1",                             file = file.Pred6, append = TRUE,sep="\n")   #n of spaghetti
    cat(".false.",                       file = file.Pred6, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.",                       file = file.Pred6, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1",                            file = file.Pred6, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'qst_maxpost.spag'",            file = file.Pred6, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.",                        file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.",                        file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY) 
    cat("'qst_Maxpost.env'",             file = file.Pred6, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".true.",                        file = file.Pred6, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." ,                      file = file.Pred6, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    ###################################################################
    file.Pred8 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_qs_ParamU.txt")
    cat("'Linear\\qs_t.txt'",            file = file.Pred8, sep="\n")
    cat(nobs.t,                          file = file.Pred8, append = TRUE, sep="\n")
    cat("1",                             file = file.Pred8, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.",                        file = file.Pred8, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.",                       file = file.Pred8, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1",                            file = file.Pred8, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'qst_ParamU.spag'",             file = file.Pred8, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.",                        file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.",                        file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'qst_ParamU.env'",              file = file.Pred8, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".true.",                        file = file.Pred8, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." ,                      file = file.Pred8, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    ###################################################################
    file.Pred7 = paste0(dir_code,"/BaM_exe/Linear/Config_Pred_qs_TotalU.txt")
    cat("'Linear\\qs_t.txt'",            file = file.Pred7, sep="\n")
    cat(nobs.t,                          file = file.Pred7, append = TRUE, sep="\n")
    cat("1",                             file = file.Pred7, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.",                        file = file.Pred7, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".true.",                        file = file.Pred7, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1",                            file = file.Pred7, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'qst_TotalU.spag'",             file = file.Pred7, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'qst_TotalU.env'", file = file.Pred7, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred7, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  }
  #-------------------------------------------------------------------------- 
  file.remnant = paste0(dir_code,"/BaM_exe/Linear/Config_RemnantSigma.txt")
  # cat("'Constant'", file = file.remnant, sep="\n")     #! Function f used in sdev=f(Qrc) 
  # cat(1, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
  # cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
  # cat(0.1, file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
  # cat("'Uniform'", file = file.remnant, append = TRUE, sep="\n")            #! Prior distribution
  # cat(0,file =file.remnant, append = TRUE, sep=",")
  # cat(",",file =file.remnant, append = TRUE, sep=",")
  # cat(100,file =file.remnant, append = TRUE, sep="\n")
  
  cat(prior.gamma.linear[1], file = file.remnant, sep="\n")     #! Function f used in sdev=f(Qrc) 
  if (prior.gamma.linear[1] == "'Linear'"){ 
    cat(2,                   file = file.remnant, append = TRUE, sep="\n")         #! Number of parameters gamma
  } else if (prior.gamma.linear[1] == "'Constant'"){
    cat(1,                   file = file.remnant, append = TRUE, sep="\n") 
  }
  cat(prior.gamma.linear[2], file = file.remnant, append = TRUE, sep="\n")         #! Parameter Name
  cat(prior.gamma.linear[3], file = file.remnant, append = TRUE, sep="\n")         #! Initial Guess
  cat(prior.gamma.linear[4], file = file.remnant, append = TRUE, sep="\n")         #! Prior distribution
  cat(prior.gamma.linear[5], file = file.remnant, append = TRUE, sep=",")
  cat(","                  , file = file.remnant, append = TRUE, sep=",")
  cat(prior.gamma.linear[6], file = file.remnant, append = TRUE, sep="\n")
  
  if (prior.gamma.linear[1] == "'Linear'"){ 
    cat(prior.gamma.linear[7],  file = file.remnant, append = TRUE, sep="\n")      #! Parameter Name
    cat(prior.gamma.linear[8],  file = file.remnant, append = TRUE, sep="\n")      #! Initial Guess
    cat(prior.gamma.linear[9],  file = file.remnant, append = TRUE, sep="\n")      #! Prior distribution
    cat(prior.gamma.linear[10], file = file.remnant, append = TRUE, sep=",")
    cat(","                  ,  file = file.remnant, append = TRUE, sep=",")
    cat(prior.gamma.linear[11], file = file.remnant, append = TRUE, sep="\n")
  }
  
  ###################################################################   RUN OPTIONS
  file.run = paste0(dir_code,"/BaM_exe/Linear/Config_RunOptions.txt")
  if (simMCMC == TRUE) {  
    cat(".true.", file =file.run,sep="\n")                     #Do MCMC?
  } else {
    cat(".false.", file =file.run,sep="\n")                    #Do MCMC?    
  }
  cat(".true.", file =file.run, append = TRUE, sep="\n")       #Do MCMC summary?
  cat(".true.", file =file.run, append = TRUE, sep="\n")       #Do Residual diagnostics?
  if ((prediction == TRUE) | (prediction.t == TRUE)) {
    cat(".true.", file =file.run, append = TRUE, sep="\n")     #Do Predictions?
    cat(".false.", file =file.run, append = TRUE, sep="\n")    #Do Predictions?
  }
  
  ###################################################################   COOKING CONFIG
  file.cooking = paste(dir_code,"/BaM_exe/Linear/Config_Cooking.txt",sep="")
  cat("'Results_MCMC_Cooked.txt'" , file =file.cooking ,sep="\n")    #Result file
  cat(0.5, file =file.cooking, append = TRUE, sep="\n")            #Burn factor
  cat(10, file =file.cooking, append = TRUE, sep="\n")             #Nslim
  ###################################################################   RESIDUALS CONFIG
  file.residuals = paste(dir_code,"/BaM_exe/Linear/Config_Residuals.txt",sep="")
  cat("'Results_Residuals.txt'" , file =file.residuals ,sep="\n")    #Result file
  ###################################################################   SUMMARY CONFIG
  file.summary = paste(dir_code,"/BaM_exe/Linear/Config_Summary.txt",sep="")
  cat("'Results_Summary.txt'" , file =file.summary ,sep="\n")    #Result file
}




































# NO MORE USED #######
#####################################################################################################
# Analyze results of the sediment transport detection:
ST.results.analysis =  function(dir.retro.an, dir.SPD.exe, dir.sed.transp, case_study_name,
                                results.ST,
                                data.annotate.step2.sort,
                                df.limni,  t_Gaug, h_Gaug, Q_Gaug, uQ_Gaug,
                                officialShiftsTime, 
                                a.prior, st_a.prior,   
                                c.prior, st_c.prior,  
                                b.prior.STEP2    = b.prior.STEP2, 
                                st_b.prior.STEP2 = st_b.prior.STEP2,
                                Bw.prior, Cr.prior, g.prior, Bc.prior, KS.prior, S0.prior,
                                st_Bw.prior, st_Cr.prior, st_g.prior, st_Bc.prior, st_KS.prior, st_S0.prior,
                                ncontrol,
                                M, 
                                dg.prior, st_dg.prior, dl.prior, st_dl.prior,   d.Bc1.prior, st_dBc1.prior, #width changes priors
                                g1.prior, g2.prior, g1.distr.type, g2.distr.type, remnant = "Linear" ,#remnant error model parameters
                                isVar, 
                                margins.bac=c('Gaussian','LogNormal','Gaussian'), 
                                pred = FALSE,  
                                nsim, 
                                Ncycles = 500, 
                                changes.method = "cumsum", 
                                ncolumns=8, rowmax=20, 
                                FinalColors =colo,
                                global.change = TRUE,  local.change = TRUE,  width.change = FALSE, 
                                n.parvar = 2,
                                grid_RC.xlim, grid_RC.xstep, grid_RC.ylim, grid_RC.ystep, ticks_RC.y.log ) {

#####################################################################################################      
# In this section we analyse the results of the sediment transport proxy model. 
# Define the new set of shift times obtained with the sediment transport proxy model.
# Define the gaugings dataset per periods with the new set of shift times.
# be careful!!! some periods may contain no gaugings, thus some NA values (-999) 
# are added to the gaugings
# dataset for those periods.
  
# Then Apply the SPD estimation of the RCs for each period.
# Dircetories:
                  dir.create(paste0(dir.retro.an,"/STEP1_SPD"))
                  dir.SPD.exe <- paste0(dir_code,"/BaM_exe/BaRatin_SPD")
                  dir.SPD.Sed.Transp.results <- paste0(dir.retro.an,"/STEP1_SPD")
                  
  #Application of BaRatin-SPD:
                  # write.table(Gaug.sed.transp[2:5], paste0(dir.SPD.Sed.Transp.results,"/Gaugings_data_SPD.txt"),
                  #             sep="\t", row.names=FALSE , col.names = c("h","Q", "uQ", "Period"))  
                  #gaugings = rbind(gaugings,c(time=-999, h=-999, Q=-999, uQ=-999, Period=10))
                  gaug.step1 = read.table(paste0(dir.retro.an,"/STEP1_data_with_periods.txt"),header=TRUE)
                  nperiod.step1.gaug = tail(gaug.step1$Period,1)
                  
                  
                  
                  write.table(gaug.step2[1:4], paste0(dir.SPD.exe,"/Gaugings_data_SPD.txt"),
                              sep="\t", row.names=FALSE , col.names = c("h","Q", "uQ", "Period"))
                  write.table(gaug.step2[1:4], paste0(dir.retro.an,"/STEP2_SPD/Gaugings_data_SPD.txt"),
                              sep="\t", row.names=FALSE , col.names = c("h","Q", "uQ", "Period"))
                  
                  BaRatin_SPD.bac_app(dir_code           =  dir_code,
                                      dir.BaM            =  dir.exe, 
                                      dir.SPD.config     =  dir.SPD.exe, 
                                      dir.SPD.results    =  dir.SPD.Sed.Transp.results, 
                                      nperiod            =  nperiod.step2.gaug,
                                      a.prior            =  a.prior, 
                                      st_a.prior         =  st_a.prior,
                                      c.prior            =  c.prior, 
                                      st_c.prior         =  st_c.prior, 
                                      b.prior            =  b.prior.STEP2, 
                                      st_b.prior         =  st_b.prior.STEP2,
                                      Bw.prior, 
                                      Cr.prior, 
                                      g.prior, 
                                      Bc.prior, 
                                      KS.prior, 
                                      S0.prior,
                                      st_Bw.prior,
                                      st_Cr.prior, 
                                      st_g.prior, 
                                      st_Bc.prior,
                                      st_KS.prior,
                                      st_S0.prior,
                                      ncontrol, 
                                      M, #hydraulic configuration matrix and controls
                                      dg.prior, st_dg.prior, dl.prior, st_dl.prior,   #global and local changes priors
                                      d.Bc1.prior, st_dBc1.prior,   # width changes priors
                                      g1.prior, g2.prior, 
                                      g1.distr.type, g2.distr.type, 
                                      remnant           =  remnant ,#remnant error model parameters
                                      isVar             =  isVar,   # determined the parameter that are varying over time
                                      margins.bac       =  c('Gaussian','LogNormal','Gaussian'), 
                                      pred              =  pred,    # prediction experiment TRUE or FALSE
                                      nsim              =  nsim,    # number of MC samples for prior progation
                                      Ncycles           =  Ncycles,  # number of MCMC cycles for Metropolis-Hastings
                                      changes.method    =  changes.method,
                                      ncolumns          =  ncolumns, 
                                      rowmax            =  rowmax, 
                                      FinalColors       =  colo,
                                      global.change     =  global.change,
                                      local.change      =  local.change,
                                      width.change      =  width.change,
                                      n.parvar          =  n.parvar,
                                      case_study_name   =  case_study_name)
  #plotting results of BaRatin-SPD in terms of rating curves :
                  plot.SPD(dir.SPD.exe         =  dir.SPD.exe, 
                           dir.SPD.results     =  dir.SPD.Sed.Transp.results, 
                           dir.SPD.config      =  dir.SPD.exe,
                           nperiod             =  nperiod.step2.gaug, 
                           df.limni            =  df.limni,
                           FinalColors         =  colo,
                           ylim.wind           =  grid_RC.ylim, 
                           xlim.wind           =  grid_RC.xlim,
                           breaks.lin.x        =  seq(grid_RC.xlim[1],grid_RC.xlim[2],grid_RC.xstep), 
                           labels.lin.x        =  seq(grid_RC.xlim[1],grid_RC.xlim[2],grid_RC.xstep), 
                           breaks.lin.y        =  seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep), 
                           labels.lin.y        =  seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep),
                           ylim.log.wind       =  c(0.1, grid_RC.ylim.log[2]),  
                           breaks.log          =  ticks_RC.y.log,   
                           labels.log          =  ticks_RC.y.log,
                           case_study_name     =  case_study_name) 

  #compute the delta b2:
                  #boxplots of b1 and b2:
                  SPD.summary = read.table(file=paste0(dir.SPD.Sed.Transp.results, "/Results_Summary.txt"))
                  SPD.mcmc.cooked = read.table(file=paste0(dir.SPD.Sed.Transp.results, "/Results_MCMC_Cooked.txt"), header=TRUE)
                  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked[,1:nperiod.step1.gaug]; names(SPD.mcmc.cooked.b1) = seq(1, nperiod.step1.gaug,1)
                  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked.b1 %>% gather(period, b1, na.rm = FALSE, convert = FALSE)
                  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked[,(nperiod.step1.gaug + 3):(2*nperiod.step1.gaug+3)]; 
                  names(SPD.mcmc.cooked.b2) = seq(1,nperiod.step1.gaug,1)
                  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked.b2 %>% gather(period, b2, na.rm = FALSE, convert = FALSE)
                  #
                  SPD.bt = ggplot()+
                           geom_boxplot(data=SPD.mcmc.cooked.b1, aes(x= period, y= b1, color="b1"), outlier.size = 0.1, lwd =0.2)+
                           geom_boxplot(data=SPD.mcmc.cooked.b2, aes(x= period, y= b2, color="b2"), outlier.size = 0.1, lwd =0.2)+
                           stat_summary(fun.y = mean, geom="point", shape=16, size=2)+
                           scale_y_continuous(name = TeX("$b_1, b_2$"))+
                           scale_color_manual(name = "Legend", labels=c("b1", "b2") , values = c("blue", "red"))+
                           scale_x_discrete(name = "Period", limits =seq(1, nperiod.step2.gaug,1))+
                           theme_bw()
                           ggsave(SPD.bt, filename = paste0(dir.SPD.Sed.Transp.results,"/bt_boxplots.png"), 
                                  bg = "transparent", device = "png", width = 6, height =4, dpi = 400)
                  #      
                  summary.SPD <- read.table(paste0(dir.SPD.Sed.Transp.results,"/Results_Summary.txt"), header= TRUE) 
                  # this part is specific to Meyras hydraulic configuration:
                  nperiod.step2.gaug = nperiod.step1.gaug
                  
                  bt1.SPD = rbind(summary.SPD[16, 1:nperiod.step2.gaug],
                                  summary.SPD[7, 1:nperiod.step2.gaug], 
                                  summary.SPD[10, 1:nperiod.step2.gaug],
                                  summary.SPD[11, 1:nperiod.step2.gaug],
                                  summary.SPD[5, 1:nperiod.step2.gaug])
                  bt2.SPD=  rbind(summary.SPD[16, (nperiod.step2.gaug+3):(2*nperiod.step2.gaug+2)], 
                                  summary.SPD[7,(nperiod.step2.gaug+3):(2*nperiod.step2.gaug+2)], 
                                  summary.SPD[10, (nperiod.step2.gaug+3):(2*nperiod.step2.gaug+2)],
                                  summary.SPD[11, (nperiod.step2.gaug+3):(2*nperiod.step2.gaug+2)],
                                  summary.SPD[5, (nperiod.step2.gaug+3):(2*nperiod.step2.gaug+2)]
                  )
                  delta.b1.MAP =NULL; delta.b1.stdev = NULL; delta.b1.mean =NULL;
                  for (evento in 1:(nperiod.step2.gaug-1)) {
                    delta.b1.MAP[evento] = (bt1.SPD[[evento+1]][1] - bt1.SPD[[evento]][1])
                    delta.b1.mean[evento] = (bt1.SPD[[evento+1]][5] - bt1.SPD[[evento]][5])
                    delta.b1.stdev[evento] = (bt1.SPD[[evento+1]][4]^2 + bt1.SPD[[evento]][4]^2)^0.5
                  }
                  delta.b2.MAP =NULL; delta.b2.stdev = NULL; delta.b2.mean =NULL;
                  for (evento in 1:(nperiod.step2.gaug-1)) {
                    delta.b2.MAP[evento] = (bt2.SPD[[evento+1]][1] - bt2.SPD[[evento]][1])
                    delta.b2.mean[evento] = (bt2.SPD[[evento+1]][5] - bt2.SPD[[evento]][5]) 
                    delta.b2.stdev[evento] = (bt2.SPD[[evento+1]][4]^2 + bt2.SPD[[evento]][4]^2)^0.5
                  }
                  
 # # if you want to cumulate the events because no enough gaugings:
 #                  realqs.cum = 0
 #                  jkl = 0
 #                  for ( real in 1:length(finalqs.cum)) {
 # 
 #                    if(any(index.ts.without.gaugings == real)) {
 #                      realqs.cum[jkl] = realqs.cum[jkl] + finalqs.cum[real]
 #                    } else {
 #                      jkl=jkl+1
 #                      realqs.cum[jkl] = finalqs.cum[real]
 #                    }
 #                  }
                  df.rel.deltab.qscum.TOT = data.frame(deltab1 = delta.b1.mean,
                                                       stdev_deltab1 = delta.b1.stdev,
                                                       deltab2 = delta.b2.mean,
                                                       stdev_deltab2 = delta.b2.stdev,
                                                       Period = seq(1,length(delta.b1.mean),1))
                  
  #final plot:
                  # finding which estimated deltab has a corresponding qsc ??
                  #*********************************************************
                  df.rel.deltab.qscum = df.rel.deltab.qscum.TOT[which(data.annotate.step2.sort$index>0,),]
                  df.rel.deltab.qscum = cbind(data.annotate.ST, df.rel.deltab.qscum)
                  #plot:
                  db1_db2_qsc.plot = ggplot(data=df.rel.deltab.qscum) +
                                     geom_point(aes(y=deltab1, x= qsc.mean), color ="red" , size=2)+
                                     geom_errorbar(aes(ymin = deltab1 - 2*stdev_deltab1, 
                                                       ymax = deltab1 + 2*stdev_deltab1,
                                                       x = qsc.mean),
                                                       color = "red",  size = 0.2, width=20)+
                                     geom_errorbarh(aes(y = deltab1, 
                                                        xmin = qsc.mean- 2*qsc.stdev,
                                                        xmax = qsc.mean+ 2*qsc.stdev), 
                                                       color = "red",  size = 0.2, width=20)+
                                     geom_point(aes(y=deltab2, x = qsc.mean), 
                                                    color ="blue" , size=2)+
                                     geom_errorbar(aes(ymin = deltab2 - 2*stdev_deltab2, 
                                                       ymax = deltab2 + 2*stdev_deltab2,
                                                       x = qsc.mean),
                                                       color = "blue",  size = 0.2, width=20)+
                                     geom_errorbarh(aes(y = deltab2, 
                                                        xmin = qsc.mean- 2*qsc.stdev,
                                                        xmax = qsc.mean+ 2*qsc.stdev), 
                                                        color = "blue",  size = 0.2, width=20)+
                                     annotate("text", x=df.rel.deltab.qscum$qsc.mean, y=0.5, 
                                              label= df.rel.deltab.qscum$Period, color = "green", size=5)+
                                     ylab("Delta b1, Delta b2 [m]")+ xlab("Cumulative bedload, qsc[m3]")+
                                     theme_bw(base_size = 15)
                                     ggsave(db1_db2_qsc.plot, filename = paste0(dir.SPD.Sed.Transp.results,"/db1_db2_qsc.png"), 
                                            bg = "transparent", device = "png", width = 6, height =4, dpi = 400)
                  # write.table(df.rel.deltab.qscum, paste0(dir.SPD.Sed.Transp.results,"/df.rel.deltab.qscum.txt"),
                  #             sep="\t", row.names=FALSE , col.names = c("q2", "MAP", "q97", "t.adj", "index", "qsc", 
                  #                                                       "deltab1","stdev_deltab1", 
                  #                                                       "deltab2","stdev_deltab2",
                  #                                                       "Period"))  
                  return(list(df.rel.deltab.qscum))
}
  
