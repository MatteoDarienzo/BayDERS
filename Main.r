###############################################################################################
#                                       BayDERS                                               # 
###############################################################################################
#                     BAYesian Detection and Estimation of Rating Shifts                      #
###############################################################################################
#                                                                                             #   
# Author:        Matteo Darienzo                                                              #
# Institute:     INRAE Lyon-Grenoble, UR RiverLy, France                                      #
# Year:          2019                                                                         #
# version:       v.1.0.0  - last update: 09/02/2022                                           #
#                                                                                             #
#                                                                                             #
#---------------------------------------------------------------------------------------------#
#                             Main objective of the program                                   #
#---------------------------------------------------------------------------------------------#
# This program performs the automatic detection of rating shifts in retrospective based on the#
# segmentation of gaugings and on a stage-recession analysis accounting for both uncertainties#
#---------------------------------------------------------------------------------------------#
#                                   Acknowledgements                                          #
#---------------------------------------------------------------------------------------------#
# This program has been written under the supervision of:                                     #
#                            Benjamin Renard (INRAE)                                          #
#                            Jerome Le Coz (INRAE)                                            #
#                            Michel Lang (INRAE)                                              # 
#                                                                                             #
# The code BayDERS makes use of:                                                              #
#                            "DMSL" package of prof. Dmitri Kavetski                          #
#                            "BaM.exe" software developped by Benjamin Renard, Inrae          #
#                            BaRatin method (Le Coz et al., 2014)                             #
#                            BaRatin-SPD method (Mansanarez et al., 2019)                     #
#                                                                                             #
# Project funded by:                                                                          #
#                            Electricite de France (EDF)                                      #
#                            Compagnie Nationale du Rhone (CNR)                               #
#                            SCHAPI (France)                                                  #
#                            INRAE (France)                                                   #
#                                                                                             #
#                                                                                             #
#---------------------------------------------------------------------------------------------#
#                             This is the "Main" program file!                                #
#---------------------------------------------------------------------------------------------#
# Read the "readme.txt" for user info.                                                        #
# Main references:                                                                            #
# - Mansanarez et al., 2019 (multi-period RC estimation)                                      #
# - Mansanarez et al., 2016 (PhD thesis, non-unique stage-discharge relations)                #
# - Le Coz et al., 2014 (RC estimation, BaRatin method)                                       #
# - Darienzo et al., 2021 (Method for gaugings segmentation)                                  #
# - Darienzo, 2021 (PhD thesis)                                                               #
# - Renard et al. 2006 (Bayesian Modelling, multi-block Metropolis MCMC approach)             # 
#                                                                                             #
# --------------------------------------------------------------------------------------------#
# - Recursive Segmentation of gaugings                                                        #
#---------------------------------------------------------------------------------------------#
# Make sure to have all following modules files in the apposite "Modules" folder:             # 
# - module_BaRatin.R                                                                          #
# - module_gaugings_segmentation.R                                                            #
# - module_gaugings_segmentation_plots.R                                                      #
# - module_recessions_analysis.R                                                              #
# - module_recessions_analysis_plots.R                                                        #
# - module_read_input.R                                                                       #
# - module_BaM_Results.R                                                                      #
# - module_shift_generator.R                                                                  #
# - module_prior_propagation.R                                                                #
# - module_BaRatinSPD.R                                                                       #
# - module_retrospective_analysis.R                                                           #
# - module_retrospective_analysis_plots.R                                                     #
#                                                                                             #
# and all following exe files for running BaM (fortran compiled exe for Bayesian modelling,   #
# developped by Benjamin Renard) within the apposite "BaM_exe" folder:                        #
# - BaM_2exp_pool2.exe                                                                        #
# - BaM_BaRatin_2.exe                                                                         #
# - BaM_BaRatin_SPD.exe                                                                       #
# - BaM_recession_multi_model_final.exe                                                       #
# - BaM_Segmentation2.exe)                                                                    #
# ant the mainconfiguration files "Config_BaM.txt" and the subfolders containing the          # 
# configuration files (mcmc options, input data, ...) and results of the Bayesian inference   #
# specific to each application.                                                               #                   
#                                                                                             #
# --------------------------------------------------------------------------------------------#
# IMPORTANT !!!!                                                                              #
# Please, notice that a few .dll libraries sometimes are needed to successfully run the exe.  #
# In particular, be sure to have the following .dll libraries in the folder "BaM_exe":        #
# - libgcc_s_dw2-1.dll                                                                        #
# - libgcc_s_sjlj-1.dll                                                                       #
# - libgfortran-3.dll                                                                         #
# - libquadmath-0.dll                                                                         #
# - libwinpthread-1.dll                                                                       #
#                                                                                             #
# All files with input data and options in the apposite "Case_studies/..." folder             #
# - "gaugings.csv" (input file with gaugings data)                                            #
# - "limni.csv" (input file with stage record (in cm))                                        #
# - "Options_General_input.R" (input R  file with the options for the RC estimation and info) #
# - "Options_Segment_gaugings.R" (input R file with the options for the segmentation)         #
# - "Options_Recession_analysis.R" (input R file with the options for recession analysis)     #
# - "Options_BaRatinSPD.R" (input R file with the options for the BaRatin-SPD)                #
#                                                                                             #
# !!! for input data files please use the same format as in the example of Meyras !!!         #
# --------------------------------------------------------------------------------------------#
#                                       Disclaimer                                            #
# --------------------------------------------------------------------------------------------#
# Please, notice that BayDERS is an experimental software. Further analysis is required for   #
# validating the proposed tools. We are not responsible for any loss (of data, profits,       #
# business, customers, orders, other costs or disturbances) derived by their use in the       #
# operational practice.                                                                       #
# The authorized user accepts the risks associated with using the software given its nature   #
# of free software. It is reserved to expert Users (developers or professionals) having prior # 
# knowledge in computer science, hydraulics and statistics.                                   #
#                                                                                             #
# For any question or feedback please contact us at one of the following adresses:            #
# - matteo.darienzo@cimafoundation.org                                                        #
# - jerome.lecoz@inrae.fr                                                                     #
# - benjamin.renard@inrae.fr                                                                  #
###############################################################################################
















##############################################################################################
#                               INITIALIZATION    (just run the following lines)             #
##############################################################################################
# Get the main directory folder automatically:         
                   if(!is.null(dev.list())) dev.off() # Clear plots
                   cat("\014") # Clear console
                   rm(list=ls())# Clean workspace
                   
                   #Libraries used  (it automatically detects and installs missing packages):
                   # options(repos = c(CRAN = "https://cran.rstudio.com"))
                   #Packages:
                   pack = c("rstudioapi",
                     "methods", "lattice", "gridExtra", "reshape","reshape2", "ggplot2", 
                     "GGally",  "grid", "devtools", "gtable",  "extrafont",  "chron", 
                     "randomcoloR",  "plotly", "mcmc", "coda", "RColorBrewer", "ggpubr", 
                     "cowplot", "svDialogs", "tcltk", "psych",  "mosaicData", 
                     "mosaicCore", "latex2exp",  "scales", "pracma" , 
                     "segclust2d",  "changepoint",  "minpack.lm", "gganimate", 
                     "magick", #"patchwork", "egg",
                     "dplyr", "tidyr", "gifski", "png", "transformr", "lubridate" 
                   )
                   # for future versions add: "airGR", "rmarkdown"
                   
                   install_pack <- function(x){
                     for( i in x ){
                       #  require returns TRUE invisibly if it was able to load package
                       if( ! require( i , character.only = TRUE ) ){
                         #  If package was not able to be loaded then re-install
                         install.packages( i , dependencies = TRUE )
                         #  Load package after installing
                         require( i , character.only = TRUE )
                       }
                     }
                   }
                   #  Then try/install packages:
                   install_pack( pack ) 

                   # set the directory
                   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
                   dir_code <- getwd()
                   setwd(dir_code)
                   # define the directory with the modules files:
                   dir.modules = paste0(dir_code, "/Modules")
                   
                   
                   
                   ###########################################################################
                   # Define the case study name (the same name of the case study folder !!!): 
                   # select the index of the case study folder that will appear in the popup 
                   # R window:
                   all_case_studies = list.files(paste0(dir_code, "/Case_studies"))
                   case_study_name <- all_case_studies[ as.numeric(
                                       dlgInput(c("Insert the name of case study folder: \n",
                                                 paste0(seq(1, length(all_case_studies), 1), " = ",
                                                        all_case_studies)), 
                                               Sys.info()[" "])$res)]
                   # or insert the folder name manually:
                   #case_study_name = "Meyras"     
                   ##########################################################################
                   


                   
                   
                   
                   
                   
                   
                   
                   
##############################################################################################
#                            READ GENERAL INPUTS FOR THE CASE STUDY                          #
##############################################################################################           
# Read all input data and options:
# IMPORTANT:
# Please, make sure that the configuration file with the general settings "Options_General.r" 
# is correctly configured and located in the case-study folder.                   
                   source(paste0(dir.modules,"/module_read_input.r"))

                   
                   
                  
                 
              
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
##############################################################################################
#                             1) SEGMENTATION OF GAUGINGS                                    #
##############################################################################################
# Author: Matteo Darienzo, Inrae
# Last update: 15/02/2021
# Objective: detection of ruptures in the time series of residuals between gaugings and base RC
# Procedure: 
# - 1) Estimation of RC0 (rating curve of reference with quantitative uncertainty)
# - 2) Computation of residuals between gaugings and RC0
# - 3) Segmentation of the residuals time series (call of Bam-segmentation)
# - 4) Adjustment of detected shift times and determination of sub-periods
# - 5) Iterative re-estimation of RC0 and re-segmentation of residuals (steps 1-4) 
#      for each sub-period
# For a detailed description of each step see Darienzo et al., 2021.
# Please, be sure to have completely filled the input files in the case study folder:
# - "Options_Segment_gauging.R"
# - "Options_General.R"
# - "Options_BaRatinSPD.R" (for estimating multiperiod RC with the results of segmentation)
###############################################################################################
                    dir.segment.g        =  paste0(dir.case_study,"/Results/segmentation_gaugings")
                    file.options.general =  paste0(dir.case_study,"/Options_General.R")
                    file.options.segment =  paste0(dir.case_study, "/Options_Segment_gaugings.R")
                    file.options.SPD     =  paste0(dir.case_study,"/Options_BaRatinSPD.R")
                    dir.create(dir.segment.g)
                    
                    #------------------------------------------------------------------------------------------------------------------------------
                    # - The function below "gaugings.segmentation" is located in "module_gaugings_segmentation.r" and 
                    #   uses other functions located in "module_BaRatin.r" and "module_gaugings_segmentation_plots.r" 
                    #   It performs the segmentation of gaugings, based on a Bayesian approach and identifies the periods with stable RC. 
                    # - It uses the software "BaM.exe" (Renard et al., 2006) for the bayesian estimation of the segmentation model.
                    # - inputs and options:
                    #   > dir_code             =  directory of the main program file, set by default.
                    #   > dir.exe              =  directory of the BaM executable, set by default.
                    #   > dir.case_study       =  directory of the case study folder, set by default.
                    #   > dir.segment.g        =  directory where to save the results, set by default
                    #   > file.options.general =  directory of the R-file containing the input options for the RC estimation and station info.
                    #                             By default this file must be located in the casestudy folder and named "General_input.R". 
                    #   > file.options.segment =  directory of the R-file containing the input options for the segmentation. 
                    #                             By default this file must be located in the casestudy folder and named 
                    #                             "Segment_gauging_options.R". 
                    #   > name.folder.results  =  name of the folder containing the segmentation results (e.g."Recursive_with_uncertainties")
                    #   > stage.record         =  dataframe of the stage record (columns: time, stage)
                    #   > gaugings             =  dataframe with the gaugings data-set (columns: time, stage, u_stage, discharge, u_discharge)
                    #   > official.dates       =  vector of dates of RC update given by the hydrometric service. by default equal to official.dates.
                    #                             which is defined in module_read_input.r and equal to FALSE if no official dates are available.
                    #   > colors.period        =  it designs the vector with the colors (highly contrasted) to illustrate the detected periods.
                    #   > save.all.results     =  if TRUE, all intermediate results will be saved, also the mcmc convergency tests.
                    #   > plot.results.only    =  if you have already performed the segmentation and you only want to plot results then type TRUE,
                    #                             otherwise type FALSE. 
                    #------------------------------------------------------------------------------------------------------------------------------
                    gaug.segment  = gaugings.segmentation(dir_code              =  dir_code, 
                                                          dir.exe               =  dir.exe,  
                                                          dir.case_study        =  dir.case_study,
                                                          dir.segment.g         =  dir.segment.g,
                                                          file.options.general  =  file.options.general,
                                                          file.options.segment  =  file.options.segment,
                                                          stage.record          =  df.limni,         
                                                          gaugings              =  data4BaRatin,       
                                                          official.dates        =  officialShiftsTime, 
                                                          colors.period         =  colo)
  
                    
                    
                    
                    # the function below (located in "module_gaugings_segmentation.r") computes the performance evaluation of the segmentation of gaugings:
                    # it gives a list with several metrics (e.g. accuracy, rmse). See the article for a better description.
                    performance.gaug.segm =   read.results.segmentation(dir.with.gaugings                 =  dir.case_study,
                                                                        dir.segment.g                     =  dir.segment.g,
                                                                        file.options.general              =  file.options.general, 
                                                                        file.options.segment              =  file.options.segment,
                                                                        gaugings                          =  data4BaRatin,
                                                                        known.shift.times                 =  officialShiftsTime,
                                                                        detected.shift.times.filename     =  "shift_times.txt",
                                                                        pdf_detected.shift.times.filename =  "pdf_ts.txt",
                                                                        use.other.seg.method              =  FALSE,
                                                                        compute.dates                     =  TRUE,
                                                                        initial.time                      =  initial.time)
                                               
                                                        
                    #If you want you can apply BaRatin-SPD to the obtained multi-period rating curve (Mansanarez et al.,2019):
                    # Please, be sure to have filled the file with options "Options_BaRatinSPD.R"
                    BaRatin_SPD.bac_app(dir_code                 =  dir_code,
                                        dir.BaM                  =  dir.exe, 
                                        dir.segment.g            =  dir.segment.g, 
                                        file.options.general     =  file.options.general, 
                                        file.options.segment     =  file.options.segment,
                                        file.options.SPD         =  file.options.SPD,
                                        stage.record             =  df.limni)

                    
                    
                
                
                    
 
                    
                    
                    
                    
                    
                    
                  
                    
                    
                    
                    
                    
####################################################################################################
#                                   STAGE-RECESSION ANALYSIS                                       #
####################################################################################################
# Procedure:
# - 1) Extraction of all recessions curves.
# - 2) Exponential regression of the curves (bayesian method: using BaM)
# - 3) Segmentation of the asymptote time series
# - 4) Plotting results
# - 5) Possibility to apply BaRatin-SPD to the obtained periods.
#---------------------------------------------------------------------------------------------------
                    # working directory:
                    dir.create(paste0(dir.case_study,"/Results/segmentation_recessions"))
                    dir.segment.rec          = paste0(dir.case_study,"/Results/segmentation_recessions")
                    file.options.general     =  paste0(dir.case_study,"/Options_General.R")
                    file.options.SPD         =  paste0(dir.case_study,"/Options_BaRatinSPD.R")
                    file.options.recessions  =  paste0(dir.case_study,"/Options_Recession_analysis.R")
                    
                    # !!! ATTENTION !!!!
                    # if limni is too dense you can filter it by setting "limni_filter" in the Options_general.R" file:  
                    
                    # Extraction of stage-recessions (from "module_recession_analysis.R"):
                    rec.extraction = recession.selection(dir.exe               =  dir.exe,
                                                         dir.segment.rec       =  dir.segment.rec,
                                                         file.options.general  =  file.options.general,
                                                         file.options.recess   =  file.options.recessions,
                                                         stage.record          =  df.limni)
                  
                    
                    # Estimation of stage-recessions (from "module_recession_analysis.R"):
                    rec.estimation = recession.regression(dir.exe               =  dir.exe,
                                                          dir.segment.rec       =  dir.segment.rec,
                                                          file.options.general  =  paste0(dir.case_study,"/Options_General.R"),
                                                          file.options.recess   =  paste0(dir.case_study,"/Options_Recession_analysis.R"),
                                                          data.recess           =  list(rec.extraction$d.h.selected, rec.extraction$t.real.good.h),
                                                          initial.time.rec      =  rec.extraction$t.real.good.h,
                                                          which.recession       =  c(1:rec.extraction$Ncurves),
                                                          stage.record          =  df.limni)  
                    
                    
                    # Segmentation of the time series of the stage-recession specific parameters (e.g. the asymptotic stage): 
                    rec.segmentation = recession.segmentation(dir.exe               = dir.exe,
                                                              dir.segment.rec       = dir.segment.rec,
                                                              file.options.general  = paste0(dir.case_study,"/Options_General.R"),
                                                              file.options.recess   = paste0(dir.case_study,"/Options_Recession_analysis.R"),
                                                              stage.record          = df.limni,
                                                              gaugings              = data4BaRatin,
                                                              initial.time          = initial.time,
                                                              colors.period         = colo)

                    
                    #If you want you can apply BaRatin-SPD to the obtained multi-period rating curve (Mansanarez et al.,2019):
                    # Please, be sure to have filled the file with options "Options_BaRatinSPD.R"
                    BaRatin_SPD.bac_app(dir_code                 =  dir_code,
                                        dir.BaM                  =  dir.exe, 
                                        dir.segment.g            =  dir.segment.rec, 
                                        file.options.general     =  file.options.general, 
                                        file.options.segment     =  file.options.segment,
                                        file.options.SPD         =  file.options.SPD,
                                        stage.record             =  df.limni)
                    
                    
                    
                    # Performance evaluation (to do yet!!!):
                    #######################################
                    # Plot DIC and other performance criteria to compare all used models and "Chi":
                    # plot.performance.model.comparison(model.names         = rec.model, # stage-recssion models
                    #                                   model.titles        = c(paste0("M", seq(1,length(rec.model)))),
                    #                                   bt.from.gaugingsss  = read.bt,  #river bed estimation
                    #                                   pdf.ts.gaugings     = pdf.ts.results.1,
                    #                                   data.annotate.off   = data.annotate.off,
                    #                                   time.limits         = limni.time.limits,
                    #                                   grid_limni.ylim     = c(-2 , 4.5 , 1))
              
            
                  

   
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
################################################################################################################
#   RESULTS OF STEP 1 OF THE RETROSPECTIVE ANALYSIS:  Gaugings + recessions
################################################################################################################
                    # This section performs the combining of results of segmentation of gaugings and of 
                    # stage-recession analysis. You will be asked to select the folder containing the 
                    # corresponding results and then, eventually, to remove some shift times, because in common
                    # between the two methods or because suspicous or not reliable according to your expertise.
                    
                    # charge all main directories (by deafault): 
                    dir.segment.g            =  paste0(dir.case_study,"/Results/segmentation_gaugings")
                    dir.segment.rec          =  paste0(dir.case_study,"/Results/segmentation_recessions")
                    dir.segment.gaug.rec     =  paste0(dir.case_study,"/Results/gaugings_&_recessions")
                    file.options.segment     =  paste0(dir.case_study,"/Options_Segment_gaugings.R")
                    file.options.general     =  paste0(dir.case_study,"/Options_General.R")
                    file.options.SPD         =  paste0(dir.case_study,"/Options_BaRatinSPD.R")
                    file.options.recessions  =  paste0(dir.case_study,"/Options_Recession_analysis.R")
                    
                    # read shift detection results from all available tools (gauings + recessions:
                    res_gaug_recess = read_all_results_shift_detection(dir.segment.g             = dir.segment.g, 
                                                                       dir.segment.rec           = dir.segment.rec,
                                                                       dir.segment.gaug.rec      = dir.segment.gaug.rec,
                                                                       stage.record              = df.limni,         
                                                                       gaugings                  = data4BaRatin,       
                                                                       official.dates            = officialShiftsTime, 
                                                                       colors.period             = colo,
                                                                       file.options.general      = file.options.general,
                                                                       file.options.segment      = file.options.segment,
                                                                       file.options.SPD          = file.options.SPD,
                                                                       file.options.recessions   = file.options.recessions)


                    
                    
                    
                    
                    
            
                    
                    
                    
                    
                    
                    
                    
                    
                    
####################################################################################################
#                                  SEDIMENT TRANSPORT PROXY ANALYSIS                               #
####################################################################################################
# Based on Darienzo 2021 PhD thesis manuscript.  
# it uses: BaM software, Benjamin Renard (BaRatin-SPD model, Linear model)
# last update : 15/11/2021
#
# Steps of the method:
#*********************
# 1)  estimate the b(k) through BaRation-SPD for each stable period found with other methods
# 2)  estimate the b(t) through a linear interpolation during the shifting events
# 3)  estimate the water depth y(t)
# 4)  fix a d50 and a taucritic ==> ycr
# 5)  define all events m for which y(t) > ycrit(ystart, ypeak, yend) 
# 6)  calculate the bedload at each instant t using the MPM model: qs(t)
# 7)  compute the cumulative bedload for each event m ==> qsc(m)
# 8)  Among the events k (found with other methods) find the event with the lowest qsc  ==>
#     qsc* = min(qsc(k))
# 9)  determine if other events have qsc(m) >= qsc* 
# 10) if yes, a new set of shift times ts(p) is found, thus a new BaRatinSPD can be 
#     performed on the new periods.
# 11) to each known qsc(p) is then associated a deltab (b1 or/and b2).
# 12) Linear estimation of the relation between qsc and deltab
#---------------------------------------------------------------------------------------------------
                    # working directory:
                    dir.Sedim.transp            =  paste0(dir.case_study,"/Results/sediment_transport")
                    file.options.general        =  paste0(dir.case_study,"/Options_General.R")
                    file.options.segment        =  paste0(dir.case_study,"/Options_Segment_gaugings.R")
                    file.options.SPD            =  paste0(dir.case_study,"/Options_BaRatinSPD.R")
                    file.options.recessions     =  paste0(dir.case_study,"/Options_Recession_analysis.R")
                    file.options.ST             =  paste0(dir.case_study,"/Options_Sediment_transport.R")
                    
                    
                    # read the reference morphogenic events (coming from user or from recession+gaugings analysis):
                    # the user is asked to confirm or to modify the set of shift times:
                    reference.shift.times       =  read.refer.morphogenic.events(dir.exe                =  dir.exe,
                                                                                 dir.case_study         =  dir.case_study,
                                                                                 dir.Sedim.transp       =  dir.Sedim.transp,
                                                                                 file.options.general   =  file.options.general,
                                                                                 file.options.segment   =  file.options.segment,
                                                                                 file.options.recessions=  file.options.recessions,
                                                                                 file.options.ST        =  file.options.ST,
                                                                                 stage.record           =  df.limni,         
                                                                                 gaugings               =  data4BaRatin,       
                                                                                 official.dates         =  officialShiftsTime, 
                                                                                 colors.period          =  colo,
                                                                                 res_gaug_recess        =  res_gaug_recess)   
                    # Apply a SPD Baratin analysis to the reference periods just defined:
                    # this will provide the shift magnitudes of the morphogenic shifts.
                    SPD.reference.periods       = BaRatin_SPD.bac_app(dir_code                 =  dir_code,
                                                                      dir.BaM                  =  dir.exe, 
                                                                      dir.segment.g            =  dir.Sedim.transp, 
                                                                      file.options.general     =  file.options.general, 
                                                                      file.options.segment     =  file.options.segment,
                                                                      file.options.SPD         =  file.options.SPD,
                                                                      stage.record             =  df.limni)
                                        
                    
                    # Compute sediment transport: 
                    # and find other potential shift times:
                    ST.segm = sediment.transport.segmentation(dir.sed.transp       = SPD.reference.periods$dir.reference, 
                                                              dir.sed.transp.SPD   = SPD.reference.periods$dir.SPD,
                                                              res.SPD              = SPD.reference.periods$res.SPD,
                                                              file.options.general = file.options.general,
                                                              file.options.ST      = file.options.ST,
                                                              df.limni.ST          = df.limni, 
                                                              officialShiftsTime   = officialShiftsTime)  
                    
                    
                    
                    # apply a linear regression to the relation V - deltab:
                    # using deltab (from baratinSPD) and V (from cumulative sediment transport for each reference event):
                    ST.linear.estim = linear.estimation(dir_code             = dir_code, 
                                                        dir.sed.transp       = SPD.reference.periods$dir.reference, 
                                                        df.deltab.V          = ST.segm$df.rel.deltab.V.TOT,      # df.rel.deltab.qscum.TOT,   #ST.SPD[[1]],
                                                        file.options.general = file.options.general,
                                                        file.options.ST      = file.options.ST)
                                                        
                    
                  
                

                    
                    
                    
                    
                    
                    
                    
                
                  
          