#############################################################################################################
#                                   SEGMENTATION Of GAUGINGS (RECURSIVE METHOD):
#############################################################################################################
gaugings.segmentation <- function(dir_code, 
                                  dir.exe, 
                                  dir.case_study,
                                  dir.segment.g,
                                  file.options.general, 
                                  file.options.segment,
                                  stage.record,
                                  gaugings,
                                  official.dates,
                                  colors.period) {
#############################################################################################################
   # create directories for results:
   dir.segmentation      = paste0(dir_code,"/BaM_exe/Segmentation")
   # read inputs and options for computation:
   source(file.options.general)
   source(file.options.segment)
   dir.create(paste0(dir.segment.g,"/", name.folder.results))
   dir.segment.gaug      = paste0(dir.segment.g,"/", name.folder.results) # dir. with the results of gaugings segmentation
   # grid_stage  = define.grid(gaugings = gaug.synthet[[sim.index]][[1]],  # Defining the grid of stage for this data-set:
   #                           sequence = seq(-10, 50, 0.5))  # MODIFY THIS, IF NEEDED !!
   
   message("****************************************************************")
   message("                    SEGMENTATION OF GAUGINGS                    ")
   message("****************************************************************")
   

   #RC Parameter a:
   if (propagat == TRUE){
      message("- Prior a: propagating the priors of geometric/physical properties.")
      for (c in 1:ncontrols){
         if (control.type[c] == "rect.channel") {
            if (a.distr== "LogNormal"){
               # Transform priors from gaussian to lognormal:
               st_Bc.prior[c] = round(Transf_Gauss_lognorm(E = Bc.prior[c], stdev = st_Bc.prior[c])$sd, digits = 2)
               st_KS.prior[c] = round(Transf_Gauss_lognorm(E = KS.prior[c], stdev = st_KS.prior[c])$sd, digits = 2)
               st_S0.prior[c] = round(Transf_Gauss_lognorm(E = S0.prior[c], stdev = st_S0.prior[c])$sd, digits = 2)
            }
         } else if (control.type[c] == "rect.weir") {
            if (a.distr== "LogNormal"){
               # Transform priors from gaussian to lognormal:
               st_Cr.prior[c] = round(Transf_Gauss_lognorm(E = Cr.prior[c], stdev = st_Cr.prior[c])$sd, digits = 2)
               st_g.prior[c]  = round(Transf_Gauss_lognorm(E = g.prior[c],  stdev = st_g.prior[c])$sd,  digits = 2)
               st_Bw.prior[c] = round(Transf_Gauss_lognorm(E = Bw.prior[c], stdev = st_Bw.prior[c])$sd, digits = 2)
            }
         } else {
            message("**** Fatal Input error: you have selected a wrong control type!")
            message("control types available: 'rect.channel' and 'rect.weir'.")
            message("Please, check again!")
         }
      }
   } else {
      # prior for 'a' is given.
      for (c in 1:ncontrols){
         if (a.distr== "LogNormal"){
            # Transform priors from gaussian to lognormal:
            st_a.prior[c] = round(Transf_Gauss_lognorm(E = a.prior[c], stdev = st_a.prior[c])$sd, digits = 4)
         }
      }
   } 
   
   
   #Parameter b:
   for (c in 1:ncontrols){
      if (b.distr== "LogNormal"){
         # Transform priors from gaussian to lognormal:
         # function which takes gaussian mean and stdev and gives lognormal mean and stdev
         st_b.prior[c] = round(Transf_Gauss_lognorm(E = b.prior[c], stdev = st_b.prior[c])$sd, digits = 4)
      }
   }
   #stdev.var.param.initial = st_b.prior

   #Parameter c:
   for (c in 1:ncontrols){
      if (c.distr== "LogNormal"){
         # Transform priors from gaussian to lognormal:
         # function which takes gaussian mean and stdev and gives lognormal mean and stdev
         st_c.prior[c] = round(Transf_Gauss_lognorm(E = c.prior[c], stdev = st_c.prior[c])$sd, digits = 4)
      }
   }
   


   
   
   
   
   # start:
   #########################################################################################################
   if (plot.results.only == FALSE) {    # if TRUE it goes directly at the end, to read and plot the results   
   #########################################################################################################
   ################
   #Initialisation:
   ################
   # iteration indexes: "depth search tree approach"
            i = 1;  seg.period = 1;  seg.iter  = 1; level = 0;  iteration.list = list()
            end.end    = FALSE;  i_init  = 0; i_final = 0;   i_init[1] = 1; 
            i_final[1] = length(Q_Gaug); tss_tot_ns = c(1);  final.period = NULL;
            acf_resid  = NULL; Pacf_resid =NULL;  autocorr_lag_resid =NULL; 
   # segments means:
            mean.of.segments = NULL;  mu.results.df = NULL;
   #shift times:
            times.of.shift.MAP <- NULL;  t.q10 <- NULL;  t.q90 <- NULL; t.q2 <- NULL;  
            t.q97  <- NULL;  ts.all.real <- NULL; ts.all.real.2 <- NULL; ts.all.MAP <- NULL; 
            ts.all.q2 = NULL; ts.all.q10 = NULL; ts.all.q90 = NULL;  ts.all.q97 = NULL;
            ts.morpho.real = NULL; ts.morpho.MAP = NULL; ts.morpho.q2 = NULL; ts.morpho.q97 = NULL;
            tau.results.df = NULL; pdf.ts = cbind();           
            shift.results.df =  data.frame(tMAP  = double(), 
                                           treal = double(), 
                                           t2    = double(),
                                           t10   = double(), 
                                           t90   = double(),
                                           t97   = double())
   # structural error model:
            gamma1.P      = NULL; 
            gamma2.P      = NULL; 
            g1.distr      = g1.distr.type; 
            g2.distr      = g2.distr.type;
            gamma1.P[[1]] = g1.prior;  
            gamma2.P[[1]] = g2.prior;
   # plots:
            criteria.plot = NULL; ts.ggplot = NULL; gaug.plot = NULL; seg.plot = NULL; 
   # all gaugings:
            t_Gaug         = gaugings$t
            Q_Gaug         = gaugings$Q
            uQ_Gaug        = gaugings$uQ
            h_Gaug         = gaugings$h
   
             
            
            
   #*************************************************************************************************
   # "Top-down" Recursive segmentation (RC estimation and residuals segmentation at each subperiod):
   #*************************************************************************************************
   message("Options selected:")
   if (recursive == TRUE){
            print("- Top-down Recursive segmentation.")
   } else {
            print("- Single-pass segmentation.")
   }
   if (resid.uncertaint == TRUE) {
            print("- Accounting for residuals uncertainties (from gaugings + RC)")
   } else {
            print("- Neglecting residuals uncertainties!")
   }
   if (shift.time.adjustment.type == 1){
            print("- Shift time adjustment: always the Maximum A Posterior estimate.")
   } else if  (shift.time.adjustment.type == 2){
            print("- Shift time adjustment: always the largest flood peak in the 95% CI.") 
   } else {
            print("- Shift time adjustment: manual selection at each iteration") 
   }
   print(paste0("- Criterion for the choice of the number of shifts: ", criterion))
   print(paste0("- Prior for the segments means ~ N(",prior.mu.segm[1],",", prior.mu.segm[2],")"))
   print(paste0("- minimum time between two shifts  =",tmin, " days"))
   print(paste0("- minimum number of points in a segment  =",Nmin))
   print(paste0("- Maximum number of segments in the series at each iteration  =",nSmax))
   
   
   
   
   
   
   
message("
Start segmentation ...
*****************************************************************
A few information:
*****************************************************************
- Please, notice that it may take some time (e.g., hours). 
- Each iteration is composed of (in sequence):
  1) the baseline RC estimation (BaRatin method)
  2) computation of residuals between gaugings and RC
  3) segmentation of residuals with increasing number of segments
  4) choice of the optimal segmentation (e.g. minimising the BIC)
  5) adjustment of shift times (if any)
*****************************************************************
")
   
   while(end.end == FALSE) {
      #folder creation for the current iteration:
      dir.create(paste0(dir.segment.gaug,"/it",seg.iter))
      dir.create(paste0(dir.segment.gaug,"/it",seg.iter,"/BaRatin"))
      dir.create(paste0(dir.segment.gaug,"/it",seg.iter,"/Segmentation"))
      #associate the name of the iteration: e.g. 3.2.1.1 (3rd level, branch 2, branch 1, branch 1)
      # iterat.info[[seg.iter]] = data.frame(level        = level[seg.iter],
      #                                      n.digits     = level[seg.iter] + 1,
      #                                      nS.index.new     = nS.index[seg.iter],
      #                                      nS.index.old = iterat.info[[seg.period]]$nS.index)
      # name.iteration[seg.iter] = paste0("Iteration ")
      # for (lev in 1:iterat.info[[seg.iter]]$n.digits) {
      #      name.iteration[seg.iter] = paste0(name.iteration[seg.iter], ".", (lev-1))
      # }
      
      #directories saving:
      dir.seg.iter    <- paste0(dir.segment.gaug,"/it",seg.iter)
      dir.seg.gaug    <- paste0(dir.seg.iter,"/Segmentation")
      dir.seg.BaRatin <- paste0(dir.seg.iter,"/BaRatin")
      
      #gaugings of the current period "P":
      tP             = c(t_Gaug[i_init[seg.period] : i_final[seg.period]])
      QP             = c(Q_Gaug[i_init[seg.period] : i_final[seg.period]])
      uQP            = c(uQ_Gaug[i_init[seg.period] : i_final[seg.period]])
      hP             = c(h_Gaug[i_init[seg.period] : i_final[seg.period]])
      data_BaRatin.P = data.frame("h"      = hP,
                                  "Q"      = QP,
                                  "uQ"     = uQP,
                                  "Period" = 1)
      nobs.gaug.P    = length(hP)
      
      
      #limni (stage record) of the current period "P":
      if (is.null(stage.record$t_limni)== FALSE) {
         t_limni.P =0; h_limni.P =0; j =0
         message("- loading the stage record of the period:")
         pb <- txtProgressBar(min = 0,               # Minimum value of the progress bar
                              max = length(stage.record$t_limni), # Maximum value of the progress bar
                              style = 3,             # Progress bar style (also available style = 1 and style = 2)
                              width = 50,            # Progress bar width. Defaults to getOption("width")
                              char = "=")            # Character used to create the bar
         
         for (i in 1:length(stage.record$t_limni))  {
            if ((stage.record$t_limni[i] >= tP[1])&(stage.record$t_limni[i] <= tail(tP,1)))  {
               j =j+1
               t_limni.P[j] =stage.record$t_limni[i]
               h_limni.P[j] =stage.record$h_limni[i]
            }
            setTxtProgressBar(pb, i)
         }     
         close(pb)
         write(h_limni.P, file =paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/limni.txt"),
               ncolumns = 1, sep = "\n")
         nobs.limni.P = length(h_limni.P)
      }
      
      
      
      # Estimation of Rating curve RC0 for the current period (BaRatin application, Le Coz,2014):
      write("h", file =paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Gaugings_data.txt"), sep ="\n")
      dir.BaRatin.config = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2")
      write.table(data_BaRatin.P, file =paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Gaugings_data.txt"),
                  sep="\t",row.names=FALSE, col.names = TRUE)
      setwd(dir.exe)
      Hmin.P = round(min(hP), digits = 2)
      
      # RC estimation (BaRatin method):
      BaRatin_app(dir               = dir.BaRatin.config, #directory of config files for BaM
                  t_Gaug            = tP, 
                  Q_Gaug            = QP, 
                  uQ_Gaug           = uQP, 
                  h_Gaug            = hP,   # gaugings 
                  nobs.gaug         = nobs.gaug.P, 
                  nlimni            = nobs.limni.P,  # gaugings and limni
                  propagat          = propagat,    #do propagation of priors ?
                  b.distr           = b.distr , # RC priors
                  a.distr           = a.distr, # RC priors
                  c.distr           = c.distr,# RC priors
                  a.prior           = a.prior, # RC priors
                  st_a.prior        = st_a.prior, # RC priors
                  c.prior           = c.prior, # RC priors
                  st_c.prior        = st_c.prior, # RC priors
                  b.prior           = b.prior, # RC priors
                  st_b.prior        = st_b.prior,  # RC priors
                  Bw.prior, Cr.prior, g.prior,  # if propagation TRUE 
                  Bc.prior, KS.prior, S0.prior,
                  st_Bw.prior, st_Cr.prior, st_g.prior, 
                  st_Bc.prior, st_KS.prior, st_S0.prior,
                  ncontrol          = ncontrols, #hydraulic configuration
                  M                 = M, #hydraulic configuration
                  remnant.err.model = remnant.err.model, 
                  g1.prior          = gamma1.P[[seg.iter]], 
                  g2.prior          = gamma2.P[[seg.iter]],
                  g1.distr.type     = g1.distr,
                  g2.distr.type     = g2.distr,  # remnant error model priors
                  predictionRC      = predictionRC, 
                  predictionQt      = predictionQt, 
                  predictionPrior   = predictionPrior, 
                  simMCMC           = simMCMC,  
                  mcmc.prior        = mcmc.prior,        # predictions choice fro BaM
                  Ncycles           = Ncycle.baratin,    # mcmc oprtions: n.cycles
                  Hmin              = Hmin.P, 
                  Hmax              = grid_RC.xlim[2],   # grid limits for plots
                  iter              = seg.iter)          # copy the results files from BaM to the new folder:

      
      Dir.BaRatin.exe <- paste0(dir_code,"/BaM_exe/BaM_BaRatin_2")
      list.of.files   <- c(
         #paste(Dir.BaRatin.exe,"/Qt_maxpost.spag", sep=""), 
         #paste(Dir.BaRatin.exe,"/Qt_TotalU.spag", sep=""), paste(Dir.BaRatin.exe,"/Qt_TotalU.env", sep=""),
         paste0(Dir.BaRatin.exe,"/Qrc_maxpost.spag"),
         paste0(Dir.BaRatin.exe,"/Qrc_TotalU.spag"),
         paste0(Dir.BaRatin.exe,"/Qrc_TotalU.env"),
         paste0(Dir.BaRatin.exe,"/Qrc_ParamU.spag"),
         paste0(Dir.BaRatin.exe,"/Qrc_ParamU.env"),
         paste0(Dir.BaRatin.exe,"/Results_MCMC_Cooked.txt"), 
         paste0(Dir.BaRatin.exe,"/Results_Residuals.txt"),
         paste0(Dir.BaRatin.exe,"/Results_Summary.txt"),
         paste0(Dir.BaRatin.exe,"/Config_Model.txt"),
         paste0(Dir.BaRatin.exe,"/Hgrid.txt"),
         paste0(Dir.BaRatin.exe,"/Gaugings_data.txt")
         )
      for (ll in 1:length(list.of.files)) {
         file.copy(list.of.files[ll], dir.seg.BaRatin, overwrite = TRUE)
      }
      # Convergence and correlation tests (function from "module_BaM_results.R"):
      converg = Convergence.test(dir.seg  = dir.seg.BaRatin, 
                                 npar     = ncontrols*3 + 3 + ncontrols,  #(b,a,c) + (g1,g2) + (logpost) + (k of ncontrols)
                                 dir.plot = dir.seg.BaRatin)
      #Plot BaM results for the mcmc: traceplots, pdf, autocorrelation,Gelman, and RC estimation (Y-lin, Y-log, Ylog-Xlog scales):
      plot.Baratin <- plot.Results_BaM(workspace       =  dir.seg.BaRatin, 
                                       iter            =  seg.iter,
                                       plot.RC.lin     =  plot.RC.lin, 
                                       plot.RC.log     =  plot.RC.log, 
                                       plot.RC.loglog  =  plot.RC.loglog,
                                       grid_RC.ylim    =  grid_RC.ylim, 
                                       grid_RC.xlim    =  grid_RC.xlim, 
                                       grid_RC.xstep   =  grid_RC.xstep, 
                                       grid_RC.ystep   =  grid_RC.ystep,
                                       ticks_RC.y.log  =  ticks_RC.y.log, 
                                       RC.x.labels     =  RC.x.labels, 
                                       RC.y.labels     =  RC.y.labels, 
                                       u.m.Qgaug       =  u.m.Qgaug,  
                                       u.m.Hgaug       =  u.m.Hgaug,
                                       ncontrol        =  ncontrols, 
                                       hP              =  hP,  # current Period gaugings
                                       QP              =  QP, 
                                       uQP             =  uQP, 
                                       h_Gaug          =  h_Gaug, # all gaugings
                                       Q_Gaug          =  Q_Gaug, 
                                       uQ_Gaug         =  uQ_Gaug)
      
      df.gaug.P            = data.frame(hP, QP, uQP, tP)
      df.gaug.tot          = data.frame(h_Gaug, Q_Gaug, uQ_Gaug, t_Gaug)
      summary.BaRatin      = read.table(file = paste0(dir.seg.BaRatin, "/Results_summary.txt"), header=TRUE)
      
      # Update structural error model parameters gamma1 and gamma2 for the next subperiods:
      g1.distr             = "Uniform"
      g2.distr             = "Uniform"
      
      if ( type.stdev.gamma.baratin == 1) {
         gamma1.P[[seg.iter]] = c(0, summary.BaRatin[5, ncontrols*3+1], 
                                  summary.BaRatin[5, ncontrols*3+1]/2)
         gamma2.P[[seg.iter]] = c(0, summary.BaRatin[5, ncontrols*3+2], 
                                  summary.BaRatin[5, ncontrols*3+2]/2)
      } else if ( type.stdev.gamma.baratin == 2){
         gamma1.P[[seg.iter]] = c(0, summary.BaRatin[5, ncontrols*3+1] + 2*summary.BaRatin[11, ncontrols*3+1], 
                                  summary.BaRatin[5, ncontrols*3+1]/2)
         gamma2.P[[seg.iter]] = c(0, summary.BaRatin[5, ncontrols*3+2] + 2*summary.BaRatin[11, ncontrols*3+2], 
                                  summary.BaRatin[5, ncontrols*3+2]/2)

      } else if (type.stdev.gamma.baratin == 3) {
         g1.distr = g1.distr.type
         g2.distr = g2.distr.type
         gamma1.P[[seg.iter]] = g1.prior
         gamma2.P[[seg.iter]] = g2.prior
      }
      
      
      # other priors for structural error model parameters:
      #***************************************************
      # gamma1.P[[seg.iter]] = c(0, 2*summary.BaRatin[5, ncontrols*3+1], summary.BaRatin[5, ncontrols*3+1]/2)
      # gamma2.P[[seg.iter]] = c(0, 2*summary.BaRatin[5, ncontrols*3+2], summary.BaRatin[5, ncontrols*3+2]/2)
      #
      # gamma1.P[[seg.iter]] = c(log(summary.BaRatin[5, ncontrols*3+1]),
      #                          (log(1+ summary.BaRatin[11, ncontrols*3+1]/
      #                                  (summary.BaRatin[5, ncontrols*3+1]^2)))^0.5, # 0.5, #(summary.BaRatin[11, ncontrols*3+1]),
      #                          summary.BaRatin[5, ncontrols*3+1])
      # gamma2.P[[seg.iter]] = c(log(summary.BaRatin[5, ncontrols*3+2]),
      #                         (log(1+ summary.BaRatin[11, ncontrols*3+2]/
      #                                 (summary.BaRatin[5, ncontrols*3+2]^2)))^0.5, #0.5, #(summary.BaRatin[11, ncontrols*3+2]),
      #                          summary.BaRatin[5, ncontrols*3+2])
      
      
      data_baratin_with_time =  cbind(data_BaRatin.P, tP)
      write.table(data_baratin_with_time, paste0(dir.seg.BaRatin,"/data_with_time.txt"))
      # maxQt <-  read.table(file = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Qt_maxpost.spag",sep=""))
      # maxQt.data <- data.frame(t_limni, maxQt)
      # Qt <- Qtimeseries(t_limni,maxQt[,1],m)
      # Qt.df <- data.frame(t_limni, Qt)
      # ht.df <- stage.record

      
      
      
      
      #***********************************************
      #compute residuals and associated uncertainties:
      #***********************************************
      residuals   <- read.table(file=paste0(dir.seg.BaRatin,"/Results_Residuals.txt"), header = TRUE)
      resid       <- residuals[,6]
      stand.resid <- residuals[,7]
      sigma.tot   <- resid/stand.resid
      Shift.Q     <- alpha_segm.Q(tP, QP, resid, sigma.tot)
      
      
      
      
      
      
      #**************************************************************************************************************     
      # Segmentation of residuals by package R "changepoint" (multi change point) cpt.mean() or cpt.meanvar
      #**************************************************************************************************************
      # functions: cpt.mean() ==> detection of changes in the mean 
      #            cpt.meanvar() ==> detection of changes in the mean and the variance of data
      if (use.other.seg.method != FALSE) {
      ####################################
         print("Using R package changepoint (function cpt.mean) instead !")
         
         cp1  = Shift.Q$alpha
         nmax = min(length(cp1)/2+1, (nSmax-1))   # nSmax - 1 = max number of change points for BinSeg !!
         cpd1.Rpackage = cpt.mean(data            = cp1, 
                                  test.stat       = 'Normal', 
                                  method          = use.other.seg.method, 
                                  Q               = nmax, 
                                  penalty         = criterion, 
                                  class           = TRUE,
                                  param.estimates = TRUE, 
                                  minseglen       = Nmin)
         cpd1 = cpts(cpd1.Rpackage)
         cpd1 = cpd1 + 1
         param.est(cpd1.Rpackage)
         summary(cpd1.Rpackage)
         shift.times.other.method = Shift.Q$alpha_t[cpd1]
         jpeg(paste0(dir.segment.gaug,"/residuals_segm_",use.other.seg.method,".jpg"), width = 10, height = 6, units = 'in', res = 300)
         plot(cpd1.Rpackage, cpt.width=3, cpt.col='blue')
         dev.off()
         dev.set(dev.prev())
         break
      }
      

      
      
      
      
      #**************************************************************************************************************     
      # Segmentation of residuals (iterative segmentation with increasing number of segments):
      #**************************************************************************************************************
      # initialisation:   AIC = Aikake Information Criterion
      #                   BIC = Bayesian Information Criterion or SIC (Schwarz)
      #                   mBIC = modified Bayesianb Information Criterion
      #                   HQC = Hannan-Quinn information Criterion
      #                   DIC = Deviance Informatiojn Criterion
      AIC =0; AICc=0; BIC=0; mBIC = 0; HQC = 0; DIC=0;
      npar =0 ; maxpost = 0; varLogpost=0; loglikelihood = 0;  maxLikelihood = 0; varLogLikelihood=0;
      N = length(Shift.Q$alpha) #Number of observations
      
      if (resid.uncertaint==TRUE) {
         Data.segm <- data.frame("t"=Shift.Q$alpha_t, "Y"=Shift.Q$alpha, "uY"=sigma.tot)
      } else {
         Data.segm <- data.frame("t"=Shift.Q$alpha_t, "Y"=Shift.Q$alpha, "uY"=0)
      }
      write.table(Data.segm, paste0(dir_code,"/BaM_exe/Segmentation/Segm_data.txt"),sep="\t",row.names=FALSE)
      
      
      ##########
      if (N>2) {
      ##########
         end.seg = FALSE
         i = 1
         while ((i <= nSmax) & ((N-i) >=1)) {
            start_time = Sys.time()
            nS         = i
            npar[nS]   = 2*nS
            Ncycles    = Ncycle.segment
            converg    = FALSE
            print(c("n.segm = ",nS))
            
            while ((converg == FALSE) & (Ncycles < Ncycles.max)){
            setwd(dir.exe)
            dir.create(paste0(dir.seg.gaug,"/nS",nS))
            dir.nS <- paste0(dir.seg.gaug,"/nS",nS)
            if (resid.uncertaint == TRUE){
                prior.gamma.segm = prior.gamma.segm.with.U
            } else {
                prior.gamma.segm = prior.gamma.segm.without.U
            }
            # launch BaM.exe with segmentation model: 
            Segmentation_app(dir_code         = dir_code,
                             nobs             = length(Shift.Q$alpha), 
                             nS               = nS, 
                             tmin             = tmin,
                             Nmin             = Nmin,
                             tfirst           = Shift.Q$alpha_t[1],
                             tfin             = tail(Shift.Q$alpha_t,1),
                             ncycles          = Ncycles,
                             nmcmc            = Nmcmc.segment,
                             Nslim            = Nslim.segment,
                             tP               = Shift.Q$alpha_t,
                             resid.uncertaint = resid.uncertaint,
                             prior.mu         = prior.mu.segm,
                             prior.gamma      = prior.gamma.segm)
            end_time    = Sys.time()
            comput_time = end_time-start_time
            write.table(comput_time, file=paste0(dir.nS,"/comput_time_nS",nS,".csv"), sep=";")
            list.of.files.segment = c(
                                      paste0(dir.segmentation,"/Config_model.txt"),
                                      paste0(dir.segmentation,"/Segm_data.txt"),
                                      paste0(dir.segmentation,"/Results_MCMC_cooked.txt")
                                      )
                                     for (ll in 1:length(list.of.files.segment)) {
                                           file.copy(list.of.files.segment[ll], dir.nS, overwrite = TRUE)
                                     }
            mcmc.segm    <- read.table(file=paste0(dir.segmentation,"/Results_MCMC_cooked.txt"),header=TRUE)
            resid.segm   <- read.table(file=paste0(dir.segmentation,"/Results_Residuals.txt"),header=TRUE)
            summary.segm <- read.table(file=paste0(dir.segmentation,"/Results_Summary.txt"),header=TRUE)
            write.table(mcmc.segm,    file=paste0(dir.nS,"/Results_MCMC_cooked_","it",seg.iter,"_","nS",nS,".csv"), sep=";")
            write.table(resid.segm,   file=paste0(dir.nS,"/Results_Residuals_","it",seg.iter,"_","nS",nS,".csv"), sep=";")
            write.table(summary.segm, file=paste0(dir.nS,"/Results_Summary_","it",seg.iter,"_","nS",nS,".csv"), sep=";")
            
            if (save.all.results == TRUE) {
               plot.mcmc.segment(workspace = dir.nS, 
                                 seg.iter  = seg.iter, 
                                 nS        = nS)
               converg = Convergence.test.segment(dir.seg  = dir.segmentation, 
                                                  npar     = npar[nS], 
                                                  dir.plot = dir.nS)
               if (converg==FALSE) {
                   Ncycles = Ncycles + 1000
                   print(paste0("Segmentation ",nS," does not converged well !!!"))
                   print(paste0("increasing mcmc cycles ",Ncycles))
               }
            } else {
               converg=TRUE
            }
            }
            logpost           = mcmc.segm[,(2*nS+1)]
            loglikelihood     = 0;
            singlelikelihoods = 0; 
            single.prior.mu   = 0; 
            single.prior.tau  = 0;
            len.mcmc          = length(mcmc.segm[,1])
            priors.mu         = matrix(NA, nrow = length(mcmc.segm[,1]), ncol = nS)
            for (j in 1:nS) {
               priors.mu[,j]  = dnorm(mcmc.segm[,j], mean = prior.mu.segm[1], sd = prior.mu.segm[2], log = TRUE)
            }
            #
            if (nS > 1) {
               priors.tau = matrix(0, nrow =len.mcmc, ncol = nS-1)
            } else {
               priors.tau = matrix(0, nrow = len.mcmc, ncol = 1)
            }
            #
            priors.gamma = 0
            if (resid.uncertaint == FALSE) { 
               priors.gamma = dunif(mcmc.segm[,2*nS],
                                    min = prior.gamma.segm.without.U[1], 
                                    max = prior.gamma.segm.without.U[2], log = TRUE)
            } else {
               priors.gamma = dlnorm(mcmc.segm[,2*nS], 
                                     meanlog = prior.gamma.segm.with.U[1], 
                                     sdlog   = prior.gamma.segm.with.U[2], log = TRUE)
            }
            #
            logprior = 0
            for (ll in 1:len.mcmc){
               logprior[ll] = sum(priors.mu[ll,]) + sum(priors.tau[ll,])  + priors.gamma[ll]   
            }
            loglikelihood = logpost - logprior
            #
            df.mcmc.LL = data.frame(loglikelihood = loglikelihood,
                                    logprior      = logprior,
                                    logposterior  = logpost)
            write.table(df.mcmc.LL, file=paste0(dir.nS,"/likelihood_","it",seg.iter,"_","nS",nS,".csv"), sep=";")
            #
            maxpost[nS]          = max(logpost)
            maxLikelihood[nS]    = max(loglikelihood)
            Likelihood.maxpost   = loglikelihood[which.max(logpost)]
            varLogpost[nS]       = var(logpost)
            varLogLikelihood[nS] = var(loglikelihood)
            #MeanDev             =-2*mean(logpost)
            MeanDev              = -2*mean(loglikelihood)
  
            # criteria for model selection:
            #------------------------------------------------------------------------------------------
            AIC[i]  =  2*npar[nS]- 2*maxLikelihood[nS]    # ref: Aikake,1974
            BIC[i]  =  log(N)*(npar[nS]) -2*maxLikelihood[nS]  #ref: Schwarz,1978
            HQC[i]  =  2*(npar[nS])*log(log(N)) - 2*maxLikelihood[nS] # ref: Hannan-Quinn
            DIC[i]  =  MeanDev + 2*varLogLikelihood[nS]   #ref: Gelman 2004 "Bayesian data analysis"
            #AICc[i] =  AIC[i] +  2*(npar[nS])*(npar[nS]+1)/(N-npar[nS]-1) #ref: Hurvich and Tsai,1989
            #------------------------------------------------------------------------------------------
            i = i+1
         }
         
         
         
         
         
         
         
         #*******************************************************************************************************
         # Analysis of segmentation results for this period:
         #*******************************************************************************************************
         # min values:
         BICmin = which.min(BIC);
         AICmin = which.min(AIC);
         HQCmin = which.min(HQC);
         DICmin = which.min(DIC);
         #AICcmin = which.min(AICc, na.rm=TRUE);
         
         # dataframe with all criteria (for the plot):
         criteria.df = data.frame(BIC = BIC, AIC = AIC, HQC = HQC, DIC = DIC, # AICc
                                  x   = seq(1, nS, 1),  #grid: increasing Number of segments
                                  BICmin = BICmin, AICmin = AICmin, HQCmin = HQCmin, DICmin = DICmin)  #AICcmin
         criteria.plot[[seg.iter]] = model.selection.plot(criteria.df, 
                                                          dir.seg.gaug, 
                                                          seg.iter)
         
         if (criterion == "AIC") {
            nS = AICmin
         } else if (criterion == "AICc") {
            nS = AICcmin
         } else if (criterion == "BIC") {
            nS = BICmin
         } else if (criterion == "DIC") {
            nS = DICmin
         } else if (criterion == "HQC") {
            nS = HQCmin
         }
         
         print(paste0("=========> Optimal number of segments (considering the minimum ", criterion,") = ", nS))
         dir.nS.ok <- paste(dir.seg.gaug,"/nS",nS, sep="")
         
         #gamma error structural update for next subperiods:
         if (nS >1) {
            if (seg.iter >1) {
               gamma1.P = append (gamma1.P, rep(list(gamma1.P[[seg.iter]]), nS), after=seg.iter)
               gamma2.P = append (gamma2.P, rep(list(gamma2.P[[seg.iter]]), nS), after=seg.iter)
            } else {
               gamma1.P = rep(list(gamma1.P[[seg.iter]]), nS+1)
               gamma2.P = rep(list(gamma2.P[[seg.iter]]), nS+1)
            }
         }
         Residuals    <- read.table(file=paste0(dir.nS.ok,"/Results_Residuals_","it",seg.iter,"_","nS",nS,".csv"),sep=";",header=TRUE)
         mu.s         <- as.numeric(Residuals[,5])
         Results.seg  <- read.table(file=paste0(dir.nS.ok,"/Results_Summary_","it",seg.iter,"_","nS",nS,".csv"),sep=";",header=TRUE)
         mcmc.segment <- read.table(file=paste0(dir.nS.ok,"/Results_MCMC_cooked_","it",seg.iter,"_","nS",nS,".csv"),sep=";",header=TRUE)
         #initialisation of results:
         Q2.ts = NULL; Q97.ts = NULL; Q2.mu = NULL; Q97.mu = NULL;
         
         if ( nS > 1) { # if at least one changepoint:
            # change point times "tau":
            #*************************
            ts.res       <- as.numeric(c(Results.seg[16,(nS+1):(npar[nS]-1)]))  #the maxpost of all mcmc
            ts.mean      <- as.numeric(c(Results.seg[5,(nS+1):(npar[nS]-1)]))  #the maxpost of all mcmc
            ts.median    <- as.numeric(c(Results.seg[6,(nS+1):(npar[nS]-1)]))  #the maxpost of all mcmc
            ts.res.plus  <- as.numeric(c(Results.seg[16,(nS+1):(npar[nS]-1)], tail(Shift.Q$alpha_t,1)))
            ts.res.before <- as.numeric(c(Shift.Q$alpha_t[1],Results.seg[16,((nS+1):(npar[nS]-1))]))
            stdev.ts     <- as.numeric(c(Results.seg[11,(nS+1):(npar[nS]-1)]))
            for (i in 1:(nS-1)) {
               Q2.ts[i]  <- c(quantile(mcmc.segment[,(nS+i)], p = c(0.025)))
               Q97.ts[i] <- c(quantile(mcmc.segment[,(nS+i)], p = c(0.975)))
            }
            Q10.ts       <- as.numeric(c(Results.seg[7,(nS+1):(npar[nS]-1)]))
            Q90.ts       <- as.numeric(c(Results.seg[10,(nS+1):(npar[nS]-1)]))
            
            # segments mean "mu":
            #********************
            mu.res       <- as.numeric(Results.seg[16,1:nS])
            mu.mean      <- as.numeric(c(Results.seg[5,1:nS]))  #the maxpost of all mcmc
            mu.median    <- as.numeric(c(Results.seg[6,1:nS]))  #the maxpost of all mcmc
            for (j in 1:nS) {
               Q2.mu[j]  <- c(quantile(mcmc.segment[,j], p = c(0.025)))
               Q97.mu[j] <- c(quantile(mcmc.segment[,j], p = c(0.975)))
            }
            Q10.mu.res   <- as.numeric(Results.seg[7,1:(nS)])
            Q90.mu.res   <- as.numeric(Results.seg[10,1:(nS)])
            stdev.mu     <- as.numeric(Results.seg[11,1:(nS)])
            
            #structural error parameter "gamma":
            #**********************************
            gamma.MAP    <- as.numeric(Results.seg[16,npar[nS]])
            gamma.stdev  <- as.numeric(Results.seg[11,npar[nS]])
            gamma.mean   <- as.numeric(Results.seg[5,npar[nS]])
            final.period = c(final.period,1)
         
         ######################################################   
         } else { # if no change points ==> only one segment !!
         ######################################################
            ts.res      <- NULL; stdev.ts <- NULL; ts.res.before <- NULL; 
            ts.res.plus <- NULL; ts.mean <- NULL;
            ts.median   <- NULL; stdev.ts <- NULL; Q2.ts <- NULL; 
            Q97.ts      <- NULL; Q10.ts <- NULL; Q90.ts <-NULL;
            #-------------------------------------------------------
            mu.res       <- as.numeric(Results.seg[16,1])
            mu.mean      <- as.numeric(Results.seg[5,1])  #the maxpost of all mcmc         
            mu.median    <- as.numeric(Results.seg[6,1])  #the maxpost of all mcmc
            Q2.mu        <- quantile(mcmc.segment[,1], p = c(0.025))
            Q97.mu       <- quantile(mcmc.segment[,1], p = c(0.975))
            Q10.mu.res   <- as.numeric(Results.seg[7,1])
            Q90.mu.res   <- as.numeric(Results.seg[10,1])
            stdev.mu     <- as.numeric(Results.seg[11,1])
            #-------------------------------------------------------
            gamma.MAP    <- as.numeric(Results.seg[16,2])
            gamma.stdev  <- as.numeric(Results.seg[11,2])
            gamma.mean   <- as.numeric(Results.seg[5,2])
            final.period = c(final.period,0)
         }
         
         #saving to a data.frame:
         tau.results.df[[seg.iter]] =  data.frame(# change point times "tau":
                                       tau.MAP    = ts.res, 
                                       tau.q2     = Q2.ts, 
                                       tau.q10    = Q10.ts, 
                                       tau.q90    = Q90.ts, 
                                       tau.q97    = Q97.ts, 
                                       tau.stdev  = stdev.ts, 
                                       tau.mean   = ts.mean, 
                                       tau.median = ts.median)
         mu.results.df[[seg.iter]]  =  data.frame(# Segment mean "mu":
                                       mu.MAP     = mu.res, 
                                       mu.q2      = Q2.mu, 
                                       mu.q10     = Q10.mu.res,
                                       mu.q90     = Q90.mu.res, 
                                       mu.q97     = Q97.mu,
                                       mu.stdev   = stdev.mu, 
                                       mu.mean    = mu.mean,  
                                       mu.median  = mu.median)

         # saving to vectors:
         times.of.shift.MAP <- c(times.of.shift.MAP, ts.res)
         mean.of.segments   <- c(mean.of.segments,mu.res)
         t.q2   <- c(t.q2, Q2.ts)
         t.q10  <- c(t.q10, Q10.ts)
         t.q90  <- c(t.q90, Q90.ts)
         t.q97  <- c(t.q97, Q97.ts)
         t.q2   <- t.q2[c(TRUE, !t.q2[-length(t.q2)] == t.q2[-1])]
         t.q10  <- t.q10[c(TRUE, !t.q10[-length(t.q10)] == t.q10[-1])]
         t.q90  <- t.q90[c(TRUE, !t.q90[-length(t.q90)] == t.q90[-1])]
         t.q97  <- t.q97[c(TRUE, !t.q97[-length(t.q97)] == t.q97[-1])]
         times.of.shift.MAP <- times.of.shift.MAP[c(TRUE, !times.of.shift.MAP[-length(times.of.shift.MAP)] == times.of.shift.MAP[-1])]
         mean.of.segments   <-  mean.of.segments[c(TRUE, !mean.of.segments[-length( mean.of.segments)] ==  mean.of.segments[-1])]
         
         
         
         
         
         
         
         
         
         
         #******************************************************************************************
         #SHIFT TIME ADJUSTMENT: Plot stage record with segmentation results:
         #******************************************************************************************
         #Interval of the shift time:     
         interval = list(); 
         ts.real  = NULL; 
         hflood   = NULL; 
         tflood   = NULL;
         hflood2  = NULL; 
         tflood2  = NULL;
         
         if (nS >1) {  # Next analysis and plots only if at least a change point has been found !!!
         #################
            #Analysis of the stage record only if the stage record is uploaded!!!
            if (is.null(stage.record)==FALSE) {
               
               for (i in 1:length(tau.results.df[[seg.iter]]$tau.MAP)) {
                  interval[[i]] = which((t_limni >= min(tau.results.df[[seg.iter]]$tau.q2[i], 
                                                        tau.results.df[[seg.iter]]$tau.MAP[i])) &
                                        (t_limni <= max(tau.results.df[[seg.iter]]$tau.MAP[i], 
                                                        tau.results.df[[seg.iter]]$tau.q97[i])))
                  #maximum stage value (if known):
                  if (length(interval[[i]]) != 0){
                     hflood[i]     = max(h_limni[interval[[i]]])
                     tflood[i]     = t_limni[which.max(h_limni[interval[[i]]])+interval[[i]][1]]
                  } else {
                     interval[[i]] = 0
                     tflood[i]     = tau.results.df[[seg.iter]]$tau.MAP[i]
                  }
               }
               tau.results.df[[seg.iter]] = cbind(tau.results.df[[seg.iter]],
                                                  tflood = tflood)
            }
            CdT.P <- data.frame(tP, QP, hP)
            # plot segmentation results before times adjustment:
            initial.tsplot =  initial.ts.plot(CdT.P             = CdT.P, 
                                              stage.record      = stage.record, 
                                              tshift            = tau.results.df[[seg.iter]],
                                              limni.labels      = limni.labels, 
                                              grid_limni.ylim   = grid_limni.ylim,
                                              dir.seg.gaug      = dir.seg.gaug, 
                                              seg.iter          = seg.iter,
                                              t_Gaug            = t_Gaug,
                                              h_Gaug            = h_Gaug, 
                                              mcmc.segment      = mcmc.segment, 
                                              nS                = nS)
            X11()  # pupup plot of stage record with proposed segmentation
            print(initial.tsplot)
            
 
            
            
            
            
            
            #***********************************************************************
            # Identification of the TRUE shift times:
            #***********************************************************************
            # Looking inside the u95% interval of the shift times.
            #
            # Options:  1) if there is a flood, assign the shift time to the max peak !
            #           2) if the shift is due to other causes (vegetation, works, ice ...)
            #              then assign the shift time to the maxpost
            #           3) if the shift time is known then fix it.
            user.choice.ts=NULL
            i =0
            while (i < length(tau.results.df[[seg.iter]]$tau.MAP)) {
                i = i +1
                 if (shift.time.adjustment.type == 1) {    # option 1: always chose the MAP of the shift time !!!
                     ts.real[i] = tau.results.df[[seg.iter]]$tau.MAP[i] 
                 }  else if (shift.time.adjustment.type == 2) {
                     if (any(ts.real==tflood[i]) | any(ts.all.real==tflood[i])){
                         print("searching for a second flood in the interval")
                         hflood2[i] = hflood[i]
                         tflood2[i] = tflood[i] + 0.01   #to improve this in future !!!
                         # hflood2[i] = sort(h_limni[interval[[i]]], TRUE)[2]
                         # tflood2[i] = t_limni[which_nth_highest_vaalue(x=h_limni[interval[[i]]], n=2)[1]
                         #             + interval[[i]][1]]
                         ts.real[i] = tau.results.df[[seg.iter]]$tau.MAP[i]
                     } else {
                         hflood2[i] = hflood[i]
                         tflood2[i] = tflood[i]
                         ts.real[i] = tflood2[i]
                     }
                     # if (!is.null(stage.record)){
                     # if (i >1) {
                     #    if(((ts.real[i] -ts.real[i-1]) <5) & (ts.real[i] > ts.real[i-1])){
                     #       ts.real[i] = tflood2[i] + 10  # if in days !!
                     #       print(paste0("shift time too close from the previous one ==> delayed of 10 days!"))
                     #    }
                     # }
                     # }
                     ts.morpho.real = c(ts.morpho.real, ts.real[i])
                     ts.morpho.MAP  = c(ts.morpho.MAP, ts.res[i])
                     ts.morpho.q2   = c(ts.morpho.q2, Q2.ts[i])
                     ts.morpho.q97  = c(ts.morpho.q97, Q97.ts[i])
               }  else {
                  # Manual selection:
                  # pop-up window for user choice of the TRUE shift time option:
                  user.choice.ts[i] <- dlgInput(paste0("Which Adjusted shift time do you chose for ts",i," ? \n",
                                                       "1 = MAP shift time =  ", ts.res[i], " days     interval=[ ", Q2.ts[i], " ; ", Q97.ts[i], " ]  \n",
                                                       "2 = stage max (morphogenic flood at t_flood =", tflood[i], "   ;   with stage h_flood = ", hflood[i], " ) \n",
                                                       "3 = other time, e.g. earthquake, cyclon, ..."), Sys.info()[" "])$res
                  #choices for the adjustment:
                  if (user.choice.ts[i] == 1) {  # MAP
                     ts.real[i]      = ts.res[i]
                  } else if (user.choice.ts[i] == 2) {  # flood event
                     ts.real[i]      = tflood[i]
                     ts.morpho.real  = c(ts.morpho.real, ts.real[i])
                     ts.morpho.MAP   = c(ts.morpho.MAP, ts.res[i])
                     ts.morpho.q2    = c(ts.morpho.q2, Q2.ts[i])
                     ts.morpho.q97   = c(ts.morpho.q97, Q97.ts[i])
                  } else if (user.choice.ts[i] == 3) {  # other events
                     ts.real[i]      <- dlgInput("insert the date (in days) ", Sys.info()[" "])$res
                     ts.real[i]      = as.numeric(ts.real[i])
                     morpho.ask      <- dlgInput("is it related to sediment transport dynamics ? [Y/N]", Sys.info()[" "])$res
                     if (morpho.ask == "Y") {
                        ts.morpho.real = c(ts.morpho.real, ts.real[i])
                        ts.morpho.MAP  = c(ts.morpho.MAP, ts.res[i])
                        ts.morpho.q2   = c(ts.morpho.q2, Q2.ts[i])
                        ts.morpho.q97  = c(ts.morpho.q97, Q97.ts[i])
                     }
                  }    
                  # if (flood[i] < flood.threshold) {         #  ===> Posssibly yo improve with a sediment transport analysis !!!!!!!
                  #    ts.real[i] = ts.res[i] 
                  # } else {
                  #   ts.real[i] = t_limni[which(h_limni==flood[i],1)]
                  #   ts.morpho.real= c(ts.morpho.real, ts.real[i])
                  #   ts.morpho.MAP = c(ts.morpho.MAP, ts.res[i])
                  #   ts.morpho.q2 = c(ts.morpho.q2, Q2.ts[i])
                  #   ts.morpho.q97 = c(ts.morpho.q97, Q97.ts[i])
                  # }
               }
               
            }
            dev.off()  
            #########
            

            df.shift.times = data.frame(Q2.ts, Q97.ts, ts.res, ts.real)
            
            df.shift.times.plus = data.frame(ts.res.before, ts.res.plus,  
                                             Q2.mu, Q10.mu.res, Q90.mu.res, Q97.mu,  mu.res)
            write.table(mcmc.segment, paste0(dir.seg.gaug,"/mcmc_segmentation.txt"),
                        sep ="\t", row.names=FALSE)
            write.table(tau.results.df[[seg.iter]], paste0(dir.seg.gaug,"/df_tau_it", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)
            write.table(mu.results.df[[seg.iter]], paste0(dir.seg.gaug,"/df_mu_it", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)
            write.table(df.shift.times, paste0(dir.seg.gaug,"/df.shift.times_it", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)
            write.table(df.shift.times.plus, paste0(dir.seg.gaug,"/df.shift.times.plus_it", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)

         } else {
            write.table(mcmc.segment, paste0(dir.seg.gaug,"/mcmc_segmentation.txt"),
                        sep ="\t", row.names=FALSE)
            write.table(mu.results.df[[seg.iter]], paste0(dir.seg.gaug,"/df_mu_it", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)
         }
         
         
         
         #------------------------------------------------------------------------------------------
         # arrange detected shift times:
         if (!is.null(ts.real[1])) {
             ts.real        = as.numeric(ts.real)
             ts.all.real    = c(ts.all.real, ts.real)
             ts.all.real    = sort(ts.all.real)
             ts.all.real    = ts.all.real[c(TRUE, !ts.all.real[-length(ts.all.real)] == ts.all.real[-1])]
             ts.all.real.2  = c(ts.all.real.2, ts.real)
             ts.all.real.2  = sort(ts.all.real.2)
             tau.results.df[[seg.iter]] = cbind(tau.results.df[[seg.iter]], 
                                                treal = ts.real)         
             shift.results.df = rbind(shift.results.df ,
                                      tau.results.df[[seg.iter]])
         } else {
             tau.results.df[[seg.iter]] = FALSE;
         }



         
         
         
         
         
         
         
         
         #**************************************
         # Separation of gaugings per period
         #**************************************
         if (!is.null(ts.real[1])) {
            tss = 0
            for (j in 1:(nS-1)) {
               for (i in 1:(length(t_Gaug)-1)) {
                  if((  t_Gaug[i+1] >= ts.real[j]) & ((t_Gaug[i] <= ts.real[j] ))) {
                     tss[j] = i+1
                  }
               }
            }
            t.shift.plus <- tP[tss+1]
            t.shift.before <- tP[tss]
         } else {
            tss <- NULL
            t.shift.plus <- NULL
            t.shift.before <- NULL
         }
         
         
         
         
         
         #***********************************************
         #definition of the "stable" periods of gaugings:
         #***********************************************
         tss_tot    =  c(tss_tot_ns, tss)
         tss_tot    =  tss_tot[c(TRUE, !tss_tot[-length(tss_tot)] == tss_tot[-1])]
         tss_tot    =  sort(tss_tot)
         tss_tot_ns =  tss_tot
         i_init     =  0; 
         i_final    =  0;
         for (i in 1:(length(tss_tot)-1)) {
               i_init[i]  = tss_tot[i]
               i_final[i] = tss_tot[i+1]-1
         }   
         if (tss_tot[length(tss_tot)] == length(t_Gaug)) {
               i_init[length(tss_tot)]  = tss_tot[length(tss_tot)]
               i_final[length(tss_tot)] = tss_tot[length(tss_tot)]
         } else {
               i_init[length(tss_tot)]  = tss_tot[length(tss_tot)]
               i_final[length(tss_tot)] = length(t_Gaug)
         }
         
         
         
         
         
         #**********************************************************************
         #Plot residuals time series with segmentation:
         #**********************************************************************
         seg.plot[[seg.iter]] <- segm_P_plot(dir.seg.gaug     = dir.seg.gaug, 
                                             Shift.Q          = Shift.Q, 
                                             Q10.ts           = Q2.ts, 
                                             Q90.ts           = Q97.ts, 
                                             ts.res           = ts.res, 
                                             ts.res.before    = ts.res.before, 
                                             ts.res.plus      = ts.res.plus, 
                                             Q10.mu.res       = Q10.mu.res, 
                                             Q90.mu.res       = Q90.mu.res, 
                                             mu.res           = mu.res, 
                                             seg.iter         = seg.iter, 
                                             df.limni         = stage.record, 
                                             ts.real          = ts.real,
                                             resid.uncertaint = resid.uncertaint)

         
         
         
         
         #************************************
         # Plotting segmentated rating curves:
         #************************************
         color = 0; Period.P = 0;
         if (!is.null(ts.real[1])) {
            for (i in 1:length(Shift.Q$alpha)) {
               if(tP[i] <= ts.real[1]) {
                  #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
                  color[i]    = colors.period[1]
                  Period.P[i] = 1
               }
            }
            for (j in 2:nS) {
               for (i in 1:length(Shift.Q$alpha)) {
                  if ((tP[i] <= tail(Shift.Q$alpha_t,1)) & (tP[i] > ts.real[j-1])) {
                     #points(x=hP[i], y=QP[i], log ="y", col = colo[j],pch=1,lwd=4)
                     color[i]    = colors.period[j]
                     Period.P[i] = j
                  }
               }
            }
         } else {
            for (i in 1:length(Shift.Q$alpha)) {
               #points(x=hP[i], y=QP[i], log ="y", col = colo[1],pch=1,lwd=4)
               color[i]    = colors.period[1]
               Period.P[i] = 1
            }
         }
         df.RC <- data.frame(hP, QP, uQP, color, Period.P)
         CdT.P <- data.frame(hP, QP, uQP, Period.P, tP, color)
         #dev.off()
         gaug.plot[[seg.iter]] <- GaugingsShifts_P_plot(dir.seg.gaug   = dir.seg.gaug,
                                                        grid_RC.ylim   = grid_RC.ylim,
                                                        grid_RC.xlim   = grid_RC.xlim, 
                                                        grid_RC.xstep  = grid_RC.xstep, 
                                                        grid_RC.ystep  = grid_RC.ystep,
                                                        ticks_RC.y.log = ticks_RC.y.log, 
                                                        RC.x.labels    = RC.x.labels, 
                                                        RC.y.labels    = RC.y.labels,
                                                        gaugings       = gaugings, 
                                                        df.RC          = df.RC, 
                                                        seg.iter       = seg.iter)
  
         
         #*******************************************************************************
         # stage time series  with first set of shift times:
         #*******************************************************************************
         #dev.set(dev.prev())
         ts.ggplot[[seg.iter]] <- ts.plot(dir.segment.res   = dir.seg.gaug,
                                          CdT.P             = CdT.P, 
                                          all.gaugings      = gaugings,
                                          df.limni          = stage.record, 
                                          df.mu             = mu.results.df[[seg.iter]],
                                          df.shift.times    = tau.results.df[[seg.iter]],
                                          axis.labels       = limni.labels,
                                          grid_limni.ylim   = grid_limni.ylim,
                                          iteration         = seg.iter,
                                          mcmc.segment      = mcmc.segment,
                                          nS                = nS,
                                          known.shift.times = official.dates)
         if ( nS > 2) {
            for (nn in 1:length(ts.ggplot[[seg.iter]][[2]])) {
                 pdf.ts = cbind(pdf.ts, ts.ggplot[[seg.iter]][[2]][,nn])
            }
         } else if ( nS==2 ) {
                 pdf.ts = cbind(pdf.ts, ts.ggplot[[seg.iter]][[2]])
         }
         
         
         
         
         
         
         # #Figure with  example of shift time adjustment for paper Darienzo et al., 2021:
         #################################################################################
         # exemple_ts.plot(CdT.P, stage.record, Q10.ts=Q2.ts, Q90.ts=Q97.ts,
         #                 ts.res, ts.real, station.name, limni.labels,
         #                 grid_limni.ylim,  dir.seg.gaug, seg.iter,
         #                 t_Gaug, h_Gaug, pdf.ts = pdf.ts, time_period = c(1800, 2700), save_name = "Figure3")
         
         
         
         
         
         
         
         # Check autocorrelation of residuals:
         #####################################
         acf_resid[[seg.iter]]            = acf(x= Shift.Q$alpha,  plot = FALSE)
         pdf(paste0(dir.seg.gaug,"/Acf_it", seg.iter,".pdf"),  width = 10, height = 5, useDingbats=F)
         plot( acf_resid[[seg.iter]], main= paste0("Iteration ", seg.iter))
         dev.off()
         # partial autocorrelation:
         Pacf_resid[[seg.iter]]            = pacf(x= Shift.Q$alpha,  plot = FALSE)
         pdf(paste0(dir.seg.gaug,"/Partial Acf_it", seg.iter,".pdf"),  width = 10, height = 5, useDingbats=F)
         plot( Pacf_resid[[seg.iter]], main= paste0("Iteration ", seg.iter))
         dev.off()
         # with another R function:
         par(mar=c(1,1,1,1))
         autocorr_lag_resid[[seg.iter]] = autocorr.plot(x           = Shift.Q$alpha,  
                                                        auto.layout = FALSE)[[1]]
         dev.off()
         pdf(paste0(dir.seg.gaug,"/Autocorrelation_it", seg.iter,".pdf"),  width = 10, height = 5, useDingbats=F)
             autocorr.plot(x = Shift.Q$alpha,  auto.layout = FALSE)
             title(main = paste0("Iteration ", seg.iter), cex.main = 2)
         dev.off()

         

         
         # Save main results in the iteration folder:
         #############################################
         if (save.all.results ==TRUE) {
            write.table(df.gaug.P, paste0(dir.segment.gaug,"/it",seg.iter,"/BaRatin/df_gaug_it", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)
            write.table(df.gaug.tot, paste0(dir.segment.gaug,"/it",seg.iter,"/BaRatin/df_gaug_tot", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)
            
            write.table(criteria.df, paste0(dir.seg.gaug,"/BIC.df_it", seg.iter,".txt"),
                        sep ="\t", row.names=FALSE)

            write.table(Shift.Q, paste0(dir.seg.gaug,"/df_residuals_it", seg.iter,".txt"),
                     sep ="\t", row.names=FALSE)

            write.table(stage.record, paste0(dir.seg.gaug,"/df.limni.txt"),
                     sep ="\t", row.names=FALSE)
            write.table(CdT.P, paste0(dir.seg.gaug,"/df.gaug.periods.txt"),
                     sep ="\t", row.names=FALSE)
            write.table(gaugings, paste0(dir.seg.gaug,"/df.gaug.all.txt"),
                     sep ="\t", row.names=FALSE)

         }
         
         #***********************
         # end of segmentation:
         #**********************
         if ((nS ==1) & (i_final[seg.period] != (length(Q_Gaug))) & (recursive ==TRUE)){
            seg.period = seg.period +1
            seg.iter = seg.iter + 1
            # segmentation is finished for this current period!!! only one segment has been found in this period!
            # we pass to another period
         } else if ((nS != 1) & (recursive ==TRUE)) {
            seg.period = seg.period
            seg.iter = seg.iter + 1
         } else if ((nS ==1) & (i_final[seg.period] == (length(Q_Gaug)))) { #all periods have been completely segmentated ) {
            end.end = TRUE
            #final.period = 1
            #segmentation is finished. stop !!!!!
            print("segmentation finished!")
         } else if (recursive ==FALSE) {   #"single pass, only one iteration is done !!!
            
            end.end = TRUE
            print("segmentation finished!")
         }
      } else { # only one or two points in this period:
         print(paste0("No segmentation: there are only one or two points to segment for this period !!"))
         final.period = c(final.period,0)
         if (tail(tP,1) == (tail(t_Gaug,1))) {
            end.end = TRUE
            #segmentation is finished. stop !!!!!
            print("segmentation finished!")
         } else {
            seg.period = seg.period +1
            seg.iter = seg.iter + 1         
         }
      }
      
   }
            print("Segmentation ended correctly.")
   } else {
            print("You have selected the option for plotting results only!")
            print("No new segmentation is computed")
            print(paste0("reading folder ", dir.segment.gaug))
   }

   
   
   
   
   
   
   
   
   
   # FINAL RESULTS AND PLOTS:
   print("Preparing final results ...")
   #***********************************************************************************
   if (use.other.seg.method != FALSE) {
   #***********************************************************************************
      if (plot.results.only == FALSE) {
        
      tshift.interval.start = t_Gaug[cpd1-1]
      tshift.interval.end = t_Gaug[cpd1]
      tshift.interval.MAP = (t_Gaug[cpd1-1] + t_Gaug[cpd1])/2
     
      t_limni.interval =NULL; 
      h_limni.interval =NULL; 
      tshift.interval.adj =0
      if ((!is.null(stage.record$t_limni)) & (floodpeakalways == TRUE))  {
         for (interv in 1:length(tshift.interval.start)) {
               t_limni.interval[[interv]] = stage.record$t_limni[which((stage.record$t_limni >= tshift.interval.start[interv]) &
                                                                      (stage.record$t_limni <= tshift.interval.end[interv]))]
               h_limni.interval[[interv]] = stage.record$h_limni[t_limni.interval[[interv]]]
               tshift.interval.adj[interv] = t_limni.interval[[interv]][which.max(h_limni.interval[[interv]])]
         }
      } else {
         for (interv in 1:length(tshift.interval.start)) {
              tshift.interval.adj[interv] = tshift.interval.MAP[interv]
         }
      }
      shift.times.gaugings <- data.frame(t     = tshift.interval.MAP, 
                                         tadj  = tshift.interval.adj, 
                                         start = tshift.interval.start, 
                                         end   = tshift.interval.end, 
                                         cpt   = cpd1 ) 
      #shift.times.gaugings <- data.frame(t=shift.times.other.method, cpt= ch.point)
      pdf.ts =data.frame(V=runif(n=100000, min= tshift.interval.start[1],
                                            max= tshift.interval.end[1]))
      for (stg in 2:length(shift.times.gaugings$t)) {
           pdf.ts = cbind(pdf.ts,  V = runif(n=100000, min= tshift.interval.start[stg],
                                                       max= tshift.interval.end[stg]) )
      }
      write.table(shift.times.gaugings, paste0(dir.segment.gaug,"/shift_times.txt"), sep ="\t", row.names=FALSE)
      nS = length(shift.times.gaugings$cpt)+1

      c_Gaug = 0; P_Gaug = 0;
      if (!is.null(shift.times.gaugings$tadj)) {
         for (i in 1:length(t_Gaug)) {
            if(t_Gaug[i] <= shift.times.gaugings$tadj[1]) {
               #points(x=hP[i], y=QP[i], log ="y", col = colors.period[1],pch=1,lwd=4)
               c_Gaug[i] = colors.period[1]
               P_Gaug[i] = 1
            }
         }
         for (j in 2:nS) {
            for (i in 1:length(t_Gaug)) {
               if ((t_Gaug[i] <= tail(t_Gaug,1)) & (t_Gaug[i] > shift.times.gaugings$tadj[j-1])) {
                  #points(x=hP[i], y=QP[i], log ="y", col = colors.period[j],pch=1,lwd=4)
                  c_Gaug[i] = colors.period[j]
                  P_Gaug[i] = j
               }
            }
         }
      } else {
         for (i in 1:length(t_Gaug)) {
            #points(x=hP[i], y=QP[i], log ="y", col = colors.period[1],pch=1,lwd=4)
            c_Gaug[i] = colors.period[1]
            P_Gaug[i] = 1
         }
      }
      df.RC <- data.frame(h_Gaug, Q_Gaug, uQ_Gaug, P_Gaug, t_Gaug, c_Gaug)
      names(df.RC) = c("X.h.","X.Q.", "X.uQ.", "X.Period.", "X.tP.","color")
      write.table(df.RC, paste0(dir.segment.gaug,"/data_with_periods.txt"), sep ="\t", row.names=FALSE)
      write.table(pdf.ts, paste0(dir.segment.gaug,"/pdf_ts.txt"), sep ="\t", row.names=FALSE) 
      
      print("****************")
      print("   All done!    ")
      print("****************")
        return(list(df.RC, pdf.ts))
      } else {
         print("****************")
         print("   All done!    ")
         print("****************")
        return()
      }
      
      
   #************************************************************************************
   } else {
   #************************************************************************************   
      
      if (recursive ==FALSE) {
         if (plot.results.only == FALSE) {
         # Saving shift times in a file:
         times.of.shift.final <- c(times.of.shift.MAP, tail(Shift.Q$alpha_t,1))
         shift.times.gaugings <- data.frame(tMAP  = times.of.shift.MAP, 
                                            treal = ts.all.real, 
                                            t2    = t.q2,
                                            t10   = t.q10, 
                                            t90   = t.q90,
                                            t97   = t.q97)

         
         #write.table(shift.times.gaugings, paste0(dir.segment.gaug,"/shift_times.txt"), sep ="\t", row.names=FALSE)
         write.table(shift.times.gaugings, paste0(dir.segment.gaug,"/shift_times.txt"), sep ="\t", row.names=FALSE)         
         names(CdT.P) = c("X.h.","X.Q.", "X.uQ.", "X.Period.", "X.tP.","color")
         write.table(CdT.P, paste0(dir.segment.gaug,"/data_with_periods.txt"), sep ="\t", row.names=FALSE)
         write.table(pdf.ts, paste0(dir.segment.gaug,"/pdf_ts.txt"), sep ="\t", row.names=FALSE) 
         print("****************")
         print("   All done!    ")
         print("****************")
         return(list(CdT.P, pdf.ts))
         } else {
            print("****************")
            print("   All done!    ")
            print("****************")
           return()
         }
         
      } else {
         if (plot.results.only == FALSE) {
         print(c(paste0("Number of detected stable periods =",  seg.period), 
                 paste0("Number of iterations =",        seg.iter), 
                 paste0("Tot number of used gaugings =", i_final[seg.period])))
         # n <- 60
         #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
         #colo <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
         #colo <- c("grey","red","pink","green","blue","orange","black")
         env <- NULL
         env.par <- NULL
         RCmaxpost <- NULL
         gaug <- NULL
         period.index = 0
         maxpost.data <- NULL
         Data4SPD <- data.frame(NULL)
         for (i in 1:(length(final.period))) {
            if(final.period[i] == 0) {
               period.index    = period.index + 1
               dir.i          <- paste0(dir.segment.gaug,"/it",i,"/BaRatin")
               file.X         <- paste0(dir.i,"/Hgrid.txt")
               file.Y         <- paste0(dir.i,"/Qrc_TotalU.env")
               file.Ypar      <- paste0(dir.i,"/Qrc_ParamU.env")
               file.model     <- paste0(dir.i,"/Config_Model.txt")
               file.maxpost   <- paste0(dir.i,"/Qrc_maxpost.spag")
               file.summary   <- paste0(dir.i,"/Results_Summary.txt")
               gaug[[i]]      <- read.table(paste0(dir.i,"/data_with_time.txt"),fileEncoding="UTF-8", quote="",
                                            sep="", dec=".", header=TRUE);
               env[[i]]       <- Results_env.tot(file.Y, file.X, file.model, colors.period[i])
               env.par[[i]]   <- Results_env.par(file.Ypar, file.X, file.model, colors.period[i])
               RCmaxpost[[i]] <- data.frame(read.table(file = file.X), read.table(file = file.maxpost))
               maxpost        <- read.table(file = file.maxpost)
               maxpost.par    <- read.table(file = file.summary)
               Hg             <- read.table(file = file.X)
               gamma1         <- maxpost.par[16,ncontrols*3+1]; 
               gamma2         <- maxpost.par[16,ncontrols*3+2]
               k1             <- maxpost.par[16,1]; 
               str.err        =  0; 
               str.err        =  gamma1 + gamma2*maxpost
               maxpost.data[[i]] <- data.frame(x  = Hg[which(Hg >= k1),1], 
                                               y  = maxpost[which(Hg>= k1),1],
                                               y2 = str.err[which(Hg>= k1),1])
               gaug[[i]]$X.Period. = period.index
               gaug[[i]]$color = colors.period[period.index]
               Data4SPD = rbind( Data4SPD, gaug[[i]])
               
            } else {
               env[[i]]          <- NULL
               env.par[[i]]      <- NULL
               RCmaxpost[[i]]    <- NULL
               gaug[[i]]         <- NULL
               maxpost.data[[i]] <- NULL
            }
            # assign(paste("env_",i, sep = ""), env)
            # assign(paste("env.par_",i, sep = ""), env.par)
            # assign(paste("RCmaxpost_",i, sep = ""), RCmaxpost)
         }
         
         # Saving shift times in a file:
         times.of.shift.final <- c(times.of.shift.MAP, tail(Shift.Q$alpha_t,1))
         shift.times.gaugings <- data.frame(tMAP  = times.of.shift.MAP,
                                            treal = ts.all.real,
                                            t2    = t.q2, 
                                            t10   = t.q10, 
                                            t90   = t.q90, 
                                            t97   = t.q97)
         shift.results.df = shift.results.df[order(shift.results.df$treal),]
         shift.times.gaugings2 <- data.frame(tMAP  = shift.results.df$tau.MAP, 
                                             treal = shift.results.df$treal, 
                                             t2    = shift.results.df$tau.q2,
                                             t10   = shift.results.df$tau.q10, 
                                             t90   = shift.results.df$tau.q90,
                                             t97   = shift.results.df$tau.q97)
         
         write.table(shift.times.gaugings2, paste0(dir.segment.gaug,"/shift_times.txt"), sep ="\t", row.names=FALSE)
         morpho.shift.times <- data.frame(tMAP  = ts.morpho.MAP, 
                                          t2    = ts.morpho.q2,
                                          treal = ts.morpho.real, 
                                          t97   = ts.morpho.q97)
         # in the case some shift times are identified as due to morphogenic floods:
         write.table(morpho.shift.times, paste0(dir.segment.gaug,"/morpho_shift_times.txt"), sep ="\t", row.names=FALSE)
         nperiods = period.index
 
         
         
         
         #******************************************************************************
         # plot gaugings per periods:
         #******************************************************************************
         #colo = c("pink","orange","black", "green", "red", "blue")
         #colo = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
         # library(RColorBrewer)
         # n <- 60
         # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
         # colo = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
         print("Plotting RC with periods.")
         final.RC.plot(gaug, 
                       RCmaxpost, 
                       colorGauging  = colors.period, 
                       dir           = dir.segment.gaug,
                       grid_RC.ylim, 
                       grid_RC.xlim, 
                       grid_RC.xstep, 
                       grid_RC.ystep,
                       ticks_RC.y.log, 
                       RC.x.labels, 
                       RC.y.labels, 
                       grid_RC.ylim.log)
         
         
         #Figure for paper:
         #################
         # Figure.article(ts1= ts.ggplot[[1]], ts2= ts.ggplot[[2]], ts3= ts.ggplot[[3]],
         #                gp1 =gaug.plot[[1]], gp2= gaug.plot[[2]], gp3= gaug.plot[[3]]
         #                ) 
         # ggsave(plot.all, filename = paste0(dir.case_study,"/Results/segmentation_gaugings/it",seg.iter,"/all_plots.png"),
         #        bg = "transparent", width = 10, height =10, dpi = 400)
         # leg.fin <- legend.final( dir.case_study, Q10.ts, Q90.ts, ts.res, ts.res.before, ts.res.plus,
         #                          Hmin_grid, Hmax_grid, seg.iter, t_Gaug, h_Gaug,
         #                          Shift.Q, BIC.df, stage.record, CdT.P,
         #                          Q10.mu.res, Q90.mu.res, mu.res,
         #                          plot.Baratin[[6]] )
         
         
         
         #return the dates of shift found and segmentated gaugings
         nperiods.by.gaugings <- period.index
         #gaugings.per.period <- data.frame(h= h_Gaug, Q= Q_Gaug, uQ= uQ_Gaug, Period = period.index)
         write.table(Data4SPD[1:4], paste0(dir.segment.gaug,"/Gaugings_data_SPD.txt"),
                     sep="\t", row.names=FALSE , col.names = c("h","Q", "uQ", "Period"))  
         write.table(Data4SPD, paste0(dir.segment.gaug,"/data_with_periods.txt"), sep ="\t", row.names=FALSE)
         write.table(pdf.ts, paste0(dir.segment.gaug,"/pdf_ts.txt"), sep ="\t", row.names=FALSE) 
         
      }
         
         
      }
     
     
      # read and plot the results of the segmentation of gaugings:
      ###########################################################################################################
            g.s.results.1                    =  read.table(paste0(dir.segment.gaug,"/data_with_periods.txt"), sep ="\t", header =TRUE)
            nperiod.from.segm.gaugings.1     =  tail(g.s.results.1[[4]], 1)
            if (nperiod.from.segm.gaugings.1 >1){
                pdf.ts.results.1             =  read.table(paste0(dir.segment.gaug,"/pdf_ts.txt"), sep ="\t", header =TRUE) 
                pdf.ts                       =  pdf.ts.results.1
                shift.times.gaugings.1       =  read.table(file = paste0(dir.segment.gaug,"/shift_times.txt"), header = TRUE)
                data.annotate.gaug.1         =  data.frame(q2     =  shift.times.gaugings.1[,3],  
                                                           MAP    =  shift.times.gaugings.1[,1], 
                                                           q97    =  shift.times.gaugings.1[,6],
                                                           t.adj  =  shift.times.gaugings.1[,2])
            } else {
                data.annotate.gaug.1   = NULL
                shift.times.gaugings.1 = NULL
                pdf.ts.results.1       = NULL
                pdf.ts                 = NULL
            }

            Data4SPD =  g.s.results.1
            if (is.null(official.dates)==FALSE) {
               data.annotate.off = data.frame(xeffect = official.dates)
            } else {
               data.annotate.off = NULL
            }
            

            
            
            
            
            ######################################
            if (limni.labels[1] != "Time [year]"){   # transform times data from numeric to dates:
            ######################################
               
               if (file_limni != FALSE) {
                     init.date = stage.record$t_limni.true[1] 
                     # transform stage record times:
                     limni.date             = as.Date(floor(stage.record$t_limni.true), origin = date_origin)
                     limni.date             = as.POSIXct(as.POSIXlt(limni.date))
                     limni.date             = limni.date  + (stage.record$t_limni.true- floor(stage.record$t_limni.true)) * 24 * 3600
                     stage.record$date      = limni.date
                     
               } else if (file_limni == FALSE) {
                     init.date = gaugings$t.true[1]      
               }
               # transform gaugings times:
               g.s.results.1$date     = as.Date(floor(g.s.results.1$X.tP. + init.date), origin = date_origin)
               g.s.results.1$date     = as.POSIXct(as.POSIXlt(g.s.results.1$date))
               g.s.results.1$date     = g.s.results.1$date  + (g.s.results.1$X.tP. + init.date  - floor(g.s.results.1$X.tP. + init.date)) * 24 * 3600

                # transform shift times:
                if (!is.null(data.annotate.gaug.1)){
                   dates.1 = 0
                   for (dat in 1:length(data.annotate.gaug.1$t.adj)) {
                      # dates.1[dat]     = data.annotate.gaug.adjust.1$t.adj[dat]  +   Gaugings$X.3[1]
                      dates.1[dat]       = data.annotate.gaug.1$t.adj[dat] + init.date #Gaugings[1,ncol(Gaugings)]
                   }
                   RealDates.1            = as.Date(dates.1,  origin = date_origin)
                   RealDates.1.new        = strptime(as.character(RealDates.1), "%Y-%m-%d" )
                   shift.times.date       = format( RealDates.1.new, "%d/%m/%Y")
                } else {
                   shift.times.date = NULL
                }
                # pdf.shift.times.date = pdf.ts.results.1 
                # dates.pdf = pdf.ts.results.1
                # for (dat in 1:length(pdf.ts.results.1)) {
                #   dates.pdf[,dat] =  pdf.ts.results.1[,dat] + init.date #Gaugings[1,ncol(Gaugings)]
                #   dates.pdf[,dat] =  as.Date(dates.pdf[,dat],  origin = date_origin)
                #   pdf.shift.times.date[,dat] = strptime(as.character(dates.pdf[,dat]), "%Y-%m-%d" )
                #   pdf.shift.times.date[,dat] = format( dates.pdf[,dat], "%d/%m/%Y")
                # }
                
            } else {
              shift.times.date = NULL
              #pdf.shift.times.date = NULL
            }
            
            
            
            print("Plotting stage record with periods.")
            # plot of the time series with detected shift times.
            if(plot_dates == TRUE){
               shiftdates_limni = shift.times.date
            } else {
               shiftdates_limni = NULL
            }
            plot.time.shifts.gaugings(dir                       = dir.segment.gaug, 
                                      g.s                       = g.s.results.1, 
                                      data.annotate.off         = data.annotate.off, 
                                      data.annotate.gaug        = data.annotate.gaug.1,
                                      t_Gaug                    = t_Gaug, 
                                      h_Gaug                    = h_Gaug, 
                                      df.limni                  = stage.record, 
                                      limni.labels              = limni.labels,
                                      grid_limni.ylim           = grid_limni.ylim,
                                      plot.shift.times.on.limni = FALSE,
                                      pdf.ts                    = pdf.ts.results.1,
                                      dates                     = shiftdates_limni) 
          
            # plot the gaugings with different colors per period:
            print("Plotting RC with periods.")
            final.RC.plot.2(gaug              = g.s.results.1, 
                            dir               = dir.segment.gaug,
                            grid_RC.ylim      = grid_RC.ylim, 
                            grid_RC.xlim      = grid_RC.xlim, 
                            grid_RC.xstep     = grid_RC.xstep,
                            grid_RC.ystep     = grid_RC.ystep,
                            ticks_RC.y.log    = ticks_RC.y.log,
                            RC.x.labels       = RC.x.labels, 
                            RC.y.labels       = RC.y.labels, 
                            grid_RC.ylim.log  = grid_RC.ylim.log,
                            colors.period     = colors.period)
            
            # #plot autocorrelation for all iterations in unique plot (check the independancy of residuals):
            ##############################################################################################
            # name.iterat = c("0", 
            #                 "1.1",
            #                 "2.1.1", "2.1.2",    
            #                 "1.2", "2.2.1", "2.2.2",  
            #                 "1.3", "2.3.1", "2.3.2",
            #                 "3.3.2.1", "3.3.2.2", "3.3.2.3")
            # #ci.type = "ma"
            # 
            # autocorr.lag1 = c()
            # pdf(paste0(dir.segment.gaug,"/autocorrelation_all_iterations.pdf"), width = 8, height = length(acf_resid), useDingbats=F)
            # par(mfrow = c(ceil(length(acf_resid)/2), 2))   #c(nrows, ncol)
            # for (acf_i in 1:(length(acf_resid))) {
            #   par(mar = c(3, 3, 3, 3))  # Set the margin on all sides 
            #   autocorr.lag1 = c(autocorr.lag1, acf_resid[[acf_i]]$acf[2]) 
            #   if (!is.null(acf_resid[[acf_i]])){
            #       plot(acf_resid[[acf_i]], main = paste0("Iteration ", name.iterat[acf_i],": Acf lag1 =", round(acf_resid[[acf_i]]$acf[2], digits=3)), cex.main = 4)
            #   } else {
            #       plot(0, ylab = '', xlab = "lag", pch ='',
            #            main = paste0("Iteration ", name.iterat[acf_i]), cex.main = 1.3) #xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
            #       text(0, "Less than three residuals",   cex = 2, col="gray") 
            #   }
            # }
            # dev.off()
            # 
            # pdf(paste0(dir.segment.gaug,"/Partial_autocorrelation_all_iterations.pdf"), width = 8, height = length(acf_resid), useDingbats=F)
            # par(mfrow = c(ceil(length(acf_resid)/2), 2))   #c(nrows, ncol)
            # for (acf_i in 1:(length(Pacf_resid))) {
            #   par(mar = c(3, 3, 3, 3))  # Set the margin on all sides 
            #   if (!is.null(Pacf_resid[[acf_i]])){
            #     plot(Pacf_resid[[acf_i]],  main = paste0("Iteration ", name.iterat[acf_i]), cex.main = 4)
            #   } else {
            #     plot(0, ylab = '', xlab = "lag", pch ='',
            #          main = paste0("Iteration ", name.iterat[acf_i]), cex.main = 1.3) #xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
            #     text(0, "Less than three residuals",   cex = 2, col="gray") 
            #   }
            # }
            # dev.off()
            
            
            
            
            
            # plot Figure for article Darienzo et al., 2020 with the results of 3 crucial iterations!!!
            ############################################################################################
            # plot.figure.article(dir               =  dir.segment.gaug,
            #                     iteration.indexes =  c(1,5,8),   #indexes of 3 selected iterations
            #                     df.limni          =  df.limni,
            #                     grid_RC.ylim      =  grid_RC.ylim,
            #                     grid_RC.xlim      =  c(-0.8, 2.2),
            #                     grid_RC.xstep     =  0.4,
            #                     grid_RC.ystep     =  grid_RC.ystep,
            #                     grid_RC.ylim.log  =  c(0.08, 500),
            #                     ticks_RC.y.log    =  ticks_RC.y.log,
            #                     RC.x.labels       =  RC.x.labels,
            #                     RC.y.labels       =  RC.y.labels,
            #                     u.m.Qgaug         =  u.m.Qgaug,
            #                     u.m.Hgaug         =  u.m.Hgaug,
            #                     ncontrol          =  ncontrols,
            #                     plot.Q.log        =  TRUE)

            
            # plot a figure for each iteration, with RC, criteria, residuals:     NOT DONE YET !!!
            ######################################################################################
            # plot.all.iterations(dir               =  dir.segment.gaug,
            #                     df.limni          =  df.limni,
            #                     grid_RC.ylim      =  grid_RC.ylim,  
            #                     grid_RC.xlim      =  c(-0.8, 2.2), 
            #                     grid_RC.xstep     =  0.4, 
            #                     grid_RC.ystep     =  grid_RC.ystep,
            #                     grid_RC.ylim.log  =  c(0.08, 500),
            #                     ticks_RC.y.log    =  ticks_RC.y.log, 
            #                     RC.x.labels       =  RC.x.labels,
            #                     RC.y.labels       =  RC.y.labels, 
            #                     u.m.Qgaug         =  u.m.Qgaug, 
            #                     u.m.Hgaug         =  u.m.Hgaug,
            #                     ncontrol          =  ncontrols,
            #                     plot.Q.log        = TRUE)  
            
            print("****************")
            print("   All done!    ")
            print("****************")
         #**********************************
         return(list(Data4SPD, pdf.ts))
         #**********************************
         # FINISH !!!
      }
}











                                                                                                                                                                                                                                                                                                                                                                                            


























#############################################################################################################
log.posterior <- function(param, x.data, y.data, uy.data, nS) {
#############################################################################################################
   pred = 0; singlelikelihoods=0; sd_lin = 0; sd_remn=0; eps=0; single.prior.mu=0; single.prior.tau=0;
   npar= length(param)
   n.mu = npar/2
   N = length(x.data)
   if (nS>1) {
      n.tau = n.mu -1
      tau =  param[(n.mu+1):(npar-1)]
   }
   mu = param[1:n.mu]
   for (i in 1:N) { # N data
      if (nS==1) {
         pred[i] = mu[1]
      } else {
         for (j in 1:n.tau) {   
            if (x.data[i] < tau[j]) {
               pred[i] = j
               
               
            }
         }
         if (x.data[i] > tau[n.tau]) {
            pred[i] = nS
         }
      }
      # Calculation of the St Dev of the Remnant error using a linear relantionship
      sd_remn[i] = param[npar]
      eps[i]= (sd_remn[i]^2 + (uy.data[i])^2)^0.5
      # Calculation of the likelihood for the each coupled data h,Q
      singlelikelihoods = dnorm(y.data[i], mean = pred[i], sd = eps[i], log = TRUE)
   }
   print(pred)
   LH = sum(singlelikelihoods)
   print(LH)
   #*****************************
   for (j in 1:n.mu) {
      single.prior.mu[j] = dunif(x=1, min = -1e+30, max = 1e+30, log = TRUE)
   }
   if (nS>1) {
      for (j in 1:n.tau) {
         single.prior.tau[j] = dunif(x = 1, min = 0, max = 1e+30, log = TRUE)
      }
   }
   single.prior.gamma = dunif(x = 1, min = -1e+30, max = 1e+30, log = TRUE)
   apriori = sum(single.prior.mu) + sum(single.prior.tau)  +single.prior.gamma
   print(apriori)
   #*****************************
   Log.poster = LH + apriori
   print(Log.poster)
   return(Log.poster)  
}















################################################
start.values.time = function(nS, tfirst, tfin) {
################################################  
    tstart = 0
    for (i in 1:(nS-1)) {
      tstart[i] = (tfirst+(tfin-tfirst)/nS*i)
    }
    REG.INTERV = rep(FALSE, nS-1)
    #while(any(REG.INTERV) ==FALSE) {
    if (nS >2) {
      for (ii in 1:(N-2)) {
        for (i in 1:(length(tstart)-1)) {
          if ((tP[ii] <= tstart[i]) &  (tP[ii+1] >= tstart[i+1])) {
            tstart[i+1]= tP[ii+1] + 0.001
            
          } else {
            REG.INTERV[i] = TRUE
          }
        }
      }
      if (tP[N-1] <= tstart[length(tstart)-1]) {
        REG.INTERV[length(tstart)-1] = FALSE
        tstart[length(tstart)-1]= tP[N-1] - 0.001
      }
    }
  return(tstart)
}
  





















###########################################################################################
Segmentation_app <- function(dir_code,
                             nobs, 
                             nS, tmin, Nmin, tfirst, tfin, 
                             ncycles, nmcmc, Nslim, 
                             tP, 
                             resid.uncertaint, 
                             prior.mu,
                             prior.gamma.segm) {
###########################################################################################
   Segmentation_config(dir_code, nobs, nS, tmin, Nmin, tfirst, tfin, ncycles, 
                       nmcmc, Nslim, tP, resid.uncertaint, prior.mu, prior.gamma.segm)
   message("Applying Segmentation - BaM !!!  Wait ... "); flush.console()
   # system2(paste(dir_code,"/BaM_exe/BaM_Segmentation.exe",sep=""),stdout =NULL, stderr = NULL);
   system2(paste(dir_code,"/BaM_exe/BaM_Segmentation2.exe",sep=""),stdout =NULL, stderr = NULL); 
}















####################################################################################################
Segmentation_config <- function(dir_code, 
                                nobs, nS, tmin, Nmin, tfirst, tfin, 
                                ncycles, nmcmc, Nslim, 
                                tP, 
                                resid.uncertaint,
                                prior.mu,
                                prior.gamma.segm) {
####################################################################################################
   npar = nS + nS - 1
   N = length(tP)
   file.name1 = paste(dir_code,"/BaM_exe/Segmentation/Config_Model.txt",sep="")
   cat('"Segmentation2"', file =file.name1,sep="\n")
   cat(1, file = file.name1, append = TRUE,sep="\n")
   cat(1, file =file.name1, append = TRUE,sep="\n")
   cat(npar, file =file.name1, append = TRUE,sep="\n")
   for (i in 1:nS) {
      cat(paste0('"k',i,'"'), file =file.name1, append = TRUE,sep="\n")
      cat(prior.mu[1], file =file.name1, append = TRUE,sep="\n")
      cat('"Gaussian"', file = file.name1, append = TRUE, sep="\n")            #! Prior distribution
      cat(prior.mu[1], file =file.name1, append = TRUE, sep=",")
      cat(",",file =file.name1, append = TRUE, sep=",")
      cat(prior.mu[2], file =file.name1, append = TRUE, sep="\n")
      # cat("'FlatPrior'", file =file.name1, append = TRUE,sep="\n")
      # cat(" ", file =file.name1, append = TRUE,sep="\n")
   }
   #********************************************
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
         cat(tstart[i], file =file.name1, append = TRUE,sep="\n")
         cat("'FlatPrior+'", file =file.name1, append = TRUE,sep="\n")
         cat(" ", file =file.name1, append = TRUE,sep="\n")
         # cat('"Gaussian"', file = file.name1, append = TRUE, sep="\n")            #! Prior distribution
         # cat((tfirst+(tfin-tfirst)/nS*i), file =file.name1, append = TRUE, sep=",")
         # cat(",",file =file.name1, append = TRUE, sep=",")
         # cat((tfin-tfirst)/4, file =file.name1, append = TRUE, sep="\n")
      }
   }
   #####################################################
   file.name2 = paste(dir_code,"/BaM_exe/Segmentation/Config_Data.txt",sep="")
   cat("'Segmentation\\Segm_data.txt'", file =file.name2,sep="\n")
   cat(1, file = file.name2, append = TRUE,sep="\n")
   cat(nobs, file = file.name2, append = TRUE,sep="\n") 
   cat(3, file =file.name2, append = TRUE,sep="\n")
   cat(1, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(2, file =file.name2, append = TRUE,sep="\n")
   cat(3, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   #####################################################
   file_BaM <- paste(dir_code,"/BaM_exe/Config_BaM.txt",sep="")
   #creation of Config_BaM.txt
   cat('"Segmentation/"', file =file_BaM , sep="\n", append = FALSE)
   cat('"Config_RunOptions.txt"', file = file_BaM , sep="\n", append = TRUE)   
   cat('"Config_Model.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_xtra.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_Data.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_RemnantSigma.txt"', file = file_BaM , sep="\n", append = TRUE)                                      
   cat('"Config_MCMC.txt"', file = file_BaM , sep="\n", append = TRUE)                                           
   cat('"Config_Cooking.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_Summary.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_Residuals.txt"', file = file_BaM , sep="\n", append = TRUE)
   #cat('"Config_Pred_Master.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('" "', file = file_BaM , sep="\n", append = TRUE)
   #####################################################
   # file extra:
   # file.xtra = paste(dir_code,"/BaM_exe/Segmentation/Config_xtra.txt",sep="")
   # cat(nS, file =file.xtra,sep="\n")
   # cat(tmin, file = file.xtra, append = TRUE,sep="\n")
   file.xtra = paste(dir_code,"/BaM_exe/Segmentation/Config_xtra.txt",sep="")
   cat(nS, file =file.xtra,sep="\n")
   cat(tmin, file = file.xtra, append = TRUE, sep="\n") # tmin
   cat(Nmin, file = file.xtra, append = TRUE, sep="\n")   # Nmin
   cat(1, file = file.xtra, append = TRUE, sep="\n")   # option Nmin (1) or tmin(2)
   ###################################################################   RUN OPTIONS
   file.run = paste(dir_code,"/BaM_exe/Segmentation/Config_RunOptions.txt",sep="")
   cat(".true.", file =file.run,sep="\n")                    #Do MCMC?
   cat(".true.", file =file.run, append = TRUE, sep="\n")    #Do MCMC summary?
   cat(".true.", file =file.run, append = TRUE, sep="\n")    #Do Residual diagnostics?
   cat(".false.", file =file.run, append = TRUE, sep="\n")    #Do Predictions?
   ###################################################################   RESIDUALS CONFIG
   file.residuals = paste(dir_code,"/BaM_exe/Segmentation/Config_Residuals.txt",sep="")
   cat("Results_Residuals.txt" , file =file.residuals ,sep="\n")    #Result file
   ###################################################################   SUMMARY CONFIG
   file.summary = paste(dir_code,"/BaM_exe/Segmentation/Config_Summary.txt",sep="")
   cat("Results_Summary.txt" , file =file.summary ,sep="\n")    #Result file
   ###################################################################   COOKING CONFIG
   file.cooking = paste(dir_code,"/BaM_exe/Segmentation/Config_Cooking.txt",sep="")
   cat("Results_MCMC_Cooked.txt" , file =file.cooking ,sep="\n")    #Result file
   cat(0.5, file =file.cooking, append = TRUE, sep="\n")            #Burn factor
   cat(Nslim, file =file.cooking, append = TRUE, sep="\n")             #Nslim
   ###################################################################   REMNANT ERROR CONFIG
   file.remnant = paste0(dir_code,"/BaM_exe/Segmentation/Config_RemnantSigma.txt")
   cat("'Constant'", file = file.remnant, sep="\n")                       #! Function f used in sdev=f(Qrc)
   cat(1, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
   cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
   if (resid.uncertaint==FALSE) {
      cat(0.1,                      file = file.remnant, append = TRUE, sep="\n")    #! Initial Guess
      cat('Uniform',                file = file.remnant, append = TRUE, sep="\n")    #! Prior distribution
      cat(prior.gamma.segm[1],      file = file.remnant, append = TRUE, sep=",")
      cat(",",                      file = file.remnant, append = TRUE, sep=",")
      cat(prior.gamma.segm[2],      file = file.remnant, append = TRUE, sep="\n")
   } else {
      cat(exp(prior.gamma.segm[1]), file = file.remnant, append = TRUE, sep="\n")    #! Initial Guess
      cat('"LogNormal"',            file = file.remnant, append = TRUE, sep="\n")   #! Prior distribution
      cat(prior.gamma.segm[1],      file = file.remnant, append = TRUE, sep=",")
      cat(",",                      file = file.remnant, append = TRUE, sep=",")
      cat(prior.gamma.segm[2],      file = file.remnant, append = TRUE, sep="\n")
      # it should be fixed to zero, but option not available in BaM!
   }
   ###################################################################   MCMC
   file.mcmc = paste(dir_code,"/BaM_exe/Segmentation/Config_MCMC.txt",sep="")
   cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
   cat(nmcmc,    file = file.mcmc, append = TRUE,sep="\n")    # Nadapt
   cat(ncycles,  file = file.mcmc, append = TRUE,sep="\n")    # Ncycles
   cat(0.1,      file = file.mcmc, append = TRUE,sep="\n")    # minMoveRate
   cat(0.5,      file = file.mcmc, append = TRUE,sep="\n")    # maxMoveRate
   cat(0.9,      file = file.mcmc, append = TRUE,sep="\n")    # DownMult
   cat(1.1,      file = file.mcmc, append = TRUE,sep="\n")    # UpMult
   cat(0,        file = file.mcmc, append = TRUE,sep="\n")    # mode for init jump distr
   cat("****",   file = file.mcmc, append = TRUE,sep="\n")
   cat(0.1,      file = file.mcmc, append = TRUE,sep="\n")    # MultFact
   cat(0.1,      file = file.mcmc, append = TRUE, sep=",")    # RC MultiFact
   cat(0.1,      file = file.mcmc, append = TRUE, sep=",")    
   cat(0.1,      file = file.mcmc, append = TRUE,sep="\n")
   cat(0.1,      file = file.mcmc, append = TRUE, sep=",")    # Remnant MultiFact
   cat(0.1,      file = file.mcmc, append = TRUE,sep="\n")
}







   
   
   
   
   
   


##########################################################################
which_nth_highest_vaalue <- function(x, n) {
##########################################################################
   ux <- unique(x)
   nux <- length(ux)
   which(x == sort(ux, partial = nux - n + 1)[nux - n + 1])
} 






   









#************************************************************************************
Segmentation_app_1MCMC <- function(nobs,nS, tmin, tfirst, tfin,  theta_mean) {
#************************************************************************************
   Segmentation_config_1MCMC(nobs,nS,tmin,tfirst, tfin,  theta_mean)
   #message("Applying Segmentation - BaM !!!  Wait ... "); flush.console()
   system2(paste(dir_code,"/BaM_exe/BaM_Segmentation2.exe",sep=""),stdout =NULL, stderr = NULL);
}










#*******************************************************************************
Segmentation_config_1MCMC <- function(nobs,nS,tmin,tfirst,tfin, theta_mean) {
   #*******************************************************************************
   npar = 2*nS - 1
   file.name1 = paste(dir_code,"/BaM_exe/Segmentation/Config_Model.txt",sep="")
   cat('"Segmentation2"', file =file.name1,sep="\n")
   cat(1, file = file.name1, append = TRUE,sep="\n")
   cat(1, file =file.name1, append = TRUE,sep="\n")
   cat(npar, file =file.name1, append = TRUE,sep="\n")
   
   for (i in 1:nS) {
      cat(paste('"k',i,'"',sep=""), file =file.name1, append = TRUE,sep="\n")
      cat(theta_mean[i], file =file.name1, append = TRUE,sep="\n")
      cat("'FlatPrior'", file =file.name1, append = TRUE,sep="\n")
      cat(" ", file =file.name1, append = TRUE,sep="\n")
   }
   if (nS>1) {
      for (i in 1:(nS-1)) {
         cat(paste('"tau',i,'"',sep=""), file =file.name1, append = TRUE,sep="\n")
         cat(theta_mean[i+nS] , file =file.name1, append = TRUE,sep="\n")
         cat("'FlatPrior+'", file =file.name1, append = TRUE,sep="\n")
         cat(" ", file =file.name1, append = TRUE,sep="\n")
      }
   }
   #####################################################
   file.name2 = paste(dir_code,"/BaM_exe/Segmentation/Config_Data.txt",sep="")
   cat("'Segmentation\\Segm_data.txt'", file =file.name2,sep="\n")
   cat(1, file = file.name2, append = TRUE,sep="\n")
   cat(nobs, file = file.name2, append = TRUE,sep="\n") 
   cat(3, file =file.name2, append = TRUE,sep="\n")
   cat(1, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(2, file =file.name2, append = TRUE,sep="\n")
   cat(3, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   cat(0, file =file.name2, append = TRUE,sep="\n")
   #####################################################
   file_BaM <- paste(dir_code,"/BaM_exe/Config_BaM.txt",sep="")
   #creation of Config_BaM.txt
   cat('"Segmentation/"', file =file_BaM , sep="\n", append = FALSE)
   cat('"Config_RunOptions.txt"', file = file_BaM , sep="\n", append = TRUE)   
   cat('"Config_Model.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_xtra.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_Data.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_RemnantSigma.txt"', file = file_BaM , sep="\n", append = TRUE)                                      
   cat('"Config_MCMC.txt"', file = file_BaM , sep="\n", append = TRUE)                                           
   cat('"Config_Cooking.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_Summary.txt"', file = file_BaM , sep="\n", append = TRUE)
   cat('"Config_Residuals.txt"', file = file_BaM , sep="\n", append = TRUE)
   #cat('"Config_Pred_Master.txt"', file = file_BaM , sep="\n", append = TRUE)
   #####################################################
   cat('" "', file = file_BaM , sep="\n", append = TRUE)
   file.xtra = paste(dir_code,"/BaM_exe/Segmentation/Config_xtra.txt",sep="")
   cat(nS, file =file.xtra,sep="\n")
   cat(tmin, file = file.xtra, append = TRUE,sep="\n")
   ###################################################################   RUN OPTIONS
   file.run = paste(dir_code,"/BaM_exe/Segmentation/Config_RunOptions.txt",sep="")
   cat(".true.", file =file.run,sep="\n")                    #Do MCMC?
   cat(".true.", file =file.run, append = TRUE, sep="\n")    #Do MCMC summary?
   cat(".true.", file =file.run, append = TRUE, sep="\n")    #Do Residual diagnostics?
   cat(".false.", file =file.run, append = TRUE, sep="\n")    #Do Predictions?
   ###################################################################   RESIDUALS CONFIG
   file.residuals = paste(dir_code,"/BaM_exe/Segmentation/Config_Residuals.txt",sep="")
   cat("Results_Residuals.txt" , file =file.residuals ,sep="\n")    #Result file
   ###################################################################   SUMMARY CONFIG
   file.summary = paste(dir_code,"/BaM_exe/Segmentation/Config_Summary.txt",sep="")
   cat("Results_Summary.txt" , file =file.summary ,sep="\n")    #Result file
   ###################################################################   COOKING CONFIG
   file.cooking = paste(dir_code,"/BaM_exe/Segmentation/Config_Cooking.txt",sep="")
   cat("Results_MCMC_Cooked.txt" , file =file.cooking ,sep="\n")    #Result file
   cat(0.5, file =file.cooking, append = TRUE, sep="\n")            #Burn factor
   cat(10, file =file.cooking, append = TRUE, sep="\n")             #Nslim
   ###################################################################   REMNANT ERROR CONFIG
   file.remnant = paste(dir_code,"/BaM_exe/Segmentation/Config_RemnantSigma.txt",sep="")
   cat("'Constant'", file = file.remnant, sep="\n")                       #! Function f used in sdev=f(Qrc)
   cat(1, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
   cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
   cat(theta_mean[2*nS], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
   cat('Uniform', file = file.remnant, append = TRUE, sep="\n")            #! Prior distribution
   cat(0,file =file.remnant, append = TRUE, sep=",")
   cat(",",file =file.remnant, append = TRUE, sep=",")
   cat(1000,file =file.remnant, append = TRUE, sep="\n")
   ###################################################################   MCMC
   file.mcmc = paste(dir_code,"/BaM_exe/Segmentation/Config_MCMC.txt",sep="")
   cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
   cat(1, file = file.mcmc, append = TRUE,sep="\n")   #Nadapt
   cat(1, file = file.mcmc, append = TRUE,sep="\n")  #Ncycles
   cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #minMoveRate
   cat(0.5, file =file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
   cat(0.9, file =file.mcmc, append = TRUE,sep="\n")    #DownMult
   cat(1.1, file =file.mcmc, append = TRUE,sep="\n")    #UpMult
   cat(0, file =file.mcmc, append = TRUE,sep="\n")      #mode for init jump distr
   cat("****", file =file.mcmc, append = TRUE,sep="\n")
   cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #MultFact
   cat(0.1,file =file.mcmc, append = TRUE, sep=",")     #RC MultiFact
   cat(0.1,file =file.mcmc, append = TRUE, sep=",")    
   cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
   cat(0.1,file =file.mcmc, append = TRUE, sep=",")      #Remnant MultiFact
   cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
}






























#****************************************************************
alpha_segm <- function(t_Gaug, h_Gaug, Qt, time, Q_Gaug, h_CM, uh) {
   #****************************************************************
   alpha = 0
   alpha_t = 0
   indexx =0
   # for (i in 1:length(Q_Gaug)) {
   #    alpha[i] = (h_Gaug[i] - h_CM[i])/(uh[i])
   # }
   for (i in 1:length(Q_Gaug)) {
      #if(h_Gaug[i] <= 0.3) {
      alpha[i] = (h_Gaug[i] - h_CM[i])#/uh[i]
      alpha_t[i] = t_Gaug[i]
      indexx[i] = i
   }
   # else {
   #     alpha[i] = NA
   #     alpha_t[i] =  NA
   #     indexx[i] = NA
   #    }
   # }
   alpha <- alpha[!is.na(alpha)]
   alpha_t <- alpha_t[!is.na(alpha_t)]
   indexx <- indexx[!is.na(indexx)]
   Shift <- data.frame(alpha_t, alpha, uh, indexx)
   return(Shift)
}




























#****************************************************************
alpha_segm.Q <- function(t_Gaug, Q_Gaug, resid, sigma.tot) {
   #**************************************************************** 
   alpha = 0
   alpha_t = 0
   indexx =0
   for (i in 1:length(Q_Gaug)) {
      #if(h_Gaug[i] <= 0.3) {
      alpha[i] = resid[i]#/(resid[i]+Q_Gaug[i])
      alpha_t[i] = t_Gaug[i]
      indexx[i] = i
   }
   alpha <- alpha[!is.na(alpha)]
   alpha_t <- alpha_t[!is.na(alpha_t)]
   indexx <- indexx[!is.na(indexx)]
   Shift.Q <- data.frame(alpha_t, Q_Gaug, alpha, indexx, sigma.tot)
   return(Shift.Q)
}
















####################################################################################################
define.grid = function(gaugings,
                       sequence.grid,
                       df.limni) {
#################################################################################################### 
  # this function defines the stage grid of the RC/limni plot on the basis of the min/max gauged stage:
  grid.h       =  seq(sequence.grid[1], sequence.grid[2], sequence.grid[3]);
  
  if (is.null(df.limni)){
  min_grid.h   =  NULL;    #searching for minimum stage:
  for (i in 1:length(grid.h)) {
    if ((grid.h[i+1] > min(gaugings$h)) &  (grid.h[i] < min(gaugings$h))) {
      min_grid.h   =  grid.h[i]
    }
  }
  #
  max_grid.h   =  NULL   # searching for maximium stage:
  for (i in 1:length(grid.h)) {
    if((grid.h[i+1] >= max(gaugings$h)) & (grid.h[i] <= max(gaugings$h))) {
      max_grid.h   =  grid.h[i+1]
    }
  }
  ########
  } else {
  ########
    min_grid.h   =  NULL;    #searching for minimum stage:
    for (i in 1:length(grid.h)) {
      if ((grid.h[i+1] > min(df.limni$h_limni)) &  (grid.h[i] < min(df.limni$h_limni))) {
        min_grid.h   =  grid.h[i]
      }
    }
    #
    max_grid.h   =  NULL   # searching for maximium stage:
    for (i in 1:length(grid.h)) {
      if((grid.h[i+1] >= max(df.limni$h_limni)) & (grid.h[i] <= max(df.limni$h_limni))) {
        max_grid.h   =  grid.h[i+1]
      }
    }
  }
  
  
  start.y.legend       =  min_grid.h  # change this for each case study !!!!!!!!   <<===
  grid_limni.ylim.sim  =  c(min_grid.h, max_grid.h,  sequence.grid[3]) 
  
  
  return(list(min_grid.h           =  min_grid.h,
              max_grid.h           =  max_grid.h,
              start.y.legend       =  start.y.legend,
              grid_limni.ylim.sim  =  grid_limni.ylim.sim))
}


                          
                          
                          
                          
                    














#****************************************************************
segmentation <- function(h_Gaug, Q_Gaug, alpha, alpha_t, indexx) {
   #****************************************************************
   #ChangePoint segmentation method:
   shift_seg <- cpt.meanvar(alpha,Q=2,test.stat="Normal",penalty="SIC",method='BinSeg')
   ts = c(cpts(shift_seg))
   param.est(shift_seg)
   #plotting:
   #---------
   plot(shift_seg,cpt.width=3,cpt.col='blue',xlab = "index", type="o", pch = 19,
        ylab ="Alpha",
        main ="Segmentation of alpha time series")
   abline(h=0, col = "red", lty = 3)
   shift_time = alpha_t[cpts(shift_seg)+1]
   return(ts)
}















##########################################################################################
plot_segm <- function(h_Gaug, Q_Gaug, indexx, ts) {
##########################################################################################
   colo <- c("red","yellow","green","blue","black","grey","pink")
   plot(h_Gaug[1:indexx[ts[1]]],log(Q_Gaug[1:indexx[ts[1]]]), xlim =c(-1,1.3),
        col = colo[1],pch=1,lwd=4, ylab = "Discharge Q [m3/s]", xlab="Stage h [m]",
        main = "Stage-Discharge relation of Ardèche River at Meyras")
   abline(h=0, col="red", lty =3)
   
   if (length(ts) > 1) {
      for (i in 2:length(ts)) {
         points(h_Gaug[indexx[ts[i-1]+1]:indexx[ts[i]]],
                log(Q_Gaug[indexx[ts[i-1]+1]:indexx[ts[i]]]),
                col = colo[i],pch=1,lwd=4)
      }
      points(h_Gaug[indexx[ts[length(ts)]]:length(h_Gaug)],
             log(Q_Gaug[indexx[ts[length(ts)]]:length(h_Gaug)]),
             col = "orange",pch=1,lwd=4)
   } else if (length(ts) == 1) {
      points(h_Gaug[1:indexx[ts[1]]],
             log(Q_Gaug[1:indexx[ts[1]]]),col = colo[i],pch=1,lwd=4)
      
      points(h_Gaug[indexx[ts[1]+1]:length(h_Gaug)],
             log(Q_Gaug[indexx[ts[1]+1]:length(h_Gaug)]),
             col = "orange",pch=1,lwd=4)
   }
   #-------------------------
   # plot(t_Gaug,Q_Gaug, pch = 19, xlab= "Time [days]", ylab = "Streamflow Q = [m3/S]",
   #      main = "Streamflow time series of Ardèche River at Meyras station", sub = "Shift time detection",
   #      ylim = c(0,300))  #,xlim =c(350,380))
   # #lines(time,Qt, col="blue", lwd = 0.2)
   # abline(v=ts1, col="red", lwd ="3")
}























# fonction that computes the accuracy of the segmentation:
############################################################################################
accuracy.segmentation = function(time, ts, nshift.found, nshift.true, t2, t97, tMAP , dir) {
############################################################################################ 
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
   for (j in 1:nshift.true) {
      k =0; 
      for (s in 1:nshift.found) {
         if ((t2[s] <= ts[j]) & (t97[s] >= ts[j])) {
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
   TN= nobs.synthet - FP - FN - TP  #==================> TN number of true non cp classified as non-cp 
   if (tail(ts,1) > tail(time,1)) {
      FN=FN +1
   } 
   print("Assigning a type to each data:")
   print(c("TN=",TN, "FP=", FP, "FN=", FN, "TP=",TP))
   points = data.frame(time = time, 
                       y    = rep(0, nobs.synthet), 
                       type = type.point)
   # plot results of the data points classification into TN, FP, FN, TP:
   plot.accuracy = ggplot()+
                   scale_x_continuous(name="Time [year]",  expand = c(0.1,0.1)) +
                   scale_y_continuous(name=element_blank(), limits = c(-1,1), expand = c(0,0)) +
                   geom_point(aes(x= ts, y=rep(0, nshift.true)), size = 4, shape =4, stroke=2, col ="black" )+
                   geom_vline(xintercept = tMAP, col="blue")+
                   geom_rect(mapping= aes(xmin= t2, xmax=t97 ,ymin=-Inf, ymax=Inf), fill="blue", alpha=0.1)+
                   annotate(geom="text", x=time, y=rep(-0.1, nobs.synthet), label=type.point, color="red", size=3)+
                   geom_point(data = points, aes(x=time, y=y, color =type), size=3, shape =16)+
                   theme_bw(base_size = 25)
   ggsave(plot.accuracy, filename =paste0(dir,"/accuracy_points.png"),  device = "png", width = 24, height =8,
          dpi = 400, units = "in")
   
   # calculation of Accuracy:
   accuracy = (TP+TN)/(TP+FP+FN+TN)
   print(c("Accuracy=",accuracy))
   
   
   return(list(accuracy    =  accuracy, 
               tMAP.TP     =  tMAP.TP, 
               tTRUE.TP    =  tTRUE.TP,
               point.type  =  data.frame(TN=TN, FP=FP, FN=FN, TP=TP)))
}






























#######################################################################################################################
read.results.segmentation = function(dir.with.gaugings,
                                     dir.segment.g,
                                     file.options.general,
                                     file.options.segment,
                                     gaugings,
                                     known.shift.times,
                                     detected.shift.times.filename,
                                     pdf_detected.shift.times.filename,
                                     use.other.seg.method,
                                     compute.dates,
                                     initial.time){
#######################################################################################################################
  message("*******************************************"); flush.console()  
  message("Performance of the segmentation.  Wait ... "); flush.console()
  message("*******************************************"); flush.console()
  
  # read results and input data:
  #*****************************
  source(file.options.general)
  source(file.options.segment)
  dir.create(paste0(dir.segment.g,"/", name.folder.results))
  dir.with.segmentation.results  = paste0(dir.segment.g,"/", name.folder.results) # dir. with the results of gaugings segmentation
  Ngaug                          =  length(gaugings$Q)
  gaugings.segment               =  read.table(paste0(dir.with.segmentation.results, "/data_with_periods.txt"), sep="", header=TRUE)
  nperiods                       =  tail(gaugings.segment[,4], 1)
     
  if (tail(names(gaugings),1)  ==  "t"){
      gaugings$time = gaugings$t 
  }
  if (!is.null(known.shift.times)) {
      nshift.true            =  length(known.shift.times)
      data.annotate.off      =  data.frame(xeffect = known.shift.times)
  } else {
      data.annotate.off      = NULL
      nshift.true            = NULL
  }
  
  
  #read results of segmentation (if at least one shift is detected):
  #*****************************************************************
  if (file.exists(paste0(dir.with.segmentation.results, "/", detected.shift.times.filename)) & nperiods >1) {  #only if there are detected shift times !!
    detected.shift.times        =  read.table(paste0(dir.with.segmentation.results, "/", detected.shift.times.filename), sep="", header=TRUE)
    pdf.detected.shift.times    =  read.table(paste0(dir.with.segmentation.results, "/", pdf_detected.shift.times.filename), sep="", header=TRUE)  
    nshift.detected             =  length(detected.shift.times$tMAP)
    data.annotate.segm.gaug     =  data.frame(q2    = detected.shift.times$t2,  
                                              MAP   = detected.shift.times$tMAP, 
                                              q97   = detected.shift.times$t97,
                                              t.adj = detected.shift.times$treal)
    
    if (!is.null(known.shift.times)) {
    if (use.other.seg.method == FALSE){
    #####################################

        accuracy.sim                =  accuracy.segmentation(time         = gaugings$t, 
                                                             ts           = known.shift.times, 
                                                             nshift.found = nshift.detected,
                                                             nshift.true  = nshift.true,
                                                             t2           = detected.shift.times$t2, 
                                                             t97          = detected.shift.times$t97,
                                                             tMAP         = detected.shift.times$treal,
                                                             dir          = dir.with.segmentation.results)
        print("Computing RMSE between detected shift times and the known ones:")
        rootmse.sim = rmserr(x  = accuracy.sim$tMAP.TP,
                             y  = accuracy.sim$tTRUE.TP, 
                             summary = TRUE)
        oversegm    =  accuracy.sim$point.type$TP/(accuracy.sim$point.type$FP + accuracy.sim$point.type$TP)   # Precision P
        undersegm   =  accuracy.sim$point.type$TP/(accuracy.sim$point.type$FN + accuracy.sim$point.type$TP)   # Sensitivity S
        accu.sim    =  accuracy.sim$accuracy; 
        rmse.sim    =  rootmse.sim$rmse;
        DF.sim      =  (nshift.detected- nshift.true)/nshift.true
    ##########
    } else {   # with R package (no uncertainty available on the change point !!!)
    ##########  
        nshift.detected      =  length(detected.shift.times$t)
        data.annotate.segm.gaug = data.frame(tstart = detected.shift.times$start,
                                             t      = detected.shift.times$t,
                                             tend   = detected.shift.times$end,
                                             t.adj  = detected.shift.times$tadj,
                                             t2     = detected.shift.times$start,
                                             MAP    = detected.shift.times$t,
                                             t97    = detected.shift.times$end)
        diffff        = nshift.true - nshift.detected
        over          = 0; 
        under         = 0;
        trueshifttime = known.shift.times
        if (diffff > 0) {  #undersegmentation
            under          = diffff
            foundshifttime = detected.shift.times$t
            for (h in 1:under) {
                for (g in 1:length(nshift.detected)) {
                    elimin         = which.max(abs(detected.shift.times$t[g] - trueshifttime))
                    trueshifttime  = trueshifttime[trueshifttime != trueshifttime[elimin]]
                }
            }
        } else if (diffff < 0) { #oversegmentation
            over           = - diffff
            foundshifttime = detected.shift.times$t
            for (hh in 1:over) {
                for (d in 1:length(nshift.true)) {
                    elimin          =  which.max(abs(foundshifttime - known.shift.times[d]))
                    foundshifttime  =  foundshifttime[foundshifttime != foundshifttime[elimin]]
                }
            }
      
        } else if (diffff == 0) {
            trueshifttime  = known.shift.times
            foundshifttime = detected.shift.times$t
        }
        oversegm      = length(foundshifttime)/(over + length(foundshifttime))    #Precision P  : over =FP
        undersegm     = length(foundshifttime)/(under + length(foundshifttime))  #Sensitivity S : under=FN
        accuracy.sim  = (length(gaugings$Q) - over - under) / length(gaugings$Q)
        print(c("Accuracy=",accuracy.sim))
        accu.sim      = accuracy.sim
        print("Computing RMSE between detected shift times and the known ones:")
        rootmse.sim   = rmserr(x = trueshifttime,
                               y = foundshifttime,
                               summary = TRUE)
        rmse.sim      = rootmse.sim$rmse;
        DF.sim        =  (nshift.detected- nshift.true)/nshift.true
    }
    
    # plot results of segmentation:
    ##############################"
    # plot.time.shifts.gaugings(dir                = dir.validat.exp1.sim.test1, 
    #                           g.s                = g.s1[[1]],
    #                           data.annotate.off  = data.annotate.off,
    #                           data.annotate.gaug =  data.annotate.gaug.1, 
    #                           data.annotate.gaug.adjust = data.annotate.gaug.adjust.1, 
    #                           colo, 
    #                           t_Gaug =gaug.synthet[[1]]$time , h_Gaug =gaug.synthet[[1]]$h, 
    #                           df.limni=NULL, limni.labels, start.y.legend, 
    #                           grid_limni.ylim = grid_limni.ylim.sim,
    #                           pdf.ts =  pdf.ts.1, dates = NULL)
    #########
  } else {
    #########  
    
    nshift.detected        = 0
    accuracy.sim           = NULL
    rootmse.sim            = NULL
    oversegm               = 1  # P, precision
    undersegm              = 0  # S, sensitivity
    rmse.sim               = 0
    if (!is.null(known.shift.times)) {
      accu.sim             = (Ngaug - nshift.true)/Ngaug
      DF.sim               = (nshift.detected- nshift.true)/nshift.true
    } else {
      accu.sim             = NULL    
      DF.sim               = NULL
    }
    print("The performance cannot be evaluated since there are no official/true shift times to be compared with!")
  }
    
    
    
    
    if (compute.dates == TRUE) {
       dates.1 = 0
       for (dat in 1:length(data.annotate.segm.gaug$t.adj)) {
          dates.1[dat] <-  data.annotate.segm.gaug$t.adj[dat] + initial.time
       }
       RealDates.1        <- as.Date(dates.1,  origin = c(date_origin) )  #"1899-12-30") # DEFINE the origin date = 
       RealDates.1.new    <- strptime(as.character(RealDates.1), "%Y-%m-%d" )
       RealDates.1.newnew <- format( RealDates.1.new, "%d/%m/%Y")
       
    } else {
       RealDates.1.newnew = NULL
       
    }
  } else {
     
     print("No shift detected !")
     accuracy.sim             = NULL
     accu.sim                 = NULL
     rootmse.sim              = NULL
     rmse.sim                 = NULL
     oversegm                 = NULL
     undersegm                = NULL
     DF.sim                   = NULL
     data.annotate.segm.gaug  = NULL
     pdf.detected.shift.times = NULL
     RealDates.1.newnew       = NULL
  }
  

  
  
  
  
  print("****************")
  print("   All done!    ")
  print("****************")
  
  return(list(accuracy.summary         = accuracy.sim,
              accuracy                 = accu.sim,
              rmse.summary             = rootmse.sim,
              rmse                     = rmse.sim,
              precision                = oversegm,
              sensitivity              = undersegm,
              detec.factor             = DF.sim,
              data.annotate.segm.gaug  = data.annotate.segm.gaug,
              pdf.detected.shift.times = pdf.detected.shift.times,
              data.annotate.known      = data.annotate.off,
              gaugings                 = gaugings,
              gaugings.segment         = gaugings.segment,
              real.dates.shifts        = RealDates.1.newnew
  ))
}




















