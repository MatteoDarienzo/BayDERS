############################################################################################################################
#BaRatin - SPD (Mansanarez et al., 2019)                   
BaRatin_SPD.bac_app <- function(dir_code,
                                dir.BaM, 
                                dir.segment.g,
                                file.options.general,
                                file.options.segment,
                                file.options.SPD,
                                stage.record){
#                                gaug.periods.df){
############################################################################################################################
      # upload inputs and options:
                message("################################################################")
                message("                 BaRatin-SPD (Mansanarez et al., 2019)          ")
                message("################################################################")
                message("                  Stage-Period-Discharge model                  ")
               
  
                source(file.options.general)
                source(file.options.SPD)
                #source(file.options.segment)
                
                
                
                ########################
                nsim                   = 10000 # number of MC samples for prior progation
                user.choice.folder.SPD = 1
                folder = dir.segment.g
                while (user.choice.folder.SPD >0){
                  list_folders = list.dirs(path = folder, recursive=FALSE)
                  if (length(list_folders) > 0) {
                        user.choice.folder.SPD  =dlgInput(c("Looking for file 'data_with_periods.txt':",
                                                            "************************************",
                                                            "move to folder: [index]", " ",
                                                            paste0(0, " = ", folder, " (stay in same folder)"),
                                                            " ",
                                                            paste0(seq(1, length(list_folders), 1), " = ",
                                                            list_folders)),  Sys.info()[" "])$res
                        if (user.choice.folder.SPD != "0"){
                          folder = list_folders[as.numeric(user.choice.folder.SPD)]
                        }
                  } else {
                        print("*****No folders. Please check the path to 'data_with_periods.txt'!")
                        break
                  }
                } 

                print(paste0("Opening folder: ", folder))
                list_selected_results_SPD   = list.files(path = folder, recursive=FALSE)
                print(list_selected_results_SPD)
                dir.segment.SPD= folder
                  
                if (length(list_selected_results_SPD) > 0) {
                    if (any(list_selected_results_SPD == "data_with_periods.txt")) {
                      g.s.results.1                = read.table(paste0(dir.segment.SPD,"/data_with_periods.txt"), sep ="\t", header =TRUE) 
                      #print(g.s.results.1) 
                      pdf.ts.results.1             = read.table(paste0(dir.segment.SPD,"/pdf_ts.txt"), sep ="\t", header =TRUE) 
                      nperiod.from.segm.gaugings.1 = tail(g.s.results.1[[4]],1)
                      shift.times.gaugings.1       = read.table(file =paste0(dir.segment.SPD,"/shift_times.txt"), header = TRUE)
                      
                      #names(shift.times.gaugings.1) = c("t2", "t10", "tMAP", "t90","t97", "t.adj")
                      names(shift.times.gaugings.1) = c("tMAP",	"t.adj", "t2", "t10", "t90", "t97")
                      data.annotate.gaug.1 <- data.frame( q2   = shift.times.gaugings.1$t2,
                                                          q10  = shift.times.gaugings.1$t10,
                                                          MAP  = shift.times.gaugings.1$tMAP,
                                                          q90  = shift.times.gaugings.1$t90,
                                                          q97  = shift.times.gaugings.1$t97,
                                                          t.adj= shift.times.gaugings.1$t.adj)
                      data.annotate.gaug.1.sort = data.annotate.gaug.1[  with( data.annotate.gaug.1, order(t.adj)),]
                      data.annotate.gaug.1.sort = cbind(data.annotate.gaug.1.sort,  t.adj = shift.times.gaugings.1$t.adj)
                      pdf.ts.results.1.sort     = pdf.ts.results.1[,  with( data.annotate.gaug.1, order(t.adj)) ]
                      
                    } else {
                      err1 = paste0("******** ERROR: No file 'data_with_periods.txt' in selected folder ", dir.segment.SPD, " Please check!")
                      return(err1)
                    }
                }
                
                
                
                
                
                
                # dir.create(paste0(dir.segment.g,"/", name.folder.results))
                dir.segm.results= dir.segment.SPD    #paste0(dir.segment.g,"/", name.folder.results) # dir. with the results of gaugings segmentation
                dir.create(paste0(dir.segm.results,"/", name.folder.results.SPD))
                dir.SPD.results = paste0(dir.segm.results,"/", name.folder.results.SPD)
                dir.SPD.config  = paste0(dir.BaM,"/BaRatin_SPD")
                message(paste0("Description:
-  it estimates the multi-period RC (through a Bayesian 
   approach) in a unique model where some parameters 
   are in common with all periods and other parameters 
   are period-specific.
    
You need to fill the options file: '../", basename(file.options.SPD), "':"))
message(paste0("-  define the parameters that are unstable. 
-  specify priors for the global/local changes between periods.
-  specify priors for all parameters (as for BaRatin method).

For this application you are:
-  reading gaugings and periods from (segmentation results): \n  '../",  basename(dir.segm.results),"'
-  saving SPD results in: \n  '../", basename(dir.SPD.results),"'
-  using 'module_BaRatinSPD.r' and 'module_prior_propagation.r' 

*****************************************************************
                        Configuration ...
*****************************************************************"))
                
      # save the gaugings dataframes with the periods in the BaM configuration folder:
                gaug.periods.df = read.table(paste0(dir.segm.results,"/data_with_periods.txt"), header= TRUE)
                write.table(gaug.periods.df[,1:4], 
                            paste0(dir.SPD.config,"/Gaugings_data_SPD.txt"), 
                            sep ="\t", row.names=FALSE, col.names = c("h","Q","uQ","Period"))
      # Read the gaugings data:
                data4BaRatinSPD <- read.table(paste0(dir.SPD.config,"/Gaugings_data_SPD.txt"), header= TRUE)
                n.periods = tail(data4BaRatinSPD$Period, 1) 
      
      # Priors for incremental global changes:
                message("You have selected:")
                if (global.change == TRUE) {
                  if (n.periods >1){
                    message("- global changes (vertical traslation of the whole section):")
                    message(paste0("  ==> prior of dg ~ N(", dg.prior[1], ",",dg.prior[2], ")."))
                    d.g <- matrix(rnorm(n    = nsim*(n.periods -1), 
                                        mean = dg.prior[1], 
                                        sd   = dg.prior[2]), 
                                  nrow = nsim,
                                  ncol = n.periods - 1)
                  } else {
                    message("- No global changes: one period only")
                    d.g = NULL
                  }
                } else {
                  d.g = NULL
                  message("- No global changes.")
                }
                
      ############
      # Priors for incremental local changes:
                if (local.change == TRUE) {
                  if (n.periods >1){
                    message("- local changes (vertical/horiz traslation of the lowflow riffle):")
                    message(paste0("  ==> prior of dl ~ N(", dl.prior[1], ",",dl.prior[2], ")."))
                    d.l <- matrix(rnorm(n    = nsim*(n.periods -1), 
                                        mean = dl.prior[1], 
                                        sd   = dl.prior[2]),
                                  nrow = nsim, 
                                  ncol = n.periods -1)
                  } else {
                    d.l = NULL
                    message("- No local changes: one period only")
                  }
                } else {
                  d.l =NULL
                  message("- No local changes.")
                }
                
                
      ############
      # Priors for incremental global changes:
                if (width.change == TRUE) {
                  if (n.periods >1){
                    message("- changes of width (transversal variations of the control):")
                    message(paste0("  ==> prior of dB ~ N(", dB.prior[1], ",",dB.prior[2], ")."))
                    d.B <- matrix(rlnorm(n       = nsim*(n.periods-1),
                                         meanlog = dB.prior[1], 
                                         sd      = dB.prior[2]),
                                  nrow = nsim,
                                  ncol = n.periods -1)
                  } else {
                    d.B = NULL
                    message("- No changes of the channel width: one period only")
                  }
                } else {
                  d.B = NULL
                  message("- No changes of the channel width.")
                }
                
                
                message(paste0("- Number of hydraulic controls (param. b,a,c): ", ncontrols)) 
                message(paste0("- Dynamic parameters: ",param.var))
                message(paste0("- Number of mcmc cycles: ",Ncycles))
                

                
                
                
      # Define priors on "physical" parameters:
                names.bac = c('b','a','c')#RC Parameter a:
                if (propagat == TRUE){
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
                stdev.var.param.initial = st_b.prior
                
                #Parameter c:
                for (c in 1:ncontrols){
                  if (c.distr== "LogNormal"){
                    # Transform priors from gaussian to lognormal:
                    # function which takes gaussian mean and stdev and gives lognormal mean and stdev
                    st_c.prior[c] = round(Transf_Gauss_lognorm(E = c.prior[c], stdev = st_c.prior[c])$sd, digits = 4)
                  }
                }
                b=NULL; c=NULL; Bw=NULL; Cr=NULL; g=NULL; Bc=NULL; KS=NULL; S0=NULL;  # Priors are defined through Monte Carlo samples

                for (i in 1:ncontrols) {
                    b[[i]] <- rnorm(n=nsim, mean=b.prior[i], sd=st_b.prior[i]) # offset (stage for which Q=0 for that control)
                    c[[i]] <- rnorm(n=nsim, mean=c.prior[i], sd=st_c.prior[i]) # exponent
                    
                    if (  propagat == TRUE){
                    #IN CASE OF BARATIN :
                    if (control.type[i] == "rect.weir") {
                          Bw[[i]] <- rlnorm(n=nsim, meanlog=log(Bw.prior[i]), sdlog=st_b.prior[i])  # weir width
                          Cr[[i]] <- rlnorm(n=nsim, meanlog=log(Cr.prior[i]), sdlog=st_Cr.prior[i]) # Discharge coefficient
                          g[[i]]  <- rlnorm(n=nsim, meanlog=log(g.prior[i]),  sdlog=st_g.prior[i])  # gravity
                    } else if (control.type[i] == "rect.channel") {
                          Bc[[i]] <- rlnorm(n=nsim, meanlog=log(Bc.prior[i]), sdlog=st_Bc.prior[i]) # channel width
                          KS[[i]] <- rlnorm(n=nsim, meanlog=log(KS.prior[i]), sdlog=st_KS.prior[i]) # Strickler coefficient
                          S0[[i]] <- rlnorm(n=nsim, meanlog=log(S0.prior[i]), sdlog=st_S0.prior[i]) # slope
                    }
                    } else {
                          a[[i]] <- rnorm(nsim, mean=a.prior[i], sd=st_a.prior[i]) # prior of parameter a
                      
                    }
                }
                
                
                if (n.periods ==1){
                    message("- Estimation of one RC only (one period) !")
                } else {
                    message(paste0("- Number of periods (RCs): ",n.periods))
                }
                
                
                #d.b <- matrix(rnorm(nsim*(nperiod-1), mean=d.b1.prior, sd=st_db1.prior), nrow=nsim, ncol=nperiod-1)
                d.b1 = NULL
                #d.a <- matrix(rnorm(nsim*(nperiod-1), mean=d.a1.prior, sd=st_d.a1.prior), nrow=nsim, ncol=nperiod-1)
                d.a1 = NULL
  
  
      # Starting point for all parameters:
                message("- Defining the starting values of all parameters...")
                start <- list()
                for (i in 1:ncontrols) {
                    if (control.type[i] == "rect.weir") {
                         start[[paste0("b",i)]]  = b.prior[i] 
                         start[[paste0("c",i)]]  = c.prior[i]
                         start[[paste0("Bw",i)]] = Bw.prior[i]
                         start[[paste0("Cr",i)]] = Cr.prior[i]
                         start[[paste0("g",i)]]  = g.prior[i]
                    } else if (control.type[i] == "rect.channel") {
                         start[[paste0("b",i)]]  =  b.prior[i] 
                         start[[paste0("c",i)]]  = c.prior[i]
                         start[[paste0("Bc",i)]] = Bc.prior[i]
                         start[[paste0("KS",i)]] = KS.prior[i]
                         start[[paste0("S0",i)]] = S0.prior[i]
                    } else if (control.type[i] == "other.function") {
                         # TO DO !!!!!!!!
                         ###################
                         # start[[paste0("b",i)]]  =  b.prior[i] 
                         # start[[paste0("c",i)]]  = c.prior[i]
                         # start[[paste0("a",i)]]  = a.prior[i]
                    }
                }
      # changes parameters:
                for (i in 1:(n.periods-1)) {
                  if (global.change == TRUE) {
                       start[[paste0("dg", i)]] = dg.prior[1]
                  } 
                  if (local.change == TRUE) {
                       start[[paste0("dl", i)]] = dl.prior[1]
                  }
                  if (width.change == TRUE) {
                       start[paste0("dB", i)]   = dB.prior[1]
                  }
                }
                #print(unlist(start))
                
                list.prior.RC = list(b  = b, 
                                     Bw = Bw, Cr=Cr, g=g,   # section control
                                     Bc = Bc, KS=KS, S0=S0, # channel control
                                     c  = c)
                list.prior.changes = list(d.b1  = d.b1,
                                          d.B   = d.B,
                                          d.g   = d.g,
                                          d.l   = d.l)
                
      # Perform Monte-Carlo propagation   (call function from "Prior_Propagation.R")
                message("- Propagation of priors for BaRatin-SPD ...")
                MC = propagation(list.prior.RC       =  list.prior.RC,
                                 list.prior.changes  =  list.prior.changes,
                                 start               =  start,    # starting vector
                                 changes.method      =  changes.method, #choice between "cumsum" and "cumprod" methods ofr incremental changes   
                                 n.parvar            =  n.parvar,
                                 isVar               =  isVar,
                                 global.change       =  global.change,
                                 local.change        =  local.change,
                                 width.change        =  width.change,
                                 ncontrols           =  ncontrols)
            
      # Fit prior distribution on MCMC samples:
                margins = c();
                names   = c();
                k       = 0;    # define marginal prior distributions and parameter names
                margins.bac  = c(b.distr, a.distr, c.distr) 
                for(i in 1:ncontrols){
                  for(j in 1:3){
                    k=k+1
                    if(isVar[k]){
                         margins = c(margins, rep(margins.bac[j], MC$nperiod))
                         names   = c(names, paste0(names.bac[j],i,'_',1:MC$nperiod))
                    }else{
                         margins = c(margins, margins.bac[j])
                         names   = c(names, paste0(names.bac[j], i))
                    }
                  }
                }
      # fit multivariate prior distribution :  call function from "Prior_Propagation.R"
                prior = fit(sim       = MC$sim, 
                            margins   = margins, 
                            names     = names, 
                            ncol      = ncolumns,
                            rowmax    = rowmax, 
                            file.save = dir.SPD.results)  
  
      # write config files from "Prior_Propagation.R" module
                message("- Configuration of files for BaM in '../BaM_exe/BaRatin_SPD'")
                colperiod = 4     # the column where period will be stored in BaM data file
                writeConfigFilesSPD(prior     = prior, 
                                    start     = MC$start2, 
                                    ncontrol  = ncontrols,
                                    nperiod   = MC$nperiod, 
                                    colperiod = colperiod, 
                                    isVar     = isVar, 
                                    model     = 'Config_Model.txt', 
                                    correl    = 'PriorCorrelation.txt',
                                    names.bac = c('b','a','c'), 
                                    directory = dir.SPD.config)
                nobs = length(data4BaRatinSPD$Q)
                BaRatin_SPD_config(dir.BaM           = dir.BaM, 
                                   dir.SPD.config    = dir.SPD.config, 
                                   pred              = FALSE, 
                                   nobs              = nobs, 
                                   M                 = M, 
                                   remnant           = remnant.err.model, 
                                   g1.prior          = g1.prior,
                                   g2.prior          = g2.prior,
                                   g1.distr.type     = g1.distr.type,
                                   g2.distr.type     = g2.distr.type,
                                   Ncycles           = Ncycles)
     
     if (plot.results.only == FALSE) {
     # Exection of BaM.exe:
                setwd(dir.BaM)
                message("
*****************************************************************
                        Launching BaM  ...
*****************************************************************")
                system2(paste0(dir.BaM,"/BaM_2exp_pool2.exe"))
                
     # save results:
                # specification for each parameter
                k=0; m=1; fname = c(); val=list();  names.bac = c('b','a','c');
                for(i in 1:ncontrols){ # loop on each hydraulic control
                  for(j in 1:length(names.bac)){ # loop on each b-a-c parameter of the control
                    k=k+1
                    # parname
                    pname = paste0(names.bac[j], i)
                    val   = c(val, pname)
                    if(isVar[k]){ # period-specific parameter
                        fname = c(fname, paste0(dir.SPD.config,"/Config_", pname, "_VAR.txt"))
                    }
                  }
                }  
                print("- Reading results files from BaM ...")
                list.of.files.spd   <- c(
                  paste0(dir.SPD.config,"/Config_ControlMatrix.txt"),
                  paste0(dir.SPD.config,"/PriorCorrelation.txt"),
                  paste0(dir.SPD.config,"/Results_MCMC_Cooked.txt"), 
                  paste0(dir.SPD.config,"/Results_Residuals.txt"),
                  paste0(dir.SPD.config,"/Results_Summary.txt"),
                  paste0(dir.SPD.config,"/Config_Model.txt"),
                  paste0(dir.SPD.config,"/Gaugings_data_SPD.txt")
                )
                
                for (ll in 1:length(list.of.files.spd)) {
                     file.copy(list.of.files.spd[ll], dir.SPD.results, overwrite = TRUE)
                }    
                for (ll in 1:length(fname)) {
                     file.copy(fname[ll], dir.SPD.results, overwrite = TRUE)
                }    
                
                   
     } else {
        message("No computation --> Reading and plotting the results from existing past simulations! ") 
     }     
      
      
                
                
      ######################################################################################################
      # Plotting results of BaRatin-SPD in terms of rating curves :
      SPD = plot.SPD(dir.BaM              = dir.BaM, 
                     dir.SPD.results      = dir.SPD.results,
                     dir.SPD.config       = dir.SPD.config,
                     nperiod              = MC$nperiod,
                     stage.record         = stage.record,
                     gaug.periods.df      = gaug.periods.df,
                     file.options.general = file.options.general,
                     file.options.SPD     = file.options.SPD,
                     data.annotate.gaug.1 = data.annotate.gaug.1.sort,
                     activation.stages    = show.activation.stages) 

      if (!is.null(SPD$err)){
         return(list(SPD$err))
        
      } else {
                
         message("****************")
         message("   All done!    ")
         message("****************")
         return(list(dir.reference = dir.segm.results,
                     dir.SPD       = dir.SPD.results,
                     res.SPD       = SPD))
      }
}  
 


  
























  
  








###########################################################################################
BaRatin_SPD_config <- function(dir.BaM, 
                               dir.SPD.config,
                               pred, 
                               nobs, 
                               M, 
                               remnant,
                               g1.prior, g2.prior,    
                               g1.distr.type, g2.distr.type,                              
                               Ncycles) {                                               
###########################################################################################
      # creation of Config_BaM.txt:
                file.BaM <- paste0(dir.BaM,"/Config_BaM.txt")
                cat('"BaRatin_SPD/"',             file = file.BaM , sep="\n", append = FALSE)
                cat('"Config_RunOptions.txt"',    file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_Model.txt"',         file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_ControlMatrix.txt"', file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_Data.txt"',          file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_RemnantSigma.txt"',  file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_MCMC.txt"',          file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_Cooking.txt"',       file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_Summary.txt"',       file = file.BaM , sep="\n", append = TRUE)
                cat('"Config_Residuals.txt"',     file = file.BaM , sep="\n", append = TRUE)
                if (pred == TRUE) {
                  cat('"Config_Pred_Master.txt"', file = file.BaM , sep="\n", append = TRUE)
                } else {
                  cat('""',                       file = file.BaM , sep="\n", append = TRUE)
                }
      # creation of Config_Data.txt
                file.name2 = paste0(dir.SPD.config,"/Config_Data.txt")
                cat("'BaRatin_SPD\\Gaugings_data_SPD.txt'", file =file.name2, sep="\n")
                cat(1, file = file.name2, append = TRUE,sep="\n")
                cat(nobs, file = file.name2, append = TRUE,sep="\n")
                cat(4, file =file.name2, append = TRUE,sep="\n")
                cat(1, file =file.name2, append = TRUE,sep="\n")
                cat(0, file =file.name2, append = TRUE,sep="\n")
                cat(0, file =file.name2, append = TRUE,sep="\n")
                cat(0, file =file.name2, append = TRUE,sep="\n")
                cat(2, file =file.name2, append = TRUE,sep="\n")
                cat(3, file =file.name2, append = TRUE,sep="\n")
                cat(0, file =file.name2, append = TRUE,sep="\n")
                cat(0, file =file.name2, append = TRUE,sep="\n")
      # creation of Config_MCMC.txt
                file.mcmc = paste(dir.SPD.config,"/Config_MCMC.txt",sep="")
                cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
                cat(100,     file = file.mcmc, append = TRUE,sep="\n")   #Nadapt
                cat(Ncycles, file = file.mcmc, append = TRUE,sep="\n")  #Ncycles
                cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")    #minMoveRate
                cat(0.5,     file = file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
                cat(0.9,     file = file.mcmc, append = TRUE,sep="\n")    #DownMult
                cat(1.1,     file = file.mcmc, append = TRUE,sep="\n")    #UpMult
                cat(0,       file = file.mcmc, append = TRUE,sep="\n")      #mode for init jump distr
                cat("****",  file = file.mcmc, append = TRUE,sep="\n") 
                cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")    #MultFact
                cat(0.1,     file = file.mcmc, append = TRUE, sep=",")     #RC MultiFact
                cat(0.1,     file = file.mcmc, append = TRUE, sep=",")     
                cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")
                cat(0.1,     file = file.mcmc, append = TRUE, sep=",")      #Remnant MultiFact
                cat(0.1,     file = file.mcmc, append = TRUE,sep="\n")
      # creation of Config_ControlMatrix.txt
                file.matrix = paste0(dir.SPD.config,"/Config_ControlMatrix.txt")
                write.table(M, file = file.matrix, row.names = FALSE, col.names = FALSE)   #M control matrix
                cat("10.",     file = file.matrix,  append = TRUE, sep="\n")  #hmax
      # Creation of Config_RemnantSigma.txt
                file.remnant = paste0(dir.SPD.config,"/Config_RemnantSigma.txt")
                if (remnant == "Linear") {
                    cat("'Linear'",    file = file.remnant, sep="\n")                                   #! Function f used in sdev=f(Qrc) 
                    cat(2,             file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
                    cat("gamma1",      file = file.remnant, append = TRUE, sep="\n")                    #! Parameter Name
                    cat(g1.prior[3],   file = file.remnant, append = TRUE, sep="\n")                    #! Initial Guess
                    cat(g1.distr.type, file = file.remnant, append = TRUE, sep="\n")                    #! Prior distribution
                    cat(g1.prior[1],   file = file.remnant, append = TRUE, sep=",")
                    cat(",",           file = file.remnant, append = TRUE, sep=",")
                    cat(g1.prior[2],   file = file.remnant, append = TRUE, sep="\n")
                  
                    cat("gamma2",      file = file.remnant, append = TRUE, sep="\n")                    #! Parameter Name
                    cat(g2.prior[3],   file = file.remnant, append = TRUE, sep="\n")                    #! Initial Guess
                    cat(g2.distr.type, file = file.remnant, append = TRUE, sep="\n")                    #! Initial Guess
                    cat(g2.prior[1],   file = file.remnant, append = TRUE, sep=",")
                    cat(",",           file = file.remnant, append = TRUE, sep=",")
                    cat(g2.prior[2],   file = file.remnant, append = TRUE, sep="\n")
                  
                } else if (remnant == "Constant") {
                    cat("'Constant'",  file = file.remnant, sep="\n")                                   #! Function f used in sdev=f(Qrc) 
                    cat(1,             file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
                    cat("gamma1",      file = file.remnant, append = TRUE, sep="\n")                    #! Parameter Name
                    cat(g1.prior[3],   file = file.remnant, append = TRUE, sep="\n")                    #! Initial Guess
                    cat(g1.distr.type, file = file.remnant, append = TRUE, sep="\n")                    #! Prior distribution
                    cat(g1.prior[1],   file = file.remnant, append = TRUE, sep=",")
                    cat(",",           file = file.remnant, append = TRUE, sep=",")
                    cat(g1.prior[2],   file = file.remnant, append = TRUE, sep="\n")
                }
      # Creation of Config_RunOptions.txt:
                file.run = paste(dir.SPD.config,"/Config_RunOptions.txt",sep="")
                cat(".true.",        file = file.run, sep="\n")                             # Do MCMC?
                cat(".true.",        file = file.run, append = TRUE, sep="\n")              # Do MCMC summary?
                cat(".true.",        file = file.run, append = TRUE, sep="\n")              # Do Residual diagnostics?
                if (pred == TRUE) {
                     cat(".true.",   file = file.run, append = TRUE, sep="\n")         # Do Predictions?
                } else {
                     cat(".false.",  file = file.run, append = TRUE, sep="\n")        # Do Predictions?
                }
     # creation of Config_Residuals.txt
                file.residuals = paste(dir.SPD.config,"/Config_Residuals.txt",sep="")
                cat('"Results_Residuals.txt"',   file = file.residuals, sep="\n")     # Result file
     # creation of Config_Summary.txt
                file.summary = paste(dir.SPD.config,"/Config_Summary.txt", sep="")
                cat('"Results_Summary.txt"',     file = file.summary, sep="\n")         # Result file
     # creation of Config_Cooking.txt
                file.cooking = paste(dir.SPD.config,"/Config_Cooking.txt", sep="")
                cat('"Results_MCMC_Cooked.txt"', file = file.cooking, sep="\n")     # Result file
                cat(0.5,                         file = file.cooking, append = TRUE, sep="\n")               # Burn factor
                cat(10,                          file = file.cooking, append = TRUE, sep="\n")                # Nslim
}
    
  
  



























  



















############################################################################################################
# Plotting multiperiod rating curves in linear and log scales:                                              
plot.SPD <- function(dir.BaM, 
                     dir.SPD.results, 
                     dir.SPD.config,
                     nperiod, 
                     stage.record,
                     gaug.periods.df,           
                     file.options.general,
                     file.options.SPD,
                     data.annotate.gaug.1,
                     activation.stages) {                                              
############################################################################################################
     # Settings:
                source(file.options.general)
                source(file.options.SPD)
                message("- Reading results of BaRatin-SPD. Wait ... "); flush.console()
     # Data loading:
                if (op_cooked){ # if cooked mcmc results:
                   if (file.exists(paste0(dir.SPD.results,"/Results_MCMC_Cooked.txt"))){
                       data.MCMC.cooked = as.matrix(read.table(paste0(dir.SPD.results,"/Results_MCMC_Cooked.txt"), header=TRUE,dec=".", sep="")) # cooked mcmc
                   } else {
                       errr = paste0("********* ERROR: file '", paste0(dir.SPD.results,"/Results_MCMC_Cooked.txt"), "' does not exist !! Please check path")
                       return(list(err = errr))
                   }
                ########
                } else {
                   if (file.exists(paste0(dir.SPD.results,"/Results_MCMC.txt"))){
                       data.MCMC.raw     = as.matrix(read.table(paste0(dir.SPD.results,"/Results_MCMC.txt"), header=TRUE,dec=".", sep="")) # Raw MCMC
                       nsim.tmp          = length(data.MCMC.raw[,1])
                       data.MCMC.cooked  = data.MCMC.raw[seq(nburn*nsim.tmp, nsim.tmp, nslim),]
                       rm("data.MCMC.raw")
                  } else {
                       errr = paste0("********* ERROR: file '", paste0(dir.SPD.results,"/Results_MCMC.txt"), "' does not exist !! Please check path")
                       break
                       return(list(err = errr))
                  }
                }

                logPost.ncol      = which(colnames(data.MCMC.cooked) == "LogPost")
                
                if (file.exists(paste0(dir.SPD.results,"/Results_Summary.txt"))){
                   data.MCMC.MaxPost = as.numeric(read.table(paste0(dir.SPD.results, "/Results_Summary.txt"),row.names=1, dec=".",sep="", skip = 16))
                } else {
                   errr2 = paste0("********* ERROR: file '", dir.SPD.results,"/Results_Summary.txt", "' does not exist !! Please check path")
                   break
                   return(list(err = errr2))
                }
                
                nsample           = length(data.MCMC.cooked[,1])
                min.grid          = min(data.MCMC.cooked[,1]) 
                hgrid             = seq(xlim.wind[1] - 1, xlim.wind[2] + 1,  0.01) 
                # hgrid           = seq(-4,  xlim.wind[2] + 1,  0.001) 
  
                
                
                
     # Propagation: Convert MCMC into rating curve results
     # Rating curve with discharge error propagation:
                ##########################################################################################
                RC_controls = function(theta, h, M, ncontrols, op=1) {
                ##########################################################################################
                               Q        = 0*h
                               mask     = list();
                               #stop     = FALSE
                               if  (ncontrols >1) {
                               for (nc in 1:(ncontrols - 1)) {
                                   # if ((abs(theta[ncontrols*3 + 3 + nc] - theta[ncontrols*3 + 3 + nc+1]) 
                                   #      <= 0.01)|(any(h <= theta[ncontrols*3 + 3 +1]))) {
                                   #   
                                   #    break 
                                   #    stop= TRUE
                                   # } else {
                                      mask[[nc]] = which(((h >  theta[ncontrols*3 + 3 + nc]) +
                                                          (h <= theta[ncontrols*3 + 3 + nc+1])) == 2)
                                  # }
                               }
                               }
                               
                               #if (stop==FALSE){
                                   mask[[ncontrols]] = which(h > theta[ncontrols*3 + 3 + ncontrols])
                                   for (segm in 1:ncontrols){  # for each segment of stage h
                                      Q[mask[[segm]]] = 0 
                                      for (ccc in 1:ncontrols){  # for each control 
                                         if (M[segm, ccc] ==1){
                                            Q[mask[[segm]]] =   Q[mask[[segm]]]  +
                                            theta[3*ccc - 1]*(h[mask[[segm]]] - theta[3*ccc -2])^theta[3*ccc]
                                         }
                                      }
                                   }   
                                   # 
                                   if(op == 1) {
                                       resQ = sapply(Q, 
                                               function(Q,theta){Q + rnorm(1, 
                                                                           mean = 0, 
                                                                           sd = theta[1]+ theta[2]*Q)}, 
                                               theta = c(theta[3*ncontrols +1], theta[3*ncontrols +2]))
                                   } else {
                                       resQ = Q
                                   }
                               # } else {
                               #     resQ =NA
                               # }
                return(resQ)
                }
                
                
                
                
                
                
     # Rating curve for maximum posterior:   
                ###########################################################################################
                RC_controls_mp = function(theta, h, M, ncontrols, op=1){ 
                ###########################################################################################
                                 Q        = 0*h
                                 mask     = NULL
                                 if  (ncontrols >1) {
                                   for (nc in 1:(ncontrols - 1)) {
                                       mask[[nc]] = which(((h >  theta[ncontrols*3 + 2 + nc]) +
                                                           (h <= theta[ncontrols*3 + 2 + nc+1])) == 2)
                                   }
                                 }
                                 mask[[ncontrols]] = which(h> theta[ (ncontrols*3 + 2 + ncontrols)])
                                 for (segm in 1:ncontrols){  # for each segment of stage h
                                    Q[mask[[segm]]] = 0 
                                   for (ccc in 1:ncontrols){  # for each control 
                                     if (M[segm, ccc] ==1){
                                       Q[mask[[segm]]] =   Q[mask[[segm]]]  + 
                                       theta[3*ccc - 1]*(h[mask[[segm]]] - theta[3*ccc -2])^theta[3*ccc]
                                     }
                                   }
                                 }  
        
                  return(Q)
                }
                
     # Initialisation:
                MCMC.save    =  matrix(NA, nrow = nperiod*nsample, ncol = 3*ncontrols +3+ ncontrols + 1)  
                ### b1 a1 c1 b2 a2 c2 b3 a3 c3 gamma1 gamma2 logpost k1 k2 k3  (nper)  # 1  2  3  4  5  6  7  8  9 10 11 12  13(14)
                MaxPost.save =  matrix(NA, nrow = nperiod,  ncol = 3*ncontrols +2+ ncontrols + 1)  
                ### b1 a1 c1 b2 a2 c2 b3 a3 c3 gamma1 gamma2  k1 k2 k3  (nper)  # 1  2  3  4  5  6  7  8  9  10  11 12
                # assign each column of the mcmc results file to the matrixes MCMC.save and MaxPost.save:
                for(i in 1:nperiod){
                  # For all mcmc:
                  col.num = c(); 
                  #----------------------------------------------------------
                  m=1
                  if (isVar[1] == TRUE){   # b1
                      col.num[m] = i 
                  } else {
                      col.num[m] = 1
                  }
                  #----------------------------------------------------------
                  m=2
                  if (isVar[2] == TRUE){    # a1
                    if (isVar[1] == TRUE){
                       col.num[m] = nperiod + i
                    } else {
                       col.num[m] = col.num[m-1] + i
                    }
                  } else {
                    if (isVar[1] == TRUE){
                       col.num[m] = nperiod + 1
                    } else {
                       col.num[m] = col.num[m-1] + 1 
                    }
                  }
                  #----------------------------------------------------------
                  m=3
                  if (isVar[3] == TRUE){   #c1
                    if (isVar[2] == TRUE){
                        col.num[m] = col.num[m-1] - i + nperiod + i
                    } else {
                        col.num[m] = col.num[m-1] + i
                    }
                  } else {
                    if (isVar[2] == TRUE){
                        col.num[m] =  col.num[m-1] - i + nperiod + 1
                    } else {
                        col.num[m] = col.num[m-1] + 1
                    }
                  }
                  #----------------------------------------------------------
                  if (ncontrols >1){
                  for (ccc in 2:ncontrols) {
                    #if (any(isVar[((ccc-1)*3+1):(ccc*3)] == TRUE)){
                    for (par in 1:3){
                       m= m + 1
                       if (isVar[m] == TRUE){   
                         if (isVar[m-1] == TRUE){
                              col.num[m] = col.num[m-1] - i + nperiod + i
                         } else {
                              col.num[m] = col.num[m-1] + i
                         }
                       } else {
                         if (isVar[m-1] == TRUE){
                              col.num[m] = col.num[m-1] - i + nperiod + 1
                         } else {
                              col.num[m] = col.num[m-1] + 1
                         }
                       }
                    }
                  }
                  }
                  m = m+1; col.num[m] =  col.num[m-1] +1   # gamma1 
                  m = m+1; col.num[m] =  col.num[m-1] +1   # gamma2
                  m = m+1; col.num[m] =  col.num[m-1] +1   # logpost
                  
                  for (ccc in 1:ncontrols){
                     m = m+1                
                     if (any(isVar== TRUE)){
                         col.num[m] = col.num[3*ncontrols +3] + ccc + (i-1)*(ncontrols*2)
                     } else {
                         col.num[m] = col.num[3*ncontrols +3] + ccc
                     }
                  }
                  
                  
                  
                  col.num.mp = c(); 
                  #---------------------------------------------------------
                  m=1
                  if (isVar[1] == TRUE){  # b1
                    col.num.mp[m] = i 
                  } else {
                    col.num.mp[m] = 1
                  }
                  #---------------------------------------------------------
                  m=2
                  if (isVar[2] == TRUE){   #a1
                    if (isVar[1] == TRUE){
                      col.num.mp[m] = nperiod + i
                    } else {
                      col.num.mp[m] = col.num.mp[m-1] + i
                    }
                  } else {
                    if (isVar[1] == TRUE){
                      col.num.mp[m] = nperiod + 1
                    } else {
                      col.num.mp[m] = col.num.mp[m-1] + 1 
                    }
                  }
                  #--------------------------------------------------------
                  m=3
                  if (isVar[3] == TRUE){   #c1
                    if (isVar[2] == TRUE){
                      col.num.mp[m] = col.num.mp[m-1] - i + nperiod + i
                    } else {
                      col.num.mp[m] = col.num.mp[m-1] + i
                    }
                  } else {
                    if (isVar[2] == TRUE){
                      col.num.mp[m] =  col.num.mp[m-1] - i + nperiod + 1
                    } else {
                      col.num.mp[m] = col.num.mp[m-1] + 1
                    }
                  }
                  #--------------------------------------------------------
                  if (ncontrols >1) {
                  for (ccc in 2:ncontrols) {
                    #if (any(isVar[((ccc-1)*3+1):(ccc*3)] == TRUE)){
                    for (par in 1:3){
                      m= m + 1
                      if (isVar[m] == TRUE){   
                        if (isVar[m-1] == TRUE){
                          col.num.mp[m] = col.num.mp[m-1] - i + nperiod + i
                        } else {
                          col.num.mp[m] = col.num.mp[m-1] + i
                        }
                      } else {
                        if (isVar[m-1] == TRUE){
                          col.num.mp[m] = col.num.mp[m-1] - i + nperiod + 1
                        } else {
                          col.num.mp[m] = col.num.mp[m-1] + 1
                        }
                      }
                    }
                  }
                  }
                  m = m+1; col.num.mp[m] =  col.num.mp[m-1] +1  # gamma1 
                  m = m+1; col.num.mp[m] =  col.num.mp[m-1] +1  # gamma2
                  for (ccc in 1:ncontrols){
                    m = m+1                
                    if (any(isVar== TRUE)){
                      col.num.mp[m] = col.num.mp[3*ncontrols +2] + ccc + (i-1)*(ncontrols*2)
                    } else {
                      col.num.mp[m] = col.num.mp[3*ncontrols +2] + ccc
                    }
                  }
                  ##############################################################################################
                  MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample))
                  MaxPost.save[i,]                          = c(data.MCMC.MaxPost[col.num.mp],i)
                  ##############################################################################################
                }
                
                
                

                
                
                ##############################
                # APPLY THE PROPAGATION TO RC:
                ##############################
                message("- Applying propagation to the RC. Wait ... "); flush.console()
                RC.Post    =  apply(MCMC.save,     MARGIN = 1,  RC_controls,     h = hgrid,  M = M, ncontrols = ncontrols) # Apply RC for hgrid with structural error:
                RC.MaxPost =  apply(MaxPost.save,  MARGIN = 1,  RC_controls_mp,  h = hgrid,  M = M, ncontrols = ncontrols) # Maximum posterior computing
                data.gaug  =  read.table(paste0(dir.SPD.results,"/Gaugings_data_SPD.txt"), header=TRUE, dec=".", sep="") # Gaugings loading
                #data.gaug =  data.gaug[-which(data.gaug$Q==-9999),]
                  
              
     # Quantiles for Figure RC: 
                message("- Plotting the multi-period RC (both linear and log scales. "); flush.console()
                message("Wait ... "); flush.console()
                List.RC.quants = list(NULL)
                for(i in 1:nperiod){
                    data.tmp            = apply(RC.Post[,(nsample*(i-1)+1):(nsample*i)], MARGIN=1, quantile, probs=c(0.025,0.975), na.rm=TRUE)
                    data.tmp            = apply(data.tmp, MARGIN=c(1,2), function(x){ifelse(x<0,0,x)})
                    List.RC.quants[[i]] = data.frame(cbind(hgrid, t(data.tmp),  RC.MaxPost[,i]))
                    colnames(List.RC.quants[[i]]) = c("h", "inf", "sup", "maxpost")
                }
                inter.per = seq(1,nperiod,1)   # inter.per=c(25:35)
                nRC       = length(inter.per)
     # Palette for plot:
                # palette.per  = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] #brewer.pal(nRC,"Spectral")
                # palette.per  = colorRampPalette(c("red","orange","yellow","green","blue","grey","purple")) 
                # palette.per = distinctColorPalette(k = nperiod, altCol = FALSE, runTsne = FALSE)
                # palette.per = c("red","blue","green","orange","brown3","lightblue","black","pink", "yellow","chartreuse", "brown", "violet", "grey","blue4",
                #                 "cyan", "blueviolet", "cyan4", "darkgreen")
                ##########
                max.number.periods   = 500   # maximum number of colors for the periods and plots !!!!
                colo         = color.generation(n = max.number.periods)  
                palette.per  = unique(gaug.periods.df$color)
                palette.per  = palette.per[1:nperiod]
                data.plot.RC = data.frame(do.call("rbind", List.RC.quants[inter.per]),
                                          Period = rep(inter.per, each = length(hgrid)))
                inter.null   = which(data.plot.RC$sup==0)
     # Initialise the plot:
                vect.break   = seq(1,nperiod,1) #seq(1,9,1)
                vect.emp     = rep("",22)
                breaks.wind  = c(vect.break*0.01, vect.break*0.1, vect.break, vect.break*10, vect.break*100, vect.break*1000)
                labels.wind  = c(0.01, vect.emp, 0.1, vect.emp, 1 ,vect.emp, 10, vect.emp, 100, vect.emp, 1000, vect.emp)
                
                qpos.num     = function(x.int){inter = which(x.int==data.gaug$period); return(inter)}
                inter.gaug   = unlist(sapply(inter.per, qpos.num), recursive = TRUE, use.names = TRUE)
     # data for the plot:
                data.obs         = data.gaug[inter.gaug,]
                data.obs$period  = factor(data.obs$period)
                if (length(inter.null) >0){
                  data.RC          = data.plot.RC[-inter.null,]
                } else {
                  data.RC          = data.plot.RC
                }
                data.RC$Period   = factor(data.RC$Period)
                data.gaug        = read.table(paste0(dir.SPD.results,"/Gaugings_data_SPD.txt"),header=TRUE,dec=".",sep="") # Gaugings loading
                #data.gaug       = read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD_official.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
                data.gaug$Period = as.factor(data.gaug$Period)
                # save results of the RCs for each stable period in one file with Maxpost, Q2.5, Q97.5 total
                write.table(data.RC, paste0(dir.SPD.results,"/RC_SPD_env.txt"), sep ="\t", row.names=FALSE)
                
     # RC ggplot::
                plot.RC= ggplot(data.RC) 
                         if (activation.stages == TRUE){
                            for (ccc in 1:ncontrols){
                               if (ccc==1){
                                  plot.RC = plot.RC +
                                  geom_vline(xintercept = MaxPost.save[inter.per, ncontrols*3 +2 + ccc], 
                                             colour = palette.per, size=0.6, linetype ="dashed")
                               } else if (ccc==2){
                                  plot.RC = plot.RC +
                                  geom_vline(xintercept = MaxPost.save[inter.per, ncontrols*3 +2 + ccc], 
                                             colour = palette.per, size=0.6)
                               } else if (ccc==3){
                                  # plot.RC = plot.RC +
                                  # geom_vline(xintercept = MaxPost.save[inter.per, ncontrols*3 +2 + ccc], 
                                  #            colour = palette.per, size=1)
                               }
                            }
                         }
                         plot.RC = plot.RC + 
                         geom_smooth(aes(x=h, y=maxpost, ymax=sup, ymin=inf, group=Period, fill=Period), size=1, stat='identity', alpha = 0.1) +
                         geom_path(aes(x=h, y=maxpost, group=Period, colour=Period), na.rm=TRUE, size=1)+
                         ### Gaugings
                         geom_linerange(aes(x=h, ymax=Q+2*uQ, ymin=Q-2*uQ, colour=Period), data=data.gaug, size=0.5)+
                         geom_point(aes(x=h,y=Q,colour=Period),data=data.gaug, na.rm=na.rm, shape=16, size=size.gaugings)+
                         ### Labels
                         xlab(expression(paste("Stage h [m]",sep="")))+
                         ylab(expression(paste("Discharge Q [",m^3,".",s^-1,"]",sep="")))+
                         labs(colour = "Period")+
                         scale_colour_manual(values = palette.per)+
                         scale_fill_manual(values = palette.per)+
                         #scale_y_continuous(breaks=breaks.wind,labels=labels.wind)+
                         coord_cartesian(ylim=ylim.wind, xlim=xlim.wind)+
                         scale_x_continuous(breaks=breaks.lin.x, labels=labels.lin.x, expand=c(0,0))+
                         scale_y_continuous(breaks=breaks.lin.y, labels=labels.lin.y, expand=c(0,0))+
                         ### Theme
                         theme_bw(base_size=20)+
                         theme( axis.text        = element_text(size=15)
                               ,axis.title       = element_text(size=20,face="bold")
                               ,panel.grid.major = element_blank()  # element_line(size=1.2)
                               ,panel.grid.minor = element_blank()  # element_line(size=0.8)
                               ,legend.text      = element_text(size=20)
                               ,legend.title     = element_text(size=30)
                               ,plot.margin      = unit(c(0.5,0.5,0.5,0.5),"cm")
                               ,legend.key.size  = unit(1.5, "cm")
                               ,legend.position  = "none")
                         ggsave(plot.RC, filename = paste0(dir.SPD.results,"/RC_SPD.png"), 
                                device = "png", width = 14, height =8, dpi = 400, units = "in")
     # RC plot Logarithmic scale:
                plot.RClog = plot.RC + 
                             coord_cartesian(ylim = ylim.log.wind,
                                             xlim = xlim.wind) +
                             scale_y_log10(breaks = breaks.log, labels =labels.log, expand=c(0,0))+
                             annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black", 
                                                 size = 1, linetype = 1)
                             # ggsave(plot.RClog, filename=paste0(dir.SPD.results,"/RClog_SPD.png"), 
                             # device = "png", width = 14, height =8, dpi = 400, units = "in")
                             pdf(paste0(dir.SPD.results,"/RClog_SPD.pdf"), 14, 8 ,useDingbats=F)
                             print(plot.RClog)
                             dev.off()
                        
                             
                             
                                  
     #save results:
                # list.files.pool <- c( paste0(dir.SPD.exe,"/Results_MCMC_Cooked.txt"),
                #                       paste0(dir.SPD.exe,"/Results_Residuals.txt"),
                #                       paste0(dir.SPD.exe,"/Results_Summary.txt"),
                #                       paste0(dir.SPD.exe,"/Config_Model.txt"),
                #                       paste0(dir.SPD.exe,"/Gaugings_data_SPD.txt"))
                #                       #paste(dir.BaM.rec,"/tgrid.txt", sep=""),paste(dir.BaM.rec,"/ht_Maxpost.spag", sep=""),
                #                       #paste(dir.BaM.rec,"/ht_ParamU.spag", sep=""),paste(dir.BaM.rec,"/ht_ParamU.env", sep=""),
                #                       #paste(dir.BaM.rec,"/ht_TotalU.env", sep=""),paste(dir.BaM.rec,"/ht_TotalU.spag", sep=""))
                #                       for (i in 1:length(list.files.pool)) {
                #                         file.copy(list.files.pool[i], dir.SPD.results ,overwrite = TRUE)
                #                       }
                setwd(dir.SPD.results)
                Results_summary.SPD   <- read.table(paste(dir.SPD.results,"/Results_Summary.txt",sep=""), header =TRUE)
                Results_mcmc.SPD      <- read.table(paste(dir.SPD.results,"/Results_MCMC_Cooked.txt",sep=""), header =TRUE)
                Results_residuals.SPD <- read.table(paste(dir.SPD.results,"/Results_Residuals.txt",sep=""), header =TRUE)
                
                
                
                
                
     #Save MCMC plots:
                message("- Plotting mcmc of SPD !!!  Wait ... "); flush.console()
                setwd(dir.SPD.results)
                vertical.length = (ncontrols*3+3+ ncontrols*2)*nperiod
                mcmc = MCMCplot(doLogPost = T,
                                doPar     = T,
                                doDPar    = T, 
                                MCMCfile  = "Results_MCMC_Cooked.txt" , 
                                type      = "trace",  #="trace", # "histogram","density","scatterplot"
                                xlab      = '',
                                ylab      = '',
                                ncol      = 1, 
                                prior     = NULL,
                                burn      = 0, 
                                slim      = 1,
                                theme_size= 15)
                
                mcmc2 = MCMCplot(doLogPost = T,
                                 doPar     = T,
                                 doDPar    = T, 
                                 MCMCfile  = "Results_MCMC_Cooked.txt" , 
                                 type      = 'density', 
                                 prior     = NULL,
                                 xlab      = '',
                                 ylab      = '',
                                 ncol      = 1, 
                                 burn      = 0, 
                                 slim      = 1,
                                 theme_size= 15)  
                # ggsave(mcmc2, filename =paste(workspace,"/mcmc2_it",iter,".png", sep=""), 
                # bg = "transparent", width = 12, height =vertical.length, dpi = 100)
                pdf(paste0(dir.SPD.results,"/mcmc_SPD.pdf"), 17, vertical.length, useDingbats=F)
                print(plot_grid(mcmc, mcmc2,
                                nrow=1, ncol = 2, 
                                labels = c("Trace plots", "Density plots"),
                                label_size = 20))
                dev.off()
                
                
                
                
                
                
                
                
                
     # Boxplots of parameters b (b1, b2, ...):               
                message("- Plotting the boxplots of the posterior RC parameter 'b' for each period. Wait ... "); flush.console()
                ################################################################################################################
                # SPD.summary               = read.table(file=paste0(dir.SPD.results, "/Results_Summary.txt"))
                # SPD.mcmc.cooked           = read.table(file=paste0(dir.SPD.results, "/Results_MCMC_Cooked.txt"), header=TRUE)

                SPD.mcmc.cooked.b  = SPD.maxpost.b = SPD.box.b = dftemp =list();
                i = 1;  col.num = c();
                #---------------------------
                m=1
                if (isVar[1] == TRUE){  # b1
                  col.num[m] = i
                } else {
                  col.num[m] = 1
                }
                SPD.mcmc.cooked.b[[1]]        = Results_mcmc.SPD[,  col.num[m]:(col.num[m]+ nperiod - 1)];
                SPD.maxpost.b[[1]]            = Results_summary.SPD[16,  col.num[m]:(col.num[m]+ nperiod - 1)]
                SPD.box.b[[1]]                = unlist(data.frame(t(c(quantile( SPD.mcmc.cooked.b[[1]][,1], p = c(0.025, 0.5, 0.975)),
                                                               mean    =  mean(SPD.mcmc.cooked.b[[1]][,1]), 
                                                               stdev   =  std( SPD.mcmc.cooked.b[[1]][,1]), 
                                                               maxpost =  SPD.maxpost.b[[1]][1]))))
                names(SPD.box.b[[1]]) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
                if (nperiod ==1) { 
                  message("******* ERROR: only period, no SPD analysis !")
                }
                for (i in 2:nperiod){
                     dftemp[[1]] = unlist(data.frame(t(c(quantile(SPD.mcmc.cooked.b[[1]][,i], p = c(0.025, 0.5, 0.975)), 
                                                             mean    = mean(SPD.mcmc.cooked.b[[1]][,i]), 
                                                             stdev   = std(SPD.mcmc.cooked.b[[1]][,i]), 
                                                             maxpost = SPD.maxpost.b[[1]][i]))))
                     names(  dftemp[[1]]) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
                     SPD.box.b[[1]]   = rbind(SPD.box.b[[1]],  dftemp[[1]])
                }
                names(SPD.mcmc.cooked.b[[1]]) = seq(1, nperiod, 1)
                SPD.mcmc.cooked.b[[1]]        = SPD.mcmc.cooked.b[[1]] %>% gather(period, b, na.rm = FALSE, convert = FALSE)
                SPD.mcmc.cooked.b[[1]]        = cbind(SPD.mcmc.cooked.b[[1]], label =rep(paste0("b_",1), length(SPD.mcmc.cooked.b[[1]]$period)))


                #---------------------------------------------
                m=2
                if (isVar[2] == TRUE){   #a1
                  if (isVar[1] == TRUE){
                    col.num[m] = nperiod + i
                  } else {
                    col.num[m] = col.num[m-1] + i
                  }
                } else {
                  if (isVar[1] == TRUE){
                    col.num[m] = nperiod + 1
                  } else {
                    col.num[m] = col.num[m-1] + 1
                  }
                }
                #---------------------------------------------
                m=3
                if (isVar[3] == TRUE){   #c1
                  if (isVar[2] == TRUE){
                    col.num[m] = col.num[m-1] - i + nperiod + i
                  } else {
                    col.num[m] = col.num[m-1] + i
                  }
                } else {
                  if (isVar[2] == TRUE){
                    col.num[m] = col.num[m-1] - i + nperiod + 1
                  } else {
                    col.num[m] = col.num[m-1] + 1
                  }
                }


                #---------------------------------------------------
                if (ncontrols >1){
                  for (ccc in 2:ncontrols) {
                    for (par in 1:3){
                      m= m + 1
                      if (isVar[m] == TRUE){
                        if (isVar[m-1] == TRUE){
                          col.num[m] = col.num[m-1] - 1 + nperiod + 1
                        } else {
                          col.num[m] = col.num[m-1] + 1
                        }
                      } else {
                        if (isVar[m-1] == TRUE){
                          col.num[m] = col.num[m-1] - 1 + nperiod + 1
                        } else {
                          col.num[m] = col.num[m-1] + 1
                        }
                      }
                      
                      if (par ==1){ # "b"
                        if (isVar[m] == TRUE) {
                          SPD.mcmc.cooked.b[[ccc]]        = Results_mcmc.SPD[,  col.num[m]:(col.num[m]+ nperiod - 1)];
                          SPD.maxpost.b[[ccc]]            = Results_summary.SPD[16,  col.num[m]:(col.num[m]+ nperiod - 1)]
                          SPD.box.b[[ccc]]                = unlist(data.frame(t(c(quantile( SPD.mcmc.cooked.b[[ccc]][,1], p = c(0.025, 0.5, 0.975)),
                                                                         mean    =  mean(SPD.mcmc.cooked.b[[ccc]][,1]), 
                                                                         stdev   =  std( SPD.mcmc.cooked.b[[ccc]][,1]), 
                                                                         maxpost =  SPD.maxpost.b[[ccc]][1]))))
                          names(SPD.box.b[[ccc]]) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
                          for (i in 2:nperiod){
                            dftemp[[ccc]] = unlist(data.frame(t(c(quantile(SPD.mcmc.cooked.b[[ccc]][,i], p = c(0.025, 0.5, 0.975)), 
                                                         mean    = mean(SPD.mcmc.cooked.b[[ccc]][,i]), 
                                                         stdev   = std(SPD.mcmc.cooked.b[[ccc]][,i]), 
                                                         maxpost = SPD.maxpost.b[[ccc]][i]))))
                            names(  dftemp[[ccc]]) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
                            SPD.box.b[[ccc]]   = rbind(SPD.box.b[[ccc]],  dftemp[[ccc]])
                          }
                          names(SPD.mcmc.cooked.b[[ccc]]) = seq(1, nperiod, 1)
                          SPD.mcmc.cooked.b[[ccc]]        = SPD.mcmc.cooked.b[[ccc]] %>% gather(period, b, na.rm = FALSE, convert = FALSE)
                          SPD.mcmc.cooked.b[[ccc]]        = cbind(SPD.mcmc.cooked.b[[ccc]],label = rep(paste0("b_",ccc), length(SPD.mcmc.cooked.b[[ccc]]$period)))
                        } else {
                          
                          SPD.mcmc.cooked.b[[ccc]]        = data.frame(cbind(replicate(nperiod, Results_mcmc.SPD[,  col.num[m]])))
                          SPD.maxpost.b[[ccc]]            = Results_summary.SPD[16,  col.num[m]]
                          SPD.box.b[[ccc]]                = data.frame(t(c(quantile( SPD.mcmc.cooked.b[[ccc]][,1], p = c(0.025, 0.5, 0.975)),
                                                                           mean    =  mean(SPD.mcmc.cooked.b[[ccc]][,1]), 
                                                                           stdev   =  std( SPD.mcmc.cooked.b[[ccc]][,1]), 
                                                                           maxpost =  SPD.maxpost.b[[ccc]][1])))
                          names(SPD.box.b[[ccc]]) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
                          names(SPD.mcmc.cooked.b[[ccc]]) = seq(1, nperiod, 1)
                          SPD.mcmc.cooked.b[[ccc]]        = SPD.mcmc.cooked.b[[ccc]] %>% gather(period, b, na.rm = FALSE, convert = FALSE)
                          SPD.mcmc.cooked.b[[ccc]]        = cbind(SPD.mcmc.cooked.b[[ccc]],label = rep(paste0("b_",ccc), length(SPD.mcmc.cooked.b[[ccc]]$period)))
                        }
                      }
                    }
                  }
                }

                
                ## Assign COLORS OF b per control:
                if ((ncontrols <=1)) {
                  colorss = c("black")
                } else if (ncontrols ==2){
                  colorss = c("black", "gray70")
                } else if (ncontrols > 4) {
                  colorss = distinctColorPalette(k = ncontrols, altCol = FALSE, runTsne = FALSE)
                } else {
                  colorss = brewer.pal(n=ncontrols, "Greys")
                }
                labelss = c(paste0("b_",seq(1:ncontrols)))
                
                
                
                # PLOTTING b statistics over periods:
                #####################################
                SPD.bt = ggplot()
                for (ccc in 1:ncontrols){ 
                  # save results of parameter b each control
                  btt = SPD.box.b[[ccc]]
                  write.table(SPD.box.b[[ccc]], paste0(dir.SPD.results,"/bt", ccc, "_df.txt"),  sep ="\t", row.names=FALSE)
      
                  SPD.bt = SPD.bt +
                    geom_boxplot(data = SPD.mcmc.cooked.b[[ccc]],
                                 aes(x= period, y= b, color = label),
                                 outlier.size = 0.1, lwd =0.2)
                  write.table( SPD.mcmc.cooked.b[[ccc]],  paste0(dir.SPD.results,"/b", ccc, "_per_period.csv"),  sep ="\t", row.names=FALSE)
                }
                SPD.bt = SPD.bt +
                  stat_summary(fun = mean, geom="point", shape=16, size=2)+
                  scale_y_continuous(name = TeX("$b$"))+
                  scale_color_manual(name = "Legend",
                                     labels = labelss,
                                     values = colorss)+
                  scale_x_discrete(name = "Period", limits = factor(seq(1,nperiod, 1)) )+
                  theme_bw()+
                  theme( panel.grid.major = element_blank()
                         ,panel.grid.minor = element_blank())
                ggsave(SPD.bt, filename = paste0(dir.SPD.results,"/b_boxplots.png"),
                       bg = "transparent", device = "png", width = 6, height =4, dpi = 400)

                
                
                # PLOTTING deltab (shifts) statistics over periods:
                ###################################################
                message("- Plotting the boxplots of shifts ('deltab') for each period and control 'deltab_boxplots.png'. Wait ... "); flush.console()
                delta.b.SPD =list()
                for (ccc in 1:ncontrols){
                  if (nperiod > 1){
                    delta.b.SPD[[ccc]] = data.frame(shifttime = "1",
                                                    deltab   = SPD.mcmc.cooked.b[[ccc]]$b[which(SPD.mcmc.cooked.b[[ccc]]$period== 2)] -
                                                               SPD.mcmc.cooked.b[[ccc]]$b[which(SPD.mcmc.cooked.b[[ccc]]$period== 1)],
                                                    label    = rep(paste0("Deltab_", 1)))
                    for (tss in 2:(nperiod-1)){
                      delta.b.SPD[[ccc]] = rbind(delta.b.SPD[[ccc]],
                                                 data.frame(
                                                   shifttime = as.character(tss),
                                                   deltab    = SPD.mcmc.cooked.b[[ccc]]$b[which(SPD.mcmc.cooked.b[[ccc]]$period== (tss+1))] -
                                                               SPD.mcmc.cooked.b[[ccc]]$b[which(SPD.mcmc.cooked.b[[ccc]]$period== tss)],
                                                   label     = rep(paste0("Deltab_", ccc))))
                    }
                  } else {
                    delta.b.SPD[[ccc]] = NULL
                  }
                }
                labelss = c(paste0("Deltab_",seq(1:ncontrols)))
                
                # boxplot of deltab:
                SPD.delta.bt = ggplot()
                for (ccc in 1:ncontrols){
                  SPD.delta.bt = SPD.delta.bt +
                    geom_boxplot(data = delta.b.SPD[[ccc]],
                                 aes(x= shifttime, y= deltab, color = label),
                                 outlier.size = 0.1, lwd =0.2)
                  write.table( delta.b.SPD[[ccc]],  paste0(dir.SPD.results,"/Deltab", ccc, "_per_period.csv"),  sep ="\t", row.names=FALSE)

                }
                SPD.delta.bt = SPD.delta.bt +
                  stat_summary(fun = mean, geom="point", shape=16, size=2)+
                  scale_y_continuous(name = TeX("$Shifts \\; \\Delta b \\;[m]$"))+
                  scale_color_manual(name = "Legend",
                                     labels = labelss,
                                     values = colorss)+
                  scale_x_discrete(name = "Rating shift index", limits = factor(seq(1,nperiod-1, 1)) )+
                  theme_bw()+
                  theme( panel.grid.major = element_blank()
                         ,panel.grid.minor = element_blank())
                ggsave(SPD.delta.bt, filename = paste0(dir.SPD.results,"/deltab_boxplots.png"),
                       bg = "transparent", device = "png", width = 6, height =4, dpi = 400)

                
                
                ####################################################################################### 
                # PLot of b1(t) with limni
                ##########################
                if (length(SPD.box.b[[1]])>6){
                  message("- Plotting b1 estimates against stage record 'stage_record_with_b1.pdf'...")
                  names(gaug.periods.df) = c("h","Q","uQ","Period","t","color")               
                
                  #  Read shift times:
                  if (length(data.annotate.gaug.1$t.adj)> 1){
                  shifts.all = data.annotate.gaug.1$t.adj
                  t.shifts.before = c(0, data.annotate.gaug.1$t.adj)
                  if (!is.null(df.limni)){
                    t.shifts.plus = c(data.annotate.gaug.1$t.adj, tail(stage.record$t_limni,1))
                  } else{  
                    t.shifts.plus = c(data.annotate.gaug.1$t.adj, tail(gaug.periods.df$t,1))
                  }
                  }

                  bt1.df = cbind(SPD.box.b[[1]],  data.frame(t.shifts.before = t.shifts.before,   
                                                             t.shifts.plus = t.shifts.plus) )       
                  bt.plot = ggplot()+
                  geom_segment(data = bt1.df, mapping= aes(x    = t.shifts.before, 
                                                         y    = maxpost, 
                                                         xend = t.shifts.plus, 
                                                         yend = maxpost), color = unique(gaug.periods.df$color),  size = 1.1) +
                  geom_rect(data = bt1.df , mapping = aes(xmin  = t.shifts.before, 
                                                        xmax  = t.shifts.plus,
                                                        ymin  = `2.5%`, 
                                                        ymax  = `97.5%`), 
                                                        fill  = unique(gaug.periods.df$color), alpha=0.3) 
                  
                  if (!is.null(df.limni)){
                    bt.plot =  bt.plot + geom_line(data = df.limni, aes(x = t_limni, y= h_limni), color="gray40", size =0.4)
                  }
                  if (!is.null(gaug.periods.df)){
                    bt.plot =  bt.plot +
                    geom_point(aes(x=gaug.periods.df$t, y = gaug.periods.df$h), color= gaug.periods.df$color, size=4)
                  }
                  
                  bt.plot =  bt.plot +
                  scale_y_continuous(name = expression("Stage h [m]"), expand = c(0,0)) + 
                  scale_x_continuous(name = expression("Time [days]"), expand = c(0,0)) +
                  coord_cartesian(clip = 'off')+
                  theme_bw(base_size=20)+
                  theme(axis.text       = element_text(size=15)
                       ,axis.title      = element_text(size=20,face="bold")
                       ,legend.text     = element_text(size=20)
                       ,legend.title    = element_text(size=30)
                       ,legend.key.size = unit(1.5, "cm")
                       ,legend.position ="none"
                       ,panel.grid      = element_blank())+
                
                  geom_segment(aes(x     = data.annotate.gaug.1$t.adj ,   
                                   y     = bt1.df$maxpost[- length(bt1.df$maxpost)] ,
                                   xend  = data.annotate.gaug.1$t.adj, 
                                   yend  = bt1.df$maxpost[-1]),
                                   arrow = arrow(length = unit(0.1, "cm"), ends = "both"), size =0.7, color= "black")
                  pdf(paste0(dir.SPD.results,"/stage_record_with_b1.pdf"), 16, 9 ,useDingbats=F)
                  print(bt.plot)
                  dev.off()
                
                }
                ########################################################################################
                # perform hydrograph for each stable period
                # if (do.final.hydrograph == TRUE){
                #      TO DO YET !!!!
                #   
                # }
                message("
############################################
#              All done !                  #  
############################################")
                
                return(list(err= NULL))       
} 


































#########################################################################################################
bt.SPD <- function(nperiods,                                                                            #
                   df.limni, dir.SPD.results,                                                           #
                   t.shift.for.b,                                                                       #
                   h_G, t_G, color_G,                                                                   #
                   times.uncert,                                                                        #
                   officialShiftsTime,                                                                  #
                   ylimits) {                                                                           #
#########################################################################################################
  #Reading and plotting b1 and b2:
  #*******************************
  # Boxplots of b1 and b2:
  SPD.summary     = read.table(file=paste0(dir.SPD.results, "/Results_Summary.txt"))
  SPD.mcmc.cooked = read.table(file=paste0(dir.SPD.results, "/Results_MCMC_Cooked.txt"), header=TRUE)
  ##################
  # b1:
  SPD.mcmc.cooked.b1        = SPD.mcmc.cooked[,1:nperiods]; 
  names(SPD.mcmc.cooked.b1) = seq(1,nperiods,1)
  SPD.mcmc.cooked.b1        = SPD.mcmc.cooked.b1 %>% gather(period, b1, na.rm = FALSE, convert = FALSE)
  
  bt1.mcmc    = SPD.mcmc.cooked[,1:nperiods]
  bt1.summary = SPD.summary[,1:nperiods]
  bt1.df      = data.frame(t(c(quantile(bt1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                           mean    =  mean(bt1.mcmc[,1]), 
                           stdev   =  std(bt1.mcmc[,1]), 
                           maxpost =  bt1.summary[16,1])))
  names(bt1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
  for (i in 2:nperiods){
    bt1.df = rbind(bt1.df, t(c(quantile(bt1.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                               mean    = mean(bt1.mcmc[,i]), 
                               stdev   = std(bt1.mcmc[,i]), 
                               maxpost = bt1.summary[16,i]) ))
  }
  write.table(bt1.df, paste0(dir.SPD.results,"/bt1_df.txt"), 
              sep ="\t", row.names=FALSE)
  
  
  
  ##################
  SPD.mcmc.cooked.b2        = SPD.mcmc.cooked[,(nperiods+3 ):(2*nperiods+2)]; 
  names(SPD.mcmc.cooked.b2) = seq(1,nperiods,1)
  SPD.mcmc.cooked.b2        = SPD.mcmc.cooked.b2 %>% gather(period, b2, na.rm = FALSE, convert = FALSE)
  
  bt2.mcmc    = SPD.mcmc.cooked[,(nperiods+3 ):(2*nperiods+2)]; 
  bt2.summary = SPD.summary[,(nperiods+3 ):(2*nperiods+2)]; 
  bt2.df      = data.frame( t(c(quantile(bt2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                            mean     =    mean(bt2.mcmc[,1]), 
                            stdev    =   std(bt2.mcmc[,1]), 
                            maxpost  = bt2.summary[16,1])))
  names(bt2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
  for (i in 2:nperiods){
    bt2.df = rbind(bt2.df, t(c(quantile(bt2.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                               mean    = mean(bt2.mcmc[,i]), 
                               stdev   = std(bt2.mcmc[,i]), 
                               maxpost = bt2.summary[16,i]) ))
  }
  write.table(bt2.df, paste0(dir.SPD.results,"/bt2_df.txt"), 
              sep ="\t", row.names=FALSE)
  
  #boxplots of b1 and b2 parameters for each period:
  bt.boxplot   = ggplot()+
                 geom_boxplot(data=SPD.mcmc.cooked.b1, aes(x= period, y= b1, color="b1"), outlier.size = 0.1, lwd =0.2)+
                 geom_boxplot(data=SPD.mcmc.cooked.b2, aes(x= period, y= b2, color="b2"), outlier.size = 0.1, lwd =0.2)+
                 stat_summary(fun = mean, geom="point", shape=16, size=2)+
                 scale_y_continuous(name = TeX("$b_1, b_2$"))+
                 scale_color_manual(name = "Legend", labels=c("b1", "b2") , values = c("blue", "red"))+
                 scale_x_discrete(name = "Period", limits = seq(1,nperiods, 1) )+
                 theme_bw()
                 ggsave(bt.boxplot, filename = paste0(dir.SPD.results,"/bt_boxplots.png"), 
                        bg = "transparent", device = "png", width = 6, height =4, dpi = 400)
                 # df$asDate = as.Date(dtm[,1], "%d.%m.%Y")   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DATE FORMAT !!!
               
                 
                 
                 
   # Read shift times:
                shifts.all = t.shift.for.b$t.adj
                if (times.uncert ==TRUE) {   # /!\ this has to be changed !!!!!!!!!!!!!!
                     t.shifts.before = c(0, t.shift.for.b$t.adj)
                     if (!is.null(df.limni)){
                        t.shifts.plus = c(t.shift.for.b$t.adj, tail(df.limni$t_limni,1))
                     } else{  
                        t.shifts.plus = c(t.shift.for.b$t.adj, tail(t_G,1))
                     }
                     bt2.df = cbind(bt2.df, 
                                    t.shifts.before = t.shifts.before , 
                                    t.shifts.plus   = t.shifts.plus)

                     bt1.df = cbind(bt1.df, 
                                    t.shifts.before =    t.shifts.before , 
                                    t.shifts.plus   =    t.shifts.plus)
                } else {
                      t.shifts        = shifts[-length(shifts)]
                      t.shifts.before = c(0,t.shifts)
                      if (!is.null(df.limni)){
                        t.shifts.plus   = sort(c(t.shifts, tail(df.limni$t_limni,1)))
                      } else{  
                        t.shifts.plus   = sort(c(t.shifts, tail(t_G,1)))
                      }
                      
                }
                
                
     # PLot of b(t) with limni
     #########################
     bt.plot = ggplot()
                          if (times.uncert ==TRUE) {
                               #geom_line(data = df.bt, aes(x=t.shifts.before, y =bt.MAP), color = "red", size = 1 )+
                               bt.plot = bt.plot +
                               # geom_rect(mapping= aes(xmin= shifts$t2, xmax=shifts$t97 ,ymin=-Inf, ymax=Inf), 
                               #   fill="blue", alpha=0.1) +
                               #geom_vline(aes(xintercept = shifts$tMAP), color = "blue", lwd =0.5, linetype = "solid")+
                                 
                               #geom_vline(aes(xintercept = bt2.df$t.shifts.before[-1]), color = "red", lwd =0.3, linetype = "dashed")+
                               # geom_segment(data = bt2.df, mapping= aes(x =t.shifts.before , 
                               #                                          y = maxpost, 
                               #                                          xend = t.shifts.plus, 
                               #                                          yend = maxpost),color = "red", size = 1) +
                               # geom_rect(data = bt2.df , mapping = aes(xmin= t.shifts.before, 
                               #                                         xmax=t.shifts.plus,
                               #                                         ymin= `2.5%`, 
                               #                                         ymax= `97.5%`), 
                               #           fill="red", alpha=0.3) +
                               geom_segment(data = bt1.df, mapping= aes(x    = t.shifts.before, 
                                                                        y    = maxpost, 
                                                                        xend = t.shifts.plus, 
                                                                        yend = maxpost), color = unique(color_G),  size = 1.2) 
                               # geom_rect(data = bt1.df , mapping = aes(xmin= t.shifts.before, 
                               #                                           xmax=t.shifts.plus,
                               #                                           ymin= `2.5%`, 
                               #                                           ymax= `97.5%`), 
                               #           fill=color_G, alpha=0.3) +
                               
                               if (!is.null(df.limni)){
                                   bt.plot =  bt.plot + geom_line(data = df.limni, aes(x = t_limni, y= h_limni), color="gray40", size =0.3)
                               }
                               bt.plot =  bt.plot +
                               geom_point(aes(x=t_G, y = h_G), color= color_G, size=4) +
                               scale_y_continuous(name = expression("Stage h [m]"), expand = c(0,0)) + 
                               scale_x_continuous(name = expression("Time [days]"), expand = c(0,0))+
                               coord_cartesian(clip = 'off')+
                               theme_bw(base_size=20)+
                               theme(axis.text=element_text(size=15),
                                     axis.title=element_text(size=20,face="bold")
                                     ,legend.text=element_text(size=20),legend.title=element_text(size=30)
                                     ,legend.key.size=unit(1.5, "cm"),legend.position="none"
                                     ,panel.grid = element_blank())+
                               geom_segment(aes(x    = t.shift.for.b$t.adj ,   
                                                y    = bt1.df$maxpost[- length(bt1.df$maxpost)] ,
                                                xend = t.shift.for.b$t.adj, 
                                                yend = bt1.df$maxpost[-1]),
                                            arrow = arrow(length = unit(0.1, "cm"), ends = "both"), size =0.7, color= "black")
                               # geom_point(aes(x = t.shift.for.b$t.adj, 
                               #                y = df.limni$h_limni[index.hmax]),
                               #              size = 10, color="red", pch=21, fill="NA")
                          } else {
                          #########
                               if (!is.null(df.limni)){
                                   bt.plot =  bt.plot + geom_line(data = df.limni, aes(x = t_limni, y= h_limni), color="gray40", size =0.3)
                                 }
                               bt.plot = bt.plot +                                 
                               geom_vline(aes(xintercept = shifts.all), color = "blue", lwd =0.7, linetype = "solid") +
                               geom_segment(mapping = aes(x =t.shifts.before , y = bt.MAP, 
                                                         xend = t.shifts.plus, yend = bt.MAP), color = "red", size = 1) +
                               geom_rect(mapping = aes(xmin= t.shifts.before, xmax=t.shifts.plus, ymin=bt.q10,
                                         ymax= bt.q90), fill="red", alpha=0.3) +
                               scale_y_continuous(name = expression("Stage h [m]"), limits =ylimits, expand = c(0,0)) + 
                               scale_x_continuous(name = expression("Time [days]"), expand = c(0,0))+
                               coord_cartesian(clip = 'off')+
                               geom_point(aes(x=t_G, y = h_G), color=color_G, size=4)+
                               theme_bw(base_size=20)+
                               theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold")
                                     ,legend.text=element_text(size=20),legend.title=element_text(size=30)
                                     ,legend.key.size=unit(1.5, "cm"),legend.position="none"
                                     ,panel.grid = element_blank())
                          }
                          # ggsave(bt.plot, filename = paste0(dir.SPD.results,"/bt.png"), 
                          #        width = 16, height =9, dpi = 400)
                          pdf(paste0(dir.SPD.results,"/bt.pdf"), 16, 9 ,useDingbats=F)
                          print(bt.plot)
                          dev.off()
                          
      # return the parameter b statistics
      return(list(bt1.df        = bt1.df, 
                  bt2.df        = bt2.df,
                  t.shift.for.b = t.shift.for.b))
}
  
  





































