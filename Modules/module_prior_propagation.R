######################################################################################################################################
propagation <- function(list.prior.RC,
                        list.prior.changes,
                        start,   # starting values vector 
                        changes.method,
                        n.parvar,
                        isVar, 
                        global.change,  
                        local.change,  
                        width.change,
                        ncontrols) {
######################################################################################################################################
  # Initialisation:  
  nsim      =  length(list.prior.RC$b[[1]])
  print(nsim)
  chang     =  which(!unlist(lapply(list.prior.changes, is.null)))[1]
  if (!is.na(chang)) {
     nperiod   =  NCOL(list.prior.changes[[chang]]) + 1   
  } else {
     nperiod = 1
  }
  # initialise parameters that will be computed by propagation
  a = NULL
  for (i in 1:ncontrols) {
      a[[i]] = vector(mode   = 'double', length = nsim)
  }
  nparvar = n.parvar
  all_parvar = NULL; 
  for (j in 1:nparvar) {
      all_parvar[[j]] = matrix(NA, nsim, nperiod)
  }
  

  
  # Monte-Carlo propagation:
  ##############################################################################
  for (i in 1:nsim){
  ##############################################################################
    # parameters a for all simulations:
    varr=0
    for (n in 1:ncontrols) {
    # b for all periods
      if (isVar[n*3- 2] == T){
          varr = varr +1
          all_parvar[[varr]][i,1] = list.prior.RC$b[[n]][i]
          if (nperiod >1) {
          if (n==1){ #  for lowest control only:
             if (!is.null(list.prior.changes$d.l)){
                 all_parvar[[varr]][i,2:nperiod] = list.prior.RC$b[[n]][i] - cumsum(list.prior.changes$d.g[i,] + list.prior.changes$d.l[i,])
             } else {
                all_parvar[[varr]][i,2:nperiod]  = list.prior.RC$b[[n]][i] - cumsum(list.prior.changes$d.g[i,])
             }
          } else {
             all_parvar[[varr]][i,2:nperiod] = list.prior.RC$b[[n]][i] - cumsum(list.prior.changes$d.g[i,])
          }
          }
      }
      if (isVar[n*3- 1] == T){
        varr = varr +1
        # then param "a" is dynamic!!
        if (control.type[n] == "rect.weir") {
          all_parvar[[varr]][i,1] = list.prior.RC$Cr[[n]][i]*list.prior.RC$Bw[[n]][i]*sqrt(2*list.prior.RC$g[[n]][i])
          if (nperiod >1) {
          all_parvar[[varr]][i,2:nperiod] =  list.prior.RC$Cr[[n]][i]*list.prior.RC$Bw[[n]][i]*sqrt(2*list.prior.RC$g[[n]][i])/cumprod(list.prior.changes$d.B[i,])
          }
        } else if (control.type[n] == "rect.channel") {
          all_parvar[[varr]][i,1] =  list.prior.RC$KS[[n]][i]*list.prior.RC$Bc[[n]][i]*sqrt(list.prior.RC$S0[[n]][i])
          if (nperiod >1) {
          all_parvar[[varr]][i,2:nperiod] =  list.prior.RC$KS[[n]][i]*list.prior.RC$Bc[[n]][i]*sqrt(list.prior.RC$S0[[n]][i])/cumprod(list.prior.changes$d.B[i,])
          }
        }
      } else {
        if (control.type[n] == "rect.weir") {
          a[[n]][i] = list.prior.RC$Cr[[n]][i]*list.prior.RC$Bw[[n]][i]*sqrt(2*list.prior.RC$g[[n]][i])
        } else if (control.type[n] == "rect.channel") {
          a[[n]][i] = list.prior.RC$KS[[n]][i]*list.prior.RC$Bc[[n]][i]*sqrt(list.prior.RC$S0[[n]][i])
        }
      }
    }
    ###############
    # Instead, parameter c is always assumed stable over periods !!!
    ###############
}
  
    
    
    
    
    
  
  
  
  
  
  
  
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  all_parvar.s = NULL
  a.s          = NULL
  varr         = 0
  for (n in 1:ncontrols) {
      # b for all periods:
      if (isVar[n*3- 2] == T){
        varr=varr+1
        all_parvar.s[[varr]] = vector(mode='double',length=nperiod)
        all_parvar.s[[varr]][1] = start[paste0("b",n)]            # eval(parse(text=paste("b",1,sep="")))
        if (nperiod >1) {
        if (global.change == TRUE) {
           if (local.change == TRUE) {
               all_parvar.s[[varr]][2:nperiod] = unlist(start[paste0("b",n)])  - cumsum(unlist(start[paste0("dg",n)]) + unlist(start[paste0("dl",n)]))
           } else {
               all_parvar.s[[varr]][2:nperiod] = unlist(start[paste0("b",n)])  - cumsum(unlist(start[paste0("dg",n)]))
           }
        } else {
           if (local.change == TRUE) {
               all_parvar.s[[varr]][2:nperiod] = unlist(start[paste0("b",n)]) - cumsum(unlist(start[paste0("dl",n)]))
           } else {
               all_parvar.s[[varr]][2:nperiod] = start[paste0("b",n)]
           }
        }
        }
      }
    
      # a for all periods:
      if (isVar[n*3- 1] == T){
         varr=varr+1
         # param "a" is dynamic!!
         all_parvar.s[[varr]] = vector(mode='double', length=nperiod)
         if (control.type[n] == "rect.weir") {
           all_parvar.s[[varr]][1]         = unlist(start[[paste0("Cr",n)]]) *unlist(start[[paste0("Bw",n)]])* sqrt(2*unlist(start[[paste0("g",n)]]))
           if (nperiod >1) {
           all_parvar.s[[varr]][2:nperiod] = unlist(start[[paste0("Cr",n)]]) *unlist(start[[paste0("Bw",n)]])* sqrt(2*unlist(start[[paste0("g",n)]]))/cumprod(unlist(start[[paste0("dB",n)]]))
           }
         } else if (control.type[n] == "rect.channel") {
           all_parvar.s[[varr]][1]         = unlist(start[[paste0("KS",n)]]) *unlist(start[[paste0("Bc",n)]])* sqrt(unlist(start[[paste0("S0",n)]]))
           if (nperiod >1) {
           all_parvar.s[[varr]][2:nperiod] = unlist(start[[paste0("KS",n)]]) *unlist(start[[paste0("Bc",n)]])* sqrt(unlist(start[[paste0("S0",n)]]))/cumprod(unlist(start[[paste0("dB",n)]]))
           }
         }
      } else {
        if (control.type[n] == "rect.weir") {
          a.s[[n]] = unlist(start[[paste0("Cr",n)]]) *unlist(start[[paste0("Bw",n)]])* sqrt(2*unlist(start[[paste0("g",n)]]))
        } else if (control.type[n] == "rect.channel") {
          a.s[[n]] = unlist(start[[paste0("KS",n)]]) *unlist(start[[paste0("Bc",n)]])* sqrt(unlist(start[[paste0("S0",n)]]))
        }
      }
  }
    
  
  
  
  
  
  
  
     
     
  sim    = list();
  start2 = list();
  varr   = 0
  for (n in 1:ncontrols){       
     if (isVar[n + (n-1)*2] == TRUE) {
         varr = varr+1
         sim[[paste0("b",n)]]    =  all_parvar[[varr]]
         start2[[paste0("b",n)]] =  all_parvar.s[[varr]]
     } else {
         sim[[paste0("b",n)]]    =  list.prior.RC$b[[n]]
         start2[[paste0("b",n)]] =  start[paste0("b",n)]
     }
     
     if (isVar[n +1 + (n-1)*2] == TRUE) {
         varr = varr+1
         sim[[paste0("a",n)]]    =  all_parvar[[varr]]
         start2[[paste0("a",n)]] =  all_parvar.s[[varr]]
     } else {
         sim[[paste0("a",n)]]    =  a[[n]]
         start2[[paste0("a",n)]] =  a.s[n]
     }   
     
     if (isVar[n + 2+ (n-1)*2] == TRUE) {
         varr = varr+1
         sim[[paste0("c",n)]]    =  all_parvar[[varr]]
         start2[[paste0("c",n)]] =  all_parvar.s[[varr]]
     } else {
         sim[[paste0("c",n)]]    =  list.prior.RC$c[[n]]
         start2[[paste0("c",n)]] =  start[paste0("c",n)]
     }   
  }
  
  # return results:
  ###############################################################################
    return(list(sim     = sim,
                start2  = start2,
                nperiod = nperiod))
  ###############################################################################
}

  




































######################################################################################
fit<-function(sim, # Monte Carlo simulations produced by function propagate_XXX$sim
              margins, #=rep('Gaussian',NCOL(sim)), # marginal distributions
              names, #=paste('par',1:NCOL(sim),sep=''), # parameter names
              ncol, rowmax, # graphical parameters (here 7*4 parameters shown per figure)
              file.save
){
########################################################################################
  #^* PURPOSE: fit a multivariate prior distribution to Monte-Carlo samples
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 16/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [numeric matrix] sim, Monte Carlo simulations (nsim*npar)
  #^*    2. [character vector] margins, marginal distributions (npar)
  #^*    3. [character vector] names, parameter names (npar)
  #^*    4. [integer] ncol, rowmax: graphical parameters controlling the number of subplots per figure
  #^* OUT
  #^*    A list containing:
  #^*    1. [character vector] margins, marginal distributions (same as the one used as input)
  #^*    2. [list of vectors] priorpar, prior parameters for each marginal prior    
  #^*    3. [numeric vector] start, starting points (e.g. prior mean)    
  #^*    4. [numeric matrix] cor, prior correlation (in Gaussian space)    
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* COMMENTS: only normal and lognormal margins available for the moment
  #^******************************************************************************
  #^* TO DO: can probably be generalized to any marginal distribution
  #^******************************************************************************
  
  # check sizes match and if so, initialize
  sim.df = unlist(sim[[1]])
  nsim   = 10000
  for (ll in 2:length(sim)) {
    sim.df = cbind(sim.df, unlist(sim[[ll]]))
  }
  
  n = ncol(sim.df)
  if(length(margins)!=n){
      message('Fatal:size mismatch [sim,margins]'); return(NA)
  }
  priorpar  = vector("list",n)
  transform = matrix(NA, nsim, n)
  
  # start computations for each margin:
  for(i in 1:n){
    # marginal prior parameters:
    priorpar[[i]]= switch(margins[i],
                          'Gaussian'  = {c(mean(sim.df[,i]),       sd(sim.df[,i]))},
                          'LogNormal' = {c(mean(log(sim.df[,i])),  sd(log(sim.df[,i])))},
                          NA)
    # Transform simulations into Gaussian space to estimate the correlation of the Gaussian copula:
    transform[,i]= switch(margins[i],
                          'Gaussian'  = {qnorm(pnorm(sim.df[,i],    mean = priorpar[[i]][1],     sd = priorpar[[i]][2]))},
                          'LogNormal' = {qnorm(plnorm(sim.df[,i],   meanlog = priorpar[[i]][1],  sdlog = priorpar[[i]][2]))},
                          NA)
  }
  # prior correlation:
  corel=cor(transform)
  # Verify that marginal fit is acceptable graphically:
  nrow=min(ceiling(n/ncol),rowmax)
  
  
  
  
  pdf(paste0(file.save,"/priors_density.pdf"), 13, 16, useDingbats = F)
  for(i in 1:n){
    if(i%%(nrow*ncol)==1) {par(mfrow=c(nrow, ncol))}
    x = seq(min(sim.df[,i]), 
            max(sim.df[,i]),100)
    y = switch(margins[i],
              'Gaussian'  = {dnorm(x,  mean    = priorpar[[i]][1], sd    = priorpar[[i]][2])},
              'LogNormal' = {dlnorm(x, meanlog = priorpar[[i]][1], sdlog = priorpar[[i]][2])},
             NA)
    hist(sim.df[,i], breaks=50, freq=F, main=NULL, xlab=names[i], ylab='')
      
    lines(x,y,col='red')
    }
  dev.off()

  
  
  # show correlation matrix
  pdf(paste0(file.save,"/matrix_correlation.pdf"), 6, 6, useDingbats = F)
  image(1:n, 1:n, corel)
  dev.off()

  
  
  #####################################################################
  return(list(margins  = margins, 
              priorpar = priorpar, 
              cor      = corel))
  #####################################################################
}




































###############################################################################################
writeConfigFilesSPD <- function(prior, # prior object produced by function 'fit'
                           start, # starting vector produced by function propagate_XXX$start
                           ncontrol, # number of hydraulic controls
                           nperiod,
                           colperiod, # number of periods and column where it is stored in data file
                           isVar, # flag varying (period-specific) parameters
                           model, # Model configuration file
                           correl, #prior correlation file
                           names.bac, # base name of the 3 parameters for one control
                           directory
){
################################################################################################
  #^* PURPOSE: write configuration files used by BaM
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 17/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [list] prior, object produced by function 'fit'
  #^*    2. [numeric vector] start, starting vector produced by function propagate_XXX$start
  #^*    3. [integer] ncontrol, number of hydraulic controls
  #^*    4. [integer] nperiod, number of periods 
  #^*    5. [integer] colperiod, column where period is stored in data file
  #^*    6. [logical vector] isVar, flag varying (period-specific) parameters
  #^*    7. [character] model, Model configuration file
  #^*    8. [character] correl, prior correlation file
  #^*    9. [character vector]  base name of the 3 parameters for one control
  #^* OUT
  #^*    Nothing - just write files
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  
  #---------------------------------------
  # Config_Model file
  #---------------------------------------
  # Start assembling values to be written in file and associated comments
  comment=c(
    'Model ID',
    'nX: number of input variables',
    'nY: number of output variables',
    'nPar: number of parameters theta'
  )
  val = list('BaRatinBAC', # 'Model ID',
             1,            # 'nX: number of input variables',
             1,            # 'nY: number of output variables',
             3*ncontrol)   # 'nPar: number of parameters theta'
  # specification for each parameter
  k=0; m=1
  for(i in 1:ncontrol){ # loop on each hydraulic control
    for(j in 1:length(names.bac)){ # loop on each b-a-c parameter of the control
      k=k+1
      # comments
      comment=c(comment,
                'Parameter name',
                'Initial guess',
                'Prior distribution',
                'Prior parameters')
      # parname
      pname = paste0(names.bac[j], i)
      val   = c(val, pname)
      if(!isVar[k]){ # Static parameter
          # starting point
          val=c(val, start[[paste0(names.bac[j], i)]])
          # prior distribution
           val=c(val, prior$margins[m])
          # prior parameters
          val=c(val, list(prior$priorpar[[m]]))
          m=m+1      
          
      } else { # period-specific parameter
          # starting point
          val=c(val, -666.666) # unused, can use a dummy number
          # prior distribution: use 'VAR'
          val=c(val,'VAR')
          # Name of the configuration file for this varying parameter
          fname = paste0('Config_', pname, '_VAR.txt')
          val   = c(val, fname)
          # Write this additional configuration file
          comment0 = c(
          'Number of periods',
          'Column where period is written in data file'
          )
          val0 = list(nperiod, colperiod)
          for(ii in 1:nperiod){ # one parameter block per period
             # comments
             comment0=c(comment0,
                       'Parameter name',
                       'Initial guess',
                       'Prior distribution',
                       'Prior parameters')
             # par name
             val0=c(val0, paste0(pname, '_', ii))
             # starting point
             val0=c(val0, start[[paste0(names.bac[j], i)]][ii])
             # prior distribution
             val0=c(val0, prior$margins[m])
             # prior parameters
             val0=c(val0, list(prior$priorpar[[m]]))
             m=m+1
        }
        writeConfig(val0, comment0, directory, fname)
      }       
    }
  }  
  
  writeConfig(val,
              comment,
              directory,
              model)
  
  #---------------------------------------
  # Prior correlations
  #---------------------------------------
  # add 2 lines/columns, uncorrelated with all others, for structural uncertainty parameters gamma1 and gamma2
  n=NCOL(prior$cor)
  foo=matrix(0,n+2,n+2)
  foo[1:n,1:n]=prior$cor
  foo[n+1,n+1]=1
  foo[n+2,n+2]=1
  # write resulting matrix
  write(foo,file=paste(directory,"/",correl,sep=""),ncolumns=n+2)
}




























##################################################################################
writeConfig<-function(val,comment,dir,fname,addQuote=T){
##################################################################################
  #^* PURPOSE: Generic function to write a configuration file
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 12/10/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [list] val, list of values to be written in the configuration file 
  #^*    2. [string vector] comment, comments (same length as val) 
  #^*    3. [string] dir, destination directory 
  #^*    4. [string] fname, file name 
  #^*    5. [logical] addQuote, add double quotes to any string encountered in val? 
  #^* OUT
  #^*    1. nothing  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* TO DO: 
  #^******************************************************************************
  #^* COMMENTS: Overwrites any pre-existing file. Directory created if needed
  #^******************************************************************************
  
  n=length(val)
  txt=vector("character",n)
  for (i in 1:n){
    foo=val[[i]]
    # add double quotes to strings
    if(is.character(foo) & addQuote){foo=addQuotes(foo)}
    # transform R logicals to Fortran logicals
    if(is.logical(foo)){foo=paste(ifelse(foo,'.true.','.false.'))}
    # stitch values and comment
    if(is.null(comment[i])){
      txt[i]=paste(foo,collapse=',')
    } else {
      txt[i]=paste(paste(foo,collapse=','),'!',comment[i])
    }
  }
  
  if(!dir.exists(dir)){dir.create(dir,recursive=T)}
  file=file.path(dir,fname)
  
  write.table(matrix(txt, ncol = 1, nrow = length(txt)), file = file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}




























##################################################################################
addQuotes<-function(txt){
##################################################################################
  #^* PURPOSE: add double quotes to a string or to each element of a string vector
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 12/10/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [string vector] txt
  #^* OUT
  #^*    1. [string vector] double-quoted string vector 
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* TO DO: 
  #^******************************************************************************
  #^* COMMENTS: 
  #^******************************************************************************
  
  n=length(txt)
  out=txt
  for(i in 1:n){
    out[i]=paste('"',txt[i],'"',sep='')
  }
  return(out)
}






