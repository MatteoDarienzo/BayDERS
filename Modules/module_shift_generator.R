generation.datasets <- function(dir.general,       # dir. where the general input options are
                                dir.sim,           # dir. where the generated data-set will be saved
                                i) {
################################################################################################
  # generate the synthtetic case study "i":
  for (syntheticc in 1:max.try.sim) {
    print(syntheticc)
    iii <- try(generateShift(directory    = dir,      #directories
                             t0           = t0, 
                             tmax         = tmax, 
                             M            = M, # configuration matrix
                             RCpar        = RCpar, #Rating curve parameters(b,a,c).
                             ncontrols    = ncontrols,
                             rate.shift   = rate.shift,   #below one RC shift every 5 years on average.
                             deltab       = deltab,  #Shifts magnitude sdev
                             uQ.LF        = uQ.LF, 
                             uQ.HF        = uQ.HF, #Relative gauging uncertainty in Qobs (standard deviation)
                             rate.gauging = rate.gauging, # Gauging rate - below 5 gaugings per year on average
                             meanQ        = meanQ,
                             stdevQ       = stdevQ))
    if (inherits(iii, "try-error")) {
      print("error in building the RC")
      next
    }
    
    if (length(iii[[2]]) > n.shifts.max) {    #maximum number of true shift times
      print(c("Too many shifts", length(iii[[2]]))) 
      next
    } else if ( abs(iii[[3]][1] - iii[[3]][length(iii[[3]])]) < deltab.min ) {   # minimum overall shift magnitude
      print(c("shifts too small", abs(iii[[3]][1] - iii[[3]][length(iii[[3]])])))
      next
    } else if (length(iii[[1]]$Q) > n.gaug.max ) { #maximum number of gaugings
      print("Too many gaugings")
      next
    } else {
      print("RC correctly defined") 
      break
    }
  }
  
}   





  
  





  
  





#################################################################################################
generateShift <- function(
                           #directories
                           dir.general,
                           dir.sim,
                           #sim
                           i) {
#################################################################################################
  # read input options:
  dir.modules = paste0(dir_code, "/Modules")
  source(paste0(dir.modules,"/module_read_input.r"))  # read all modules and input options !!
  
  

  
  
  
  # "Flow duration curve", aka quantile function:
  #**********************************************
  # Used to generate the streamflow of the gauging.
  # Here a lognormal distribution is used, could use something else but I'm not sure it's 
  # important here.
  # However might good to fit the parameters so that it looks like one of the case studies 
  # in the paper.
  Qquantile<-function(p){qlnorm(p=p, meanlog=log(meanQ),sdlog=stdevQ)}
  ########################################################################################

  # Function used to generate a probability between 0 and 1, that will be 
  # transformed into discharge using Qquantile above:
  #*****************************************************************
  #generateP<-function(){runif(n=1)}
  # NOTE: the uniform distribution above correspond to 
  # no gauging strategy (gaugings are made 'randomly').
  # We may mimic gauging strategies (favoring low/high flows) by 
  # using a beta, rather than a uniform, distribution.
  # I'm not sure it's a good idea to go into this sort of refinement, 
  # but if you want to try:
  
  # with beta distrib:
  # generateP<-function(){rbeta(1,0.9,0.65)} # favor high flows but still gauge everywhere
  # generateP<-function(){rbeta(1,0.9,0.1)} # strongly favor high flows
  ###################################################################################
  generateP <- function(){rbeta(1, 0.1, 0.9)} # strongly favor low flows
  ###################################################################################
  # generate gauging times
  time.gauging = generateTimes(t0,tmax,rate.gauging)
  n            = length(time.gauging)
  
  # generate shifts (times + amplitude)
  time.shift   = generateTimes(t0,
                               tmax,
                               rate.shift)
  np           = length(time.shift)
  csum         = cumsum(rnorm(n    = np+1,
                              mean = 0,
                              sd   = deltab)) #shift on parameter b
  b1s=RCpar[1]+ csum
  if (ncontrols >1) {
    b2s = RCpar[4] +csum
  }
  
  # generate gaugings
  gaugings=data.frame(time=numeric(),
                      period=integer(),
                      h=numeric(),
                      Q=numeric())
  
  for (i in 1:n){
    gtime=time.gauging[i]
    # determine RC period associated with gauging
    period=1
    for(j in 1:length(time.shift)){
      if(gtime<=time.shift[j]){break} else {period=period+1}
    }
    # generate Qtrue
    Qtrue=Qquantile(generateP())
    # compute H from inverse RC
    param=RCpar; 
    param[1]= b1s[period] # use the b1 of the relevant period
    
    if (ncontrols ==1) {
      # inverse function of RC:
      h=RCinv1(Q= Qtrue, theta=param)
      
    } else {
      param[4]= b2s[period] # use the b2 of the relevant period
      
      ktrans = 0; 
      Qtrans =0;
      ktrans[1] = param[1] ; Qtrans[1] = 0
      ktrans[3] = param[7] ; Qtrans[3] = param[5]*(ktrans[3]-param[4])^param[6]
      ktransold = param[1]
      toll = 1 ; it =1
      while (toll > 0.00001) { 
        ktrans[2] = (param[2]/param[5])^(1/param[6])*(ktransold - param[1])^(param[3]/param[6]) + param[4]
        toll = abs(ktrans[2]-ktransold)
        ktransold = ktrans[2]
        it = it +1 
      }
      Qtrans[2] = param[2]*(ktrans[2]-param[1])^param[3]
      h = RCinv3(Q= Qtrue, theta=param, Qtrans, ktrans, ncontrols = ncontrols)
    }
    
    Qobs = 0; uQ.BaM = 0;
    
    # add noise to Qtrue to generate Qobs
    for (disc in 1:length(Qtrue)) {
      if (Qtrue[disc] <= 100) {
        Qobs[disc]=Qtrue[disc]+rnorm(n=1,mean=0,sd=(Qtrue[disc]*uQ.LF/100))
        uQ.BaM[disc] = uQ.LF*Qobs[disc]/100
      } else {
        Qobs[disc]=Qtrue[disc]+rnorm(n=1,mean=0,sd=(Qtrue[disc]*uQ.HF/100))
        uQ.BaM[disc] = uQ.HF*Qobs[disc]/100
      }
    }
    # add to gaugings dataset
    gaugings=rbind(gaugings,
                   data.frame(time=gtime,period=period,h=h,Q=Qobs, uQ=uQ.BaM))
  }
  
  
  # save gaugings in a file:
  write.table(gaugings, file = paste0(directory,"/synthet_gaugings.csv"), col.names = c("time", "period", "h", "Q", "uQ"), sep=";")
  # save gaugings in a file:
  write.table(time.shift, file = paste0(directory,"/synthet_shifts.csv"), sep=";")
  
  
  
  
  # plot gaugings:
  require(ggplot2)
  g=ggplot(gaugings,aes(x=h, y=Q, fill=as.factor(period)))
  g=g+ geom_errorbar(aes(x=gaugings$h, ymin =gaugings$Q-2*gaugings$uQ , 
                         ymax =gaugings$Q+2*gaugings$uQ), width=.04, size = 0.3)
  g=g+geom_point(pch=21,size=3.5) + 
    scale_fill_brewer(palette="Spectral")
  g=g+scale_y_log10(limits=c(min(gaugings$Q-2*gaugings$uQ), max(gaugings$Q+2*gaugings$uQ))) + 
    #Qquantile(0.00000000000001),Qquantile(0.99999999999999))) +
    labs(fill="Period")+
    ### Labels
    xlab(expression(paste("Stage h [m]",sep="")))+
    ylab(expression(paste("Discharge Q [",m^3,".",s^-1,"]",sep=""))) +
    theme_bw(base_size=20)+ 
    theme( axis.text        = element_text(size=20)
          ,axis.title       = element_text(size=30,face="bold")
          ,panel.grid.major = element_line(size=1.2)
          ,panel.grid.minor = element_line(size=0.8)
          ,legend.text      = element_text(size=20)
          ,legend.title     = element_text(size=30)
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "right")
  ggsave(g, filename=paste0(directory,"/RClog_synthetic.png"), device = "png", width = 16, height =8,
         dpi = 400, units = "in")
  
  
  
  
  return(list(gaugings, time.shift, b1s ))
}






########################################################################
RCinv1 <-function(Q, theta){    # theta = b1,a1,c1, ... 
  ########################################################################
  
  h=theta[1]+(Q/theta[2])^(1/theta[3])
  return(h)
}




########################################################################
RCinv3 <-function(Q, theta, Qtrans,ktrans,  ncontrols ){    # theta = b1,a1,c1, ... 
  ########################################################################
  # require(spuRs)
  # inverse rating curve:
  # complex RC: three controls
  if ((Q> Qtrans[1]) & (Q <= Qtrans[2])) {
    h=theta[1]+(Q/theta[2])^(1/theta[3])
    
  } else if ((Q > Qtrans[2]) & (Q <= Qtrans[3])) {
    h=theta[4]+(Q/theta[5])^(1/theta[6])
    
  } else if (Q > Qtrans[3]) {
    #print(c("k3=",ktrans[3], "Q3=",Qtrans[3] , "Q=",Q ))
    hold = ktrans[3]
    toll = 10 ; it =1
    while ((toll > 0.001) & (it<100)) { 
      #print(paste0("hold=",hold))
      h = theta[4] + (Q/theta[5] - theta[8]/theta[5]*( hold - theta[7])^(theta[9]))^(1/theta[6])
      #print(c("h=",h))
      if ((h <= ktrans[3]) | (is.nan(h)== TRUE) ) {
        h = ktrans[3] + 0.0001
      }
      
      toll = abs(h-hold)
      hold = h
      it = it +1 
      
    }
  } 
  #print(h)
  
  # h = fixedpoint(
  #     ftn = ftn1, 
  #     x0 = theta[4]+((Qtrans[3]+10)/theta[5])^(1/theta[6]), 
  #     tol = 1e-6,
  #     Q=Q, theta=theta)
  
  # h = uniroot(fun.NR, 
  #             interval=c(theta[4]+((Qtrans[3])/theta[5])^(1/theta[6]),
  #                      10),
  #             theta =theta, Q=Q,
  #             tol = 0.0001)
  return(h)
}






###################################################
# Function for the third control implicit problem :
###################################################
ftn1 <- function(h, Q, theta) return(
  h = theta[4] + (Q/theta[5] - theta[8]/theta[5]*(h-theta[7])^theta[9])^(1/theta[6])
)


fun.NR <- function(h, theta, Q) {
  theta[5]*(h-theta[4])^theta[6] + theta[8]*(h-theta[7])^theta[9] - Q
}








##########################################################
generateTimes<-function(t0,tmax,rate){
  ##########################################################
  # generate times between t0 and tmax with a given rate:
  times=c()
  titi=t0
  while(titi<=tmax){
    titi=titi+rexp(n=1,rate=rate)
    times=c(times,titi)
  }
  # remove times after tmax
  times=times[times<=tmax] 
  return(times)
}





#####################################
accuracy <- function(TP,FP,FN, TN) {
  #####################################
  accu = (TP+TN)/(TP+FP+FN+TN)
  return(accu)
}













#--------------------------------------------------
plot.generate.example = function(gaugings, title) {
  #---------------------------------------------------
  gaugings$title = title
  g=ggplot(data=gaugings, aes(x=h, y=Q,
                              fill=as.factor(period)))
  g=g+ geom_errorbar(aes(x=gaugings$h,
                         ymin =gaugings$Q-2*gaugings$uQ ,
                         ymax =gaugings$Q+2*gaugings$uQ), width=0.1, size = 0.3)
  g=g+geom_point(pch=21,size=5) + scale_fill_brewer(palette="Spectral")
  g=g+scale_y_log10(limits=c(min(gaugings$Q-2*gaugings$uQ),
                             max(gaugings$Q+2*gaugings$uQ))) + #Qquantile(0.00000000000001),Qquantile(0.99999999999999))) +
    labs(fill="Period")+
    ### Labels
    xlab(expression(paste("Stage h [m]",sep="")))+
    ylab(expression(paste("Discharge Q [",m^3,".",s^-1,"]",sep=""))) 
  # theme(strip.background =element_rect(fill="lightblue") )
  g=g + facet_grid(.~ title) +
    theme_light(base_size=20)+
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="right"
          ,plot.margin= unit(c(2, 0.5, 0.5, 0.5),"cm")
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,strip.text.x = element_text(
            size = 30, color = "black")
          ,strip.background = element_rect(
            color="gray70", fill="gray85", linetype="solid", size = 0.5))
}





















# plot of examples of datasets generated (Figure for paper Darienzo et al., 2020):
######################################################################################################
plot.examples.synthetic = function(dir,
                                   name.gauging.csv.file,
                                   indexes.sim.to.plot,
                                   titles.each.subplot) {
######################################################################################################
  gaug.synthet.example1 =   read.table(paste0(dir,"/",9,"/",  name.gauging.csv.file, ".csv"), sep= ";")
  gaug.synthet.example2 =   read.table(paste0(dir,"/",28,"/", name.gauging.csv.file, ".csv"), sep= ";")
  gaug.synthet.example3 =   read.table(paste0(dir,"/",38,"/", name.gauging.csv.file, ".csv"), sep= ";")
  gaug.synthet.example4 =   read.table(paste0(dir,"/",56,"/", name.gauging.csv.file, ".csv"), sep= ";")
  gaug.synthet.example5 =   read.table(paste0(dir,"/",85,"/", name.gauging.csv.file, ".csv"), sep= ";")
  gaug.synthet.example6 =   read.table(paste0(dir,"/",91,"/", name.gauging.csv.file, ".csv"), sep= ";")
  #
  plot1 = plot.generate.example(gaugings = gaug.synthet.example1, title = titles.each.subplot[1])
  plot2 = plot.generate.example(gaugings = gaug.synthet.example2, title = titles.each.subplot[2])
  plot3 = plot.generate.example(gaugings = gaug.synthet.example3, title = titles.each.subplot[3])
  plot4 = plot.generate.example(gaugings = gaug.synthet.example4, title = titles.each.subplot[4])
  plot5 = plot.generate.example(gaugings = gaug.synthet.example5, title = titles.each.subplot[5])
  plot6 = plot.generate.example(gaugings = gaug.synthet.example6, title = titles.each.subplot[6])
  #
  plotgrid.examples = plot_grid( plot1, plot2, plot3, plot4, plot5, plot6,
                                 ncol = 2, 
                                 nrow = 3,
                                 rel_widths = c(1,1), 
                                 rel_heights = c(1,1))
  #
  ggsave(plotgrid.examples, filename=paste0(dir,"/RCexamples.png"),
         device = "png", width = 20, height =25,
         dpi = 500, units = "in")
  #
  pdf(paste0(dir,"/Figure6.pdf"), 20, 27 ,useDingbats=F)
  print(plotgrid.examples)
  dev.off()
  #----------------------------------------------------------------------------------------------------
}

