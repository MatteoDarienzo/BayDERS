###############################################################################################
#                ==>    Insert all general inputs for the computation   <==                   #
###############################################################################################
#Hydrometric station general information:
#****************************************
station.name         = "Ardeche River at Meyras, France"  # Insert general info of the station (for plots)   
data.period          = "11/2001 - 6/2019"                 # Insert the study period (for plots)


# Gaugings:
#**********
file_gaugings        = "Meyras_gaugings_1.csv" # name of the csv file with the gaugings (put FALSE if not available!)
u.m.Hgaug            = "m"                     # options: "m", "cm", "mm" or "ft"
u.m.Qgaug            = "m^3*s^-1"              # options: "m^3.s-1" or "ft^3.s-1"
QGaug.col            = 7                       # column position in the file ...
uQGaug.col           = 8                       # column position in the file ...
hGaug.col            = 6                       # column position in the file ...
tGaug.col            = 13                      # or c(1,2,3,4,5,6) if a date format in multiple columns.
tGaug.format         = "%d/%m/%Y %H:%M"        # insert "numeric" or the date format like "%d/%m/%Y %H:%M".  
uQ.absolute          = FALSE                   # TRUE (if st dev of uncertainty) or FALSE (if % of Q)
gaugings_filter      = 1                       # consider a gauging every "gaugings_filter" gaugings. 

#settings for RC plot:
grid_RC.xlim         = c(-1,3)                    # grid for RC plots (stage h limits)
grid_RC.ylim         = c(0,200)                   # grid for RC plots (discharge Q limits)
grid_RC.xstep        = 1                          # grid for RC plots 
grid_RC.ystep        = 20                         # grid for RC plots
grid_RC.ylim.log     = c(0.01, 1000)              # grid for RC plots in y-log scale
ticks_RC.y.log       = c(0.01, 0.1,1,10,100,1000) # grid for RC-log plots
RC.y.labels          = "Discharge Q"              # RC:  y label (for plots)
RC.x.labels          = "Stage h"                  # RC:  x label (for plots)





# Stage record:
#**************
file_limni           = "Meyras_limni_cm.csv"              # name of the data .csv file with the stage record (put FALSE if not available!)
hLimni.col           = 2                                  # column position in the file ...
tLimni.col           = 3                                  # column position in the file ...
tLimni.format        = "%d/%m/%Y %H:%M"                   # insert "numeric" or the date format like "%d/%m/%Y %H:%M".  
limni_filter         = 1                                  # consider a stage value every "limni_filter" values. 
u.m.limni            = "cm"                               # ["m" / "cm" / "mm" / "ft"] unity of the stage record
grid_limni.ylim      = c(-1,4,1)                          # grid for Limni plots
limni.labels         = c("Time [days]", "Stage h [m]")    # limni: x and y labels (for plots)
date_origin          = "1899-12-30"                       # origin date to transform dates from numeric (by default is "1899-12-30").




 
# Official dates of rating shifts (or RC update):
#***********************************************
official.shift.times = "DateShiftsMeyras.csv"  # name of the csv file with official shift dates (put FALSE if not available!)
tOfficial.col        = 2                       # column position in the file ...
tOfficial.format     = "%d/%m/%Y"              # insert "numeric" or the date format like "%d/%m/%Y %H:%M".  











###################################################################################################
# Hydraulic configuration for BaRatin:  RC ==>  Q(h)  =  a* ( h - b ) ^ c  
###################################################################################################
# Where:
# - Q is the discharge,
# - h is the water level (also called stage)
# - a,b,c are the RC parameters for which priors are needed !!!
# Hydraulic configuration:
ncontrols         = 3   #number of controls
M                 = matrix(0, ncontrols, ncontrols)  #matrix of controls:
M[1,]             = c(1,0,0) # control section (rectangular weir in critical condition)
M[2,]             = c(0,1,0) # main control channel (rectangular wide channel in uniform condition)
M[3,]             = c(0,1,1) # main control channel + 2 floodways (rectangular channels)
control.type      = 0  #initialisation
control.type[1]   = "rect.weir"    
control.type[2]   = "rect.channel"
control.type[3]   = "rect.channel"




# PRIORS: (parameterization b,a,c !!!) for each control   
#------------------------------------------------------
#Parameter a:
#------------
propagat   = TRUE #[FALSE/TRUE] TRUE = propagate the priors of geom/phys properties for parameter "a"
a.prior    = c()  # initialize the prior for "a" (mean). If you already have it insert the prior.
st_a.prior = c()  # initialize the prior for "a" (st.dev). If you already have it insert the prior 
a.distr    = "LogNormal"  # Prior distribution of RC parameter a ("LogNormal"  suggested).
# Insert mean (Cr.prior) and stdev (st_Cr.prior) of Gaussian 
# distribution for each control (in the correct order!).
# or insert NA if "no control"  
#**************************************************************
#if rectangular weir: a = Cr*Bw*sqrt(2*g)
Cr.prior  = c(0.4, NA, NA);       st_Cr.prior = c(0.1, NA, NA)      # flux coefficient
g.prior   = c(9.81, NA, NA);      st_g.prior  = c(0.05, NA, NA)     # gravity
Bw.prior  = c(8, NA, NA);         st_Bw.prior = c(2, NA, NA)        # weir width
#if rectangular channel: a = Ks*S0*Bc
Bc.prior  = c(NA, 15, 30);        st_Bc.prior = c(NA, 3, 6)         # channel width
KS.prior  = c(NA, 25, 15);        st_KS.prior = c(NA, 3, 3)         # Strickler coeff.
S0.prior  = c(NA, 0.005, 0.005);  st_S0.prior = c(NA, 0.003, 0.003) # Longitudinal slope



#Parameter b:
#------------
b.distr = "Gaussian"  #"Fix" or "Gaussian" or "Uniform" or "Flat"
b.prior  = c(-0.6, -0.6, 1.2);    st_b.prior = c(0.5, 0.5, 0.2)


#Parameter c:
#------------
c.distr = "Gaussian"  #"Fix" or "Gaussian" or "Uniform" or "Flat"
c.prior = c(1.5, 1.67, 1.67);     st_c.prior= c(0.025, 0.025, 0.025)


#Remnant error model:
#--------------------
#LINEAR: err = g1+g2*Q   or    CONSTANT: err = g1):
remnant.err.model = "Linear"    #"Linear" or "Constant"
g1.prior          = c(0, 1000, 0.1)      #c(min,max, starting point)
g2.prior          = c(0, 100, 0.1)         #c(min,max, starting point) or c(mean, stdev, starting point) 
#if model = "Constant" then put FALSE !.
g1.distr.type     = "Uniform"              #  "Uniform.
g2.distr.type     = "Uniform"            # "Lognormal" or "Uniform
