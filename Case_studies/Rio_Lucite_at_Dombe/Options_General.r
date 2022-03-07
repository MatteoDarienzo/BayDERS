###############################################################################################
#                ==>    Insert all inputs for the computation           <==                   #
###############################################################################################
#Hydrometric data of the station:
#********************************
# General description:
station.name         <- "Rio Lucite at Dombe (Mozambique)"  # Insert general info of the station (for plots)   
data.period          <- "7/1/1956 - 6/30/2017"                 # Insert the study period (for plots)


# Gaugings:
file_gaugings        <- FALSE #"Rio_Lucite_Dombe_gaugings.csv"  # name of the csv file with gaugings (put FALSE if not available!)
u.m.Hgaug            <- "m"              # options: "m", "cm", "mm" or "ft"
u.m.Qgaug            <- "m^3*s^-1"       # options: "m^3.s-1" or "ft^3.s-1"
u.m.time             <- "days"           # options: "days" or date format "c(Y,M,d,h,m,s)"
QGaug.col            <- 3                # column position in the file ...
hGaug.col            <- 2                # column position in the file ...
uQGaug.col           <- 4                # column position in the file ...
tGaug.col            <- 1                # or c(1,2,3,4,5,6) if a date format in multiple columns
tGaug.format        = "%d/%m/%Y"  # insert "numeric" or the date format like "%d/%m/%Y %H:%M".  
uQ.absolute          <- FALSE            # TRUE (if st dev of uncertainty) or FALSE (if % of Q)
gaugings_filter      <- 5                # consider a gauging every "gaugings_filter" gaugings. 

#settings for RC plot:
grid_RC.xlim         = c(0, 10)          # grid for RC plots
grid_RC.ylim         = c(0, 1200)         # grid for RC plots
grid_RC.xstep        = 2                # grid for RC plots
grid_RC.ystep        = 200               # grid for RC plots
grid_RC.ylim.log     = c(0.1, 10000)     # grid for RC-log plots
ticks_RC.y.log       = c(0.1, 1, 10, 100, 1000, 10000)    # grid for RC-log plots
RC.y.labels          = "Discharge Q"    # RC:  y labels (for plots)
RC.x.labels          = "Stage h"        # RC:  y labels (for plots)
    




# Stage record:
file_limni           <- "Rio_Lucite_Dombe_limni.csv"     # name of the csv file with the stage record (put FALSE if not available!)
u.m.limni            <- "m"             # options: "m", "cm", "mm" or "ft"
hLimni.col           <- 2                # column position in the file ...
tLimni.col           <- 1                # column position in the file ...
tLimni.format        = "%d/%m/%Y"                  # insert "numeric" or the date format like "%d/%m/%Y %H:%M".  
limni_filter         <- 1                # consider a stage every "stage_filter" values
grid_limni.ylim      = c(0,10,2)         # grid for Limni plots
limni.labels         = c("Time [days]", "Stage h [m]")  #limni: x and y labels (for plots)
limni.save           = "Stage_record"   # Name of the study Limni (for plots)
date_origin          = "1899-12-30"      # to transform dates from numeric (by default is "1899-12-30").





# Official dates of rating shifts (or RC update):
#***********************************************
official.shift.times <- "Rio_Lucite_Dombe_official_shift_times.csv" # name of the csv file with official shift dates (put FALSE if not available!)
tOfficial.col        = 2                       # column position in the file .
tOfficial.format     = "%d/%m/%Y"                         # insert "numeric" or the date format like "%d/%m/%Y %H:%M".  






# Hydraulic configuration for BaRatin:      ==>     Q(h)  =  a* ( h - b ) ^ c     <==
#*************************************************************************************
# Where:
# - Q is the discharge,
# - h is the water level (also called stage)
# - a,b,c are the RC parameters for which priors are needed !!!

# Hydraulic configuration:
ncontrols         = 1         #number of controls
M                 = matrix(0,ncontrols, ncontrols)  #matrix of controls:
M[1,]             = c(1) # control section (rectangular weir in critical condition)
#M[2,]             = c(1,1) # main control channel (rectangular wide channel in uniform condition)
control.type      = 0  #initialisation
control.type[1]   = "rect.channel"
#control.type[2]   = "rect.channel"





# PRIORS: (parameterization b,a,c !!!) for each control   
#------------------------------------------------------
#Parameter a:
#------------
propagat   = TRUE #[FALSE/TRUE] TRUE = propagate the priors of geom/phys properties for parameter "a"
a.prior    = c(38.96)  # initialize the prior for "a" (mean). If you already have it insert the prior.
st_a.prior = c(7)  # initialize the prior for "a" (st.dev). If you already have it insert the prior 
a.distr    = "LogNormal"  # Prior distribution of RC parameter a ("LogNormal"  suggested).
# Insert mean (Cr.prior) and stdev (st_Cr.prior) of Gaussian 
# distribution for each control (in the correct order!).
# or insert NA if "no control"  
#**************************************************************
#if rectangular weir: a = Cr*Bw*sqrt(2*g)
Cr.prior  = c(NA);         st_Cr.prior = c(NA)      # flux coefficient
g.prior   = c(NA);         st_g.prior  = c(NA)     # gravity
Bw.prior  = c(NA);         st_Bw.prior = c(NA)        # weir width
#if rectangular channel: a = Ks*S0^0.5*Bc
Bc.prior  = c(25);        st_Bc.prior = c(5)         # channel width
KS.prior  = c(20);        st_KS.prior = c(5)         # Strickler coeff.
S0.prior  = c(0.005);        st_S0.prior = c(0.001) # Longitudinal slope



#Parameter b:
#------------
b.distr = "Gaussian"  #"Fix" or "Gaussian" or "Uniform" or "Flat"
b.prior  = c(0.536);    st_b.prior = c(0.2)


#Parameter c:
#------------
c.distr = "Gaussian"  #"Fix" or "Gaussian" or "Uniform" or "Flat"
c.prior = c( 1.67);     st_c.prior= c(0.025)


#Remnant error model:
#--------------------
#LINEAR: err = g1+g2*Q   or    CONSTANT: err = g1):
remnant.err.model = "Linear"    #"Linear" or "Constant"
g1.prior          = c(0, 1000, 0.1)      #c(min,max, starting point)
g2.prior          = c(0, 100, 0.1)         #c(min,max, starting point) or c(mean, stdev, starting point) 
#if model = "Constant" then put FALSE !.
g1.distr.type     = "Uniform"              #  "Uniform.
g2.distr.type     = "Uniform"            # "Lognormal" or "Uniform



