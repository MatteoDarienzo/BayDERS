# General options:
################################################################
name.folder.results        = "mu500_BIC_closestflood"
prior.mu.segm              = c(0, 500)          # prior for segments mean "mu" = N(mean, stdev): stdev must be of the Order of gaugings residuals
tmin                       = 0                  # mininum distance (in time) between two consecutif change points (it depends on time units!). By default =0
nSmax                      = 5                  # Maximum number of segments in the series at each iteration
Nmin                       = 1                  # minimum number of data in a segment (>=1, by default = 1) 
criterion                  = "BIC"              # Criterion for the choice of the model ("AIC","BIC", "DIC", "HQC")
shift.time.adjustment.type = 2                  # 1 = select always the MAP; 2= select always the largest closest flood; 3 = will ask to manually insert a value 
plot.results.only          = FALSE              # [FALSE/TRUE] if TRUE no computation is performed. only previous results are plotted.
plot_dates                 = TRUE               # [TRUE/FALSE] plot or not the date of each detected shift time (if any) on the stage record.
st_b.prior                 = st_b.prior         # Initial stdev of b: e.g., c(1.5, 1.5) or leave st_b.prior.
deltat_peaks               = 1000               # minimum time lag between major floods (to find all major peaks and adjust shift times), e.g., 100 or 1000 days









# Other BaM options (mcmc, ...): By default
################################################
predictionRC         = TRUE
predictionQt         = FALSE 
predictionPrior      = FALSE
simMCMC              = TRUE
mcmc.prior           = 10000
Ncycle.segment       = 1500
Ncycles.max          = 2000    # must be > Ncycle.segment !
Nmcmc.segment        = 100
Nslim.segment        = 100
Ncycle.baratin       = 500
nsim                 = 10000   # number of simulations for the prior propagation










# Advanced options  - DO NOT TOUCH if you are not sure!    /!\ /!\
##################################################################
save.all.results           =  TRUE
resid.uncertaint           =  TRUE
recursive                  =  TRUE
use.other.seg.method       =  FALSE
prior.gamma.segm.without.U = c(0, 100)  # c(min, max): uniform distrib. when resid.uncertaint = FALSE
prior.gamma.segm.with.U    = c(log(0.001), 0.001) # c(meanlog, sdlog): LogNormal distrib. when resid.uncertaint =TRUE
type.stdev.gamma.baratin   = 1  # type of stdev for the RC structural error:  1  (~U(0, mean gamma)) or 2  (~U(0, mean gamma + 2*stdev gamma)) 






