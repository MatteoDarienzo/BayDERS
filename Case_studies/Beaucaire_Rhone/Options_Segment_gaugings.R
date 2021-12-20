# General options:
################################################################
name.folder.results        = "mu500_BIC_MAP"
prior.mu.segm              = c(0, 500)          # prior for segments mean "mu" = N(mean, stdev): stdev must be of the Order of gaugings residuals
tmin                       = 0                  # mininum distance (in time) between two consecutif change points (it depends on time units!). By default =0
nSmax                      = 5                  # Maximum number of segments in the series at each iteration
Nmin                       = 1                  # minimum number of data in a segment (>=1, by default = 1) 
criterion                  = "BIC"              # Criterion for the choice of the model ("AIC","BIC", "DIC", "HQC")
shift.time.adjustment.type = 1                  # 1 = select always the MAP; 2= select always the largest closest flood; 3 = will ask to manually insert a value 
plot.results.only          = T              # [FALSE/TRUE] if TRUE no computation is performed. only previous results are plotted.
plot_dates                 = TRUE               # [TRUE/FALSE] plot or not the date of each detected shift time (if any) on the stage record.
st_b.prior                 = st_b.prior         # Initial stdev of b: e.g., c(1.5, 1.5) or leave st_b.prior.









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



# df.limni          = df.limni 
# t_Gaug            = t_Gaug
# h_Gaug            = h_Gaug
# Q_Gaug            = Q_Gaug 
# uQ_Gaug           = uQ_Gaug
# u.m.Qgaug         = u.m.Qgaug 
# u.m.Hgaug         = u.m.Hgaug
# ncontrols          = ncontrols
# M                 = M 
# b.distr           = b.distr 
# a.distr           = a.distr 
# c.distr           = c.distr 
# # a.prior           = a.prior 
# # st_a.prior        = st_a.prior 
# c.prior           = c.prior 
# st_c.prior        = st_c.prior 
# Bw.prior          = Bw.prior
# Cr.prior          = Cr.prior
# g.prior           = g.prior
# Bc.prior          = Bc.prior
# KS.prior          = KS.prior
# S0.prior          = S0.prior
# st_Bw.prior       = st_Bw.prior
# st_Cr.prior       = st_Cr.prior
# st_g.prior        = st_g.prior 
# st_Bc.prior       = st_Bc.prior
# st_KS.prior       = st_KS.prior
# st_S0.prior       = st_S0.prior # Priors for RC model
# remnant.err.model = remnant.err.model 
# g1.prior          = g1.prior 
# g2.prior          = g2.prior
# g1.distr.type     = g1.distr.type
# g2.distr.type     = g2.distr.type  







