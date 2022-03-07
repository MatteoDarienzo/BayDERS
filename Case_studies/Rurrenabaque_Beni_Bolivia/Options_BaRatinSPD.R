###############################################################
# BaRatin SPD model (Mansanarez et al.,2019):
# General options:
################################################################
name.folder.results.SPD  =  "SPD_test1"
global.change            = TRUE                          # [TRUE/FALSE] this will consider a change of all dynamic RC parameter "b" (vertical traslation of the entire channel). 
local.change             = FALSE                         # [TRUE/FALSE] this will consider a change of RC parameter "b1" (translation of the lowest control only). 
width.change             = FALSE                         # [TRUE/FALSE] this will consider a change of RC parameter "a" (due to channel width variation). 
n.parvar                 = 1                             # number of dynamic parameters
isVar                    = c(T,F,F, F,F,F)               # [T/F] Define which RC parameters are varying (b1,a1,c1, b2,a2,c2, ...)
dg.prior                 = c(0, 0.25)                    # mean and stdev for incremental global change priors
dl.prior                 = c(0, 0.25)                    # mean and stdev for incremental local change priors
dB.prior                 = c(1, 0.2)                     # mean and stdev for channel width change priors
b.prior                  = c(-0.5, 4)                       # mean for parameters "b" for the oldest period 
st_b.prior               = c(0.5, 0.2)              # stdev for parameters "b" for the oldest period
param.var                = c("b1")                   # names of parameters that vary (e.g. = c("b1", "b2) )
plot.results.only        = FALSE                         # [TRUE/FALSE] if TRUE it plot the results only without computation.
do.final.hydrograph      = FALSE                         # [TRUE/FALSE] option not supported yet!
show.activation.stages   = FALSE                         # [TRUE/FALSE] If you want to plot the activation stages of the lowest controls (max two controls). 
Ncycles                  = 500                           # Number of MCMC cycles for Metropolis-Hastings for BaM.
ylim.log.wind            = c( grid_RC.ylim.log[1], 30000)  # define Q limits for plots in log scale (e.g. grid_RC.ylim.log).









# Advanced (do not change):
###################################################################
# settings fro plot of prior correlation (default options):
ncolumns                = 8
rowmax                  = 20 
# Settings for the final multiple RC plot (default options):
op_cooked               = TRUE   # [TRUE/FALSE] if you want to consider the mcmc cooked or not.
nburn                   = 0.5    # burn fraction of mcmc
nslim                   = 10     # slim of mcmc for statistics summary
na.rm                   = TRUE
size.gaugings           = 4
ylim.wind               = grid_RC.ylim 
xlim.wind               = grid_RC.xlim
breaks.lin.x            = seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)
labels.lin.x            = seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)
breaks.lin.y            = seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep)
labels.lin.y            = seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep)
breaks.log              = ticks_RC.y.log
labels.log              = ticks_RC.y.log