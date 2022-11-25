#########################################################################################################################
#                                    Options for the recession analysis:                                                #
#########################################################################################################################
#####################################
# 1) Recession extraction:
#####################################
name.folder.results.recession = "chi_30_tgood_20"      # string with the name of the folder with results.
# name.folder.results.recession = "tburn0.1_grad-100_chi50_tgood50_dtmax50_uh0.5_glin_N10_dtmin0_fil1_mcmc5000_errold_shift150"      # string with the name of the folder with results.
uh.rec                 = 0.5                               # assign a stdev to the stage obs error [in cm], by default is 0.
tburn.rec              = 0.1                               # discard the first n days of recession, by default is 0.
Nmin.rec               = 10                                # min number of data in a recession curve
tgood                  = 20                                # min length of the recession in days
delta.t.min            = 0                                 # min days between two recess data
delta.t.max            = 20                                # max days between two recess data
chi                    = 30                                # max stage rise between two recess data
gradient.max           = -100                              # max gradient dh/dt for the recession 







#####################################
# 2) Recession estimation:
#####################################
estim.plot.results.only  = TRUE                             # [TRUE/FALSE] put TRUE if you only to plot results previously obtained.
prior.param.rec          = NULL                          # [NULL] initialise Priors for the recession model: e.g., h = a1*e^(-b1) + a2*e^(-b2) + a3*e^(-b3) + a4
rec.model = NULL

# ### "1expWithAsympt" ==> h(t) = a1(k)*exp(-b1*t) + a2(k)

# rec.model[[1]]                = c("1expWithAsympt")           # Select stage-recession models (also more than one) (see the list at the bottom):
# prior.param.rec[[1]]  = c(  0,                            1000,     "'Uniform'",     500,                        "var",
#                             -log((0.5))+log(log(2)) ,     1,        "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                             -150                       ,  50 ,      "'Uniform'",     0 ,                      "var")


rec.model[[1]]                = c("2expWithAsympt")           # Select stage-recession models (also more than one) (see the list at the bottom):
prior.param.rec[[1]]  = c(  0,                            1000,       "'Uniform'",     100,                       "var",
                            -log((0.5))+log(log(2)) ,       1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
                            0,                             100,       "'Uniform'",     50,                         "static",
                            -log(100)+log(log(2)) ,       0.5,      "'LogNormal'",   exp(-log(100)+log(log(2))),  "static",
                            -150                  ,        50,      "'Uniform'",      0 ,                         "var")


# rec.model[[2]]                = c("2expWithAsympt_bis")       # Select stage-recession models (also more than one) (see the list at the bottom):
# prior.param.rec[[2]]  = c(  0,                             1000,       "'Uniform'",     100,                      "var",
#                             -log((0.5))+log(log(2)) ,      1,          "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                             0,                             100,       "'Uniform'",     50,                         "var",
#                             -log(100)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(100)+log(log(2))),  "static",
#                             -150                       ,   50 ,      "'Uniform'",      0 ,                         "var")
# 
# 
# 
# rec.model[[3]]                = c("3expWithAsympt")                                                                          # Select stage-recession models (also more than one) (see the list at the bottom):
# prior.param.rec[[3]]  =c( 0,                            1000,       "'Uniform'",      100,                      "var",       # prior: [in cm for param a and days-1 for param b]
#                           -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",     exp(-log(0.5)+log(log(2))), "static",     # min or mean, max or stdev, distribution, starting value, "static"/"var";
#                           0,                            500,       "'Uniform'",       50,                         "static",
#                           -log(50)+log(log(2)) ,        1,      "'LogNormal'",     exp(-log(50)+log(log(2))),  "static",
#                           0,                            100,      "'Uniform'",        10,                         "static",
#                           -log(100)+log(log(2)) ,       0.5,       "'LogNormal'",  exp(-log(100)+log(log(2))),  "static",
#                           -150                 ,        50 ,      "'Uniform'",     0 ,                         "var")        #  this last line (asymptote prior) is specific to your case study !!!!
# 
# 
# 
# rec.model[[4]]                = c("3expWithAsympt_bis")
# prior.param.rec[[4]]  =c( 0,                            1000,     "'Uniform'",     200,                      "var",
#                           -log((0.5))+log(log(2)) ,     1,        "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                           0,                            500,      "'Uniform'",     50,                         "var",
#                           -log(50)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(50)+log(log(2))),  "static",
#                           0,                            100,      "'Uniform'",     10,                         "static",
#                           -log(100)+log(log(2)) ,       0.5,      "'LogNormal'",   exp(-log(100)+log(log(2))),  "static",
#                           -150                       ,  50 ,      "'Uniform'",        0 ,                       "var")
# 


# # ### "expexp"   ==> h(t) = a1(k)*exp(exp(-b1*t^eta)) + a2(k)
# rec.model[[6]]                = c("expexp")
# prior.param.rec[[6]]  = c( 0,                            1000,     "'Uniform'",      200,                       "var",
#                            -log((0.5))+log(log(2)),      1,        "'LogNormal'",    exp(-log(0.5)+log(log(2))), "static",
#                            -log((0.5))+log(log(2)),      1,        "'LogNormal'",    exp(-log(0.5)+log(log(2))), "static",
#                            -150                   ,      50,       "'Uniform'",      0,                          "var")
# 
# 
# 
# # ### "expexp bis"
# rec.model[[7]]                = c("expexp_bis")
# prior.param.rec[[7]]  = c( 0,                            1000,     "'Uniform'",     200,                        "var",
#                            -log((0.5))+log(log(2)) ,     1,        "'LogNormal'",   exp(-log(0.5)+log(log(2))), "var",
#                            -log((0.5))+log(log(2)) ,     1,        "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -150                    ,     50 ,      "'Uniform'",     0,                          "var")
# 
# 
# 
# rec.model[[8]]                = c("hyperb_bis")
# prior.param.rec[[8]]  = c( 0,                            1000,      "'Uniform'",     200,                        "var",
#                            -log((0.5))+log(log(2)) ,     1,         "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -log((0.5))+log(log(2)) ,     1,         "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -150                    ,     50,        "'Uniform'",     0,                          "var")
# 
# 
# rec.model[[9]]                = c("hyperb_bis")
# prior.param.rec[[9]]  = c( 0,                            1000,      "'Uniform'",     200,                        "var",
#                             -log((0.5))+log(log(2)) ,    1,         "'LogNormal'",   exp(-log(0.5)+log(log(2))), "var",
#                             -log((0.5))+log(log(2)) ,    1,         "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                             -150                    ,    50,        "'Uniform'",     0,                          "var")
# 


# Remnant error model (err = gamma1 + gamma2*h   or  err = gamma1):
gamma.model.rec          = "linear"                     # Structural error model: "constant" or "linear".    
prior.gamma.rec          = c(0, 100,   "'Uniform'", 10,     # Structural error prior gamma1 (intercept):  min or mean, max or stdev, distribution, starting value
                             0, 10,   "'Uniform'", 1)      # Structural error prior gamma2 (slope):  min or mean, max or stdev, distribution, starting value
Ncycles.mcmc.rec         = 500                            # number of mcmc cycles during the bayesian inference for recessions estimation
nmcmc.rec                = 100                            # number of mcmc per cycle
nslim.rec                = 10                             # number of slimmed mcmc for the summary statistics
jump.pos.rec             = 1.1                            # parameter positive jump rate
jump.neg.rec             = 0.9                            # parameter negative jump rate  
nburn.rec                = 0.5                            # fraction of initial mcmc burned.










#####################################
# 3) Recession segmentation:
#####################################
prior.mu.rec.segment           = c("Uniform", -100, 100, 0)        # c(distribution,  min or mean,  max or stdev, stanrting value)
gamma.prior.rec                = c("Uniform", 0, 10, 1)           # c(distribution,  min or mean,  max or stdev,  starting value)
nSmax.rec                      = 10                                # Maximum number of change points that can be detected 
criterion.rec                  = "BIC"                            # crtierion for the selection of optimal number of change points ("DIC", "BIC")
seg.plot.results.only          = FALSE                            # [FALSE/TRUE] if TRUE no computation is performed. only results are plotted.
Ncycle.rec.segment             = 1500   #1500                          # number of cycles mcmc.
Nmcmc.rec.segment              = 100    #100                          # number of mcmc per cycle.
Nslim.rec.segment              = 100                              # for summary statistics, we consider one mcmc each "Nslim.rec.segment".
Nburn.rec.segment              = 0.5                              # fraction of mcmc burned!
tmin.rec                       = 1                                # minimum distance between shifts
shift.time.adjustment.type.rec = 2                                # 1 = select always the MAP; 2= select always the largest closest flood peak; 3 = will ask to manually insert a value 







