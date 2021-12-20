#########################################################################################################################
#                                    Options for the recession analysis:                                                #
#########################################################################################################################
#####################################
# 1) Recession extraction:
#####################################
name.folder.results.recession = "Test_4"                 # string with the name of the folder with results.
uh.rec                 = 0                               # assign a stdev to the stage obs error [in cm], by default is 0.
tburn.rec              = 0.1                             # discard the first n days of recession, by default is 0.
Nmin.rec               = 5                               # min number of data in a recession curve
tgood                  = 15                              # min length of the recession in days
delta.t.min            = 0                               # min days between two recess data
delta.t.max            = 15                              # max days between two recess data
chi                    = 20                              # max stage rise between two recess data
gradient.max           = -1                              # max gradient dh/dt for the recession 



#####################################
# 2) Recession estimation:
#####################################
estim.plot.results.only  = FALSE                         # [TRUE/FALSE] put TRUE if you only to plot results previously obtained.
rec.model                = c("3expWithAsympt")           # Select stage-recession models (also more than one) (see the list at the bottom):
prior.param.rec          = NULL                          # [NULL] initialise Priors for the recession model: e.g., h = a1*e^(-b1) + a2*e^(-b2) + a3*e^(-b3) + a4
prior.param.rec[[1]]  =c( 0,                            1000,       "'Uniform'",     100,                      "var",  # prior: [in cm for param a and days-1 for param b]
                          -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static", # min or mean, max or stdev, distribution, starting value, "static"/"var";
                          0,                            500,       "'Uniform'",     50,                         "static",
                          -log(50)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(50)+log(log(2))),  "static",
                          0,                            100,      "'Uniform'",      10,                         "static",
                          -log(80)+log(log(2)) ,        0.2,       "'LogNormal'",  exp(-log(80)+log(log(2))),  "static",
                          -150                 ,        50 ,      "'Uniform'",     0 ,                         "var")  #  this last line (asymptote prior) is specific to your case study !!!!
prior.gamma.rec          = c(0, 10, "'Uniform'",  1,      # Structural error prior gamma1:  min or mean, max or stdev, distribution, starting value
                             0, 10,  "'Uniform'", 0.1)    # Structural error prior gamma2:  min or mean, max or stdev, distribution, starting value
Ncycles.mcmc.rec         = 100                            # number of mcmc cycles during the bayesian inference for recessions estimation
nmcmc.rec                = 100                            # number of mcmc per cycle
nslim.rec                = 50                             # number of slimmed mcmc for the summary statistics
jump.pos.rec             = 1.1                            # parameter positive jump rate
jump.neg.rec             = 0.9                            # parameter negative jump rate  
nburn.rec                = 0.5                            # fraction of initial mcmc burned.



#####################################
# 3) Recession segmentation:
#####################################
prior.mu.rec.segment           = c("Uniform", -150, 50, 0)        # c(distribution,  min or mean,  max or stdev, stanrting value)
gamma.prior.rec                = c("Uniform", 0, 50, 0.1)         # c(distribution,  min or mean,  max or stdev,  starting value)
nSmax.rec                      = 6                                # Maximum number of change points that can be detected 
criterion.rec                  = "BIC"                            # crtierion for the selection of optimal number of change points ("DIC", "BIC")
seg.plot.results.only          = FALSE                            # [FALSE/TRUE] if TRUE no computation is performed. only results are plotted.
Ncycle.rec.segment             = 1500                             # number of cycles mcmc.
Nmcmc.rec.segment              = 100                              # number of mcmc per cycle.
Nslim.rec.segment              = 100                              # for summary statistics, we consider one mcmc each "Nslim.rec.segment".
Nburn.rec.segment              = 0.5                              # fraction of mcmc burned!
tmin.rec                       = 1                                # minimum distance between shifts
shift.time.adjustment.type.rec = 1                                # 1 = select always the MAP; 2= select always the largest closest flood peak; 3 = will ask to manually insert a value 









#############################################################################################################################
# Add other models if you want to do a comparison!  (see Darienzo et al., 2021)
# The following models are available (index (k) indicates the recession specific parameters):
#############################################################################################################################
# rec.model  = c("1expWithAsympt", "2expWithAsympt", "3expWithAsympt", "expexp", "hyperb", "Coutagne", 
#                "2expWithAsympt_bis", "3expWithAsympt_bis", "expexp_bis", "hyperb_bis", "Coutagne_bis", 
#                "2expWithAsympt_rel")
#
##########
# PRIORS:
##########
# ### "1expWithAsympt" ==> h(t) = a1(k)*exp(-b1*t) + a2(k)
# prior.param.rec[[1]]  = c(  0,                            1000,     "'Uniform'",     500,                        "var",
#                             -log((0.5))+log(log(2)) ,     1,        "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                             -10000                       , 10000 ,  "'Uniform'",     1000 ,                      "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "2expWithAsympt" ==> h(t) = a1(k)*exp(-b1*t) + a2*exp(-b2*t) + a3(k)
# prior.param.rec[[2]]  = c(  0,                            1000,       "'Uniform'",     200,                      "var",
#                             -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                             0,                            500,       "'Uniform'",     50,                         "static",
#                             -log(80)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(50)+log(log(2))),  "static",
#                             -10000                       , 10000 ,      "'Uniform'",  1000 ,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### 3expWithAsympt   ==> h(t) = a1(k)*exp(-b1*t) + a2*exp(-b2*t) + a3*exp(-b3*t) + a4(k)
# prior.param.rec[[3]]  =c( 0,                            1000,       "'Uniform'",     200,                      "var",
#                           -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                           0,                            200,       "'Uniform'",     50,                         "static",
#                           -log(50)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(50)+log(log(2))),  "static",
#                           0,                            100,      "'Uniform'",      10,                         "static",
#                           -log(80)+log(log(2)) ,       0.5,       "'LogNormal'",  exp(-log(80)+log(log(2))),  "static",
#                           -10000                       , 10000 ,      "'Uniform'",  1000 ,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "expexp"   ==> h(t) = a1(k)*exp(exp(-b1*t^eta)) + a2(k)
# prior.param.rec[[4]]  = c( 0,                            10000,       "'Uniform'",     200,            "var",
#                            -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -10000                       , 10000 ,      "'Uniform'",  1000 ,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "hyperb"
# prior.param.rec[[5]]  = c( 0,                            1000,       "'Uniform'",     200,                      "var",
#                            -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -10000                       , 10000 ,      "'Uniform'",  1000,                         "var")
# 
#-----------------------------------------------------------------------------------------------------------------------------
# ### "Coutagne"
# prior.param.rec[[6]]  = c(   0,                            1000,     "'Uniform'",     200,                      "var",
#                              -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                              -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                              -10000                       , 10000 ,      "'Uniform'",  1000,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# # ### "Boussinesq"
# # prior.param.rec[[7]]  = c(   0,                            1000,     "'Uniform'",     200,                      "var",
# #                              -log((0.5))+log(log(2)) ,     1,       "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
# #                              -10000                     , 10000 ,      "'Uniform'",  1000 ,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "2expWithAsympt bis"
# prior.param.rec[[7]]  = c(  0,                            1000,       "'Uniform'",     200,                      "var",
#                             -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                             0,                            500,       "'Uniform'",     50,                         "var",
#                             -log(80)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(50)+log(log(2))),  "static",
#                             -10000                       , 10000 ,      "'Uniform'",  1000 ,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "3expWithAsympt bis"
# prior.param.rec[[8]]  =c( 0,                            1000,       "'Uniform'",     200,                      "var",
#                           -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                           0,                            500,       "'Uniform'",     50,                         "var",
#                           -log(50)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(50)+log(log(2))),  "static",
#                           0,                            100,      "'Uniform'",      10,                         "static",
#                           -log(80)+log(log(2)) ,       0.5,       "'LogNormal'",  exp(-log(80)+log(log(2))),  "static",
#                           -10000                       , 10000 ,      "'Uniform'",  1000 ,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "expexp bis"
# prior.param.rec[[9]]  = c( 0,                            10000,       "'Uniform'",     200,            "var",
#                            -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "var",
#                            -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                            -10000                       , 10000 ,      "'Uniform'",  1000 ,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "hyperb bis" 
# prior.param.rec[[10]]  = c( 0,                            1000,       "'Uniform'",     200,                      "var",
#                             -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "var",
#                             -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                             -10000                       , 10000 ,      "'Uniform'",  1000,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "Coutagne bis"
# prior.param.rec[[11]]  = c(   0,                            1000,     "'Uniform'",     200,                      "var",
#                               -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "var",
#                               -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                               -10000                       , 10000 ,      "'Uniform'",  1000,                         "var")
#-----------------------------------------------------------------------------------------------------------------------------
# ### "2expWithAsympt rel"
# prior.param.rec[[12]]  = c(  0,                            1000,       "'Uniform'",     200,                      "var",
#                              -log((0.5))+log(log(2)) ,     1,      "'LogNormal'",   exp(-log(0.5)+log(log(2))), "static",
#                              0,                            100,       "'Uniform'",     0.1,                         "static",
#                              -log(80)+log(log(2)) ,        0.5,      "'LogNormal'",   exp(-log(50)+log(log(2))),  "static",
#                              -10000                       , 10000 ,      "'Uniform'",  1000# ,                         "var")
###############################################################################################################################






