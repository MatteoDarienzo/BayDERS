#########################################################################################################################
#                                    Options for the sediment transport analysis:                                       #
#########################################################################################################################
name.folder.results.ST    = "ST_test1"
file_name_ref_shift_times = "shift_times.txt"  # filename with the reference morphogenic shift times (***put FALSE if none***)
control.main.channel      = 2                       # [index] indicate which hydraulic control represents the main channel on which the average river bed is computed. 
Bc.st                     = Bc.prior[2] # average with of the river stretch characterized by sediment transport.
taucritic                 = 0.047  # by default 0.047 
d50.st                    = 0.06   # sediment characteristic diameter [m]
S0.st                     = 0.001  # average longitudinal slope of the river stretch characterized by sediment transport.
s.st                      = 2.65   # s = rho_s (density of sediments= 2650 kg/m3)/rho_w (density of water =1000 kg/m3)
alpha.sed.transp          = 8      # first coefficient of the MPM model  (8, taken from literature )
beta.sed.transp           = 1.5    # coefficient exponent in the MPM model (1.5, taken from literature )
tmax_between_events       = 1000   # max days between two peaks (to be considered part of the same morphogenic event)
deltat_peaks              = 1000   # max days between two peaks to chose the reference peaks



# options for plots:
ylimits              = c(-3, 9, 3)            # [min, max, step]  limits for stage record h(t)  y-axis (grid_limni.ylim). 
V.limits             = c(0, 500000, 100000)   # [mn, max, step]   limits for V(t) y-axis  , cumulative sediment volume [m^3]




# Options for linear estimation of relation V-deltab:
#****************************************************
model.type           = "null"    # options: "null" (if you want the mean =0) or "proportional" (if you want the mean != 0),
Nmcmc.st             = 100       # number of mcmc per cycle in BaM
Ncycles.st           = 1500      # number of mcmc cycles (tot number of mcmc = Nmcmc.st *Ncycles.st)  
# model V-deltab ==> deltab = a1 + a2*V = 0 
# in BaM
# DO NOT CHANGE !!!!
prior.param.propor    = c("'a1'", 0, "'Uniform'", -1e-13, 1e-13,    # param, start value, distribut, min, max (or mean, stdev) for a1 and a2
                          "'a2'", 0, "'Uniform'", -1e-13, 1e-13)    # priors for model 'null' -> deltab = 0 , thus a1~~0, a2~0.
#priors for error model proportional -> err= gamma1 + gamma2*V
prior.gamma.linear = c("'Linear'",
                       "'gamma1'", 1e-13, "'Uniform'", 0, 1e-13,      # intercept!     ~0 !!
                       "'gamma2'", 1e+8,  "'Uniform'", 0, 1e+13)      # slope!         










# other possible priors:
########################

# prior.param.propor    = c("'a1'", 0, "'Uniform'", -0.000000000001, 0.000000000001,
#                           "'a2'", 0, "'Uniform'", -0.000000000001, 0.000000000001)
# prior.gamma.linear = c("'Linear'",
#                        "'gamma1'", 0.00000000001, "'Uniform'", 0, 0.000000001,
#                        "'gamma2'", 100, "'Uniform'", 0, 1000000)   



# prior.param.propor = c("'a1'", 0, "'Uniform'", -0.000000000001, 0.000000000001,
#                        "'a2'", 0, "'Gaussian'", 0, 0.01)
# prior.gamma.linear = c("'Linear'",
#                        "'gamma1'", 0.0000000000001, "'Uniform'", 0, 0.00000000001,
#                        "'gamma2'", 1, "'Uniform'", 0, 1000000000)      # deltab = 0 + gamma2*V



# or maybe:
# prior.gamma.linear = c("'Linear'",
#                        "'gamma1'", 0.000000001, "'Uniform'", 0, 0.0000001, 
#                        "'gamma2'", 1, "'Uniform'", 0, 100000))









# uncomment this to test your choice of phi:
#********************************************
# d.norm         = rnorm(100000, mean= 0.1,   sd = 0.025)
# S0.norm        = rnorm(100000, mean= 0.005, sd = 0.001)
#phi.st = 10 ;  d50.st = 0.05;  S0.st = 0.005;
#phi.st = 20 ;  d50.st = 0.05;  S0.st = 0.0025;
#phi.st = 30 ;  d50.st = 0.10;  S0.st = 0.00333;
#phi.st = 15 ;  d50.st = 0.06;  S0.st = 0.004;
#phi.norm = d.norm/S0.norm
#phi.norm.stat = c(mean(phi.norm),  2*sd(phi.norm)) 
#plot(density(phi.norm))
