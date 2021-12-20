#########################################################################################################################
#                                    Options for the sediment transport analysis:                                       #
#########################################################################################################################
name.folder.results.ST    = "TS_1"
file_name_ref_shift_times = "shift_times.txt"  # filename with the reference morphogenic shift times (***put FALSE if none***)



Bc.st                = Bc.prior[2] # average with of the river stretch characterized by sediment transport.
taucritic            = 0.047  # by default 0.047 
d50.st               = 0.06   # sediment characteristic diameter
S0.st                = 0.001  # average longitudinal slope of the river stretch characterized by sediment transport.
s.st                 = 2.65
alpha.sed.transp     = 8      # coefficient  MPM model
beta.sed.transp      = 1.5    # coefficient MPM model
tmax_between_events  = 1000   # max days between two peaks (to be considered part of the same morphogenic event)
phi.st               = d50.st/S0.st  # for thesis, dont change!

# options for plots:
ylimits              = c(-3, 9, 3)  # grid_limni.ylim,
tlimits              = c(0, 57000)
V.limits             = c(0, 500000, 1000)
V.for.text           = 500000


# Options for BaM (linear estimation of relation V-deltab):
Nmcmc.st             = 100
Ncycles.st           = 1500


# deltab = 0 + gamma2*V
prior.param.propor    = c("'a1'", 0, "'Uniform'", -0.000000000001, 0.000000000001,
                          "'a2'", 0, "'Uniform'", -0.000000000001, 0.000000000001)
prior.gamma.linear = c("'Linear'",
                       "'gamma1'", 0.00000000001, "'Uniform'", 0, 0.000000001,
                       "'gamma2'", 100, "'Uniform'", 0, 1000000)   




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
