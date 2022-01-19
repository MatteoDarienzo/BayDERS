#########################################################################################################################
#                                    Options for the sediment transport analysis:                                       #
#########################################################################################################################
name.folder.results.ST = "Test_1"

d50      = 0.05;       # characteristic sediment diameter [m].
S0.ST    = 0.005;      # average longitudinal slope of the river stretch characterized by sediment transport.
Bc.st    = Bc.prior[2] # average with of the river stretch characterized by sediment transport.
d.norm   = rnorm(100000, mean= 0.1,   sd = 0.025)
S0.norm  = rnorm(100000, mean= 0.005, sd = 0.001)

phi.st   = 20     # phi = d50/S0
d50.st   = 0.1    # sediment characteristic diameter
S0.st    = 0.005  # average longitudinal slope of the river stretch characterized by sediment transport.






#phi.st = 10 ;  d50.st = 0.05;  S0.st = 0.005;
#phi.st = 20 ;  d50.st = 0.05;  S0.st = 0.0025;
#phi.st = 30 ;  d50.st = 0.10;  S0.st = 0.00333;
#phi.st = 15 ;  d50.st = 0.06;  S0.st = 0.004;








# uncomment this to test your choice of phi:
#********************************************
#phi.norm = d.norm/S0.norm
#phi.norm.stat = c(mean(phi.norm),  2*sd(phi.norm)) 
#plot(density(phi.norm))
