# om.ctl
# dval 		dlb 	dub  	phz  prior_type  p1    p2    	parameter
  6.21		0.00    10.0    1    1           6.21  1.00 	#log_bo
  6.21		0.00    10.0    1    1           6.21  1.00 	#log_b1
  0.80		0.20    1.00    1    3           1.01  1.01 	#h
  0.85		0.00    1.00    2    3           347.  56.2 	#s
  0.00		0.00    1.00   -2   -1           0.00  1.00 	#gamma
  9.65		0.00    10.0    2    4           1.26  0.00439 	#log_sigma
  1.83		0.00    10.0    2    4           3.69  0.56 	#log_tau
  0.00     -5.00    5.00    1   -1          -5.00  5.00 	#rec_devs

## The prior for s is based on M=0.15, so E(x) = exp(-0.15), Sig2 = (0.15*CV)^2 where
## the CV is arbirarilty assumed at 0.02: p1 = 347.3694, p2 = 56.2162