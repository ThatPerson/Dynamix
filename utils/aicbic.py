## https://en.wikipedia.org/wiki/Maximum_likelihood_estimation - equation under "Continuous distribution, continuous parameter space" - normal distribution PDF.

# https://www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html

import sys
import numpy as np
# arguments are python aicbic.py FOLDER N_PARAMS
# then we read FOLDER/backcalc_N.dat from N=1-56
# for each, 

folder = sys.argv[1]
params = int(sys.argv[2])
with open(sys.argv[3], "w") as q:
	for i in range(1, 57):
		x = np.loadtxt("%s/backcalc_%d.dat" % (folder, i))
		#print(x)
		if (np.size(x) == 0):
			q.write("%d, -1, -1\n" % (i))
			continue

		calc_R = x[:, 1]
		exp_R = x[:, 2]
		err_R = x[:, 3]

		

		n_data = np.size(calc_R)




		#exp_R[err_R == -1] = calc_R[err_R == -1]
		#err_R[err_R == -1

		sigma = err_R / 2.
		## errors are 2sigma (2 standard deviations). So sigma is err_R/2

		expFun = -np.power(calc_R - exp_R, 2.) / (2 * np.power(sigma, 2.))
		divisor = np.sqrt(2 * np.pi * np.power(sigma, 2.))

		errF = np.exp(expFun) / divisor
		
		LerrF = np.log(errF)
		
		
		for qw in range(0, np.size(calc_R)):
			if (calc_R[qw] < 0): # no fit
				LerrF[qw] = 0
				n_data = n_data - 1
		
		#LerrF[LerrF < 1E-308] = 0 # not entirely sure how valid this is
		#print(LerrF)
		logv = np.sum(LerrF)
		
		AIC = 2 * params - 2 * logv
		BIC = params * np.log(n_data) - 2 * logv
		#print("%d, %f, %f\n"% (i, AIC, BIC))
		q.write("%d, %f, %f\n"% (i, AIC, BIC))
