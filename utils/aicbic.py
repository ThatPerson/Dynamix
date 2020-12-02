## https://en.wikipedia.org/wiki/Maximum_likelihood_estimation - equation under "Continuous distribution, continuous parameter space" - normal distribution PDF.

# https://www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html

import sys
import numpy as np
# arguments are python aicbic.py FOLDER modeltype OUTPUTFILE
# then we read FOLDER/backcalc_N.dat from N=1-56
# for each, 

folder = sys.argv[1]
model = sys.argv[2]

adj = 0
params = -10
if (model[0] == "V"):
	model = model[1:]
	adj = 3

if (model == "SMF"):
	params = 2
elif (model == "EMF"):
	params = 3
elif (model == "EMFT"):
	params = 5
elif (model == "SMFT"):
	params = 3
elif (model == "DEMF"):
	params = 4
elif (model == "DEMFT"):
	params = 6
elif (model == "GAF"):
	params = 8
elif (model == "GAFT"):
	params = 10
elif (model == "EGAF"):
	params = 6
elif (model == "EGAFT"):
	params = 8

params = params + adj
print("Params: %d" % (params))
if (params <= 0):
	print("Incorrect model specifier")
	exit()



#params = int(sys.argv[2])
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


		# from 'Solid-State NMR Provides Evidence for Small-Amplitude Slow Domain Motions in a Multispanning Transmembrane α-Helical Protein'

		chisq = 0

		for i in range(0, len(calc_R)):
			c = np.power((exp_R[i] - calc_R[i]), 2.) / np.power(err_R[i], 2.)
			chisq += c

		df = len(calc_R) - params - 1
		if (df <= 0 ):
			AIC = 1e9
			BIC = 1e9
			AICc = 1e9
		else:
			AIC = chisq + 2 * params
			BIC = chisq + params * np.log(len(calc_R))
			AICc = AIC + ((2 * params * (params + 1)) / df)


		#print("%d, %f, %f\n"% (i, AIC, BIC))
		q.write("%d, %f, %f, %f\n"% (i, AIC, BIC, AICc))
