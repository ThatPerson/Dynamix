### todo: add in order parameter calculation

### note that you can't directly compare the chisq outputs from Dynamix because different models include different numbers of order paremters.
### (eg SMF and EMF only fit to S2NH as this is the dominant one. In principle they could be fit to the others too but this leads to worse fits).
### as the order parameter contribution is x10 it follows that this can lead to a huge difference. Instead you should use the output of this
### (which only fits to relaxation data note, not order parameters - need to add this!)


# https://www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html

import sys
import numpy as np
import models
# arguments are python aicbic.py FOLDER modeltype OUTPUTFILE
# then we read FOLDER/backcalc_N.dat from N=1-56
# for each, 

#arguments are python combine_aicbic.py smf smft demf demft... etc.


modellist = sys.argv[1:]
modelbins = np.zeros((len(modellist), 3))
num_mods = len(sys.argv) - 1

results = np.zeros((56, num_mods, 3))

with open('600.attr', 'w') as f:
	c = 0
	f.write("attribute: chisq7\n")
	f.write("match mode: 1-to-1\nrecipient: residues\n")
	for model in sys.argv[1:]:
		adj = 0
		params = -10
		modeln = model

		if (model in models.mds):
			params = models.mds[model]['n']
		else:
			exit(-1)
		print("%s Params: %d" % (model, params))
		if (params <= 0):
			print("Incorrect model specifier")
			exit()

		for i in range(1, 57):
			x = np.loadtxt("%s/backcalc_%d.dat" % (model, i))
			if (np.size(x) == 0):
				results[i-1, c, :] = [-1, -1, -1]
				continue

			# selective for 600 R1p	
			calc_R = x[7:14, 1]
			exp_R = x[7:14, 2]
			err_R = x[7:14, 3]
			n_data = np.size(calc_R)
			sigma = err_R / 2.


			# from 'Solid-State NMR Provides Evidence for Small-Amplitude Slow Domain Motions in a Multispanning Transmembrane α-Helical Protein'

			chisq = 0

			for il in range(0, len(calc_R)):
				ctmp = np.power((exp_R[il] - calc_R[il]), 2.) / np.power(err_R[il], 2.)
				chisq += ctmp

			df = 1 #len(calc_R) - params - 1
			print(df)
			print(chisq)
			AICc = 0
			if (df <= 0 or chisq == 0):
				AIC = 1e9
				BIC = 1e9
				AICc = 1e9
				results[i-1, c, :] = [AIC, BIC, AICc]
				continue

			AIC = chisq + 2 * params
			BIC = chisq + params * np.log(len(calc_R))
			#AICc = AIC + ((2 * params * (params + 1)) / df)


			if (AIC < 0 or AIC > 1000000):
				AIC = 1e9 
			if (BIC < 0 or BIC > 1000000):
				BIC = 1e9
			#if (AICc < 0 or AICc > 1000000):
			#	AICc = 1e9


			results[i-1, c, :] = [AIC, BIC, AICc]

			f.write("\t#1:%d\t%f\n" % (i, BIC))
		c = c + 1


