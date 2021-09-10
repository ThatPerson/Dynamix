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

s2weights = [1, 0, 1, 1]  ## how much extra weight to put on S2 parameters?
 # eg in Lamley I believe the weight of these was increased (as it is for fitting) but I don't know if this is also done for the analysis.

modellist = []
for i in sys.argv[1:]:
	if (i[0] == '-'):
		sel = i[1:]
		if (len(sel) != 4):
			print("sel parameter should be NH, CH, CN, CC")
			exit(-1)
		s2weights = [int(p) for p in list(i)]
	modellist.append(i)
#modellist = sys.argv[1:]
modelbins = np.zeros((len(modellist), 3))
num_mods = len(sys.argv) - 1

results = np.zeros((56, num_mods, 3))

c = 0
for model in sys.argv[1:]:
	adj = 0
	params = -10
	modeln = model
	'''while (modeln[0] == 'v' or modeln[0] == 'r'):
		if (modeln[0] == 'v'):
			adj += 3
		else:
			adj += 2 # as only fitting 15N for now.
		modeln = modeln[1:]


	if (modeln == "smf"):
		params = 2
	elif (modeln == "emf"):
		params = 3
	elif (modeln == "emft"):
		params = 5
	elif (modeln == "smft"):
		params = 3
	elif (modeln == "demf"):
		params = 4
	elif (modeln == "demft"):
		params = 6
	elif (modeln == "gaf"):
		params = 8
	elif (modeln == "gaft"):
		params = 10
	elif (modeln == "egaf"):
		params = 6
	elif (modeln == "egaft"):
		params = 8


	params = params + adj'''
	if (model in models.mds):
		params = models.mds[model]['n']
	else:
		exit(-1)
	print("%s Params: %d" % (model, params))
	if (params <= 0):
		print("Incorrect model specifier")
		exit()

	orderparams = np.loadtxt("%s/orderparams.dat" % (model))	

	for i in range(1, 57):
		x = np.loadtxt("%s/backcalc_%d.dat" % (model, i))
		if (np.size(x) == 0):
			results[i-1, c, :] = [-1, -1, -1]
			continue
		calc_R = x[:, 1]
		ign = 0
		
		exp_R = x[:, 3]
		
		cr = exp_R[exp_R > 0]
		if (len(cr) < 5):
			ign = 1
		
		err_R = x[:, 4]
		n_data = np.size(calc_R)
		sigma = err_R / 2.


		# from 'Solid-State NMR Provides Evidence for Small-Amplitude Slow Domain Motions in a Multispanning Transmembrane α-Helical Protein'

		chisq = 0

		for il in range(0, len(calc_R)):
			ctmp = np.power((exp_R[il] - calc_R[il]), 2.) / np.power(err_R[il], 2.)
			chisq += ctmp

		# order parameters
		## for now, for GAF and GAFT models I'm using S2NH, S2CH and S2CN.
		## in order to allow for direct comparison with the other models
		## note that their contribution from S2NH is tripled - I don't know how
		## statistically valid this is? But as we're not using reduced chisq
		## the altnerative will cause chisq to be much bigger for them.



		ops = orderparams[orderparams[:, 0] == i, :]
		#print(ops)
		ops_s = np.shape(ops)
		#print(ops_s)
		max_op = 4
		
		for inc in range(0, max_op): 
			calc = ops[0, 1 + (inc*3)]
			exp = ops[0, 2 + (inc*3)]
			if (exp <= 0):
				ign = 1
			err = ops[0, 3 + (inc*3)]
			ctmp = np.power(exp - calc, 2.) / np.power(err, 2.)
			chisq += s2weight[inc] * ctmp

		N_meas = sum(s2weight) + len(calc_R)
		df = N_meas - params - 1
		if (df <= 0 or chisq == 0 or ign == 1):
			AIC = 1e9
			BIC = 1e9
			AICc = 1e9
			results[i-1, c, :] = [AIC, BIC, AICc]
			continue

		AIC = chisq + 2 * params
		BIC = chisq + params * np.log(N_meas)
		AICc = AIC + ((2 * params * (params + 1)) / df)

		## errors are 2sigma (2 standard deviations). So sigma is err_R/2
		#expFun = -np.power(calc_R - exp_R, 2.) / (2 * np.power(sigma, 2.))
		#divisor = np.sqrt(2 * np.pi * np.power(sigma, 2.))

		#errF = np.exp(expFun) / divisor
		
		#LerrF = np.log(errF)
		
		
		#for qw in range(0, np.size(calc_R)):
		#	if (calc_R[qw] < 0): # no fit
		#		LerrF[qw] = 0
		#		n_data = n_data - 1
		
		#LerrF[LerrF < 1E-308] = 0 # not entirely sure how valid this is
		#print(LerrF)
		#logv = np.sum(LerrF)
		
		#AIC = 2 * params - 2 * logv
		#BIC = params * np.log(n_data) - 2 * logv

		#adjusted AIC
		#if (n_data - params - 1 > 0):
		#	corr_f = (2 * params * (params + 1)) / (n_data - params - 1)
		#	AICc = AIC + corr_f
		#else:
		#	AIC = 1e9#
		#	BIC = 1e9
		#	AICc = 1e9

		if (AIC < 0 or AIC > 1000000):
			AIC = 1e9 
		if (BIC < 0 or BIC > 1000000):
			BIC = 1e9
		if (AICc < 0 or AICc > 1000000):
			AICc = 1e9

		results[i-1, c, :] = [AIC, BIC, AICc]
	c = c + 1

with open("AIC.csv", "w") as aic, open('BIC.csv', 'w') as bic, open('AICc.csv', 'w') as aicc:
	for models in sys.argv[1:]:
		aic.write(", %s" % (models))
		bic.write(", %s" % (models))
		aicc.write(", %s" % (models))

	aic.write("\n")
	bic.write("\n")
	aicc.write("\n")
	for i in range(0, 56):
		aic.write("%d" % (i+1))
		bic.write("%d" % (i+1))
		aicc.write("%d" % (i+1))
		for k in range(0, num_mods):
			aic.write(", %f" % (results[i, k, 0]))
			bic.write(", %f" % (results[i, k, 1]))
			aicc.write(", %f" % (results[i, k, 2]))

		aicmin = np.argmin(results[i, :, 0])
		bicmin = np.argmin(results[i, :, 1])
		aiccmin = np.argmin(results[i, :, 2])

		if (-1 not in results[i, :, 0] and sum(results[i, :, 0]) != num_mods * 1e9):
			aic.write(", %s\n" % (modellist[aicmin]))
			modelbins[aicmin, 0] = modelbins[aicmin, 0] + 1
		else:
			aic.write(", \n")

		if (-1 not in results[i, :, 1] and sum(results[i, :, 1]) != num_mods * 1e9):
			bic.write(", %s\n" % (modellist[bicmin]))
			modelbins[bicmin, 1] = modelbins[bicmin, 1] + 1
		else:
			bic.write(", \n")

		if (-1 not in results[i, :, 2] and sum(results[i, :, 2]) != num_mods * 1e9):
			aicc.write(", %s\n" % (modellist[aiccmin]))
			modelbins[aiccmin, 2] = modelbins[aiccmin, 2] + 1
		else:
			aicc.write(", \n")

print("Model Choices")
print("=============")

print("\t\tAIC\tBIC\tAICc ")
for i in range(0, len(modellist)):
	print("\t%s:\t%d\t%d\t%d" % (modellist[i], modelbins[i, 0], modelbins[i, 1], modelbins[i, 2]))




#params = int(sys.argv[2])
'''with open(sys.argv[3], "w") as q:
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

		#adjusted AIC
		if (n_data - params - 1 > 0):
			corr_f = (2 * params * (params + 1)) / (n_data - params - 1)
			AICc = AIC + corr_f
		else:
			AIC = -1
			BIC = -1
			AICc = -1

		#print("%d, %f, %f\n"% (i, AIC, BIC))
		q.write("%d, %f, %f, %f\n"% (i, AIC, BIC, AICc))'''
