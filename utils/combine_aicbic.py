## https://en.wikipedia.org/wiki/Maximum_likelihood_estimation - equation under "Continuous distribution, continuous parameter space" - normal distribution PDF.

# https://www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html

import sys
import numpy as np
# arguments are python aicbic.py FOLDER modeltype OUTPUTFILE
# then we read FOLDER/backcalc_N.dat from N=1-56
# for each, 

#arguments are python combine_aicbic.py smf smft demf demft... etc.


modellist = sys.argv[1:]
modelbins = np.zeros((len(modellist), 3))
num_mods = len(sys.argv) - 1

results = np.zeros((56, num_mods, 3))

c = 0
for model in sys.argv[1:]:
	adj = 0
	params = -10
	modeln = model
	while (modeln[0] == 'v' or modeln[0] == 'r'):
		if (modeln[0] == 'v'):
			adj += 3
		else:
			adj += 3 # as only fitting 15N for now.
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


	params = params + adj
	print("%s Params: %d" % (model, params))
	if (params <= 0):
		print("Incorrect model specifier")
		exit()

	for i in range(1, 57):
		x = np.loadtxt("%s/backcalc_%d.dat" % (model, i))
		if (np.size(x) == 0):
			results[i-1, c, :] = [-1, -1, -1]
			continue
		calc_R = x[:, 1]
		exp_R = x[:, 2]
		err_R = x[:, 3]
		n_data = np.size(calc_R)
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

		if (AIC < 0 or AIC > 100000):
			AIC = -1
		if (BIC < 0 or BIC > 100000):
			BIC = -1
		if (AICc < 0 or AICc > 100000):
			AICc = -1

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

		if (-1 not in results[i, :, 0]):
			aic.write(", %s\n" % (modellist[aicmin]))
			modelbins[aicmin, 0] = modelbins[aicmin, 0] + 1
		else:
			aic.write(", \n")

		if (-1 not in results[i, :, 1]):
			bic.write(", %s\n" % (modellist[bicmin]))
			modelbins[bicmin, 1] = modelbins[bicmin, 1] + 1
		else:
			bic.write(", \n")

		if (-1 not in results[i, :, 2]):
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
