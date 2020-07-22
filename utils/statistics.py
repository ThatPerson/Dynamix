import numpy as np
import scipy.stats as sp
import sys
import matplotlib.pyplot as plt
if (len(sys.argv)<=2):
	print("No directory")
	exit()
	
direc = sys.argv[1]
n_resid = int(sys.argv[2])

n_pars = 1
for i in (3, len(sys.argv)-1):
	n_pars = n_pars * int(sys.argv[i])
	
#print("N pars %d" % (n_pars))


for i in range(1, n_resid+1):
	k = np.loadtxt("%s/backcalc_%d.dat" % (direc, i))

	calc = k[k[:, 1] > 0, 1]
	real = k[k[:, 1] > 0, 2]


	#resid = y - y_hat
	#sse = sum(resid**2)


	resid = real - calc
	sse = sum(resid ** 2)

	#print(sse)


	# Not actually AIC, BIC. I don't have a Likelihood function, just SSE.
	# A large SSE would indicate a poor model, a small SSE a good model.
	# So I've inverted the sign so that a large SSE disfavours the model (eg higher IC)

	AIC = 2 * n_pars + 2 * np.log(sse)
	#print("AIC: %f" % (AIC))

	BIC = n_pars * np.log(len(calc)) + 2 * np.log(sse)
	#print("BIC: %f" % (BIC))
	if (AIC < 0 or BIC < 0):
		AIC = -1
		BIC = -1
	print("%d, %f, %f, %f" % (i, AIC, BIC, sse))
