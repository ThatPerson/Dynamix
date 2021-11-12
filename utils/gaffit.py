import numpy as np
import argparse as ap
import gafpy as GAF
import scipy.optimize as so

parser = ap.ArgumentParser()
parser.add_argument('-S2', type=str, default=None, help='File containing order parameter.')
parser.add_argument('-column', type=int, default=1, help='Column in file to read order parameter from (starts at 1).')
parser.add_argument('-columnerrors', type=int, default=-1, help='Column in file to read order parameter errors from (starts at 1).')
parser.add_argument('-orientation', type=str, default=None, help='File containing orientations of order parameters.')
parser.add_argument('-outputS2', type=str, default="gaf.S2all", help='File to output backcalculated order parameters into')
parser.add_argument('-outputS2sum', type=str, default="gaf.S2", help='File to output backcalculated order parameters into')
parser.add_argument('-outputfit', type=str, default="gaf.fit", help='File to output fit parameters into')
parser.add_argument('-outputfitstats', type=str, default="gaf.stats", help='File to output fit statistics into')
parser.add_argument('-model', type=str, default="GAF", help='Model to fit (GAF/iso/AIMF)')
parser.add_argument('-mc', type=int, default=50, help='Monte Carlo iterations')
args = parser.parse_args()
if (args.S2 == None or args.orientation == None):
	parser.print_help()
	exit(-1)

def calc_GAF(theta, phi, sigA, sigB, sigG, alpha, beta, gamma):
	thetap, phip = GAF.rotate_orients(theta, phi, alpha, beta, gamma)
	S2 = GAF.GAF(sigA, sigB, sigG, phip, thetap, phip, thetap)
	return S2

def objf_GAF(x, S2exp, S2error, theta, phi):
	sigA, sigB, sigG, alpha, beta, gamma = x
	S2calc = calc_GAF(theta, phi, sigA, sigB, sigG, alpha, beta, gamma)
	#print(np.shape(S2exp))
#	print(np.shape(S2calc))
#	print(np.shape(S2error))
	chisq = np.sqrt(np.nansum(np.power(S2exp - S2calc, 2) / np.power(S2error, 2)))
	#print(chisq)
	return chisq

def objf_ISO(x, S2exp, S2error, theta, phi):
	S2calc = x
	return np.sqrt(np.nansum(np.power(S2exp - S2calc, 2) / np.power(S2error, 2)))
	
orients = np.loadtxt(args.orientation)
res = orients[:, 0]
theta = orients[:, 1]
phi = orients[:, 2]
if (args.columnerrors == -1):
	S2 = np.zeros((len(res, 2))) + 0.05
	d = np.loadtxt(args.S2, usecols=[args.column - 1])
	#print(d)
	S2[:, 0] = d
else:
	S2 = np.loadtxt(args.S2, usecols=[args.column - 1, args.columnerrors - 1])

S2[S2[:, 0] < 0, :] = [np.nan, np.nan]

if (args.model == "iso"):
	x0 = [0.9]
	print("Running main fit...")
	results = so.minimize(objf_ISO, x0, args=(S2[:, 0], S2[:, 1], theta, phi))
	S2bc = np.zeros((len(res), args.mc))
	S2bcn = np.zeros((len(res))) + results['x'][0]
	
	parnames = ["S2"]
	angular = [False]
	
	fititeration = np.zeros((args.mc, len(results['x'])))
	for iteration in range(0, args.mc):
		print("Running Monte Carlo %d" % (iteration))
		S2curr = np.random.normal(S2bcn, S2[:, 1]/2) # assuming S2[:, 1] is 2 standard deviations
		fit = so.minimize(objf_ISO, results['x'], args=(S2curr, S2[:, 1], theta, phi))
		fititeration[iteration, :] = fit['x']
		S2bc[:, iteration] += fit['x'][0]
	
elif (args.model == "GAF"):
	x0 = np.random.random((6)) * 0.15
	print("Running main fit...")
	results = so.minimize(objf_GAF, x0, args=(S2[:, 0], S2[:, 1], theta, phi), method="Powell", bounds=((0, None), (0, None), (0, None), (None, None), (None, None), (None, None)))	
	S2bc = np.zeros((len(res), args.mc))
	S2bcn = calc_GAF(theta, phi, results['x'][0], results['x'][1], results['x'][2], results['x'][3], results['x'][4], results['x'][5])
	fititeration = np.zeros((args.mc, len(results['x'])))
	parnames = ["sigA", "sigB", "sigG", "alpha", "beta", "gamma"]
	angular = [True, True, True, True, True, True]
	for iteration in range(0, args.mc):
		print("Running Monte Carlo %d" % (iteration))
		S2curr = np.random.normal(S2bcn, S2[:, 1]/2) # assuming S2[:, 1] is 2 standard deviations
		fit = so.minimize(objf_GAF, results['x'], args=(S2curr, S2[:, 1], theta, phi), bounds=((0, None), (0, None), (0, None), (None, None), (None, None), (None, None)))
		fititeration[iteration, :] = fit['x']
		S2bc[:, iteration] = calc_GAF(theta, phi, fit['x'][0], fit['x'][1], fit['x'][2], fit['x'][3], fit['x'][4], fit['x'][5])
	#S2bc = calc_GAF(theta, phi, 0.1, 0.2, 0, 0.4, 0.6, 0.1)	

elif (args.model == "GAFM"):
	x0 = np.random.random((6)) * 0.15
	print("Running main fit...")
	results = so.minimize(objf_GAF, x0, args=(S2[:, 0], S2[:, 1], theta, phi), method="Powell", bounds=((0, None), (0, None), (0, None), (None, None), (None, None), (-0.001, 0.001)))	
	S2bc = np.zeros((len(res), args.mc))
	S2bcn = calc_GAF(theta, phi, results['x'][0], results['x'][1], results['x'][2], results['x'][3], results['x'][4], results['x'][5])
	fititeration = np.zeros((args.mc, len(results['x'])))
	parnames = ["sigA", "sigB", "sigG", "alpha", "beta", "gamma"]
	angular = [True, True, True, True, True, True]
	for iteration in range(0, args.mc):
		print("Running Monte Carlo %d" % (iteration))
		S2curr = np.random.normal(S2bcn, S2[:, 1]/2) # assuming S2[:, 1] is 2 standard deviations
		fit = so.minimize(objf_GAF, results['x'], args=(S2curr, S2[:, 1], theta, phi), bounds=((0, None), (0, None), (0, None), (None, None), (None, None), (-0.001, 0.001)))
		fititeration[iteration, :] = fit['x']
		S2bc[:, iteration] = calc_GAF(theta, phi, fit['x'][0], fit['x'][1], fit['x'][2], fit['x'][3], fit['x'][4], fit['x'][5])
	#S2bc = calc_GAF(theta, phi, 0.1, 0.2, 0, 0.4, 0.6, 0.1)

data = np.zeros((len(S2bc), 3))
data[:, 0] = res
data[:, 1] = np.mean(S2bc, axis=1)
data[:, 2] = 2 * np.std(S2bc, axis=1)
np.savetxt(args.outputS2sum, data)
np.savetxt(args.outputS2, S2bc)
np.savetxt(args.outputfit, fititeration)

fit_mean = np.mean(fititeration, axis=0)
fit_std = np.std(fititeration, axis=0)

with open(args.outputfitstats, "w") as f:
	for mean, std, name, ang in zip(fit_mean, fit_std, parnames, angular):
		f.write("%s: %0.5f +- %0.5f " % (name, mean, 2*std))
		if (ang == True):
			f.write("(%0.2f +- %0.2f)" % (mean * (180./np.pi), 2 * std * (180./np.pi)))
		f.write("\n")
	
print(results)
