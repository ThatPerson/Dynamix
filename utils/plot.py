import numpy as np
import sys
import models
import matplotlib.pyplot as plt

# takes arguemtns -f<folder> -m<model> -e



errors = False
model = ""
folder = ""
folders = []
for i in sys.argv:

	if (i[:len('-m')] == '-m'):
		model = (i[2:])
	elif (i[:len('-e')] == '-e'):
		errors = True
	else:
		folders.append(i)
	

if (len(folders) == 0):
	print("Please give folder name")
	exit(-1)
model_notset = False
if (model == ''):
	model_notset = True

def run_plot(folder):	
	if (model_notset == True):
		model = folder
		print("Model %s" % (model))
	if (model not in models.mds):
		print("Model %s not defined" % (model))
		return 1

	mod = models.mds[model]
	print(mod)

	xdim = (mod['n'] + 1)/2
	xdim = int(xdim) + 1
	xdim = mod['n'] + 1
	print(xdim)

	fig,axs = plt.subplots(xdim, 1, figsize=(7., mod['n'] * 2), dpi=80)

	if (errors == True):
		x = np.loadtxt("%s/errors.dat" % (folder), delimiter='\t')
		print(np.shape(x))
		print(mod['n'])
	else:
		y = np.loadtxt("%s/final.dat" % (folder), delimiter='\t')
		k = np.shape(y)
		x = np.zeros((k[0], k[1] + mod['n']))
		x[:, 0:3] = y[:, 0:3]
		print(np.shape(x))
		for counter in range(0, mod['n']):
			x[:, 3 + (counter * 2) ] = y[:, 3 + counter]
		print(x) 
		print(np.shape(x))

	x[x == -1] = np.nan

	# now for plotting
	xdat = x[:, 0]
	ydat = x[:, 2] 
	axs[0].plot(xdat, ydat, 'k,')
	axs[0].set_ylim(bottom=0)
	axs[0].set_title("Chisq")

	xd = 1
	yd = 1

	for pl in range(0, mod['n']):
		print("%d, %d"% (xd, yd))
		xdat = x[:, 0]
		#if (errors == True):
		ydat = x[:, 3 + (pl * 2)]
		yerr = x[:, 4 + (pl * 2)]
	
		#else:

		#	ydat = x[:, 3 + pl]
		#	print("%s : %d : %f" % (mod['p'][pl], 3 + pl, np.mean(ydat)))
		#	yerr = np.zeros(np.shape(ydat))
			
		if ("tau" in mod['p'][pl]):
			ydat = ydat * 1e-9
			yerr = yerr * 1e-9
		if (errors == True):
			axs[xd].errorbar(xdat, ydat, yerr=yerr, fmt='k,')
		else:
			axs[xd].plot(xdat, ydat, 'k+')
		
		axs[xd].set_title("%s" % (mod['p'][pl]))
		axs[xd].set(xlabel='residue', ylabel=mod['p'][pl])
		if ("tau" in mod['p'][pl]):
			axs[xd].set_yscale("log")
			if ('taus' in mod['p'][pl]):
				# temperature dependent
				axs[xd].set_ylim(bottom=1e-9, top=1e-7)
			elif ('tauf' in mod['p'][pl]):
				# temperature dependent
				axs[xd].set_ylim(bottom=1e-13, top=1e-9)
			elif ('tau' in mod['p'][pl]):
				axs[xd].set_ylim(bottom=1e-13, top=1e-7)
			'''if ('taus' in mod['p'][pl] and 't' in model):
				# temperature dependent
				axs[xd].set_ylim(bottom=1e-18, top=1e-8)
			elif ('tauf' in mod['p'][pl] and 't' in model):
				# temperature dependent
				axs[xd].set_ylim(bottom=1e-18, top=1e-8)
			elif ('taus' in mod['p'][pl] and 't' not in model):
				# temperature independent
				axs[xd].set_ylim(bottom=1e-9, top=1e-6)
			elif ('tauf' in mod['p'][pl] and 't' not in model):
				# temperature independent
				axs[xd].set_ylim(bottom=1e-12, top=1e-9)
			elif ('tau' in mod['p'][pl] and 't' in model):
				# temperature dependent
				axs[xd].set_ylim(bottom=1e-12, top=1e-8)
			elif ('tau' in mod['p'][pl] and 't' not in model):
				# temperature independent
				axs[xd].set_ylim(bottom=1e-9, top=1e-6)'''

		if ("Ea" in mod['p'][pl]):
			axs[xd].set_ylim(bottom=0, top = 60000)
		if ('S2' in mod['p'][pl] and "dS2" not in mod['p'][pl]):
			axs[xd].set_ylim(top=1)
			kl = axs[xd].get_ylim()
			if (kl[0] < 0.2):
				axs[xd].set_ylim(bottom=0.5)
		

		if ('sA' in mod['p'][pl] or 'sB' in mod['p'][pl] or 'sG' in mod['p'][pl]):
			axs[xd].set_ylim(bottom=0)
		if (mod['p'][pl] in ['alph', 'beta', 'gamm']):
			axs[xd].set_ylim(bottom=0)
		if (mod['p'][pl] in ['papbN', 'papbC']):
			axs[xd].set_ylim(bottom=0, top=2*np.nanmedian(ydat))
		if (mod['p'][pl] in ['kex']):
			axs[xd].set_ylim(bottom=0, top=2*np.nanmedian(ydat))
		if ('600del' in mod['p'][pl]):
			print("T Delta")
			print(ydat)
			plw = np.nanstd(ydat)
			axs[xd].set_ylim(-100, 100)
			axs[xd].set_xlim(0, 56)
			
		#yd = yd + 1
		#if (yd > 1):
		#	xd = xd + 1
		#	yd = 0
		xd = xd+1

		#fprintf(fp, "%d\t%f\t%f\t%f", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);

	fig.tight_layout(h_pad=1)
	#plt.show(block=True)

	plt.savefig('%s_params.pdf' % (folder), dpi=300, facecolor='w', edgecolor='w',
	        orientation='portrait', papertype='a4', format='pdf',
	        transparent=False, bbox_inches=None, pad_inches=0.5,
	        frameon=None, metadata=None)

for f in folders:
	run_plot(f)
