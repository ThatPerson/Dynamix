## https://en.wikipedia.org/wiki/Maximum_likelihood_estimation - equation under "Continuous distribution, continuous parameter space" - normal distribution PDF.

# https://www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html

import sys
import numpy as np
import matplotlib.pyplot as plt
# arguments are python aicbic.py FOLDER modeltype OUTPUTFILE
# then we read FOLDER/backcalc_N.dat from N=1-56
# for each, 

def run_plot(folder):

	# TODO: read in final.dat length
	n = 0 #int(sys.argv[2])
	with open('%s/final.dat' % (folder), 'r') as f:
		for l in f:
			n = n + 1
	print(n)
	for i in range(1, n+1):
		try:
			k = np.loadtxt("%s/backcalc_%d.dat" % (folder, i))
			break
		except OSError:
			continue
	lk = (np.shape(k))

	n_relax = lk[0]

	data = np.zeros((n+1, n_relax, lk[1]))

	fig,axs = plt.subplots(int(np.floor(n_relax / 3.)) + 1, 3, figsize=(8, n_relax/1.2), dpi=80)



	curr_x = 0
	curr_y = 0

	for i in range(1, n+1):
		try:
			data[i, :, :] = np.loadtxt("%s/backcalc_%d.dat" % (folder, i))
		except OSError:
			continue

	x = np.arange(1, n+2)

	types = ['15N R1', '15N R1p', '13C R1', '13C R1p']

	for i in range(0, n_relax):
		print(np.shape(x))
		print(np.shape(data[:, i, 2]))

		R = data[:, i, 3]
		Rerr = data[:, i, 4]
		Rcalc = data[:, i, 1]
		Rcalcerr = data[:, i, 2]
		R = R[Rcalc > 0]
		Rerr = Rerr[Rcalc > 0]
		Rcalc = Rcalc[Rcalc > 0]

		axs[curr_x,curr_y].errorbar(x, data[:, i, 3], yerr=data[:, i, 4], fmt='k,')
		axs[curr_x,curr_y].errorbar(x, data[:, i, 1], yerr=data[:,i, 2], fmt='b,')
		axs[curr_x,curr_y].set_ylim(bottom=0)



		field = data[20, i, 5]
		wr = data[20, i, 6]/1000.
		w1 = data[20, i, 7]/1000.
		T = data[20, i, 8]
		print(data[20, i, :])
		typ =  int(data[20,i, 9])

		if (typ == 0 or typ == 2):
			krt = axs[curr_x, curr_y].get_ylim()
			if (krt[1] > 0.4):
				axs[curr_x, curr_y].set_ylim(top = 0.4)
		if (typ == 1 or typ == 3):
			krt = axs[curr_x, curr_y].get_ylim()
			if (krt[1] > 6):
				axs[curr_x, curr_y].set_ylim(top = 20)
			

		#fprintf(fp, "\t%f\t%f\t%f\t%f\t%d", resid->relaxation[i].field, resid->relaxation[i].wr, resid->relaxation[i].w1, resid->relaxation[i].T, resid->relaxation[i].type); 
		#	fprintf(fp, "\n");
		print("%s (%d MHz, %0.1f kHz, %0.1f kHz, %0.1f K)" % (types[typ], field, wr, w1,T))
		axs[curr_x,curr_y].set_title("%s (%d MHz, %0.1f kHz, %0.1f kHz, %0.1f K)" % (types[typ], field, wr, w1,T), fontsize=5)

		curr_y = curr_y + 1
		if (curr_y == 3):
			curr_y = 0
			curr_x = curr_x + 1

		#fprintf(fp, "%d\t%f\t%f\t%f", i, (calc_R<0?-1.:calc_R), resid->relaxation[i].R, resid->relaxation[i].Rerror);



	for ax in axs.flat:
	    ax.set(xlabel='residue', ylabel='rate / s^-1')

	fig.tight_layout(h_pad=1)
	#plt.show(block=True)

	plt.savefig('%s_relax.pdf' % (folder), dpi=300, facecolor='w', edgecolor='w',
	        orientation='portrait', papertype='a4', format='pdf',
	        transparent=False, bbox_inches=None, pad_inches=0.5,
	        frameon=None, metadata=None)

for i in sys.argv[1:]:
	run_plot(i)
