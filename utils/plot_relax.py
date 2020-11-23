## https://en.wikipedia.org/wiki/Maximum_likelihood_estimation - equation under "Continuous distribution, continuous parameter space" - normal distribution PDF.

# https://www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html

import sys
import numpy as np
import matplotlib.pyplot as plt
# arguments are python aicbic.py FOLDER modeltype OUTPUTFILE
# then we read FOLDER/backcalc_N.dat from N=1-56
# for each, 

folder = sys.argv[1]
n = int(sys.argv[2])

k = np.loadtxt("%s/backcalc_2.dat" % (folder))
lk = (np.shape(k))

n_relax = lk[0]

data = np.zeros((n+1, n_relax, lk[1]))

fig,axs = plt.subplots(int(np.floor(n_relax / 3.)) + 1, 3, figsize=(12, 30), dpi=80)



curr_x = 0
curr_y = 0

for i in range(1, n+1):
	data[i, :, :] = np.loadtxt("%s/backcalc_%d.dat" % (folder, i))

x = np.arange(1, n+2)

types = ['15N R1', '15N R1p', '13C R1', '13C R1p']

for i in range(0, n_relax):
	print(np.shape(x))
	print(np.shape(data[:, i, 2]))

	R = data[:, i, 2]
	Rerr = data[:, i, 3]
	Rcalc = data[:, i, 1]
	R = R[Rcalc > 0]
	Rerr = Rerr[Rcalc > 0]
	Rcalc = Rcalc[Rcalc > 0]

	axs[curr_x,curr_y].errorbar(x, data[:, i, 2], yerr=data[:, i, 3], fmt='k,')
	axs[curr_x,curr_y].plot(x, data[:, i, 1], 'b,')
	axs[curr_x,curr_y].set_ylim(bottom=0)

	field = data[20, i, 4]
	wr = data[20, i, 5]/1000.
	w1 = data[20, i, 6]/1000.
	T = data[20, i, 7]
	print(data[20, i, :])
	typ =  int(data[20,i, 8])



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

plt.savefig('%s_backcalc.eps' % (folder), dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='a4', format='eps',
        transparent=False, bbox_inches=None, pad_inches=0.5,
        frameon=None, metadata=None)


