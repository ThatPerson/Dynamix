import numpy as np
import matplotlib.pyplot as plt
import sys

try:
	fn = sys.argv[1]
	ofn = sys.argv[2]
except IndexError:
	print("python gaffit_plot.py <fn> <ofn>")
	exit(-1)

data = np.loadtxt(fn) * (180 / np.pi)
k = np.shape(data)
nsite = k[1]
fig, axes = plt.subplots(nsite, nsite)
fig.set_size_inches((12, 12))

columns = ["$\\sigma_{\\alpha}$ / $^{\circ}$", "$\\sigma_{\\beta}$ / $^{\circ}$", "$\\sigma_{\\gamma}$ / $^{\circ}$", "$\\alpha$ / $^{\circ}$", "$\\beta$ / $^{\circ}$", "$\\gamma$ / $^{\circ}$"]
sl = [[0, 10], [0,10], [0,10], [-25, 25], [-25, 25], [-25, 25]]

for sitex in range(0, nsite):
	axes[sitex, 0].set_ylabel(columns[sitex])
	axes[-1, sitex].set_xlabel(columns[sitex])
	
	
	for sitey in range(0, nsite):
		if (sitex != sitey): axes[sitex, sitey].plot(data[:, sitey], data[:, sitex], marker='o', markerfacecolor="white", markeredgecolor="red", linestyle='none')
		axes[sitex, sitey].set_xlim(sl[sitex])
		if (sitey < 3):
			axes[sitex, sitey].set_xlim([0,10])
		else:
			axes[sitex, sitey].set_xlim([-25, 25])
		
		if (sitex < 3):
			axes[sitex, sitey].set_ylim([0,10])
		else:
			axes[sitex, sitey].set_ylim([-25, 25])
		axes[sitex, sitey].set_xticks([])
		axes[sitex, sitey].set_yticks([])
for sitex in range(0, nsite):
	if (sitex < 3):
		axes[sitex, 0].set_yticks([2.5, 7.5])
		axes[-1, sitex].set_xticks([2.5, 7.5])
	else:
		axes[sitex, 0].set_yticks([-20, 0, 20])
		axes[-1, sitex].set_xticks([-20, 0, 20])
		
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.savefig(ofn, bbox_inches="tight", dpi=300)
