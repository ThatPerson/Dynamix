import numpy as np
import sys
import matplotlib.pyplot as plt
import models

try:
	ic_file = sys.argv[1]
	output_file = sys.argv[2]
	ic_type = sys.argv[3]
except IndexError:
	print("python plot_aicbic.py <IC file> <output file> <IC type>" )
	exit(-1)
	
with open(ic_file, "r") as f:
	for l in f:
		k = l.split(",")
		mls = k[1:]
		mls = [p.strip() for p in mls]
		for i in mls:
			if (i not in models.mds):
				print("Unknown model %s" % (i))
		break

mdict = {}
mdict["smf"] = "MF"
mdict["smft"] = "MF / T"
mdict["bgaf"] = "GAF"
mdict["bgaft"] = "GAF / T"
mdict["vbgaf"] = "GAF / V"
mdict["vbgaft"] = "GAF / V + T"
mdict["demf"] = "MF + MF"
mdict["demft"] = "MF + MF / T"
mdict["egaf"] = "GAF + MF"
mdict["egaft"] = "GAF + MF / T"
mdict["vegaf"] = "GAF + MF / V"
mdict["vegaft"] = "GAF + MF / V + T"
mdict["gaf"] = "GAF + GAF"
mdict["gaft"] = "GAF + GAF / T"
mdict["vgaf"] = "GAF + GAF / V"
mdict["vgaft"] = "GAF + GAF / V + T"

mls_pars = [models.mds[p]["n"] for p in mls]
mls_ticks = [(mdict[m] if (m in mdict) else m) for m in mls]
print(mls_pars)

#plt.figure(figsize=(10,8))
#plt.gcf().set_size_inches(10, 8)
#fig, axes = plt.subplots(5,2, sharex='col')
#fig.set_size_inches(10, 8)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['font.sans-serif'] = 'Cantarell'
font = {'fontname': 'Cantarell'}

plt.figure(figsize=(10,8))
plt.gcf().set_size_inches(10, 8)
fig, ax = plt.subplots(1,1, sharex='col')
fig.set_size_inches(10, 8)

data = np.loadtxt(ic_file, skiprows=1, delimiter=",", usecols=np.arange(1, 1+len(mls)))
data[data > 1e5] = np.nan
mask = ~np.isnan(data)
filtered_data = [d[m] for d, m in zip(data.T, mask.T)]
print(data)

k = np.shape(data)
X = np.arange(1, k[1]+1)
for i in range(0, k[0]):
	if (np.nan not in data[i, :]):
		plt.plot(X, data[i, :], 'k--', alpha=.1)
		

plt.boxplot(filtered_data, showfliers = False)
#plt.violinplot(filtered_data, showmedians = True, showmeans=True, showextrema=True)
for i in range(0, len(mls)):
	print(i)
	y = filtered_data[i]
	x = np.random.normal(i+1, 0.04, size=len(y))
	fmt = '.'
	color = 'tab:blue'
	print(mls_ticks[i])
	if ("V" in mls_ticks[i]):
		color = 'tab:orange'
	if ("T" in mls_ticks[i]):
		color = 'tab:green'
	if ("V + T" in mls_ticks[i]):
		color = 'tab:red'
	plt.plot(x, y, fmt, color=color, alpha=1)



mlst = []

for i in range(0, len(mls)):
	mlst.append("%s (%d)" % (mls_ticks[i].upper(), mls_pars[i]))

ax.set_xticks(np.arange(1, len(mls)+1))
#ax.set_yscale("log")
plt.xticks(rotation=90)
ax.set_xticklabels(mlst)
plt.ylim([0, 4000])
ax.set_xlabel("Model (Parameters)")
ax.set_ylabel("%s" % (ic_type))
#ax.set_yscale("log")
plt.savefig(output_file, bbox_inches="tight")
