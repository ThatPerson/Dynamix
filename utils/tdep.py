import numpy as np
import sys
import models
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

fn = sys.argv[1]

x = np.loadtxt("%s/final.dat" % (fn), delimiter='\t')

mod = models.mds[fn]['p']

#0,1,2,3

Xc = x[:, 0]
Yc = np.transpose(np.arange(1, 501))
slow_t = np.zeros((56, 500))
fast_t = np.zeros((56, 500))
tau_t = np.zeros((56, 500))



if ('taus' in mod):
	index = mod.index('taus')
	t = x[:, 3 + index]
	if ('Eas' in mod):
		index = mod.index('Eas')
		Ea = x[:, 3 + index]
	else:
		Ea = np.zeros((56))
	for temp in range(1, 501):
		slow_t[:, temp-1] = t * np.power(10., -9.) * np.exp(Ea / (8.3145 * temp))
	print(t)
else:
	print("No slow")

if ('tauf' in mod):
	index = mod.index('tauf')
	t = x[:, 3 + index]
	if ('Eaf' in mod):
		index = mod.index('Eaf')
		Ea = x[:, 3 + index]
	else:
		Ea = np.zeros((56))
	for temp in range(1, 501):
		fast_t[:, temp-1] = t * np.power(10., -9.) * np.exp(Ea / (8.3145 * temp))
else:
	print("No fast")

if ('tau' in mod):
	index = mod.index('tau')
	t = x[:, 3 + index]
	if ('Ea' in mod):
		index = mod.index('Ea')
		Ea = x[:, 3 + index]
	else:
		Ea = np.zeros((56))
	for temp in range(1, 501):
		tau_t[:, temp-1] = t * np.power(10., -9.) * np.exp(Ea / (8.3145 * temp))
else:
	print("no tau")
print(np.shape(Xc))
print(np.shape(Yc))
print(np.shape(slow_t))
print(slow_t)
slow_t[slow_t < 0] = 0
slow_t[slow_t > 1] = 1
fast_t[fast_t < 0] = 0
fast_t[fast_t > 1] = 1
print(slow_t)
plt.figure()
plt.yscale('log')
with open('d.csv', 'w') as f:
	for i in range(0, 56):
		plt.plot(Yc, slow_t[i, :])
		plt.plot(Yc, fast_t[i, :])
		for t in range(0, 500):
			f.write("%d, %d, %0.20f, %.20f\n" % (i, t + 1, slow_t[i, t], fast_t[i, t]))
#ax.plot_wireframe(Xc, Yc, slow_t)
plt.show(block=True)


