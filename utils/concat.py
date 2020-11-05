import numpy as np
import sys

data = np.zeros((56, len(sys.argv)-1))
data[:, 0] = np.arange(1, 57)
c = 1


for i in sys.argv[1:-1]:
	x = np.loadtxt(i, delimiter=',')
	data[:, c] = x[:, 3]
	c = c + 1
print(sys.argv[-1])
	
np.savetxt(sys.argv[-1], data, fmt="%f", delimiter='\t')
