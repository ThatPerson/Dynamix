import numpy as np
import sys

try:
	bmrb_fn = sys.argv[1]
	csiso_C = sys.argv[2]
	csiso_N = sys.argv[3]
except IndexError:
	print("bmrb_to_csiso.py <bmrb fn> <csiso C> <csiso N>")
	exit(-1)
	
with open(bmrb_fn, "r") as bmrb, open(csiso_C, "w") as csC, open(csiso_N, "w") as csN:
	for l in bmrb:
		k = l.split()
		if (len(k) != 24):
			continue
		resid = int(k[5])
		atom = k[7]
		shift = float(k[10])
		shift_err = float(k[11])
		if (atom == 'C'):
			csC.write("%d %f %f\n" % (resid, shift, shift_err))
		elif (atom == 'N'):
			csN.write("%d %f %f\n" % (resid, shift, shift_err))
	  

'''
0: 1
1: .
2: 1
3: 1
4: 1
5: 1
6: MET
7: HA
8: H
9: 1
10: 4.125
11: 0.05
12: .
13: 1
14: .
15: .
16: .
17: A
18: 1
19: MET
20: HA
21: .
22: 30088
23: 1
'''
