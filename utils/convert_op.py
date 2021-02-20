## convert from old style order parameter files for demft (eg all experimental being S2NH) to new style.

#fprintf(orderparams, "%d\t", l+1);
#fprintf(orderparams, "%f\t%f\t%f\t", S2, m->residues[l].S2NH, m->residues[l].S2NHe);
#fprintf(orderparams, "%f\t%f\t%f\t", S2, m->residues[l].S2CH, m->residues[l].S2CHe);
#fprintf(orderparams, "%f\t%f\t%f\t", S2, m->residues[l].S2CN, m->residues[l].S2CNe);
#fprintf(orderparams, "%f\t%f\t%f\n", S2, m->residues[l].S2CC, m->residues[l].S2CCe);

import sys
import numpy as np
try:
	in_file = sys.argv[1]
	out_file = sys.argv[2]
	system_file = sys.argv[3]
except IndexError:
	print("python convert_op.py INPUT OUTPUT SYSTEM")
	exit(-1)
	
S2NH = -1
S2CC = -1
S2CH = -1
S2CN = -1	
NH = -1 
CC = -1 
CN = -1 
CH = -1

with open(system_file, "r") as f:
	for l in f:
		k = l.split("=")
		k = [p.strip() for p in k]
		#print(k)
		if (k[0] == "S2NH"):
			print("Opening %s" % (k[1]))
			S2NH = np.loadtxt(k[1], delimiter=" ")
			NH = 1
		elif (k[0] == "S2CH"):
			print("Opening %s" % (k[1]))
			S2CH = np.loadtxt(k[1], delimiter=" ")
			CH = 1
		elif (k[0] == "S2CN"):
			print("Opening %s" % (k[1]))
			S2CN = np.loadtxt(k[1], delimiter=" ")
			CN = 1
		elif (k[0] == "S2CC"):
			print("Opening %s" % (k[1]))
			S2CC = np.loadtxt(k[1], delimiter=" ")
			CC = 1

parms = [S2NH, S2CH, S2CN, S2CC]

with open(in_file, "r") as inp, open(out_file, "w") as oup:
	for l in inp:
		k = l.split("\t")
		st = k[0]
		v = int(k[0]) - 1 
		c = 0
		for i in range(1, len(k), 3):
			#print("%d : %d \n" % (c, v))
			st = st + "\t%f\t%f\t%f" % (float(k[i]), parms[c][v, 1], parms[c][v, 2])
			c = c + 1
			print(i)
			print(st)
		oup.write(st+"\n")
		#if (len(k) == 13):
			# then S2NH, S2CH, S2CN, S2CC
		#elif (len(k) == 10):
			# then S2NH, S2CH, S2CN
		#print(len(k))
