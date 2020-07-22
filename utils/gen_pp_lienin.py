# generates the pp_ input files using data from Lienin1998.

import sys
import numpy as np
def deg2rad(f):
    return (np.pi / 180.) * f

fn = sys.argv[1]
theta = float(sys.argv[2])
phi = float(sys.argv[3])

with open(fn, "w") as f:
	for i in range(0, 56):
		f.write("%d %f %f\n" % (i+1, deg2rad(theta), deg2rad(phi)))
