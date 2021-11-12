import numpy as np
import math

h = 6.626 * math.pow(10, -34)
hbar = h / (2 * np.pi)

rNH = 1.02
rNC = 1.33
rNCA = 1.46
rCCAp = 1.525
rCCAc = 2.49
rCN = rNC
rCH = 2.04
rCHr = 1.82
rNHr = 1.8

gamma1H = 26.75 * math.pow(10,7.)
gamma13C = 6.73 * math.pow(10,7.)
gamma15N = 2.71 * math.pow(10,7.)
gammaE = 658 * gamma1H

def calc_D(r, g1, g2):
	global hbar
	r = r * math.pow(10, -10)
	d = math.pow(10, -7.) * g1 * g2 * hbar
	return d / math.pow(r, 3)
	
print("#define D_NH\t%f" % (calc_D(rNH, gamma15N, gamma1H)))
print("#define D_NC\t%f" % (calc_D(rNC, gamma15N, gamma13C)))
print("#define D_NCA\t%f" % (calc_D(rNCA, gamma15N, gamma13C)))
print("#define D_CCAp\t%f" % (calc_D(rCCAp, gamma13C, gamma13C)))
print("#define D_CCAc\t%f" % (calc_D(rCCAc, gamma13C, gamma13C)))
print("#define D_CN\t%f" % (calc_D(rCN, gamma15N, gamma13C)))
print("#define D_CH\t%f" % (calc_D(rCH, gamma13C, gamma1H)))
print("#define D_NHr\t%f" % (calc_D(rNHr, gamma15N, gamma1H)))
print("#define D_CHr\t%f" % (calc_D(rCHr, gamma13C, gamma1H)))
print("#define D_CE\t%f" % (calc_D(1, gammaE, gamma13C)))
print("#define D_NE\t%f" % (calc_D(1, gammaE, gamma15N)))

