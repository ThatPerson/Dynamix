import numpy as np
import argparse as ap


parser = ap.ArgumentParser()
parser.add_argument('-expS2', type=str, default=None, help='File containing order parameter.')
parser.add_argument('-column', type=int, default=1, help='Column in file to read order parameter from (starts at 1).')
parser.add_argument('-columnerrors', type=int, default=-1, help='Column in file to read order parameter errors from (starts at 1).')
parser.add_argument('-calcS2', type=str, default=None, help='File containing backcalculated order parameters.')
parser.add_argument('-npars', type=int, default=1, help='Number of parameters fit.')
args = parser.parse_args()

exp = np.loadtxt(args.expS2, usecols=[args.column-1, args.columnerrors-1])
calc = np.loadtxt(args.calcS2)

exp[exp[:, 0] < 0, :] = [np.nan, np.nan]

chisq = np.nansum(np.power(exp[:, 0] - calc[:, 1], 2) / np.power(exp[:, 1], 2))

print("Chisq:  %0.2f" % (chisq))

nmeas = np.count_nonzero(~np.isnan(exp[:, 0]))
df = nmeas - args.npars - 1


rchisq = chisq / df
print("rChisq: %0.2f" % (rchisq))
AIC = chisq + 2 * args.npars
BIC = chisq + args.npars * np.log(df)
AICc = AIC + ((2 * args.npars * (args.npars + 1)) / df)

print("AIC:  %0.2f" % (AIC))
print("BIC:  %0.2f" % (BIC))
print("AICc: %0.2f" % (AICc))
