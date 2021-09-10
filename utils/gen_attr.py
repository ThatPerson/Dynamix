import sys
import numpy as np
import models
# command line options are
## python3 gen_attr.py <final.dat> <model> <tag> <outputfile>
# eg python3 gen_attr.py demf/final.dat demf #1 file.attr

print(len(sys.argv))
print(sys.argv)
if (len(sys.argv) != 5):
	print("Wrong number of arguments")
	exit(-1)

fn = sys.argv[1]
model = sys.argv[2]
tag = sys.argv[3]
of = sys.argv[4]



x = np.loadtxt(fn, delimiter="\t")
print(x[:, 0])

def writer(of, k, v, n):
	with open(n+"_"+of, "w") as f:
		f.write("attribute: %s\n" % (n))
		f.write("match mode: any\nrecipient: residues\n")
		for i in range(0, len(k)):
			if (v[i] == -1):
				continue
			f.write("\t%s:%d\t%f\n" % (tag, k[i], v[i]))

if (model not in models.mds):
	print("Model doesn't exist")
	exit(-1)


writer(of, x[:, 0], x[:, 2], 'chisq')
for i in range(0, len(models.mds[model]['p'])):
	writer(of, x[:, 0], x[:, 3 + i], models.mds[model]['p'][i])

#attribute: saturationRatio
#match mode: 1-to-1
#recipient: residues
	#1.3:2	0.231700
	#1.3:4	0.135226
	#1.3:5	0.196665
