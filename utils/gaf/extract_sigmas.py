import sys
if (len(sys.argv) < 3):
	print("Need more args")
	
fn = sys.argv[1]
ty = sys.argv[2]
output = sys.argv[3]

with open(output, "w") as out:
	with open(fn, "r") as inp:
		for l in inp:
			k = l.split()
			

			if (ty == "slow"):
				sA = float(k[5]) * (180 / 3.14)
				sB = float(k[6]) * (180 / 3.14)
				sG = float(k[7]) * (180 / 3.14)
			else:
				sA = float(k[8]) * (180 / 3.14)
				sB = float(k[9]) * (180 / 3.14)
				sG = float(k[10]) * (180 / 3.14)
			if (float(k[5]) < 0):
				sA = 0
				sB = 0
				sG = 0
			out.write("%d %f %f %f\n" % (int(k[0]), sA, sB, sG))
			
			

