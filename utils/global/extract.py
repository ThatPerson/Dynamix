import sys
import numpy as np

try:
	pdb_fn = sys.argv[1]
	prefix = sys.argv[2]
	mode = sys.argv[3]
except IndexError:
	print("extract.py pdb_fn prefix")
	exit()

## I'm too lazy to use a library	
def read_pdb(fn):
	protein = []
	with open(fn, "r") as f:
		for l in f:
			if (l[:len("ATOM")] != "ATOM"):
				continue
			k = l.split()
			#print(k)
			residue = int(k[4])
			atom = k[2]
			x = float(k[5])
			y = float(k[6])
			z = float(k[7])
			while (len(protein) <= residue):
				protein.append({})
			protein[residue][atom] = {"x": x, "y": y, "z": z}
	return protein

def calc_inertial(prot):
	It = np.zeros((3, 3))
	
	mass = {"H":1, "N":15, "O":16, "C": 12}

	for k in prot:
		for at in k:
			try:
				m = mass[at[0]]
			except KeyError:
				continue
			It[0,0] = It[0,0] + m * np.power(k[at]["y"], 2.) * np.power(k[at]["z"], 2.)
			It[1,1] = It[1,1] + m * np.power(k[at]["x"], 2.) * np.power(k[at]["z"], 2.)
			It[2,2] = It[2,2] + m * np.power(k[at]["x"], 2.) * np.power(k[at]["y"], 2.)
			It[0,1] = It[0,1] - m * k[at]["x"] * k[at]["y"]
			It[0,2] = It[0,2] - m * k[at]["x"] * k[at]["z"]
			It[1,2] = It[1,2] - m * k[at]["y"] * k[at]["z"]#
	It[1,0] = It[0,1]
	It[2,1] = It[1,2]
	It[2,0] = It[0,2]
	'''[[-5216.86827896 -1962.542976   -6601.555232  ]
 [-1962.542976    -203.11043896 -2276.712192  ]
 [-6601.555232   -2276.712192   -7184.62832696]]'''
	#print(It)
	w, v = np.linalg.eig(It)
	#print(w)
	#print(v)
	
	w, v = zip(*sorted( zip(w, v)))
	#print(w)
	#print(v)
	return v

def get_vecs(prot, fr, to):
	k = []
	
	for r in range(1, len(prot)):
		m, f = fr
		n, t = to
		while (r >= len(k)):
			k.append(np.zeros((3)))
		try:
			fro = prot[r+m][f]
			too = prot[r+n][t]
		except KeyError:
			k[r] = np.zeros(3)
			continue
		
		l = [too["x"] - fro["x"], too["y"] - fro["y"], too["z"] - fro["z"]]
		k[r] = np.array(l) / np.linalg.norm(l)

	return k
	
def cross_all(a, b):
	c = []
	for i in range(0, len(a)):
		c.append(np.cross(a[i], b[i]))
	return c

def draw_vecs(a, p, fn):
	with open(fn, "w") as f:
		for i in range(0, len(a)):
			print(p[i])
			try:
				f.write(".arrow %f %f %f" % (p[i]["N"]["x"], p[i]["N"]["y"], p[i]["N"]["z"]))
				f.write(" %f %f %f\n" % (p[i]["N"]["x"] + a[i][0], p[i]["N"]["y"] + a[i][1], p[i]["N"]["z"] + a[i][2]))
			except KeyError:
				continue

def calc_orient(V, basis_set, fn, theta_div, phi_div):
	## basis set [X, Y, Z]
	with open(fn, "w") as f:
		for i in range(1, len(V)):
			bs = np.zeros((3, 3))
			bs[0, :] = basis_set[0][i]
			bs[1, :] = basis_set[1][i]
			bs[2, :] = basis_set[2][i]
			q = np.matmul(bs, V[i])
			q = q / np.linalg.norm(q)
			#print(q)
			# x = sin(theta) cos(phi)
			# y = sin(theta) sin(phi)
			# z = cos(theta)
			
			theta = np.arccos(q[2])
			#print("%d, %f, %f, %f, %f\n" % (i, q[0], np.sin(theta), q[0] / np.sin(theta), np.arccos(q[0] / np.sin(theta))))
			phid = q[0] / np.sin(theta)
			if (phid > 1):
				phid = 1
			if (phid < -1):
				phid = -1
			phi = np.arccos(phid)
			phi2= np.arcsin(q[1] / np.sin(theta))
			
			# phi2 and phi1 should be 90 degrees out. Relationship between them gives sign.
			if (np.isnan(theta)):
				f.write("%d -1 -1\n" % (i))
			else:
				f.write("%d %f %f\n" %( i, theta + theta_div, phi + phi_div))

def gen_bs(x, y, z, n):
	X = []
	Y = []
	Z = []
	for i in range(0, n):
		X.append(x)
		Y.append(y)
		Z.append(z)
	return X, Y, Z

k = read_pdb(pdb_fn)

pairs = [((-1, "C"), (0, "N"), "CN", 0, 0), 
((0, "N"), (-1, "C"), "NC", 0, 0),
((0, "N"), (-1, "C"), "CCSAxx", -136.1, 0),
((0, "N"), (-1, "C"), "CCSAyy", -46.1, 0),
((0, "N"), (-1, "C"), "CCSAzz", -48.4, +90),
((0, "N"), (0, "H"), "NH", 0, 0),
((0, "N"), (0, "H"), "NCSAxx", -68.0, 0),
((0, "N"), (0, "H"), "NCSAyy", -11.3, -90),
((0, "N"), (0, "H"), "NCSAzz", +22.0, 0),
((0, "N"), (0, "CA"), "NCA", 0, 0),
((-1, "C"), (0, "H"), "CNH", 0, 0),
((-1, "C"), (0, "CA"), "CCAc", 0, 0),
((-1, "C"), (-1, "CA"), "CCAp", 0, 0)]





CN = (get_vecs(k, (-1, "C"), (0, "N")))
NH = get_vecs(k, (0, "N"), (0, "H"))

if (mode == "local"):
	Z = get_vecs(k, (-1, "CA"), (0, "CA"))
	A = get_vecs(k, (0, "H"), (0, "N"))
	Y = cross_all(A, Z)
	X = cross_all(Z, Y)
elif (mode == "global"):
	X, Y, Z = gen_bs((1, 0, 0), (0, 1, 0), (0, 0, 1), len(k))
	
draw_vecs(X, k, "X.bild")
draw_vecs(Y, k, "Y.bild")
draw_vecs(Z, k, "Z.bild")


for i in pairs:
	f, t, n, td, pd = i
	td = td * (np.pi / 180)
	pd = pd * (np.pi / 180)
	pl = get_vecs(k, f, t)
	calc_orient(pl, [X, Y, Z], "%s_%s.csv" % (prefix, n), td, pd)
#print(k)
#ATOM      3  C   MET A   1      -0.920   0.356   3.619  1.00  8.69           C  
