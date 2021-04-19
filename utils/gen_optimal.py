## prerequisites need combine_aicbic to have been run.
## arguments. 
# 1 - IC file (eg AIC/BIC/AICc)
# 2 - slow or fast
# 3 - output.bild file

import numpy as np
import sys
import calc
import process

RYD = 8.3145

try:
	ic_file = sys.argv[1]
	mode = sys.argv[2]
	output_file = sys.argv[3]
	pdb_file = sys.argv[4]

except IndexError:
	print("Insufficient arguments. Should be python gen_optimal.py <IC file> <mode> <output>")
	exit()

v = np.loadtxt(ic_file, delimiter=",", usecols=[0], skiprows=1)
k = np.shape(v)
nres = k[0]

try:
	def_file = sys.argv[5]
	residue_inc = []
	with open(def_file, "r") as f:
		for l in f:
			residue_inc.append(float(l))
except IndexError:
	residue_inc = []


try:
	temperature = float(sys.argv[6])
except IndexError:
	temperature = 300
	print("Defaulting temperature to 300 K")

try:
	min_mag = float(sys.argv[7])
	max_mag = float(sys.argv[8])
except IndexError:
	min_mag = 0.7
	max_mag = 1.0

## load in atomic coordinates from PDB file.
pdb = process.Pdb(pdb_file, nres)
pdb.get_pp()
#traj = pt.load(pdb_file)
min_v = 1
max_v = nres

midpoints = np.zeros((nres+2, 3))
vectors = np.zeros((nres+1, 3, 3))
basis_set = np.zeros((nres+1, 3, 3))
centers = np.zeros((nres+1, 3))


for peptide_plane in range(2, nres+1):
	if ("global" in mode):
		x = [1,0,0]
		y = [0,1,0]
		z = [0,0,1]
	else:
		x = pdb.alpha(peptide_plane)
		y = pdb.beta(peptide_plane)
		z = pdb.gamma(peptide_plane)

	vectors[peptide_plane, :, 0] = x
	vectors[peptide_plane, :, 1] = y
	vectors[peptide_plane, :, 2] = z
	centers[peptide_plane, :] = pdb.center(peptide_plane)
	midpoints[peptide_plane, :] = pdb.center(peptide_plane)


def read_file(fn):
	choices = []
	models = []
	with open(fn, "r") as f:
		for l in f:
			k = l.split(",")
			try:
				resid = int(k[0])
			except ValueError:
				continue
			model = k[-1].strip()
			choices.append({"residue": resid, "model": model})
			if (model not in models and model != ""):
				models.append(model)
	return choices, models

def temp_tau(tau0, Ea, temp):
	return tau0 * np.exp(Ea / (RYD * temp))

def read_model(name, model):
	components = []
	l_tau = 100
	u_tau = -100
	u_mag = -1000
	l_mag = 1000
	with open("%s/final.dat" % (name), "r") as f:
		for l in f:
			k = l.split("\t")
			
			#print(k)
			try:
				resid = int(k[0])
				kf = [float(p) for p in k]
			except ValueError:
				continue
			# 0 - id
			# 1 - S2NH
			# 2 - chisq
			# 3 - first param
			current_res = {"residue": resid}
			if (model == "smf"):
				current_res["taus"] = kf[3]
				current_res["tauf"] = kf[3]
				current_res["S2s"] = kf[4]
				current_res["S2f"] = kf[4]
			elif (model == "smft"):
				current_res["taus"] = temp_tau(kf[3], kf[5], temperature)
				current_res["tauf"] = temp_tau(kf[3], kf[5], temperature)
				current_res["S2s"] = kf[4]
				current_res["S2f"] = kf[4]
			elif (model == "demf"):
				current_res["taus"] = kf[3]
				current_res["tauf"] = kf[5]
				current_res["S2s"] = kf[4]
				current_res["S2f"] = kf[6]
			elif (model == "demft"):
				current_res["taus"] = temp_tau(kf[3], kf[7], temperature)
				current_res["tauf"] = temp_tau(kf[5], kf[8], temperature)
				current_res["S2s"] = kf[4]
				current_res["S2f"] = kf[6]
			elif (model == "gaf" or model == "vgaf"):
				current_res["taus"] = kf[3]
				current_res["tauf"] = kf[4]
				current_res["slow"] = kf[5:8]
				current_res["fast"] = kf[8:11]
			elif (model == "gaft" or model == "vgaft"):
				current_res["taus"] = temp_tau(kf[3], kf[11], temperature)
				current_res["tauf"] = temp_tau(kf[4], kf[12], temperature)
				current_res["slow"] = kf[5:8]
				current_res["fast"] = kf[8:11]
			elif (model == "egaf" or model == "vegaf"):
				current_res["taus"] = kf[3]
				current_res["tauf"] = kf[4]
				current_res["slow"] = kf[5:8]
				current_res["S2f"] = kf[8]
			elif (model == "egaft" or model == "vegaft"):
				current_res["taus"] = temp_tau(kf[3], kf[9], temperature)
				current_res["tauf"] = temp_tau(kf[4], kf[10], temperature)
				current_res["slow"] = kf[5:8]
				current_res["S2f"] = kf[8]
			elif (model == "aimf" or model == "vaimf"):
				current_res["mags"] = kf[5:8]
				current_res["magf"] = kf[8:11]
				current_res["taus"] = kf[3]
				current_res["tauf"] = kf[4]
			elif (model == "aimft" or model == "vaimft"):
				current_res["mags"] = kf[5:8]
				current_res["magf"] = kf[8:11]
				current_res["taus"] = temp_tau(kf[3], kf[11], temperature)
				current_res["tauf"] = temp_tau(kf[4], kf[12], temperature)
				

			if (model[0] == "v"):
				current_res["orientation"] = kf[-3:]
			while (len(components) < resid+1):
				components.append({})
			components[resid] = current_res
		#	print(resid)
			if (current_res["taus"] <= 0):
				print("Oof")
				continue
			
			if (np.log(current_res["taus"]) < l_tau):
				l_tau = np.log(current_res["taus"])
			elif (np.log(current_res["tauf"]) < l_tau):
				l_tau = np.log(current_res["tauf"])
			if (np.log(current_res["taus"]) > u_tau):
				u_tau = np.log(current_res["taus"])
			elif (np.log(current_res["tauf"]) > u_tau):
				u_tau = np.log(current_res["tauf"])
			
	#print(components)
			
	return components, l_tau, u_tau
			
def generate_aimf(parms, residue, mode, min_tau, max_tau, min_mag, max_mag):
	st = ""
	tau = 0
	mags = []
	if ("slow" in mode):
		tau = np.log(parms["taus"])
		mags = parms["mags"]
	elif ("fast" in mode):
		tau = np.log(parms["tauf"])
		mags = parms["magf"]
	else:
		return ""
		
	mags = [(a - min_mag) / (max_mag - min_mag) for a in mags]
	for i in range(0, len(mags)):
		if (mags[i] < 0):
			mags[i] = 0.01
			
	x, y, z = midpoints[residue]
	#print(mags)
		
	beta = np.arcsin(-vectors[residue, 2, 0])
	gamma = np.arcsin(vectors[residue, 2, 1] / np.cos(beta))
	alpha = np.arccos(vectors[residue, 0, 0] / np.cos(beta))
	
	alpha = (180. / np.pi) * alpha
	beta = (180. / np.pi) * beta
	gamma = (180. / np.pi) * gamma
	
	#print("Rotation gamma %f about x" % (gamma))
	#print("Rotation beta %f about y" % (beta))
	#print("Rotation alpha %f about z" % (alpha))
	
	if ("orientation" in parms):
		print("Not yet implemented")
		return ""
		
	st = st + ".comment residue %d (%f, %f, %f)\n" % (residue, x, y, z)
	cspec = 1
#	cspec = (tau - min_tau) / (max_tau - min_tau) # normalise colours.
	st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".translate %f %f %f\n" % (x, y, z)
	st = st + ".scale %f %f %f\n" % (np.power(mags[1], 3.), np.power(mags[0], 3), np.power(mags[2], 3.)) # Sb along X, Sa along Y, Sc along Z
	st = st + ".rotate %f x\n" % (gamma)
	st = st + ".rotate %f y\n" % (beta)
	st = st + ".rotate %f z\n" % (alpha)
	st = st + ".sphere 0 0 0 1\n" # unit sphere
	st = st + ".pop\n.pop\n.pop\n.pop\n.pop\n"

	
	

	#
	return st
	
	
	''' *   Let Sa, Sb, Sc be the squared order parameters along the following directions;
 *     Sa - CaCa axis (eg the axis with gamma deflections about it)
 *     Sb - ~N-H axis (eg the axis with alpha deflections about it)
 *     Sc - perpendicular to plane axis (eg beta).'''


def generate_mf(parms, residue, mode, min_tau, max_tau):
	st = ""
	tau = 0
	mag = 0
	if ("slow" in mode):
		try:
			tau = np.log(parms["taus"])
			mag = parms["S2s"]
		except RuntimeWarning:
			tau = np.log(parms["tau"])
			mag = parms["S2"]
	elif ("fast" in mode):
		try:
			tau = np.log(parms["tauf"])
			mag = parms["S2f"]
		except RuntimeWarning:
			tau = np.log(parms["tau"])
			mag = parms["S2"]
	else:
		return ""
	if (tau == -1 or mag == -1):
		return ""
	x,y,z = midpoints[residue]
	#print(midpoints[residue])
	#print("(%f, %f, %f)" % (x, y, z))
	st = st + ".comment residue %d (%f, %f, %f)\n" % (residue, x, y, z)
	#cspec = (tau - min_tau) / (max_tau - min_tau) # normalise colours.
	cspec = 1
	st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".sphere %f %f %f %f\n" % (x, y, z, 2* ((1 - mag)) / 0.2)
	return st

def format_v(v):
	st = ""
	for i in v:
		st += "%f "%(i)
	return st

def generate_pp(residue):
	st = ".transparency 60\n"
	st += ".color 0.6 0.6 0.6\n"
	
	corner1 = midpoints[residue] - 1.*vectors[residue, :, 0] - 1.*vectors[residue, :, 2] 
	corner2 = midpoints[residue] + 1.*vectors[residue, :, 0] - 1.*vectors[residue, :, 2]
	corner3 = midpoints[residue] + 1.*vectors[residue, :, 0] + 1.*vectors[residue, :, 2]
	corner4 = midpoints[residue] - 1.*vectors[residue, :, 0] + 1.*vectors[residue, :, 2]
	
	st += ".polygon "+format_v(corner1) + format_v(corner2) + format_v(corner3) + format_v(corner4)+"\n"
	
	st += ".transparency 0\n"
	return st

def gen_rot_pp(residue, rot, axis):
	#def vrot(A, B, C, theta, r):
	Gamma = vectors[residue, :, 2]
	Beta = vectors[residue, :, 1]
	Alpha = vectors[residue, :, 0]
	
	if (axis == 0):
		rax = Alpha
	elif (axis == 1):
		rax = Beta
	elif (axis == 2):
		rax = Gamma
		
	Alphan, Betan, Gamman = calc.vrot(Alpha, Beta, Gamma, rot, rax)
	Alphap, Betap, Gammap = calc.vrot(Alpha, Beta, Gamma, -rot, rax)
	
	magn = 3*rot

	corner1 = midpoints[residue] - magn*Alphap - magn*Gammap
	corner2 = midpoints[residue] + magn*Alphap - magn*Gammap
	corner3 = midpoints[residue] + magn*Alphap + magn*Gammap
	corner4 = midpoints[residue] - magn*Alphap + magn*Gammap
	st = ".polygon "+format_v(corner1) + format_v(corner2) + format_v(corner3) + format_v(corner4)+"\n"
	
	corner1n = midpoints[residue] - magn*Alphan - magn*Gamman
	corner2n = midpoints[residue] + magn*Alphan - magn*Gamman
	corner3n = midpoints[residue] + magn*Alphan + magn*Gamman
	corner4n = midpoints[residue] - magn*Alphan + magn*Gamman
	st += ".polygon "+format_v(corner1n) + format_v(corner2n) + format_v(corner3n) + format_v(corner4n)+"\n"
	
	midax = midpoints[residue] + magn * rax
	midaxp = midpoints[residue] - magn * rax
	
	# faces
	if (axis == 0):
		st += ".polygon "+format_v(corner1n) + format_v(corner2n) + format_v(corner2) + format_v(corner1)+"\n"
		st += ".polygon "+format_v(corner3n) + format_v(corner4n) + format_v(corner4) + format_v(corner3)+"\n"
		
		st += ".polygon "+format_v(midax) + format_v(corner3) + format_v(corner3n) + "\n"
		st += ".polygon "+format_v(midaxp) + format_v(corner4) + format_v(corner4n) + "\n"
		st += ".polygon "+format_v(midax) + format_v(corner2) + format_v(corner2n) + "\n"
		st += ".polygon "+format_v(midaxp) + format_v(corner1) + format_v(corner1n) + "\n"
		
	elif (axis == 2):
		st += ".polygon "+format_v(corner1n) + format_v(corner4n) + format_v(corner4) + format_v(corner1)+"\n"
		st += ".polygon "+format_v(corner2n) + format_v(corner3n) + format_v(corner3) + format_v(corner2)+"\n"
		st += ".polygon "+format_v(midax) + format_v(corner4) + format_v(corner4n) + "\n"
		st += ".polygon "+format_v(midaxp) + format_v(corner2) + format_v(corner2n) + "\n"
		st += ".polygon "+format_v(midax) + format_v(corner3) + format_v(corner3n) + "\n"
		st += ".polygon "+format_v(midaxp) + format_v(corner1) + format_v(corner1n) + "\n"
	
	
	
	
	return st

def generate_gaf_new(parms, residue, mode, min_tau, max_tau):
	if (parms == {}):
		return ""
	st = ""
	tau = 0
	sigs = []
	if ("slow" in mode):
		tau = np.log(parms["taus"])
		sigs = parms["slow"]
	elif ("fast" in mode):
		tau = np.log(parms["tauf"])
		sigs = parms["fast"]
	else:
		return ""

	if ("orientation" in parms):
		print("New gen gaf doesn't support v* models")
		return ""
		
	if (tau == -1 or sigs[0] == -1):
		return ""
	
	#S = [p * 3 for p in sigs]
	S = [p * 1.96 for p in sigs]
	#st = ".transparency 50\n"
	st += ".color 1 0 0\n"
	st += gen_rot_pp(residue, S[0], 0)
	st += ".color 1 0 1\n"
	st += gen_rot_pp(residue, S[1], 1)
	st += ".color 0 0 1\n"
	st += gen_rot_pp(residue, S[2], 2)

	st += ".transparency 0\n"
	return st
	
#	".color 1 0 0"
#	".color 1 0 1"
#	".color 0 0 1"


def generate_gaf(parms, residue, mode, min_tau, max_tau):
	return generate_gaf_new(parms, residue, mode, min_tau, max_tau)
	if (parms == {}):
		return ""
	st = ""
	tau = 0
	sigs = []
	if ("slow" in mode):
		tau = np.log(parms["taus"])
		sigs = parms["slow"]
	elif ("fast" in mode):
		tau = np.log(parms["tauf"])
		sigs = parms["fast"]
	else:
		return ""
	if ("orientation" in parms):
		# then vgaf
		print("Performing rotation")
		orients = parms["orientation"]
		X = vectors[residue, :, 0]
		Y = vectors[residue, :, 1]
		Z = vectors[residue, :, 2]

		X, Y, Z = calc.undo_wigner(X, Y, Z, orients[0], orients[1], orients[2])
		vectors[residue, :, 0] = X
		vectors[residue, :, 1] = Y
		vectors[residue, :, 2] = Z
		
	if (tau == -1 or sigs[0] == -1):
		return ""
	
	S = [p * (180 / 3.14)/20. for p in sigs]
	
	if (S[0] > 1 or S[1] > 1 or S[2] > 1):
		return "" # eg deflection is so big it's no longer gaussian.
	
	alpha_start = midpoints[residue] - 3. * S[0] * vectors[residue, :, 0]
	beta_start = midpoints[residue] - 3. * S[1] * vectors[residue, :,  1]
	gamma_start = midpoints[residue] - 3. * S[2] * vectors[residue, :, 2]
	alpha_end = midpoints[residue] + 3. * S[0] * vectors[residue, :, 0]
	beta_end = midpoints[residue] + 3. * S[1] * vectors[residue, :, 1]
	gamma_end = midpoints[residue] + 3. * S[2] * vectors[residue, :, 2]
	

	x,y,z = midpoints[residue]
	#print(midpoints[residue])
	#print("(%f, %f, %f)" % (x, y, z))
	st = st + ".comment residue %d (%f, %f, %f)\n" % (residue, x, y, z)
	cspec = (tau - min_tau) / (max_tau - min_tau) # normalise colours.
	cspec = 1
	#st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".color %f 0 0\n" % (cspec)
	st = st + ".arrow %f %f %f %f %f %f %f %f\n" % (midpoints[residue, 0], midpoints[residue, 1], midpoints[residue, 2], alpha_start[0], alpha_start[1], alpha_start[2], 0.4 * S[0], 0.4 * S[0])
	st = st + ".arrow %f %f %f %f %f %f %f %f\n" % (midpoints[residue, 0], midpoints[residue, 1], midpoints[residue, 2], alpha_end[0], alpha_end[1], alpha_end[2], 0.4 * S[0], 0.4 * S[0])
	#st = st + ".cylinder %f %f %f %f %f %f %f\n" % (alpha_start[0], alpha_start[1], alpha_start[2], \
	#	alpha_end[0], alpha_end[1], alpha_end[2], 0.4 * S[0])
	#st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".color %f 0 %f\n" % (cspec, cspec)
	st = st + ".arrow %f %f %f %f %f %f %f %f\n" % (midpoints[residue, 0], midpoints[residue, 1], midpoints[residue, 2], beta_start[0], beta_start[1], beta_start[2], 0.4 * S[1], 0.4 * S[1])
	st = st + ".arrow %f %f %f %f %f %f %f %f\n" % (midpoints[residue, 0], midpoints[residue, 1], midpoints[residue, 2], beta_end[0], beta_end[1], beta_end[2], 0.4 * S[1], 0.4 * S[1])
	#st = st + ".cylinder %f %f %f %f %f %f %f\n" % (beta_start[0], beta_start[1], beta_start[2], \
	#	beta_end[0], beta_end[1], beta_end[2], 0.4 * S[1])
	#st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".color 0 0 %f\n" % (cspec)
	st = st + ".arrow %f %f %f %f %f %f %f %f\n" % (midpoints[residue, 0], midpoints[residue, 1], midpoints[residue, 2], gamma_start[0], gamma_start[1], gamma_start[2], 0.4 * S[2], 0.4 * S[2])
	st = st + ".arrow %f %f %f %f %f %f %f %f\n" % (midpoints[residue, 0], midpoints[residue, 1], midpoints[residue, 2], gamma_end[0], gamma_end[1], gamma_end[2], 0.4 * S[2], 0.4 * S[2])
	#st = st + ".cylinder %f %f %f %f %f %f %f\n" % (gamma_start[0], gamma_start[1], gamma_start[2], \
	#	gamma_end[0], gamma_end[1], gamma_end[2], 0.4 * S[2])
	return st

choices, mods = read_file(ic_file)
models = {}
max_tau = -1000
min_tau = 1000

for i in mods:
	models[i], l_tau, u_tau = read_model(i, i)
	print("%f, %f\n" % (l_tau, u_tau))
	if (l_tau < min_tau):
		min_tau = l_tau
	if (u_tau > max_tau):
		max_tau = u_tau

#print(models["smf"][2])
#print("Min: %f Max: %f" % (min_tau, max_tau))
#generate_mf(models["smf"][3], 3, "slow", min_tau, max_tau)
#generate_gaf(models["gaft"][3], 3, "slow", min_tau, max_tau) 

direct_mf = ["smf", "smft", "demf", "demft"]
direct_gaf = ["gaf", "gaft", "vgaf", "vgaft"]
slow_gaf = ["egaf", "egaft", "vegaf", "vegaft"]
aimf = ["aimf", "aimft"]

with open(output_file, "w") as of:
	for l in choices:
		if (l["residue"] not in residue_inc and len(residue_inc) != 0):
			continue
		if ("pp" in sys.argv):
			st = generate_pp(l["residue"])
			of.write(st)
		#print(l)
		if (l["model"] == ''):
			continue
		st = ""
		if (l["model"] in direct_mf):
			st = generate_mf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
		elif (l["model"] in direct_gaf):
			print("Hello")
			print(models[l["model"]][l["residue"]])
			st = generate_gaf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
		elif (l["model"] in slow_gaf):
			if ("slow" in mode):
				st = generate_gaf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
			else:
				st = generate_mf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
		elif (l["model"] in aimf): 
			st = generate_aimf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau, min_mag, max_mag)
		
		of.write(st)
		
		


		
## blocks should be sized by their order parameter or angular deflection, and coloured by their timescale.
