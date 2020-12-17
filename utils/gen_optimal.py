## prerequisites need combine_aicbic to have been run.
## arguments. 
# 1 - IC file (eg AIC/BIC/AICc)
# 2 - slow or fast
# 3 - output.bild file

import numpy as np
import sys
import pytraj as pt
import calc

RYD = 8.3145

try:
	ic_file = sys.argv[1]
	mode = sys.argv[2]
	output_file = sys.argv[3]
except IndexError:
	print("Insufficient arguments. Should be python gen_optimal.py <IC file> <mode> <output>")
	exit()

try:
	temperature = float(sys.argv[4])
except IndexError:
	temperature = 300
	print("Defaulting temperature to 300 K")

## load in atomic coordinates from PDB file.
traj = pt.load("2GI9.pdb")
min_v = 1
max_v = 56
N_coords = traj[:, ':%d-%d@N' % (min_v + 1, max_v)]
O_coords = traj[:, ':%d-%d@O' % (min_v, max_v)]
C_coords = traj[:, ':%d-%d@C' % (min_v, max_v)]
CA_coords = traj[:, ':%d-%d@CA' % (min_v, max_v)]

pp = np.zeros((58, 4, 3))

midpoints = np.zeros((58, 3))

for peptide_plane in range(2, 57):
	pp[peptide_plane, 0] = N_coords[0, peptide_plane-2]
	pp[peptide_plane, 1] = O_coords[0, peptide_plane-2]
	pp[peptide_plane, 2] = C_coords[0, peptide_plane-2]
	pp[peptide_plane, 3] = CA_coords[0, peptide_plane-2]
	midpoints[peptide_plane] = np.mean(pp[peptide_plane, :], axis=0)

pp[57, 1] = O_coords[0, 55]
pp[57, 2] = C_coords[0, 55]
pp[57, 3] = CA_coords[0, 55]

vectors = np.zeros((57, 3, 3))
basis_set = np.zeros((57, 3, 3)) 
centers = np.zeros((57, 3))


for peptide_plane in range(2, 57):
	Z = pp[peptide_plane + 1, 3] - pp[peptide_plane, 3] # CA-CA, this is Z
	D = pp[peptide_plane, 1] - pp[peptide_plane, 2] # C-O axis
	Y = np.cross(Z, D) # perpendicular to peptide plane  
	X = np.cross(Y, Z) # approximately aligned with CO axis
	# I'm not sure if these have the correct sign... (eg up or down).

	# In Lienin, the Gamma axis is shown aligned along Calpha_i, Calpha_i-1. I think this
	# is taken as the z axis, hence Calpha_i+1 - Calpha_i (my numbers are one above). 
	# Y is perpendicular to the plane, which is shown in rotation_tests.c. I think this is therefore the same as beta
	# (as I'm using the same theta and phi). Unfortunately l=2 spherical harmonics are invariant under inversion
	# so I can't think how to directly test this. X is then perpendicular and to maintain the
	#	   Z	 axis, X must be aligned along alpha.
	#	   |	 Which kind of makes sense that they would be aligned in the same sense,
	#	  / \	but I'm still skeptical.
	#	 X   Y

	x = X / np.linalg.norm(X)
	y = Y / np.linalg.norm(Y)
	z = Z / np.linalg.norm(Z)
	vectors[peptide_plane, :, 0] = x
	vectors[peptide_plane, :, 1] = y
	vectors[peptide_plane, :, 2] = z
	centers[peptide_plane, :] = (pp[peptide_plane, 3] + pp[peptide_plane + 1, 3]) / 2.


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
	with open("%s/final.dat" % (name), "r") as f:
		for l in f:
			k = l.split("\t")
			
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
			if (model[0] == "v"):
				current_res["orientation"] = kf[-3:]
			while (len(components) < resid+1):
				components.append({})
			components[resid] = current_res
			if (current_res["taus"] <= 0):
				continue
			
			if (np.log(current_res["taus"]) < l_tau):
				l_tau = np.log(current_res["taus"])
			elif (np.log(current_res["tauf"]) < l_tau):
				l_tau = np.log(current_res["tauf"])
			if (np.log(current_res["taus"]) > u_tau):
				u_tau = np.log(current_res["taus"])
			elif (np.log(current_res["tauf"]) > u_tau):
				u_tau = np.log(current_res["tauf"])
	return components, l_tau, u_tau
			


def generate_mf(parms, residue, mode, min_tau, max_tau):
	st = ""
	tau = 0
	mag = 0
	if (mode == "slow"):
		tau = np.log(parms["taus"])
		mag = parms["S2s"]
	elif (mode == "fast"):
		tau = np.log(parms["tauf"])
		mag = parms["S2f"]
	else:
		return ""
	if (tau == -1 or mag == -1):
		return ""
	x,y,z = midpoints[residue]
	#print(midpoints[residue])
	#print("(%f, %f, %f)" % (x, y, z))
	st = st + ".comment residue %d (%f, %f, %f)\n" % (residue, x, y, z)
	cspec = (tau - min_tau) / (max_tau - min_tau) # normalise colours.
	st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".sphere %f %f %f %f\n" % (x, y, z, 2* ((1 - mag)) / 0.2)
	return st

def generate_gaf(parms, residue, mode, min_tau, max_tau):
	st = ""
	tau = 0
	sigs = []
	if (mode == "slow"):
		tau = np.log(parms["taus"])
		sigs = parms["slow"]
	elif (mode == "fast"):
		tau = np.log(parms["tauf"])
		sigs = parms["fast"]
	else:
		return ""
	if ("orientation" in parms):
		# then vgaf
		orients = parms["orientation"]
		X = vectors[peptide_plane, :, 0]
		Y = vectors[peptide_plane, :, 1]
		Z = vectors[peptide_plane, :, 2]
		X, Y, Z = calc.apply_wigner(X, Y, Z, orients[0], orients[1], orients[2])
		vectors[peptide_plane, :, 0] = X
		vectors[peptide_plane, :, 1] = Y
		vectors[peptide_plane, :, 2] = Z
		
	if (tau == -1 or sigs[0] == -1):
		return ""
	
	S = [p * (180 / 3.14)/20. for p in sigs]
	
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
	st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".cylinder %f %f %f %f %f %f %f\n" % (alpha_start[0], alpha_start[1], alpha_start[2], \
		alpha_end[0], alpha_end[1], alpha_end[2], 0.4 * S[0])
	st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".cylinder %f %f %f %f %f %f %f\n" % (beta_start[0], beta_start[1], beta_start[2], \
		beta_end[0], beta_end[1], beta_end[2], 0.4 * S[1])
	st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
	st = st + ".cylinder %f %f %f %f %f %f %f\n" % (gamma_start[0], gamma_start[1], gamma_start[2], \
		gamma_end[0], gamma_end[1], gamma_end[2], 0.4 * S[0])
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
print(models["smf"][2])
print("Min: %f Max: %f" % (min_tau, max_tau))
generate_mf(models["smf"][3], 3, "slow", min_tau, max_tau)
generate_gaf(models["gaft"][3], 3, "slow", min_tau, max_tau) 

direct_mf = ["smf", "smft", "demf", "demft"]
direct_gaf = ["gaf", "gaft", "vgaf", "vgaft"]
slow_gaf = ["egaf", "egaft", "vegaf", "vegaft"]

with open(output_file, "w") as of:
	for l in choices:
		#print(l)
		if (l["model"] == ''):
			continue
		st = ""
		if (l["model"] in direct_mf):
			st = generate_mf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
		elif (l["model"] in direct_gaf):
			st = generate_gaf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
		elif (l["model"] in slow_gaf):
			if (mode == "slow"):
				st = generate_gaf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
			else:
				st = generate_mf(models[l["model"]][l["residue"]], l["residue"], mode, min_tau, max_tau)
		of.write(st)
		
## blocks should be sized by their order parameter or angular deflection, and coloured by their timescale.