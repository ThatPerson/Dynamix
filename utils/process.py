import numpy as np
import os
import models
import calc
import sys
import glob
import matplotlib.pyplot as plt
import math

def xyz(f):
	return "%f %f %f" % (f[0], f[1], f[2])

## from https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    indices = ~np.isnan(values)
    values = values[indices]
    weights = weights[indices]
    indices = ~np.isnan(weights)
    values = values[indices]
    weights = weights[indices]
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))

class Pdb:
	def get_pp(self):
		types = ["H", "N", "C", "O", "CA", "CAp"]
		H = types.index("H")
		N = types.index("N")
		C = types.index("C")
		O = types.index("O")
		CA= types.index("CA")
		CAp = types.index("CAp")
		self.vecs = np.zeros((self.nres+1, 4, 3))
		# alpha, beta, gamma, center
		for i in range(1, self.nres + 1):
			Z = self.pp[i, CA, :] - self.pp[i, CAp, :] # gamma

			D = self.pp[i, O, :] - self.pp[i, C, :]
			Y = np.cross(Z, D) # beta
			X = np.cross(Y, Z) # alpha
			if (np.linalg.norm(X) == 0 or np.linalg.norm(Y) == 0 or np.linalg.norm(Z) == 0):
				self.vecs[i, 0, :] = [0, 0, 0]
				self.vecs[i, 1, :] = [0, 0, 0]
				self.vecs[i, 2, :] = [0, 0, 0]
			else:
				self.vecs[i, 0, :] = X / np.linalg.norm(X)
				self.vecs[i, 1, :] = Y / np.linalg.norm(Y)
				self.vecs[i, 2, :] = Z / np.linalg.norm(Z)
			self.vecs[i, 3, :] = (self.pp[i, C, :] + self.pp[i, N, :])/2.

			#print("(%f, %f, %f) (%f, %f, %f), (%f, %f, %f)" % (self.vecs[i, 0, 0], self.vecs[i, 0, 1], self.vecs[i, 0, 2], self.vecs[i, 1, 0], self.vecs[i, 1, 1], self.vecs[i, 1, 2], self.vecs[i, 2, 0], self.vecs[i, 2, 1], self.vecs[i, 2, 2]))
			
	def get_vecs(self, fr, to):
		types = ["H", "N", "C", "O", "CA", "CAp"]
		if (fr not in types or to not in types):
			print("Unknown atoms %s, %s" % (fr, to))
			return -1
		fro = types.index(fr)
		too = types.index(to)
		k = np.zeros((self.nres + 1, 3))
		for i in range(1, self.nres + 1):
			froo = self.pp[i, fro, :]
			tooo = self.pp[i, too, :]
			l = tooo - froo
			k[i, :] = l / np.linalg.norm(l)
			
		return k
	
	def cross_all(self, A, B):
		if (np.shape(A) != np.shape(B)):
			print("Inconsistent shapes")
			return -1
		C = np.zeros(np.shape(A))
		k = np.shape(A)
		for i in range(0, k[0]):
			C[i, :] = np.cross(A[i, :], B[i, :])
		return C
		
	def calc_orient(self, V, basis, theta_div, phi_div):
		k = np.shape(V)
		orients = np.zeros((k[0], 2))
		for i in range(0, k[0]):
			bs = np.zeros((3, 3))
			bs[0, :] = basis[0][i, :]
			bs[1, :] = basis[1][i, :]
			bs[2, :] = basis[2][i, :]
			q = np.matmul(bs, V[i, :])
			q = q / np.linalg.norm(q)
			theta = np.arccos(q[2])
			phid = q[0] / np.sin(theta)
			phid = 1 if phid > 1 else phid
			phid = -1 if phid < -1 else phid
			phi = np.arccos(phid)
			phi2 = np.arcsin(q[1] / np.sin(theta))
			if (np.isnan(theta)):
				orients[i, :] = [-1, -1]
			else:
				orients[i, :] = [theta + theta_div, phi + phi_div]
		return orients
	
	
	def output_orient(self, fn, orie):
		print("Output %s" % (fn))
		with open(fn, "w") as f:
			for i in range(1, self.nres + 1):
				f.write("%d %f %f\n" % (i, orie[i, 0], orie[i, 1]))
	
	
	def gen_bs(self, M):
		X = np.zeros((self.nres+1, 3))
		Y = np.zeros((self.nres+1, 3))
		Z = np.zeros((self.nres+1, 3))
		for i in range(0, self.nres+1):
			X[i, :] = M[0]
			Y[i, :] = M[1]
			Z[i, :] = M[2]
		return X, Y, Z
	
	def draw_vecs(self, V, fn, atom):
		types = ["H", "N", "C", "O", "CA", "CAp"]
		if (atom not in types):
			return -1
		at = types.index(atom)
		with open(fn, "w") as f:
			for i in range(1, self.nres+1):
				try:
					atp = self.pp[i, at, :]
					vec = V[i, :]
				except KeyError:
					continue
				q = atp+vec
				if (xyz(atp+vec) == "nan nan nan"):
					continue
				f.write(".arrow %s %s\n" % (xyz(atp), xyz(atp + vec)))
				
				
	
	def gen_orients(self, prefix, M):
		pairs = [("C", "N", "CN", 0, 0),
				 ("N", "C", "NC", 0, 0),
				 ("N", "C", "CCSAxx", -136.1, 0),
				 ("N", "C", "CCSAyy", -46.1, 0),
				 ("N", "C", "CCSAzz", -48.4, +90),
				 ("N", "H", "NH", 0, 0),
				 ("N", "H", "NCSAxx", -60., 0),
				 ("N", "H", "NCSAyy", -11.3, -90.),
				 ("N", "H", "NCSAzz", +22.0, 0),
				 ("N", "CA", "NCA", 0, 0),
				 ("C", "H", "CNH", 0, 0),
				 ("C", "CA", "CCAc", 0, 0),
				 ("C", "CAp", "CCAp", 0, 0)]
		for i in pairs:
			print(i)
			f, t, n, td, pd = i
			td = td * (np.pi / 180.)
			pd = pd * (np.pi / 180.)
			pl = self.get_vecs(f, t)

			
			k = self.calc_orient(pl, M, td, pd)
			try:
				self.output_orient("%s_%s.csv" % (prefix, n), k)
			except ValueError:
				continue
	
	def gen_global(self, prefix):
		X, Y, Z = self.gen_bs([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
		self.gen_orients(prefix, [X, Y, Z])

	def gen_local(self, prefix):
		Z = self.get_vecs("CAp", "CA")
		A = self.get_vecs("H", "N")
		Y = self.cross_all(A, Z)
		X = self.cross_all(Z, Y)
		self.gen_orients(prefix, [X, Y, Z])
	
	

			
	def alpha(self, i):
		return self.vecs[i, 0, :]
	
	def beta(self, i):
		return self.vecs[i, 1, :]
	
	def gamma(self, i):
		return self.vecs[i, 2, :]
	
	def center(self, i):
		return self.vecs[i, 3, :]
		
	def read_pdb(self, fn):
		self.coords = np.zeros((self.nres+2, 5, 3))
		self.pp = np.zeros((self.nres+2, 6, 3))
		types = ["H", "N", "C", "O", "CA", "CAp"]
		# [residue, type, xyz]
		# types go [H, N, C, O]
		with open(fn, "r") as f:
			for l in f:
				if (l[:len("ATOM")] != "ATOM"):
					continue
				residue = int(l[23:27])
				x = float(l[31:39])
				y = float(l[39:47])
				z = float(l[47:55])

				atom = l[13:17].strip()
				if (atom in types):
					self.coords[residue, types.index(atom), :] = [x, y, z]
					
					if (atom == 'N' or atom == "H"):
						self.pp[residue, types.index(atom), :] = [x, y, z]
					elif (atom == "C" or atom == "O"):
						self.pp[residue + 1, types.index(atom), :] = [x, y, z]
					elif (atom == "CA"):
						self.pp[residue, types.index(atom), :] = [x, y, z]
						self.pp[residue+1, types.index("CAp"), :] = [x, y, z]
						
						
		return self.coords
		
	def __init__(self, fn="", nres=0):
		self.nres = nres
		if (fn != ""):
			self.read_pdb(fn)


class Results:
	def gen_mf(self, i, mode, pdb):
		self.get_npars()
		pars = models.mds[self.type]['p']
		con = {"slow": {"S2": "S2s", "tau": "taus"},"fast": {"S2": "S2f", "tau": "tauf"},"smf": {"S2": "S2", "tau": "tau"}}
		if (mode not in con):
			print("%s not recognized" % (mode))
			return ""
		S2p = con[mode]["S2"]
		taup = con[mode]["tau"]
		S2_ind = pars.index(S2p)
		tau_ind = pars.index(taup)
		this_res = self.final_dat[self.final_dat[:, 0] == i, :]
		#print(this_res)
		S2 = this_res[0, 3 + S2_ind]
		tau = this_res[0, 3 + tau_ind]
		#print("%f, %f" % (S2, tau))
		if (S2 < 0):
			return ""
		x, y, z = pdb.center(i)
		cspec = 1
		st = ".comment residue %d (%f, %f, %f)\n" % (i, x, y, z)
		st = st + ".color %f %f %f\n" % (0, 1. * cspec, 1.*(1-cspec))
		st = st + ".sphere %f %f %f %f\n" % (x, y, z, 2* ((1 - S2)) / 0.2)
		return st
		
	def gen_aimf(self, i, mode, pdb):
		return ".comment AIMF not implemented yet"
		
	def gen_gaf(self, i, mode, pdb):
		self.get_npars()
		alpha = pdb.alpha(i)
		beta = pdb.beta(i)
		gamma = pdb.gamma(i)
		center = pdb.center(i)
		pars = models.mds[self.type]['p']
		#print(self.final_dat)
		this_res = self.final_dat[self.final_dat[:, 0] == i, :]
		if (self.type[0] == 'v'):
			# variable orientation
			a_ind = pars.index("alph")
			b_ind = pars.index("beta")
			g_ind = pars.index("gamm")
			a = this_res[0, a_ind]
			b = this_res[0, b_ind]
			g = this_res[0, g_ind]
			alpha, beta, gamma = calc.undo_wigner(alpha, beta, gamma, a, b, g)
		
		con = {"slow": {"A": "sAs", "B": "sBs", "G": "sGs", "tau": "taus"}, "fast": {"A": "sAf", "B": "sBf", "G": "sGf", "tau": "tauf"}}
		if (mode not in con):
			print("%s not recognized" % (mode))
			return ""
		Ap = con[mode]["A"]
		Bp = con[mode]["B"]
		Gp = con[mode]["G"]
		taup = con[mode]["tau"]
		SA_ind = pars.index(Ap)
		SB_ind = pars.index(Bp)
		SG_ind = pars.index(Gp)
		tau_ind = pars.index(taup)

		#print(this_res)
		sigs = [0, 0, 0]
		sigs[0] = this_res[0, 3 + SA_ind]
		sigs[1] = this_res[0, 3 + SB_ind]
		sigs[2] = this_res[0, 3 + SG_ind]
		for q in sigs:
			if (q < 0):
				return ""
		tau = this_res[0, 3 + tau_ind]
		S = [p * (180 / 3.14) / 5. for p in sigs]
		for q in S:
			if (q > 5):
				return ""
		a_start = center - alpha * S[0]
		a_end = center + alpha * S[0]
		
		b_start = center - beta * S[1]
		b_end = center + beta * S[1]
		
		g_start = center - gamma * S[2]
		g_end = center + gamma * S[2]
		
		x, y, z = pdb.center(i)
		cspec = 1
		S = [0.1 * k for k in S]
		
		st = ".comment residue %d (%f, %f, %f)\n" % (i, x, y, z)
		st = st + ".color %f 0 0\n" % (cspec)
		st = st + ".arrow %s %s %f %f\n" % (xyz(center), xyz(a_start), S[0], S[0])
		st = st + ".arrow %s %s %f %f\n" % (xyz(center), xyz(a_end), S[0], S[0])
		st = st + ".color %f 0 %f\n" % (cspec, cspec)
		st = st + ".arrow %s %s %f %f\n" % (xyz(center), xyz(b_start), S[1], S[1])
		st = st + ".arrow %s %s %f %f\n" % (xyz(center), xyz(b_end), S[1], S[1])
		st = st + ".color 0 0 %f\n" % (cspec)
		st = st + ".arrow %s %s %f %f\n" % (xyz(center), xyz(g_start), S[2], S[2])
		st = st + ".arrow %s %s %f %f\n" % (xyz(center), xyz(g_end), S[2], S[2])
		
		return st
		
		
	
	def gen_self(self, k, pdb, fn, mode):
		self.get_npars()
		with open(fn, "w") as f:
			for i in k:
				if (mode == "fast"):
					if ("egaf" in self.type):
						f.write(self.gen_mf(i, mode, pdb))
					elif ("gaf" in self.type):
						f.write(self.gen_gaf(i, mode, pdb))
					elif ("aimf" in self.type):
						f.write(self.gen_aimf(i, mode, pdb))
					else:
						f.write(self.gen_mf(i, mode, pdb))
				elif (mode == "slow"):
					if ("gaf" in self.type):
						f.write(self.gen_gaf(i, mode, pdb))
					elif ("aimf" in self.type):
						f.write(self.gen_aimf(i, mode, pdb))
					else:
						f.write(self.gen_mf(i, mode, pdb))
				else:
					f.write(self.gen_mf(i, mode, pdb))	
			

	def combine_files(self, fn, f):
		n = 0
		with open(fn+"/%s.dat"%(f), "w") as fil:
			for i in glob.glob(fn+"/%s_proc*" % (f)):
				n = n + 1
				with open(i, "r") as q:
					for l in q:
						fil.write(l)
		return n

	def combine_output(self, fn):
		self.combine_files(fn, "final")
		self.combine_files(fn, "gaf")
		er = self.combine_files(fn, "errors")
		self.combine_files(fn, "orderparams")
		if (er == 0):
			self.errors = False
		else:
			self.errors = True
	
	def read_backcalc(self, fn, i):
		fi = "%s/backcalc_%d.dat" % (fn, i)
		if (not os.path.exists(fi)):
			return -1
		return np.loadtxt(fi)

	def read_results(self, fn = ""):
		if (fn == ""):
			fn = self.system["OUTPUT"]
		self.res = fn
		self.errors = True
		if (len(glob.glob(fn + "/errors*")) == 0):
			self.errors = False
		if (not os.path.exists(fn + "/final.dat")):
			self.combine_output(fn)
		
		self.final_dat = np.loadtxt(fn+"/final.dat")
		if (self.errors != False): self.errors_dat = np.loadtxt(fn+"/errors.dat")
		self.gaf_dat = np.loadtxt(fn+"/gaf.dat")
		self.op_dat = np.loadtxt(fn+"/orderparams.dat")
		
		self.backcalcs = []
		k = np.shape(self.final_dat)
		for i in range(1, k[0]+1):
			self.backcalcs.append(self.read_backcalc(fn, i))
		
		
	def get_npars(self):
		if ("MODEL" not in self.system):
			print("No model defined")
			return -1
		model = self.system["MODEL"].lower()
		if ("OR_VARY" in self.system):
			if (self.system["OR_VARY"] == 1):
				model = "v"+model
		self.type = model
		if (model not in models.mds):
			print("Model %s not found.\n" % (model))
			return -1
		return models.mds[model]["n"]
		
		
	def calc_IC(self):
		IC = np.zeros((len(self.backcalcs), 5))
		for i in range(0, len(self.backcalcs)):
			calc_R = self.backcalcs[i][:, 1]
			exp_R = self.backcalcs[i][:, 2]
			err_R = self.backcalcs[i][:, 3]
			
			exp_R[calc_R < 0] = 0
			err_R[calc_R < 0] = 1
			calc_R[calc_R < 0] = 0
			
			sigma = err_R / 2.
			chisq_R = np.sum(np.power((exp_R - calc_R), 2.) / np.power(err_R, 2.))
		
			max_op = 4
			ops = self.op_dat[self.op_dat[:, 0] == i+1, :]
			chisq_OP = 0
			ign = 0
			for inc in range(0, max_op):
				calc = ops[0, 1 + (inc*3)]
				exp = ops[0, 2 + (inc*3)]
				if (exp <= 0):
					ign = 1
				err = ops[0, 3 + (inc*3)]
				ctmp = np.power(exp - calc, 2.) / np.power(err, 2.)
				chisq_OP = chisq_OP + ctmp
			chisq = chisq_R + chisq_OP
			npars = self.get_npars()

			nvals = len(calc_R[calc_R != 0]) + max_op
			
			
			df = nvals - npars - 1

			if (df <= 0 or chisq == 0 or ign == 1):
				AIC = 1e9
				BIC = 1e9
				AICc = 1e9
			else:
				AIC = chisq + 2 * npars
				BIC = chisq + npars * np.log(nvals)
				AICc = AIC + ((2 * npars * (npars + 1)) / df)
			IC[i, 1:] = [AIC, BIC, AICc, chisq]
			IC[i, 0] = i + 1
		return IC
		
	def summarystats(self, units):
		self.get_npars()
		pars = models.mds[self.type]["p"]
		# data is arranged as self.errors_dat[RESIDUE, 3 + (i*2)]
		# and errors are      self.errors_dat[RESIDUE, 4 + (i*2)]
		if (self.errors == False):
			print("Only works with errors")
			return 1
			
		ret_data = {}
		for unit in units:
			# unit should be a list of residues to be included
			ret_data[unit] = {}
			print(unit)
			dat = self.errors_dat[np.isin(self.errors_dat[:, 0], units[unit]), :]
			k = np.shape(dat)
			for i in range(0, k[1]):
				print("%d: %s" % (i, ' '.join([str(pl) for pl in dat[:, i]])))
			print("\n\n\n")
			for i in range(0, len(pars)):
				vals = dat[:, 3 + (i * 2)]
				errs = dat[:, 4 + (i * 2)]
				print("%d (%d): %s" % (i, 3 + (2 * i), ' '.join([str(pl) for pl in vals])))
				if (len(vals) == 0 or len(errs) == 0):
					ret_data[unit][pars[i]] = {"val": -1, "std": 0}
					continue
				errs[errs == 0] = np.nan
				weight = 1/errs
				avg, std = weighted_avg_and_std(vals, weight)
				ret_data[unit][pars[i]] = {"val": avg, "std": std}
		return ret_data

		
	def plot(self, fn):
		self.get_npars()
		pars = models.mds[self.type]["p"]
		xdim = len(pars) + 1
		fig,axs = plt.subplots(xdim, 1, figsize=(7., xdim * 2), dpi=80)
		ic = self.calc_IC()
		axs[0].plot(ic[:, 0], ic[:, 4], 'k,')
		axs[0].set_ylim(bottom = 0)
		axs[0].set_title("Chisq")
		c = 1
		gaf = ["alph", "beta", "gamm"]
		for i in range(0, len(pars)):
			axs[c].set_title(pars[i])
			xdat = self.final_dat[:, 0]
			if (self.errors == True):
				ydat = self.errors_dat[:, 3 + (i * 2)]
				yerr = self.errors_dat[:, 4 + (i * 2)]
				axs[c].errorbar(xdat, ydat, yerr=yerr, fmt='k,')
			else:
				ydat = self.final_dat[:, 3 + i]
				yerr = np.zeros(np.shape(ydat))
				axs[c].plot(xdat, ydat, 'k+')
			axs[c].set(xlabel="residue", ylabel=pars[i])
			if ("tau" in pars[i]):
				axs[c].set_yscale("log")
			elif ("Ea" in pars[i]):
				axs[c].set_ylim(bottom=0, top=60000)
			elif ("S2" in pars[i]):
				axs[c].set_ylim(top = 1)
				kl = axs[c].get_ylim()
				if (kl[0] < 0.2):
					axs[c].set_ylim(bottom=0.5)
			elif (pars[i][0] == 's'):
				axs[c].set_ylim(bottom = 0)
			c = c + 1
		fig.tight_layout(h_pad = 1)
		plt.savefig(fn, dpi=300, facecolor='w', edgecolor='w',
	        orientation='portrait', papertype='a4', format='eps',
	        transparent=False, bbox_inches=None, pad_inches=0.5,
	        frameon=None, metadata=None)
	
	def plot_bc(self, fn):
		k = len(self.backcalcs)
		kl = np.shape(self.backcalcs[0])
		data = np.zeros((k, kl[0], kl[1]))
		for i in range(0, k):
			data[i, :, :] = self.backcalcs[i]
		plt.plot(data[:, :, 1], data[:, :, 2], 'k.')
		plt.xlim([0, 10])
		plt.ylim([0, 10])
		plt.savefig(fn)
		
	
	def plot_relax(self, fn):
		types = ['15N R1', '15N R1p', '13C R1', '13C R1p']
		k = len(self.backcalcs)
		kl = np.shape(self.backcalcs[0])
		data = np.zeros((k, kl[0], kl[1]))
		for i in range(0, k):
			data[i, :, :] = self.backcalcs[i]
		n_relax = kl[0]
		fig,axs = plt.subplots(int(np.floor(n_relax / 3.)) + 1, 3, figsize=(8, n_relax/1.2), dpi=80)
		curr_x = 0
		curr_y = 0
		x = np.arange(1, k+1)
		for i in range(0, n_relax):
			Rcalc = data[:, i, 1]
			Rexp = data[:, i, 2]
			Rerr = data[:, i, 3]
			axs[curr_x, curr_y].errorbar(x, Rexp, yerr=Rerr, fmt='k,')
			axs[curr_x, curr_y].plot(x, Rcalc, 'b,')
			axs[curr_x,curr_y].set_ylim(bottom=0)
			field = self.backcalcs[20][i, 4]
			wr = self.backcalcs[20][i, 5]/1000.
			w1 = self.backcalcs[20][i, 6]/1000.
			T = self.backcalcs[20][i, 7]
			typ = int(self.backcalcs[20][i, 8])
			if (typ in [0, 2]):
				krt = axs[curr_x, curr_y].get_ylim()
				if (krt[1] > 0.4):
					axs[curr_x, curr_y].set_ylim(top=0.4)
			else:
				krt = axs[curr_x, curr_y].get_ylim()
				if (krt[1] > 6):
					axs[curr_x, curr_y].set_ylim(top=20)
			axs[curr_x,curr_y].set_title("%s (%d MHz, %0.1f kHz, %0.1f kHz, %0.1f K)" % (types[typ], field, wr, w1,T), fontsize=5)
			axs[curr_x, curr_y].set(xlabel="residue", ylabel="rate / s^-1")
			curr_y = curr_y + 1
			if (curr_y == 3):
				curr_y = 0
				curr_x = curr_x + 1
		fig.tight_layout(h_pad=1)
		plt.savefig(fn, dpi=300, facecolor='w', edgecolor='w',
	        orientation='portrait', papertype='a4', format='eps',
	        transparent=False, bbox_inches=None, pad_inches=0.5,
	        frameon=None, metadata=None)
		print(x)
	
	def read_system(self, fn):
		self.fn = fn
		mode = 0
		self.system = {}
		self.relax_files = []

		with open(fn, "r") as f:
			for l in f:
				if (l[0] == '%'):
					continue
				if (l[:len("#RELAXATION")] == "#RELAXATION"):
					mode = 1
					continue
				if (mode == 0):
					k = [i.strip() for i in l.split()]
					if (len(k) != 3):

						continue

					try:
						self.system[k[0]] = float(k[2])
					except ValueError:
						self.system[k[0]] = k[2]
				elif (mode == 1):
					self.relax_files.append(l.strip())

				

	def __init__(self, fn=""):
		self.fn = fn
		self.errors = False
		if (fn != ""):
			self.read_system(self.fn)

def compare_IC(mods, selector, fn):
	m_ic = mods[0].calc_IC()
	m = {"aic": 1, "bic": 2, "aicc": 3, "chisq": 4}
	n = m[selector.lower()]
	k = np.shape(m_ic)
	pl = np.zeros((k[0], len(mods)+1))
	pl[:, 0] = m_ic[:, 0]
	c = 1
	
	for i in mods:
		tmp = i.calc_IC()
		pl[:, c] = tmp[:, n]
		c = c + 1
	victors = []
	with open(fn, "w") as outp:
		for i in range(0, k[0]):
			ma = np.argmin(pl[i, 1:])
			print(mods[ma].type)
			outp.write("%d\t" % (pl[i, 0]))
			for q in range(1, len(mods) + 1):
				outp.write("%f\t" % (pl[i, q]))
			outp.write("%s\n" % (mods[ma].type))
			victors.append(ma)
	
	
	print(k)
	return victors

class Simulation:
	def __init__(self):
		self.pars = {}
		self.relax = []
	
	def prompto(self, text, options):
		print(text)
		for i in range(0, len(options)):
			print("\t[%d] %s" % (i, options[i][0]))
		q = 0
		val = 0
		while (q == 0):
			k = input("> ")
			if (k == "exit"):
				exit(-1)
			try:
				val = int(k)
				ke = options[val]
				q = 1
			except IndexError:
				print("Option not recognized")
			except ValueError:
				print("Please enter integer")
		return val
		
	def promptf(self, text):
		print(text)
		q = 0
		val = 0
		while (q == 0):
			k = input("> ")
			if (k == "exit"):
				exit(-1)
			try:
				val = float(k)
				q = 1
			except ValueError:
				print("Please enter float")
		return val
	
	def prompti(self, text):
		print(text)
		q = 0
		val = 0
		while (q == 0):
			k = input("> ")
			if (k == "exit"):
				exit(-1)
			try:
				val = int(k)
				q = 1
			except ValueError:
				print("Please enter integer")
		return val
	
	
	def read_orderparams(self, typ, v):
		k = ""
		print("Please enter file containing S2 %s order parameters" % (typ))
		fn = input("> ")
		self.pars[v] = fn
		self.pars["W_%s" % (v)] = self.prompti("What weighting should I give these? [1 is relaxation]")
		
	def read_relaxation(self):
		return 1
		
	def save(self, folder):
		with open("%s/model.dx" % (folder), "w") as f:
		
			pdb = Pdb(self.pdb, self.pars["N_RESIDUES"])
			pdb.get_pp()
			if (self.pars["GLOBAL"] == 1):
				pdb.gen_global("%s/OP" % (folder))
			else:
				pdb.gen_local("%s/OP" % (folder))
			ors = ["CN", "NC", "CCSAxx", "CCSAyy", "CCSAzz", "NH", "NCSAxx", "NCSAyy", "NCSAzz", "NCA", "CNH", "CCAc", "CCAp"]
			for i in ors:
				self.pars["OR_%s"%(i)] = "%s/OP_%s.csv" % (folder, i)
				
			for i in self.pars:
				if ((i == "GLOBAL" or i == "CNCOMP") and self.pars[i] == 0):
					continue
				
				try:
					k = int(self.pars[i])
					if (float(self.pars[i]) == int(self.pars[i])):
						f.write("%s = %d\n" % (i, self.pars[i]))
					else:
						f.write("%s = %f\n" % (i, self.pars[i]))
				except ValueError:
					f.write("%s = %s\n" % (i, self.pars[i]))
	
			f.write("#RELAXATION")
	
	def prompt(self):
		model = [
					("Simple Model Free", "SMF"),
					("Extended Model Free", "EMF"), 
					("MF Fast, GAF Slow", "EGAF"),
					("Gaussian Axial Fluctuations", "GAF")
				]
		modes = [
					("Local", 0),
					("Global", 1)
				]
		yn = [
				("Yes", 1),
				("No", 0)
			]
		mod = self.prompto("Please select a model", model)
		self.pars["MODEL"] = model[mod][1]
		vt = self.prompto("Variable temperature?", yn)
		if (yn[vt][1] == 1):
			self.pars["MODEL"] = self.pars["MODEL"] + "T"
		if ("GAF" in self.pars["MODEL"]):
			vo = self.prompto("Variable orientation?", yn)
			if (yn[vo][1] == 1):
				self.pars["MODEL"] = "V" + self.pars["MODEL"]
		mod = self.prompto("Mode:", modes)
		self.pars["GLOBAL"] = mod
		self.pars["N_RESIDUES"] = self.prompti("Number of residues:")
		self.pars["N_ANNEAL_ITER"] = self.prompti("Number of annealing iterations")
		self.pars["N_NM_ITER"] = self.prompti("Number of Nelder-Mead iterations")
		self.pars["N_ERROR_ITER"] = self.prompti("Number of error iterations")
		if (self.prompto("Use default annealing parameters?", yn) == 0):
			self.pars["ANNEAL_TEMP"] = 6000
			self.pars["ANNEAL_WOBB"] = 0.05
			self.pars["ANNEAL_THERM"] = 1.1
			self.pars["ANNEAL_RESTART"] = 0.01
		else:
			self.pars["ANNEAL_TEMP"] = self.promptf("Anneal temperature")
			self.pars["ANNEAL_WOBB"] = self.promptf("Anneal wobb")
			self.pars["ANNEAL_THERM"] = self.promptf("Anneal thermostat")
			self.pars["ANNEAL_RESTART"] = self.promptf("Anneal restart probability")
		cncomp = self.prompto("Enable C/N Ratio compensation?", yn)
		if (yn[cncomp][0] == "Yes"):
			self.pars["CNCOMP"] = 1
		else:
			self.pars["CNCOMP"] = 0
		self.pars["NTHREADS"] = self.prompti("Number of threads?")

		print("Please enter PDB file")
		self.pdb = input("> ")
		self.read_orderparams("NH", "S2NH")
		self.read_orderparams("CH", "S2CH")
		self.read_orderparams("CC", "S2CC")
		self.read_orderparams("CN", "S2CN")
		self.read_relaxation()
		print("Please enter folder?")
		k = input("> ")
		if (k == "exit"):
			return -1
		self.save(k)
		

		
		
	def read_spreadsheet(self, fn, fr, to):
		return 1

if (__name__ == "__main__"):	
	k = Simulation()
	k.prompt()

	exit(-1)
	pdb = Pdb("../global/2GI9_H.pdb", 56)
	pdb.get_pp()
	pdb.gen_local("local")
	#pdb.gen_global("global")

	res = Results("../../2021system/gaft.dx")
	res.read_results("../../APMs/APM27/gaft")

	#
	#print(res.gen_mf(20, "slow", pdb))

	#def gen_self(self, k, pdb, fn, mode):

	k = np.arange(1, 56)
	res.gen_self(k, pdb, "test.bild", "slow")

	#res.plot_relax("file.eps")
	#res.plot_bc("bc.eps")
	exit(-1)
	#print(res.get_npars())
	print(res.calc_IC())
	res2 = Results("2021system/gaft.dx")
	res2.read_results()
	#print(res.get_npars())
	compare_IC([res, res2], "AIC", "test.csv")
