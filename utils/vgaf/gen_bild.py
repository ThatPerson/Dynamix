import pytraj as pt
import numpy as np
import sys
import calc

min_v = 1
max_v = 56

traj = pt.load("2QMT.pdb")
N_coords = traj[:, ':%d-%d@N' % (min_v + 1, max_v)]
O_coords = traj[:, ':%d-%d@O' % (min_v, max_v)]
C_coords = traj[:, ':%d-%d@C' % (min_v, max_v)]
CA_coords = traj[:, ':%d-%d@CA' % (min_v, max_v)]

pp = np.zeros((58, 4, 3))

for peptide_plane in range(2, 57):
	pp[peptide_plane, 0] = N_coords[0, peptide_plane-2]
	pp[peptide_plane, 1] = O_coords[0, peptide_plane-2]
	pp[peptide_plane, 2] = C_coords[0, peptide_plane-2]
	pp[peptide_plane, 3] = CA_coords[0, peptide_plane-2]

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
	#       Z     axis, X must be aligned along alpha.
	#       |     Which kind of makes sense that they would be aligned in the same sense,
	#      / \    but I'm still skeptical.
	#     X   Y

	x = X / np.linalg.norm(X)
	y = Y / np.linalg.norm(Y)
	z = Z / np.linalg.norm(Z)
	vectors[peptide_plane, :, 0] = x
	vectors[peptide_plane, :, 1] = y
	vectors[peptide_plane, :, 2] = z
	centers[peptide_plane, :] = (pp[peptide_plane, 3] + pp[peptide_plane + 1, 3]) / 2.
	print("\n\n\n")
	
variant = 1
ty = "slow"

# args are
# final.dat outputfile variant type
fn = sys.argv[1]
of = sys.argv[2]
variant = int(sys.argv[3])
ty = sys.argv[4]

with open(fn, "r") as inp:
	with open(of, "w") as out:
		for l in inp:
			k = l.split()
			residue = int(k[0])
			if (ty == "slow"):
				sA = float(k[5]) * (180 / 3.14) / 20.
				sB = float(k[6]) * (180 / 3.14) / 20.
				sG = float(k[7]) * (180 / 3.14) / 20.
			else: 
				sA = float(k[8]) * (180 / 3.14) / 20.
				sB = float(k[9]) * (180 / 3.14) / 20.
				sG = float(k[10]) * (180 / 3.14) / 20.
			if (float(k[5]) < 0):
				sA = 0
				sB = 0
				sG = 0
				continue
			if (variant == 1):
				alpha = float(k[len(k) - 3])
				beta = float(k[len(k) - 2])
				gamma = float(k[len(k) - 2])

				#vectors[peptide_plane, :, 0] = x
				#vectors[peptide_plane, :, 1] = y
				#vectors[peptide_plane, :, 2] = z

				X = vectors[peptide_plane, :, 0]
				Y = vectors[peptide_plane, :, 1]
				Z = vectors[peptide_plane, :, 2]

				X, Y, Z = apply_wigner(X, Y, Z, alpha, beta, gamma)

				vectors[peptide_plane, :, 0] = X
				vectors[peptide_plane, :, 1] = Y
				vectors[peptide_plane, :, 2] = Z

			#sA = 0.1
			#sB = 0.1
			#sG = 0.3
			alpha = [9.0]
			beta = [9.0]
			gamma = [9.0]
		
			alpha_start = centers[residue] - 3. * sA * vectors[residue, :, 0]
			beta_start = centers[residue] - 3. * sB * vectors[residue, :,  2]
			gamma_start = centers[residue] - 3. * sG * vectors[residue, :, 1]
			alpha_end = centers[residue] + 3. * sA * vectors[residue, :, 0]
			beta_end = centers[residue] + 3. * sB * vectors[residue, :, 2]
			gamma_end = centers[residue] + 3. * sG * vectors[residue, :, 1]
			
			alpha.extend(alpha_start)
			alpha.extend(alpha_end)
			beta.extend(beta_start)
			beta.extend(beta_end)
			gamma.extend(gamma_start)
			gamma.extend(gamma_end)
			
			alpha_c = [sA * 0.4, 0., 1.*sA, 1.*(1-sA), 0., 1.*sA, 1.*(1-sA)]
			beta_c = [sB * 0.4, 0., 1.*sB, 1.*(1-sB), 0., 1.*sB, 1.*(1-sB)]
			gamma_c = [sG * 0.4, 0., 1.*sG, 1.*(1-sG), 0., 1.*sG, 1.*(1-sG)]
			colors = [0.2, 1., 1., 1., 1., 1., 1.]
	
			alpha.extend(alpha_c)
			beta.extend(beta_c)
			gamma.extend(gamma_c)
			
			out.write(".comment residue %d (%f, %f, %f)\n" % (residue, sA*30., sB*30., sG*30.))
			out.write(".color %f %f %f\n" % (0, 1.*sA, 1.*(1-sA)))
			out.write(".cylinder %f %f %f %f %f %f %f\n" % (alpha_start[0], alpha_start[1], alpha_start[2],\
						alpha_end[0], alpha_end[1], alpha_end[2], 0.4 * sA))
			out.write(".color %f %f %f\n" % (0, 1.*sB, 1.*(1-sB)))
			out.write(".cylinder %f %f %f %f %f %f %f\n" % (beta_start[0], beta_start[1], beta_start[2],\
						beta_end[0], beta_end[1], beta_end[2], 0.4 * sB))
			out.write(".color %f %f %f\n" % (0, 1.*sG, 1.*(1-sG)))
			out.write(".cylinder %f %f %f %f %f %f %f\n" % (gamma_start[0], gamma_start[1], gamma_start[2],\
						gamma_end[0], gamma_end[1], gamma_end[2], 0.4 * sG))
			
			
