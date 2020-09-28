# I'm not entirely sure I have the sense of the rotations right - Lienin 1998 is a bit unclear on this point.
# I think that the axis main line is taken from Ci-1 to Ci and all rotations are anticlockwise.

import pytraj as pt
import numpy as np
import sys

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
def rotate_x(V, deg):
	rot_mat = np.zeros((3, 3))
	rot_mat[0, 0] = 1 # y, x
	rot_mat[1, 1] = np.cos(deg)
	rot_mat[2, 1] = -np.sin(deg)
	rot_mat[1, 2] = np.sin(deg)
	rot_mat[2, 2] = np.cos(deg)
	return np.matmul(rot_mat, V)

def rotate_y(V, deg):
	rot_mat = np.zeros((3, 3))
	rot_mat[1, 1] = 1 # y, x
	rot_mat[0, 0] = np.cos(deg)
	rot_mat[2, 0] = np.sin(deg)
	rot_mat[0, 2] = -np.sin(deg)
	rot_mat[2, 2] = np.cos(deg)
	return np.matmul(rot_mat, V)

def rotate_z(V, deg):
	rot_mat = np.zeros((3, 3))
	rot_mat[2, 2] = 1 # y, x
	rot_mat[0, 0] = np.cos(deg)
	rot_mat[1, 0] = -np.sin(deg)
	rot_mat[0, 1] = np.sin(deg)
	rot_mat[1, 1] = np.cos(deg)
	return np.matmul(rot_mat, V)



for peptide_plane in range(2, 57):
	A = pp[peptide_plane + 1, 3] - pp[peptide_plane, 3]
	D = pp[peptide_plane, 1] - pp[peptide_plane, 2]
	C = np.cross(A, D) # perpendicular to peptide plane
	B = np.cross(A, C) # approximately aligned with CO axis
	alpha = B / np.linalg.norm(B)
	beta = C / np.linalg.norm(C)
	gamma = A / np.linalg.norm(A)
	vectors[peptide_plane, :, 0] = alpha
	vectors[peptide_plane, :, 1] = gamma
	vectors[peptide_plane, :, 2] = beta
	
	centers[peptide_plane, :] = (pp[peptide_plane, 3] + pp[peptide_plane + 1, 3]) / 2.
	
	basis_set[peptide_plane, 0, :] = alpha
	basis_set[peptide_plane, 1, :] = gamma
	basis_set[peptide_plane, 2, :] = beta
	
	print(pp[peptide_plane + 1, 3])
	print(pp[peptide_plane, 3])
	print(A)
	print(alpha)
	
	#exit()
	print(alpha)
	print(beta)
	print(gamma)
	print(basis_set[peptide_plane])
	
	print("Vectors out of new basis set")
	print(vectors)
	vectors_bs = np.matmul(basis_set[peptide_plane], vectors[peptide_plane])
	vectors_ro = rotate_z(vectors_bs, np.pi/2)
	vectors_ro2 = rotate_y(vectors_ro, np.pi/3)
	print("Vectors before rotation;")
	print(vectors_bs)
	print("Vectors after rotation 90 about z")
	print(vectors_ro)
	#print("Vectors after rotation about 60 about y")
	#print(vectors_ro2)
	inv_basis_set = np.linalg.inv(basis_set[peptide_plane])
	vectors_ret = np.matmul(inv_basis_set, vectors_ro)
	print("Vectors back in peptide plane")
	print(vectors_ret)
	
	
	
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
				theta = -float(k[len(k) - 2])
				phi = -float(k[len(k) - 1])
				#theta = 0 * (np.pi / 180.)
				#phi = np.pi
				vectors_bs = np.matmul(basis_set[residue], vectors[residue])
				vectors_theta = rotate_z(vectors_bs, theta) # rotation about beta axis
				vectors_phi = rotate_y(vectors_theta, phi) # rotation about gamma axis
				ibs = np.linalg.inv(basis_set[residue])
				vectors[residue] = np.matmul(ibs, vectors_phi) # return.
			
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
			
			
