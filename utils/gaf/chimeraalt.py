'''
http://pymolwiki.org/index.php/bbAxis

Draws a CGO plane across the backbone atoms of neighboring amino acids

Author: Jason Vertrees, 06/2010
  Modified by Thomas Holder, 2010-2012
  Modified by Blaine Bell, 08/2011
	Modified by Ben Tatman, 09/2020

(c) 2010 Schrodinger

License: MIT
'''

from pymol import cmd

## read in fast and slow data

slow = []
fast = []
while (len(slow) < 60):
	slow.append([0, 0, 0])
	fast.append([0, 0, 0])
with open("slow.csv", "r") as f:
	for l in f:
		k = l.split()
		ids = int(k[0])
		sA = float(k[1])/30.
		sB = float(k[2])/30.
		sG = float(k[3])/30.
		slow[ids] = [sA, sB, sG]
with open("fast.csv", "r") as f:
	for l in f:
		k = l.split()
		ids = int(k[0])
		sA = float(k[1])/30.
		sB = float(k[2])/30.
		sG = float(k[3])/30.
		fast[ids] = [sA, sB, sG]
		





def bbAxis(selection='(all)', color='gray', transp=0.3, state=-1, name=None, quiet=1, val=slow):
	"""
DESCRIPTION

	Draws a plane across the backbone for a selection

ARGUMENTS

	selection = string: protein object or selection {default: (all)}

	color = string: color name or number {default: white}

	transp = float: transparency component (0.0--1.0) {default: 0.0}

	state = integer: object state, 0 for all states {default: 1}

NOTES

	You need to pass in an object or selection with at least two
	amino acids.  The plane spans CA_i, O_i, N-H_(i+1), and CA_(i+1)
	"""
	from pymol.cgo import BEGIN, TRIANGLES, COLOR, VERTEX, END, CYLINDER
	from pymol import cgo
	from chempy import cpv

	counter = 1 # counts, denoting the pp.

	# format input
	transp = float(transp)
	state, quiet = int(state), int(quiet)
	if name is None:
		name = cmd.get_unused_name("backbonePlane")

	if state < 0:
		state = cmd.get_state()
	elif state == 0:
		for state in range(1, cmd.count_states(selection) + 1):
			bbAxis(selection, color, transp, state, name, quiet)
		return

	AAs = []
	coords = dict()

	# need hydrogens on peptide nitrogen
	cmd.h_add('(%s) and n. N' % selection)

	# get the list of residue ids
	for obj in cmd.get_object_list(selection):
		
		sel = obj + " and (" + selection + ")"
		for a in cmd.get_model(sel + " and n. CA", state).atom:
			key = '/%s/%s/%s/%s' % (obj, a.segi, a.chain, a.resi)
			AAs.append(key)
			coords[key] = [a.coord, None, None]
		for a in cmd.get_model(sel + " and n. O", state).atom:
			key = '/%s/%s/%s/%s' % (obj, a.segi, a.chain, a.resi)
			if key in coords:
				coords[key][1] = a.coord
		for a in cmd.get_model(sel + " and ((n. N extend 1 and e. H) or (r. PRO and n. CD))", state).atom:
			key = '/%s/%s/%s/%s' % (obj, a.segi, a.chain, a.resi)
			if key in coords:
				coords[key][2] = a.coord

	# need at least two amino acids
	if len(AAs) <= 1:
		("ERROR: Please provide at least two amino acids, the alpha-carbon on the 2nd is needed.")
		return

	# prepare the cgo
	obj = [
		BEGIN, TRIANGLES,
	]
	

	with open("chimera_fast.bild", "w") as chimera_slow:
		for res in range(0, len(AAs) - 1):
			curIdx, nextIdx = str(AAs[res]), str(AAs[res + 1])

			# populate the position array
			pos = [coords[curIdx][0], coords[curIdx][1], coords[nextIdx][2], coords[nextIdx][0]]
			# 0 is n CA, 1 is O, 2 is N-H

			sA = val[counter][0]
			sB = val[counter][1]
			sG = val[counter][2]

			#sA = 1
			#sB = 1
			#sG = 1

			alpha = [9.0]
			beta = [9.0]
			gamma = [9.0]


			# if the data are incomplete for any residues, ignore
			if None in pos:
				if not quiet:
						(' bbAxis: peptide bond %s -> %s incomplete' % (curIdx, nextIdx))
				continue

			if cpv.distance(pos[0], pos[3]) > 4.0:
				if not quiet:
						(' bbAxis: %s and %s not adjacent' % (curIdx, nextIdx))
				continue

			center = [0, 0, 0]
			for i in range(0, len(pos)):
				for j in range(0, 3):
					center[j] = center[j] + pos[i][j]
			for j in range(0, 3):
				center[j] = center[j]/float(len(pos))
					
			beta_normal = cpv.normalize(cpv.cross_product(
				cpv.sub(pos[1], pos[0]),
				cpv.sub(pos[2], pos[0])))
				
			gamma_normal = cpv.normalize(cpv.sub(pos[3], pos[0]))
			
			alpha_normal = cpv.normalize(cpv.cross_product(gamma_normal, beta_normal))
			
			#	(pos)
			
			alpha_start = cpv.sub(center, cpv.scale(alpha_normal, 3. * sA))
			beta_start = cpv.sub(center, cpv.scale(beta_normal, 3. * sB))
			gamma_start = cpv.sub(center, cpv.scale(gamma_normal, 3 * sG))
			
			alpha.extend(alpha_start)
			beta.extend(beta_start)
			gamma.extend(gamma_start)
			
			alpha_end = cpv.add(center, cpv.scale(alpha_normal, 3. * sA))
			beta_end = cpv.add(center, cpv.scale(beta_normal, 3. * sB))
			gamma_end = cpv.add(center, cpv.scale(gamma_normal, 3. * sG))
			
			#alpha_c = [sA * 0.2, 1. * sA, 0., 0., 1. * sA, 0., 0.]
			#beta_c = [sB * 0.2, 0., 1. * sB, 0., 0., 1. * sB, 0.]
			#gamma_c = [sG * 0.2, 0., 0., 1. * sG, 0., 0., 1. * sG]
			alpha_c = [sA * 0.4, 0., 1.*sA, 1.*(1-sA), 0., 1.*sA, 1.*(1-sA)]
			beta_c = [sB * 0.4, 0., 1.*sB, 1.*(1-sB), 0., 1.*sB, 1.*(1-sB)]
			gamma_c = [sG * 0.4, 0., 1.*sG, 1.*(1-sG), 0., 1.*sG, 1.*(1-sG)]
			colors = [0.2, 1., 1., 1., 1., 1., 1.]
			
			alpha.extend(alpha_end)
			beta.extend(beta_end)
			gamma.extend(gamma_end)
			
			alpha.extend(alpha_c)
			beta.extend(beta_c)
			gamma.extend(gamma_c)
			
			chimera_slow.write(".comment residue %d (%f, %f, %f)\n" % (counter, sA*30., sB*30., sG*30.))
			chimera_slow.write(".color %f %f %f\n" % (0, 1.*sA, 1.*(1-sA)))
			chimera_slow.write(".cylinder %f %f %f %f %f %f %f\n" % (alpha_start[0], alpha_start[1], alpha_start[2],\
						alpha_end[0], alpha_end[1], alpha_end[2], 0.4 * sA))
			chimera_slow.write(".color %f %f %f\n" % (0, 1.*sB, 1.*(1-sB)))
			chimera_slow.write(".cylinder %f %f %f %f %f %f %f\n" % (beta_start[0], beta_start[1], beta_start[2],\
						beta_end[0], beta_end[1], beta_end[2], 0.4 * sB))
			chimera_slow.write(".color %f %f %f\n" % (0, 1.*sG, 1.*(1-sG)))
			chimera_slow.write(".cylinder %f %f %f %f %f %f %f\n" % (gamma_start[0], gamma_start[1], gamma_start[2],\
						gamma_end[0], gamma_end[1], gamma_end[2], 0.4 * sG))
			
			
			#alpha.append(alpha_end[0], alpha_end[1], alpha_end[2], 1., 1., 1., 1., 1., 1., 1.)
			#beta.append(beta_end[0], beta_end[1], beta_end[2], 1., 1., 1., 1., 1., 1., 1.)
			#gamma.append(gamma_end[0], gamma_end[1], gamma_end[2], 1., 1., 1., 1., 1., 1., 1.)
			
			cmd.load_cgo(alpha, "alpha%d"%(counter))
			cmd.load_cgo(beta, "beta%d"%(counter))
			cmd.load_cgo(gamma, "gamma%d"%(counter))
			counter = counter + 1
			

			normal = cpv.normalize(cpv.cross_product(
				cpv.sub(pos[1], pos[0]),
				cpv.sub(pos[2], pos[0])))

			obj.append(cgo.NORMAL)
			obj.extend(normal)

			# need to order vertices to generate correct triangles for plane
			if cpv.dot_product(cpv.sub(pos[0], pos[1]), cpv.sub(pos[2], pos[3])) < 0:
				vorder = [0, 1, 2, 2, 3, 0]
			else:
				vorder = [0, 1, 2, 3, 2, 1]

			# fill in the vertex data for the triangles;
			for i in vorder:
				obj.append(VERTEX)
				obj.extend(pos[i])

	# finish the CGO
	obj.append(END)

	# update the UI
	cmd.load_cgo(obj, name, state, zoom=0)
	cmd.set("cgo_transparency", transp, name)
	cmd.color(color, name)

#def bbAxis(selection='(all)', color='gray', transp=0.3, state=-1, name=None, quiet=1):

def slowP():
	cmd.fetch('2qmt')
	mol = "2qmt"
	obj=cmd.get_object_list(mol)[0]
   # cmd.alter(mol,"b=1.0")
	cmd.show_as("cartoon", mol) 
	#cmd.cartoon("putty", mol)
	cmd.set("cartoon_putty_radius", 0.2,obj)
	cmd.set("cartoon_color", 'white',obj)
	cmd.set("cartoon_transparency", 0.4,obj)
	cmd.bg_color('black')
	bbAxis(selection='(all)', color='white', transp=0.2, state=-1, name=None, quiet=1, val= slow)
   
def fastP():
	cmd.fetch('2qmt')
	mol = "2qmt"
	obj=cmd.get_object_list(mol)[0]
   # cmd.alter(mol,"b=1.0")
	cmd.show_as("cartoon", mol) 
	#cmd.cartoon("putty", mol)
	cmd.set("cartoon_putty_radius", 0.2,obj)
	cmd.set("cartoon_color", 'white',obj)
	cmd.set("cartoon_transparency", 0.4,obj)
	cmd.bg_color('black')
	bbAxis(selection='(all)', color='white', transp=0.2, state=-1, name=None, quiet=1, val=fast)
	 

cmd.extend("bbAxis", bbAxis)
cmd.extend("slow", slowP)
cmd.extend("fast", fastP)

# tab-completion of arguments
cmd.auto_arg[0]['bbAxis'] = cmd.auto_arg[1]['color']
cmd.auto_arg[1]['bbAxis'] = cmd.auto_arg[0]['color']

# vi:expandtab
