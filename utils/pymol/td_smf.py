from pymol import cmd
import sys
import statistics
fn = "smf/final.dat"
col = 4
if (len(sys.argv) > 2):
	fn = sys.argv[1]
	col = int(sys.argv[2])
	mi = float(sys.argv[3])
	ma = float(sys.argv[4])
	out = sys.argv[5]
	
cmd.fetch('2qmt')
mol = "2qmt"
obj=cmd.get_object_list(mol)[0]
cmd.alter(mol,"b=-1.0")
cmd.show_as("cartoon", mol)

bfacts = []
ignore = []
with open(fn, "r") as f:
	for l in f:
		k = l.split("\t")
		print(len(k))
		print(k)
		bfact = float(k[col])
		if (float(k[col]) > 0):
			bfacts.append(float(k[col]))
		else:
			ignore.append(int(k[0]))
		print("Setting residue %d to b=%f" % (int(k[0]), float(k[col])))
		cmd.alter("%s and resi %d and n. CA"%(mol, int(k[0])), "b=%f" % (float(k[col])))


cmd.set("cartoon_putty_scale_min", 0, obj)
cmd.set("cartoon_putty_scale_max", 10, obj)
cmd.set("cartoon_putty_transform", 4,obj)
cmd.set("cartoon_putty_radius", 0.2,obj)
cmd.set("cartoon_discrete_colors", 1, obj)
cmd.set("ray_trace_mode", 3)

#spectrum b, blue_white_red, minimum=20, maximum=50

cmd.spectrum("b", "blue_white_red", minimum=mi, maximum=ma)
#cmd.spectrum("b","rainbow", "%s and n. CA " %mol)
#cmd.ramp_new("count", obj, [mi, ma], "rainbow")
for i in ignore:
	cmd.alter("%s and resi %d and n. CA"%(mol, i), "b=%f" % (statistics.mean(bfacts)))
	cmd.color("black", 'resi ' + str(i))

cmd.recolor()

#cmd.cartoon("putty", mol)
cmd.zoom()
cmd.png(out, 1000, 800, dpi=150, ray=1)

