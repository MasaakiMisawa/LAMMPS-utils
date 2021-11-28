# rdf_lammps.py ver.2021-11-28
import numpy as np
import modules as mm
###############################################################################################
# Calculate radial distribution function (RDF).
# 	Input file: lammpstrj file
# 	Output file: RDF_X-Y.dat (X,Y = atom names)
###############################################################################################
def rdf_lammps():
##### Setup #####
	filename = 'test.lammpstrj' 

	initial_step = 0
	final_step   = 1
	skip_step    = 1

	atom_names      = ['Mg', 'Si', 'O']
	number_of_atoms = [384,  192,  768]
	
	dr        = 0.10      # in (angstroms)
	initial_r = 0.05      # in (angstroms)
	max_r     = 10.0      # <= 10
#################

	number_of_type = len(number_of_atoms)
	typeindex = []
	for i in range(number_of_type):
		for j in range(number_of_atoms[i]):
			typeindex.append(i)
	pair_types = 0
	for t1 in range(number_of_type):
		for t2 in range(t1, number_of_type):
			pair_types = pair_types + 1 
	fractional_coords = [[0, 0, 0] for i in range(sum(number_of_atoms))]
	real_coords = [[0, 0, 0] for i in range(sum(number_of_atoms))]
	i = 0; radius = []
	while initial_r + dr*i <= max_r:
		radius.append(initial_r + dr*i)
		i = i + 1
	rdf = [[0 for i in range(len(radius))] for j in range(pair_types)]
	rfile = open(filename, 'r')
	nsample = 0

	while True:
		dat = rfile.readline().split()
		if len(dat) == 0: break
		step = int(rfile.readline())
		target = False
		if step >= initial_step and np.mod(step - initial_step, skip_step) == 0: target = True
		rfile.readline()
		dat = int(rfile.readline())
		if int(dat) != sum(number_of_atoms): 
			print('error in number of atoms')
			return 1
		rfile.readline()
		x = np.array(rfile.readline().split(), dtype='float')
		y = np.array(rfile.readline().split(), dtype='float')
		z = np.array(rfile.readline().split(), dtype='float')
		if target: 
			box = get_lammps_box(list(x), list(y), list(z))
			cel_vectors = np.array(get_celvectors(box))
		dat = rfile.readline().split()
		i_id = 0; i_x = 2; lreal = True
		for i in range(len(dat)):
			if dat[i] == 'id': i_id = i - 2
			if dat[i] == 'x': i_x = i - 2
			if dat[i] == 'xs': 
				i_x = i - 2
				lreal = False
		for i in range(sum(number_of_atoms)):
			dat = rfile.readline().split()
			atom_id = int(dat[i_id])
			if lreal: real_coords[atom_id - 1] = [float(dat[i_x]), float(dat[i_x + 1]), float(dat[i_x + 2])]
			else: fractional_coords[atom_id - 1] = [float(dat[i_x]), float(dat[i_x + 1]), float(dat[i_x + 2])]
		if target:
			print('step: %d' %step)
			if lreal: fractional_coords = get_fractionalcoords(real_coords, cel_vectors)
			rdf_tmp = radial_DF(fractional_coords, cel_vectors, typeindex, radius)
			nsample = nsample + 1
			gr_index = -1
			for t1 in range(number_of_type):
				for t2 in range(t1, number_of_type):
					gr_index = gr_index + 1
					if sum(rdf_tmp[gr_index])==0: continue
					for i in range(len(radius)):
						rdf[gr_index][i] = rdf[gr_index][i] + rdf_tmp[gr_index][i]
		if step + 1 > final_step: break
	gr_index = -1
	for t1 in range(number_of_type):
		for t2 in range(t1, number_of_type):
			gr_index = gr_index + 1
			if sum(rdf[gr_index]) == 0: continue
			wfile = open('RDF_%s-%s.dat' %(atom_names[t1], atom_names[t2]), 'w')
			for i in range(len(rdf[gr_index])):
				rdf[gr_index][i] = rdf[gr_index][i]/nsample
			for i in range(len(radius)):
				wfile.write('%6.2f, %6.4f\n' %(radius[i], rdf[gr_index][i]))
			wfile.close()
	rfile.close()


def get_lammps_box(x,y,z):
	alpha = 90.0; beta = 90.0; gamma = 90.0; rad = 180.0/np.arccos(-1.0)
	if [len(x), len(y), len(z)] == [2, 2, 2]:
		la = x[1] - x[0]; lb = y[1] - y[0]; lc = z[1] - z[0]
		box = [la, lb, lc, alpha, beta, gamma]
		return box
	elif [len(x), len(y), len(z)] == [3, 3, 3]:
		xy = x[2]; xz = y[2]; yz = z[2]
		xlo = x[0] - min(0.0, xy, xz, xy+xz); xhi = x[1] - max(0.0, xy, xz, xy+xz)
		ylo = y[0] - min(0.0, yz); yhi = y[1] - max(0.0, yz)
		zlo = z[0]; zhi = z[1]
		lx = xhi - xlo; ly = yhi - ylo; lz = zhi - zlo
		la = lx
		lb = np.sqrt(ly*ly + xy*xy)
		lc = np.sqrt(lz*lz + xz*xz + yz*yz)
		alpha = np.arccos((xy*xz + ly*yz)/(lb*lc))*rad
		beta = np.arccos(xz/lc)*rad
		gamma = np.arccos(xy/lb)*rad
		box = [la, lb, lc, alpha, beta, gamma]
		return box
	else:
		print('error in get_lammps_box')
		return 1

def get_celvectors(lengths_angles):
	deg_rad = 180.0/np.arccos(-1.0)
	if len(lengths_angles) != 6: 
		print('error in get_celvectors'); return -1 
	ax = lengths_angles[0]
	ay = 0.0
	az = 0.0                     
	bx = lengths_angles[1]*np.cos(lengths_angles[5]/deg_rad)
	by = np.sqrt(lengths_angles[1]*lengths_angles[1] - bx*bx)
	bz = 0.0
	cx = lengths_angles[2]*np.cos(lengths_angles[4]/deg_rad)
	cy = (lengths_angles[1]*lengths_angles[2]*np.cos(lengths_angles[3]/deg_rad) - bx*cx)/by
	cz = np.sqrt(lengths_angles[2]*lengths_angles[2] - cx*cx - cy*cy)
	cel_vectors = [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]
	return cel_vectors

def get_fractionalcoords(real_coords, cel_vectors):
	inv_cel_vectors = np.linalg.inv(cel_vectors)
	number_of_atoms = len(real_coords)
	fractional_coords = [[0.0, 0.0, 0.0] for i in range(number_of_atoms)]
	for i in range(number_of_atoms):
		for j in range(3):
			fractional_coords[i][j] = np.dot(real_coords[i],inv_cel_vectors[:,j])
	return fractional_coords

def radial_DF(fractional_coords, cel_vectors, typeindex, radius):
	dr = radius[1] - radius[0]
	initial_r = radius[0]
	zeros = [0 for i in range(len(radius))]
	max_r = min(10.0, max(radius))
	multi = [1, 1, 1]
	for i in range(3):
		if cel_vectors[i][i] < max_r*2.0:
			multi[i] = int(20.0/cel_vectors[i][i]) + 1
	fractional_coords_tmp = [0, 0, 0]
	number_of_atoms = len(fractional_coords)
	for j0 in range(multi[0]):
		for j1 in range(multi[1]):
			for j2 in range(multi[2]):
				if [j0, j1, j2] == [0, 0, 0]: continue
				for i in range(number_of_atoms):
					fractional_coords_tmp[0] = fractional_coords[i][0] + j0
					fractional_coords_tmp[1] = fractional_coords[i][1] + j1
					fractional_coords_tmp[2] = fractional_coords[i][2] + j2
					fractional_coords.append([fractional_coords_tmp[0], fractional_coords_tmp[1], fractional_coords_tmp[2]])
	for i in range(len(fractional_coords)):
		for j in range(3):
			fractional_coords[i][j] = fractional_coords[i][j]/multi[j]
	for i in range(3):
		for j in range(3):
			cel_vectors[i][j] = cel_vectors[i][j]*multi[i]
	volume = np.dot(np.cross(cel_vectors[0],cel_vectors[1]), cel_vectors[2])
	number_of_types = max(typeindex) + 1
	rdf = []; n = 0
	rdf_index = [[-1 for i in range(number_of_types)] for j in range(number_of_types)]
	for i in range(number_of_types): #0-0, 0-1, 0-2, 1-1, 1-2, 2-2
		for j in range(i, number_of_types):
			rdf.append(list(zeros))
			rdf_index[i][j] = n
			n = n + 1
	for i in range(1,multi[0]*multi[1]*multi[2]):
		typeindex = typeindex + typeindex
	number_of_atoms = [0 for i in range(max(typeindex) + 1)]
	for i in range(len(typeindex)): number_of_atoms[typeindex[i]] = number_of_atoms[typeindex[i]] + 1
	fractional_distance = [0, 0, 0]; real_distance = [0, 0, 0]
	for i in range(sum(number_of_atoms)-1):
		for j in range(i+1,sum(number_of_atoms)):
			for k in range(3):
				fractional_distance[k] = fractional_coords[j][k] - fractional_coords[i][k]
				if abs(fractional_distance[k]) > 0.5: fractional_distance[k] = fractional_distance[k] - np.sign(fractional_distance[k])*1.0
				real_distance[k] = np.dot(fractional_distance[k], cel_vectors[:][k])
			rij = np.linalg.norm(real_distance)
			if rij <= max_r:
				k = int((rij + (dr/2.0) - initial_r)/dr)  
				#deno = 4.0*np.arccos(-1.0)*(radius[k]-(dr/2.0))*(radius[k]-(dr/2.0))*dr
				deno = 4.0*np.arccos(-1.0)*radius[k]*radius[k]*dr # VMD definition
				deno = deno*number_of_atoms[typeindex[i]]*number_of_atoms[typeindex[j]]/volume
				rdf[rdf_index[typeindex[i]][typeindex[j]]][k] = rdf[rdf_index[typeindex[i]][typeindex[j]]][k] + 1/deno
				rdf[rdf_index[typeindex[j]][typeindex[i]]][k] = rdf[rdf_index[typeindex[j]][typeindex[i]]][k] + 1/deno
			
	return rdf

rdf_lammps()
