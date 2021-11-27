# badf_lammps.py ver.2021-11-27
import numpy as np
###############################################################################################
# Calculate normalized bond angle distribution function (BADF).
# 	Input file: lammpstrj file
# 	Output file: BADF_X-Y-Z.dat (X,Y,Z = atom names)
###############################################################################################
def badf_lammps():
##### Setup #####
	filename = 'test.lammpstrj' 

	initial_step = 300
	final_step   = 302
	skip_step    = 2

	atom_names      = ['Mg', 'Si', 'O']
	number_of_atoms = [384,  192,  768]
	
	cutoff_distance = [[0.0 for i in range(len(number_of_atoms))] for j in range(len(number_of_atoms))]
	cutoff_distance[0][2] = 2.5 
	cutoff_distance[1][2] = 2.5
    #############################################################################################
    # cutoff_distance[X][Y]: bond-cutoff distance of atom_names[X]-atom_names[Y] in (angstroms)
    #############################################################################################

	delta_theta = 1.0   # in (degree)
#################

	rad = 180.0/np.arccos(-1.0)
	number_of_type = len(number_of_atoms)
	for j in range(number_of_type):
		for i in range(j + 1, number_of_type):
			cutoff_distance[i][j] = cutoff_distance[j][i]
	typeindex = []
	for i in range(number_of_type):
		for j in range(number_of_atoms[i]):
			typeindex.append(i)
	angle_types = 0
	for t1 in range(number_of_type):
		for t2 in range(number_of_type):
			for t3 in range(t2, number_of_type):
				angle_types = angle_types + 1 
	fractional_coords = [[0, 0, 0] for i in range(sum(number_of_atoms))]
	real_coords = [[0, 0, 0] for i in range(sum(number_of_atoms))]
	array_length = int(180/delta_theta) + 1
	angle = [i*delta_theta for i in range(array_length)]
	badf = [[0 for i in range(array_length)] for j in range(angle_types)]
	rfile = open(filename, 'r')

	while True:
#read timestep
		dat = rfile.readline().split()
		if len(dat) == 0: break
		step = int(rfile.readline())
		target = False
		if step >= initial_step and np.mod(step - initial_step, skip_step) == 0: target = True
#read number of atoms
		rfile.readline()
		dat = int(rfile.readline())
		if int(dat) != sum(number_of_atoms): 
			print('error in number of atoms')
			return 1
#read cell vectors
		rfile.readline()
		x = np.array(rfile.readline().split(), dtype='float')
		y = np.array(rfile.readline().split(), dtype='float')
		z = np.array(rfile.readline().split(), dtype='float')
		if target: 
			box = get_lammps_box(list(x), list(y), list(z))
			cel_vectors = np.array(get_celvectors(box))
#read real coordinates
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
#get fractional coordinates
		if target:
			print('step: %d' %step)
			if lreal: fractional_coords = get_fractionalcoords(real_coords, cel_vectors)
#get neighborlist
			nlist = get_neighborlist(fractional_coords, cel_vectors, typeindex, cutoff_distance)
#search bond angles
			badf_tmp = bond_angle_DF(fractional_coords, cel_vectors, nlist, typeindex, angle)
#step loop
			angle_index = -1
			for t1 in range(number_of_type):
				for t2 in range(number_of_type):
					for t3 in range(t2, number_of_type):
						angle_index = angle_index + 1
						if sum(badf_tmp[angle_index])==0: continue
						for i in range(len(angle)):
							badf[angle_index][i] = badf[angle_index][i] + int(badf_tmp[angle_index][i])
		if step + 1 > final_step: break
#dump ADF
	angle_index = -1
	for t1 in range(number_of_type):
		for t2 in range(number_of_type):
			for t3 in range(t2, number_of_type):
				angle_index = angle_index + 1
				if sum(badf[angle_index]) == 0: continue
				normal = sum(badf[angle_index])*delta_theta
				wfile = open('BADF_%s-%s-%s.dat' %(atom_names[t2], atom_names[t1], atom_names[t3]), 'w')
				for i in range(len(badf[angle_index])):
					badf[angle_index][i] = badf[angle_index][i]/normal
				for i in range(len(angle)):
					wfile.write('%6.1f, %6.4f\n' %(angle[i], badf[angle_index][i]))
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

def get_neighborlist(fractional_coords, cel_vectors, typeindex, cutoff_distances):
	number_of_atoms = len(typeindex)
	fractional_coords = tuple(fractional_coords)
	cel_vectors = tuple(cel_vectors)
	typeindex = tuple(typeindex)
	cutoff_distances = tuple(cutoff_distances)
	neighborlist = [[] for i in range(number_of_atoms)]
	fractional_distance = np.array([0.0, 0.0, 0.0])
	real_distance = np.array([0.0, 0.0, 0.0])
	cel_vectors = np.array(cel_vectors, dtype='float')
	if len(fractional_coords) != number_of_atoms:
		print('Error'); return -1 
	if cel_vectors.shape != (3, 3): 
		print('Error'); return -1 
	if np.array(cutoff_distances).shape != (max(typeindex)+1,max(typeindex)+1):
		print('Error'); return -1 
	for i in range(number_of_atoms-1):
		for j in range(i+1, number_of_atoms):
			if cutoff_distances[typeindex[i]][typeindex[j]] == 0: continue
			for k in range(3):
				fractional_distance[k] = fractional_coords[j][k] - fractional_coords[i][k]
				if abs(fractional_distance[k]) > 0.5: 
					fractional_distance[k] = fractional_distance[k] - np.sign(fractional_distance[k])*1.0
			for k in range(3):
				real_distance[k] = np.dot(fractional_distance,cel_vectors[:,k])
			if np.linalg.norm(real_distance) <= cutoff_distances[typeindex[i]][typeindex[j]]:
				neighborlist[i].append(j); neighborlist[j].append(i)
	return neighborlist

def bond_angle_DF(fractional_coords, cel_vectors, neighborlist, typeindex, angle):
	rad = 180.0/np.arccos(-1.0)
	number_of_types = max(typeindex) + 1
	number_of_atoms = len(typeindex)
	delta_theta = angle[1] - angle[0]
	angle_index = -1
	badf = [[0 for i in range(len(angle))] for j in range(number_of_types*number_of_types*number_of_types)]
	for t1 in range(number_of_types):
		for t2 in range(number_of_types):
			for t3 in range(t2,number_of_types):
				target_angle = [t2, t1, t3]
				angle_index = angle_index + 1
				for i in range(number_of_atoms):
					if typeindex[i] != target_angle[1]: continue
					if len(neighborlist[i]) < 2: continue
					for j1 in range(len(neighborlist[i])):
						for j2 in range(j1 + 1, len(neighborlist[i])):
							n1 = neighborlist[i][j1]; n2 = neighborlist[i][j2]
							if [typeindex[n1], typeindex[n2]] != [target_angle[0], target_angle[2]]:
								if[typeindex[n1], typeindex[n2]] != [target_angle[2], target_angle[0]]: continue
							in1 = np.array(fractional_coords[n1]) - np.array(fractional_coords[i])
							for k in range(3):
								if abs(in1[k]) > 0.5: in1[k] = in1[k] - np.sign(in1[k])*1.0
							for k in range(3): np.dot(in1[k], cel_vectors[:,k])
							in2 = np.array(fractional_coords[n2]) - np.array(fractional_coords[i])
							for k in range(3):
								if abs(in2[k]) > 0.5: in2[k] = in2[k] - np.sign(in2[k])*1.0
							for k in range(3): np.dot(in2[k], cel_vectors[:,k])
							theta = np.arccos(np.dot(in1,in2)/(np.linalg.norm(in1)*np.linalg.norm(in2)))*rad
							for k in range(len(angle)):
								if angle[k] - theta >= -delta_theta/2.0 and angle[k] - theta < delta_theta/2.0:
									badf[angle_index][k] = badf[angle_index][k] + 1
									break
	return badf

badf_lammps()
