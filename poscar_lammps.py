#poscar_lammps.py ver.2021.11.04
import numpy as np; import datetime

# Description #######################################################################
# This script generates a geometrical data file for LAMMPS based on VASP POSCAR file.
#   VASP POSCAR file format:
#     2nd line:    lattice scaling factor
#     3-5th lines: lattice vectors (3rd, 4th, and 5th lines = a, b and c vectors)
#     7th line:    number of atom (should be separeted according to atom types)
#     9th- lines:  atomic fractional coordinates (not "cartesian")
#####################################################################################

def make_lammpsdata():
## Setup ##
  filename = 'POSCAR.vasp'    # filename of VASP POSCAR to be read
  ntype  = 2                  # number of atom types
  nbtype = 1                  # number of bond types
  natype = 1                  # number of angle types
  ndtype = 0                  # number of dihedral types
  nitype = 0                  # number of improper types

  astyle = 'full'             # enter 'atomic', 'charge', 'molecular', or 'full' 

  q = [0 for i in range(ntype)]; mass = [0 for i in range(ntype)]           
  q[0]     = -0.82000         # charge of atom type 0（only for astyle = 'charge' or 'full'）
  q[1]     = 0.41000          # 
  mass[0]  = 16.000           # mass of atom type 0
  mass[1]  = 1.010            # 

  btype = [0 for i in range(nbtype)]
  btype[0] = [0, 1]           # component of atom types of bond type 0        

  atype = [0 for i in range(natype)]
  atype[0] = [1, 0, 1]        # component of atom types of angle type 0  

  dtype = [0 for i in range(ndtype)]
#  dtype[0] = [0, 1, 1, 1]     # component of atom types of dihedral type 0

  itype = [0 for i in range(nitype)]
#  itype[0] = [1, 1, 1, 1]     # component of atom types of improper type 0

  mul = [2, 2, 2]             # supercell to be created = mul[0]*mul[1]*mul[2] of unit cell

  bcoeff = ' 1 22.965000 1.0120000'          # sting to be written on "Bond Coeffs" term 
  acoeff = ' 1 harmonic 1.6456800 113.24000' # sting to be written on "Angle Coeffs" term
  dcoeff = ''                                # sting to be written on "Dihedral Coeffs" term
  icoeff = ''                                # sting to be written on "Improper Coeffs" term

  cutoff = [[0 for i in range(ntype)] for j in range (ntype)]
  cutoff[0][0] = 0.0          # cutoff distance of atom type 0-0 in angstroms
  cutoff[0][1] = 1.2       
  cutoff[1][1] = 0.0        

  # Note ##########################################################################
  # Neighbor list is generated based on the cutoff distances.
  # If an interatomic distance between atom i and j <= cutoff distance,        
  # the atom j is listed as a neighbor of atom i.
  # Molecule-ID, bond, angle, dihedral, and improper lists are generated
  # based on the neighbor list. 
  #################################################################################


## Prepare ##
  for i in range(ntype - 1):
    for j in range(i+1,ntype):
      cutoff[j][i] = cutoff[i][j]
  natom = [0 for i in range(ntype)]; 
  blist = []; alist = []; dlist = []; ilist = [] 

## Read data ##
  rf = open(filename, 'r')
  wf = open('lammpsdata.dat', 'w')
  for i in range(2): dat = rf.readline()
  scf = float(dat)
  # Get cell vectors
  cel = np.array([])
  for i in range(3):
    cel = np.append(cel, np.array(rf.readline().split(), dtype='float')*scf )
  cel = cel.reshape(3,3)
  # Get cell parameters
  pi = np.arccos(-1); rtd = 180/pi;
  la = np.linalg.norm(cel[0,:]); lb = np.linalg.norm(cel[1,:]); lc = np.linalg.norm(cel[2,:])
  alp = np.arccos(np.dot(cel[1,:],cel[2,:])/(lb*lc))
  bet = np.arccos(np.dot(cel[0,:],cel[2,:])/(la*lc))
  gam = np.arccos(np.dot(cel[1,:],cel[0,:])/(lb*la))
  print('unit cell: %.3f %.3f %.3f %.3f %.3f %.3f' %(la, lb, lc, alp*rtd, bet*rtd, gam*rtd))
  la = la*mul[0]; lb = lb*mul[1]; lc = lc*mul[2]
  print('supercell: %.3f %.3f %.3f %.3f %.3f %.3f' %(la, lb, lc, alp*rtd, bet*rtd, gam*rtd))
  lx = la                      # inverse relationship in https://docs.lammps.org/Howto_triclinic.html
  xy = lb*np.cos(gam)
  xz = lc*np.cos(bet)
  ly = np.sqrt(lb*lb - xy*xy)
  yz = (lb*lc*np.cos(alp) - xy*xz)/ly
  lz = np.sqrt(lc*lc - xz*xz - yz*yz)
  # Get total number of atoms
  dat = rf.readline()
  dat = rf.readline().split()
  if len(dat) != ntype: print('Warning: ntype does not much with POSCAR file!\n')
  for i in range(ntype): natom[i] = int(dat[i])
  ntot = 0
  for i in range(ntype): ntot = ntot + natom[i]
  ntot = ntot*mul[0]*mul[1]*mul[2]
  print('Number of atoms: %d' %ntot)
  dat = rf.readline()
  # Get fractional coordinates
  fc = [[] for i in range(ntot)]
  cnt = 0
  for i in range(sum(natom)):
    dat = rf.readline().split()
    for iz in range(mul[2]):
      for iy in range(mul[1]):
        for ix in range(mul[0]):
          fc[cnt] = [float(dat[0])+ix, float(dat[1])+iy, float(dat[2])+iz] 
          cnt = cnt + 1 

## Create type list ##
  tlist = []; 
  for i in range(ntype): 
    for j in range(natom[i]): 
      for k in range(mul[0]*mul[1]*mul[2]):
        tlist.append(i)

## Create neighbor list ##
  nlist = [[] for i in range(ntot)]
  for i in range(ntot):
    for j in range(ntot):
      if i == j: continue
      dfx = float(fc[j][0]) - float(fc[i][0])
      if abs(dfx) > 0.5*mul[0]: dfx = dfx - np.sign(dfx)*mul[0]
      dfy = float(fc[j][1]) - float(fc[i][1])
      if abs(dfy) > 0.5*mul[1]: dfy = dfy - np.sign(dfy)*mul[1]
      dfz = float(fc[j][2]) - float(fc[i][2])
      if abs(dfz) > 0.5*mul[2]: dfz = dfz - np.sign(dfz)*mul[2]
      dx = dfx*cel[0][0] + dfy*cel[1][0] + dfz*cel[2][0]
      dy = dfx*cel[0][1] + dfy*cel[1][1] + dfz*cel[2][1]
      dz = dfx*cel[0][2] + dfy*cel[1][2] + dfz*cel[2][2]
      df = np.sqrt(dx*dx + dy*dy + dz*dz)
      if df <= cutoff[tlist[i]][tlist[j]]: 
        nlist[i].append(j) 

## Create molecule list ##
  if astyle == 'molecular' or astyle == 'full':
    mlist = [[] for i in range(ntot)]
    nmtot = 0
    for i in range(ntot):
      i2 = 0
      if mlist[i] == []: 
        mlist[i] = nmtot
        nmtot = nmtot + 1
      for j in range(len(nlist[i])):
        if mlist[nlist[i][j]] == []: mlist[nlist[i][j]] = mlist[i]
      lflag = 1
      while True:
        if mlist[i2] == mlist[i]:
          for k in range(len(nlist[i2])):
            if mlist[nlist[i2][k]] == []:
              lflag = 0
              mlist[nlist[i2][k]] = mlist[i]
        i2 = i2 + 1
        if i2 == ntot: 
          if lflag == 1: break
          i2 = 0; lflag = 1
    print('Number of Molecules: %d' %nmtot)

## Create bond list ##
  nbtot = 0
  if nbtype != 0:
    for i in range(ntot): 
      for j in range(len(nlist[i])):
        if i < nlist[i][j]:
          for k in range(nbtype):
            if [tlist[i],tlist[nlist[i][j]]] == btype[k] or [tlist[nlist[i][j]],tlist[i]] == btype[k]: 
              nbtot = nbtot + 1
              blist.append([k, i, nlist[i][j]])

  print('Number of Bonds: %d' %nbtot)

## Create angle list ##
  natot = 0
  if natype != 0:
    for i in range(ntot): 
      if len(nlist[i]) >= 2:
        for j in range(len(nlist[i]) - 1):
          for k in range(j+1,len(nlist[i])): 
            for l in range(natype):
              if [tlist[nlist[i][j]],tlist[i],tlist[nlist[i][k]]] == atype[l] or [tlist[nlist[i][k]],tlist[i],tlist[nlist[i][j]]] == atype[l]: 
                natot = natot + 1
                alist.append([l, nlist[i][j], i, nlist[i][k]])

  print('Number of Angles: %d' %natot)

## Create dihedral list ##
  ndtot = 0
  if ndtype != 0:
    for i in range(ntot):
      for j in range(len(nlist[i])):
        for k in range(ndtype):
          if tlist[i] == dtype[k][0] and tlist[nlist[i][j]] == dtype[k][1]:
            d1 = i; d2 = nlist[i][j]
            for j2 in range(len(nlist[d2])):
              if tlist[nlist[d2][j2]] == dtype[k][2]:
                d3 = nlist[d2][j2]
                if d3 != d1:
                  for j3 in range(len(nlist[d3])):
                    if tlist[nlist[d3][j3]] == dtype[k][3]:
                      d4 = nlist[d3][j3]
                      if d4 != d2 and d1 < d4:
                        ndtot = ndtot + 1
                        dlist.append([k, i, nlist[i][j], nlist[nlist[i][j]][j2], nlist[nlist[nlist[i][j]][j2]][j3]])
  print('Number of dihedrals: %d' %ndtot)

## Create improper list ##
  nitot = 0; limp = 0;
  if nitype != 0:
    for i in range(ntot):
      if len(nlist[i]) >= 3:
        for j in range(len(nlist[i]) - 2):
          for k in range(j+1, len(nlist[i])-1):
            for l in range(k+1, len(nlist[i])):
              for m in range(nitype):
                if [tlist[i], tlist[nlist[i][j]], tlist[nlist[i][k]], tlist[nlist[i][l]]] == itype[m]: limp = 1
                elif [tlist[i], tlist[nlist[i][l]], tlist[nlist[i][j]], tlist[nlist[i][k]]] == itype[m]: limp = 1
                elif [tlist[i], tlist[nlist[i][k]], tlist[nlist[i][l]], tlist[nlist[i][j]]] == itype[m]: limp = 1
                elif [tlist[i], tlist[nlist[i][j]], tlist[nlist[i][l]], tlist[nlist[i][k]]] == itype[m]: limp = 1
                elif [tlist[i], tlist[nlist[i][k]], tlist[nlist[i][j]], tlist[nlist[i][l]]] == itype[m]: limp = 1
                elif [tlist[i], tlist[nlist[i][l]], tlist[nlist[i][k]], tlist[nlist[i][j]]] == itype[m]: limp = 1
                if limp == 1:
                  nitot = nitot + 1; limp = 0
                  ilist.append([m, i, nlist[i][j], nlist[i][k], nlist[i][l]])
  print('Number of impropers: %d' %nitot)

## Output Header ##
  dd = datetime.datetime.now()
  wf.write('LAMMPS data file generated by poscar_lammps.py (Masaaki MISAWA) on %s\n' %dd)
  wf.write(' %d atoms\n %d bonds\n %d angles\n %d dihedrals\n %d impropers\n' %(ntot, nbtot, natot, ndtot, nitot))
  wf.write(' %d atom types\n %d bond types\n %d angle types\n %d dihedral types\n %d improper types\n\n' %(ntype, nbtype, natype, ndtype, nitype))

## Output Cell ##
  wf.write(' %.6f %.6f xlo xhi\n' %(0.0, lx))
  wf.write(' %.6f %.6f ylo yhi\n' %(0.0, ly))
  wf.write(' %.6f %.6f zlo zhi\n' %(0.0, lz))
  wf.write(' %.6f %.6f %.6f xy xz yz\n' %(xy, xz, yz))  

## Output Masses ##
  wf.write('\nMasses\n\n')
  for i in range(ntype):
    wf.write(' %d %6.3f\n' %(i+1, mass[i]))

## Output Coeffs ##
  wf.write('\nBond Coeffs\n\n%s\n\n' %bcoeff)
  wf.write('\nAngle Coeffs\n\n%s\n\n' %acoeff)
  wf.write('\nDihedral Coeffs\n\n%s\n\n' %dcoeff)
  wf.write('\nImproper Coeffs\n\n%s\n\n' %icoeff)

## Output Atoms ##
  wf.write('\nAtoms\n\n')
  for i in range(ntot):
    rx = float(fc[i][0])*cel[0][0] + float(fc[i][1])*cel[1][0] + float(fc[i][2])*cel[2][0]
    ry = float(fc[i][0])*cel[0][1] + float(fc[i][1])*cel[1][1] + float(fc[i][2])*cel[2][1]
    rz = float(fc[i][0])*cel[0][2] + float(fc[i][1])*cel[1][2] + float(fc[i][2])*cel[2][2]
    if astyle == 'full':
      wf.write(' %3d %3d %2d %12.6f %12.6f %12.6f %12.6f\n' %(i+1, mlist[i]+1, tlist[i]+1, q[tlist[i]], rx, ry, rz))
    elif astyle == 'molecular':
      wf.write(' %3d %3d %2d %12.6f %12.6f %12.6f\n' %(i+1, mlist[i]+1, tlist[i]+1, rx, ry, rz))
    elif astyle == 'charge':
      wf.write(' %3d %2d %12.6f %12.6f %12.6f %12.6f\n' %(i+1, tlist[i]+1, q[tlist[i]], rx, ry, rz))
    elif astyle == 'atomic':
      wf.write(' %3d %2d %12.6f %12.6f %12.6f\n' %(i+1, tlist[i]+1, rx, ry, rz))

## Output Bonds ##
  if nbtype != 0:
    wf.write('\n\nBonds\n\n')
    for i in range(nbtot): 
      wf.write(' %3d %2d %4d %4d\n' %(i+1, blist[i][0]+1, blist[i][1]+1, blist[i][2]+1))

## Output Angles ##
  if natype != 0:
    wf.write('\n\nAngles\n\n')
    for i in range(natot): 
      wf.write(' %3d %2d %4d %4d %4d\n' %(i+1, alist[i][0]+1, alist[i][1]+1, alist[i][2]+1, alist[i][3]+1))

## Output Dihedrals ##
  if ndtype != 0:
    wf.write('\n\nDihedrals\n\n')
    for i in range(ndtot): 
      wf.write(' %3d %2d %4d %4d %4d %4d\n' %(i+1, dlist[i][0]+1, dlist[i][1]+1, dlist[i][2]+1, dlist[i][3]+1, dlist[i][4]+1))

## Output Impropers ##
  if nitype != 0:
    wf.write('\n\nImpropers\n\n')
    for i in range(nitot): 
      wf.write(' %3d %2d %4d %4d %4d %4d\n' %(i+1, ilist[i][0]+1, ilist[i][1]+1, ilist[i][2]+1, ilist[i][3]+1, ilist[i][4]+1))

  rf.close(); wf.close()
 
make_lammpsdata()

