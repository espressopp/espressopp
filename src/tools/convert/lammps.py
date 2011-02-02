# -*- coding: iso-8859-1 -*-
"""This Python module allows one to use a LAMMPS data file as the
   input to an ESPResSo++ simulation."""

def read(fin):

  f = open(fin)
  line = f.readline() # comment line
  line = f.readline() # blank line
  num_particles = int(f.readline().split()[0])
  num_bonds = int(f.readline().split()[0])
  num_angles = int(f.readline().split()[0])
  num_dihedrals = int(f.readline().split()[0])
  line = f.readline() # impropers
  line = f.readline() # blank line
  num_types = int(f.readline().split()[0]) # atom types

  # find and store size of box
  line = ''
  while not 'xlo' in line:
    line = f.readline()
  xmin, xmax = map(float, line.split()[0:2])
  ymin, ymax = map(float, f.readline().split()[0:2])
  zmin, zmax = map(float, f.readline().split()[0:2])
  Lx = xmax - xmin
  Ly = ymax - ymin
  Lz = zmax - zmin

  # find and store coordinates
  line = ''
  while not 'Atoms' in line:
    line = f.readline()
  line = f.readline()

  if(num_types == 1):
    rstart = 3
    if(num_bonds == 0): rstart = 2

    x = []
    y = []
    z = []
    for i in range(num_particles):
      rx, ry, rz = map(float, f.readline().split()[rstart:])
      x.append(rx)
      y.append(ry)
      z.append(rz)
  else:
    p_type = []
    q = []
    x = []
    y = []
    z = []
    for i in range(num_particles):
      k, rq, rx, ry, rz = map(float, f.readline().split()[2:])
      p_type.append(int(k))
      q.append(rq)
      x.append(rx)
      y.append(ry)
      z.append(rz)

  if(num_bonds != 0):
    # find and store bonds
    line = ''
    while not 'Bonds' in line:
      line = f.readline()
    line = f.readline()
    bonds = []
    for i in range(num_bonds):
      bond_id, bond_type, pid1, pid2 = map(int, f.readline().split())
      bonds.append((pid1, pid2))

  if(num_angles != 0):
    # find and store angles
    line = ''
    while not 'Angles' in line:
      line = f.readline()
    line = f.readline()
    angles = []
    for i in range(num_angles):
      angle_id, angle_type, pid1, pid2, pid3 = map(int, f.readline().split())
      angles.append((pid1, pid2, pid3))


  if(num_dihedrals != 0):
    # find and store angles
    line = ''
    while not 'Dihedrals' in line:
      line = f.readline()
    line = f.readline()
    dihedrals = []
    for i in range(num_dihedrals):
      dihedral_id, dihedral_type, pid1, pid2, pid3, pid4 = map(int, f.readline().split())
      dihedrals.append((pid1, pid2, pid3, pid4))


    f.close()


  if(num_types != 1):
    return p_type, bonds, angles, q, x, y, z, Lx, Ly, Lz

  if(num_bonds == 0 and num_angles == 0 and num_dihedrals == 0):
    return x, y, z, Lx, Ly, Lz
  elif(num_bonds != 0 and num_angles == 0 and num_dihedrals == 0):
    return bonds, x, y, z, Lx, Ly, Lz
  elif(num_bonds != 0 and num_angles != 0 and num_dihedrals == 0):
    return bonds, angles, x, y, z, Lx, Ly, Lz
  else:
    return bonds, angles, dihedrals, x, y, z, Lx, Ly, Lz
