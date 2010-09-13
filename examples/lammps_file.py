def read(fin):

  f = open(fin)
  line = f.readline()
  line = f.readline()
  num_particles = int(f.readline().split()[0])
  num_bonds = int(f.readline().split()[0])
  num_angles = int(f.readline().split()[0])

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
    f.close()

  if(num_bonds == 0 and num_angles == 0):
    return x, y, z, Lx, Ly, Lz
  elif(num_bonds != 0 and num_angles == 0):
    return bonds, x, y, z, Lx, Ly, Lz
  else:
    return bonds, angles, x, y, z, Lx, Ly, Lz
