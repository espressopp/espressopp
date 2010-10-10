"""This Python module allows one to use GROMACS data files as the
   input to an ESPResSo++ simulation."""

def read(gro_file, top_file, itp_file):

  # gro_file contains number of particles, positions, velocities and box size
  # top_file contains topology information
  # itp_file contains topology information

  f = open(gro_file)
  line = f.readline()
  num_particles = int(f.readline())

  # store coordinates
  x = []
  y = []
  z = []
  for i in range(num_particles):
    s = f.readline()[20:44]
    rx = float(s[0:8])
    ry = float(s[8:16])
    rz = float(s[16:24])
    x.append(rx)
    y.append(ry)
    z.append(rz)

  # store box size
  Lx, Ly, Lz = map(float, f.readline().split())
  f.close()

  f = open(top_file)
  # find and store number of chains
  line = ''
  while not 'molecules' in line:
    line = f.readline()
  num_chains = int(f.readline().split()[1])
  monomers = num_particles / num_chains
  f.close()

  f = open(itp_file)
  # find and store bonds
  line = ''
  while not 'bonds' in line:
    line = f.readline()

  bonds = []
  line = f.readline()
  line = f.readline()
  while(len(line) != 1):
    pid1, pid2 = map(int, line.split()[0:2])
    bonds.append((pid1, pid2))
    line = f.readline()

  # find and store angles
  line = ''
  while not 'angles' in line:
    line = f.readline()

  angles = []
  line = f.readline()
  while(len(line) != 0):
    pid1, pid2, pid3 = map(int, line.split()[0:3])
    angles.append((pid1, pid2, pid3))
    line = f.readline()
  f.close()

  # extend bonds to num_chains - 1 chains
  bonds_per_chain = len(bonds)
  for i in range(num_chains - 1):
    for j in range(bonds_per_chain):
      pid1 = bonds[j][0]
      pid2 = bonds[j][1]
      bonds.append((pid1 + (i + 1) * monomers, pid2 + (i + 1) * monomers))

  # extend angles to num_chains - 1 chains
  angles_per_chain = len(angles)
  for i in range(num_chains - 1):
    for j in range(angles_per_chain):
      pid1 = angles[j][0]
      pid2 = angles[j][1]
      pid3 = angles[j][2]
      angles.append((pid1 + (i + 1) * monomers, \
                     pid2 + (i + 1) * monomers, \
                     pid3 + (i + 1) * monomers))

  return bonds, angles, x, y, z, Lx, Ly, Lz
