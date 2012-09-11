#!/usr/bin/env python

"""This Python script generates input files for ESPResSo++, ESPResSo,
   LAMMPS and GROMACS."""

import sys

chains = 200
monomers = 200
num_particles = chains * monomers
Lx = Ly = Lz = 36.1033

###############################################################################
#                                 ESPResSo++                                  #
###############################################################################
f = open('rings.dat')
data = f.readlines()
f.close()

f = open('espressopp/espressopp_polymer_melt.start', 'w')
f.writelines(data)
f.close()
sys.stdout.write('\nESPResSo++ input file written: espressopp/espressopp_polymer_melt.start\n')


###############################################################################
#                                   LAMMPS                                    #
###############################################################################
f = open('rings.dat')
data = f.readlines()
f.close()

f = open('lammps/lammps_polymer_melt.start', 'w')
f.writelines(data)
f.close()
sys.stdout.write('LAMMPS input file written: lammps/lammps_polymer_melt.start\n')


###############################################################################
#                                  ESPResSo                                   #
###############################################################################
f = open('rings.dat')
data = f.readlines()
f.close()
data = data[20:num_particles + 20]

x = []
y = []
z = []
for d in data:
  line = d.split()[3:]
  x.append(float(line[0]))
  y.append(float(line[1]))
  z.append(float(line[2]))

bonds = []
for i in range(chains):
  for j in range(monomers):
    p = i * monomers + j
    id1 = i * monomers + j
    id2 = id1 + 1
    if(j == 0):
      bonds.append([p, 0, p + monomers - 1, 1, p + monomers - 1, p + 1])
    elif(j == monomers - 1):
      bonds.append([p, 0, p - 1, 1, p - 1, p - monomers + 1])
    else:
      bonds.append([p, 0, p - 1, 1, p - 1, p + 1])

f = open('espresso/espresso_polymer_melt.start', 'w')
f.write('{variable\n')
f.write('\t{time 0.0}\n')
f.write('\t{periodicity 1 1 1}\n')
#f.write('\t{box_l %.4f %.4f %.4f}\n' % (Lx, Ly, Lz))
f.write('\t{n_part %d}\n' % (num_particles))
f.write('}\n')

f.write('{tclvariable\n')
f.write('\t{density 0.85}\n')
f.write('}\n')

f.write('{particles {id pos type v f}\n')
for i in range(num_particles):
  f.write('\t{%d %.3f %.3f %.3f 0 0.0 0.0 0.0 0.0 0.0 0.0}\n' % (i, x[i], y[i], z[i]))
f.write('}\n')
f.write('{bonds\n')

for b in bonds:
  f.write('\t{%d { {%d %d} {%d %d %d} } }\n' % (b[0], b[1], b[2], b[3], b[4], b[5]))
f.write('}\n')
f.close()
sys.stdout.write('ESPResSo input file written: espresso/espresso_polymer_melt.start\n')


###############################################################################
#                                   GROMACS                                   #
###############################################################################
f = open('rings.dat')
data = f.readlines()
f.close()
data = data[20:num_particles + 20]

f = open('gromacs/conf.gro', 'w')
f.write('Rings\n')
f.write('40000\n')

fmt = '%5d%s%03d%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'

ct = 0
for i in range(chains):
  for j in range(1, monomers + 1):
    x = float(data[ct].split()[3])
    y = float(data[ct].split()[4])
    z = float(data[ct].split()[5])
    f.write(fmt % (i + 1, 'RING  A', j, ct + 1, x, y, z, 0.0, 0.0, 0.0))
    ct = ct + 1

s = '  %.4f %.4f %.4f\n' % (Lx, Ly, Lz)
f.write(s)
f.close()
sys.stdout.write('GROMACS input file written: gromacs/gromacs_polymer_melt.start\n\n')
