#!/usr/bin/env python

###########################################################################
#                                                                         #
#  This Python script generates an input file to be used for              #
#  LAMMPS or ESPRESSO++                                                   #
#                                                                         #
#  -> generates file data.lj                                              #
#                                                                         #
#  Input: N    (size of lattice)                                          #
#                                                                         #
#     Output file will contain N * N * N particles on a lattice           #
#                                                                         #
###########################################################################

import sys
import random
from espressopp.tools import lattice, velocities


# cubic lattice with user-defined values of N and rho
# num_particles should be a perfect cube (e.g. 25**3=15625, 32**3=32768)

particles_per_direction = 32
num_particles = particles_per_direction**3
x, y, z, Lx, Ly, Lz = lattice.create(num_particles, rho=0.8442, perfect=False)
vx, vy, vz = velocities.gaussian(T=0.6, N=num_particles, zero_momentum=True)

###############################################################################
#                                  ESPResSo++                                 #
###############################################################################
file = open("espressopp/espressopp_lennard_jones.start", "w")
file.write("LAMMPS\n")
file.write("\n")
file.write("%5d atoms\n"%num_particles)
file.write("    0 bonds\n")
file.write("    0 angles\n")
file.write("    0 dihedrals\n")
file.write("    0 impropers\n")
file.write("\n")
file.write("    1 atom types velocities\n")
file.write("    0 bond types\n")
file.write("    0 angle types\n")
file.write("    0 dihedral types\n")
file.write("    0 improper types\n")
file.write("\n")

file.write("%.4f %.4f xlo xhi\n" % (0.0, Lx))
file.write("%.4f %.4f ylo yhi\n" % (0.0, Ly))
file.write("%.4f %.4f zlo zhi\n" % (0.0, Lz))

file.write("\n")
file.write("Atoms\n")
file.write("\n")

pid = 1
for px, py, pz in zip(x, y, z):
   file.write("%d 1 %.3f %.3f %.3f\n"%(pid, px, py, pz))
   pid += 1

file.write("\n")
file.write("Velocities\n")
file.write("\n")

# write velocities
for i in range(num_particles):
  file.write("%d %.3f %.3f %.3f\n"%(i + 1, vx[i], vy[i], vz[i]))
file.close()
sys.stdout.write('\nESPResSo++ input file written: espressopp/espressopp_lennard_jones.start\n')


###############################################################################
#                                   LAMMPS                                    #
###############################################################################
file = open("lammps/lammps_lennard_jones.start", "w")
file.write("LAMMPS\n")
file.write("\n")
file.write("%5d atoms\n"%num_particles)
file.write("    0 bonds\n")
file.write("    0 angles\n")
file.write("    0 dihedrals\n")
file.write("    0 impropers\n")
file.write("\n")
file.write("    1 atom types\n")
file.write("    0 bond types\n")
file.write("    0 angle types\n")
file.write("    0 dihedral types\n")
file.write("    0 improper types\n")
file.write("\n")

file.write("%.4f %.4f xlo xhi\n" % (0.0, Lx))
file.write("%.4f %.4f ylo yhi\n" % (0.0, Ly))
file.write("%.4f %.4f zlo zhi\n" % (0.0, Lz))

file.write("\n")
file.write("Atoms\n")
file.write("\n")

pid = 1
for px, py, pz in zip(x, y, z):
   file.write("%d 1 %.3f %.3f %.3f\n"%(pid, px, py, pz))
   pid += 1

file.write("\n")
file.write("Velocities\n")
file.write("\n")

# write velocities
for i in range(num_particles):
  file.write("%d %.3f %.3f %.3f\n"%(i + 1, vx[i], vy[i], vz[i]))
file.close()
sys.stdout.write('LAMMPS input file written: lammps/lammps_lennard_jones.start\n')

###############################################################################
#                                  ESPResSo                                   #
###############################################################################
f = open('espresso/espresso_lennard_jones.start', 'w')
f.write('{variable\n')
f.write('\t{time 0.0}\n')
f.write('\t{periodicity 1 1 1}\n')
f.write('\t{box_l %.4f %.4f %.4f}\n' % (Lx, Ly, Lz))
f.write('\t{n_part %d}\n' % (num_particles))
f.write('}\n')

f.write('{tclvariable\n')
f.write('\t{density 0.8442}\n')
f.write('}\n')

f.write('{particles {id pos type v f}\n')
for i in range(num_particles):
  f.write('\t{%d %.3f %.3f %.3f 0 %.3f %.3f %.3f 0.0 0.0 0.0}\n' % (i, x[i], y[i], z[i], vx[i], vy[i], vz[i]))
f.write('}\n')
sys.stdout.write('ESPResSo input file written: espresso/espresso_lennard_jones.start\n')


###############################################################################
#                                  GROMACS                                    #
###############################################################################
f = open('gromacs/conf.gro', 'w')
f.write('Rings\n')
f.write('%d\n' % num_particles)

fmt = '%5d%s%03d%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n'

for i in range(num_particles):
  f.write(fmt % (i + 1, 'RING  A', i + 1, i + 1, x[i], y[i], z[i], vx[i], vy[i], vz[i]))

s = '  %.4f %.4f %.4f\n' % (Lx, Ly, Lz)
f.write(s)
f.close()
sys.stdout.write('GROMACS input file written: gromacs/gromacs_lennard_jones.start (conf.gro)\n\n')
