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

from espresso.tools.init_cfg import lattice


# cubic lattice with user-defined values of N and rho
# num_particles should be a perfect cube (e.g. 25**3=15625, 32**3=32768)

N = 10

num_particles = N**3
rho = 0.8442
x, y, z, Lx, Ly, Lz = lattice.create(num_particles, rho, perfect=False)

file = open("data.lj", "w")

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

file.write("%g %g xlo xhi\n"%(0.0, Lx))
file.write("%g %g ylo yhi\n"%(0.0, Ly))
file.write("%g %g zlo zhi\n"%(0.0, Lz))

file.write("\n")
file.write("Atoms \n")
file.write("\n")

pid = 1
for px, py, pz in zip(x, y, z):
   file.write("%d 1 %g %g %g\n"%(pid, px, py, pz))
   pid += 1

file.close()
