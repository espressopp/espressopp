#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

###########################################################################
#                                                                         #
#  ESPResSo++ Benchmark Python script for a Lennard Jones System          #
#                                                                         #
###########################################################################

import time
import os
import espressopp
#from espressopp.tools.chronometer import chronometer

nsteps      = 10
isteps      = 100
rc          = pow(2.0, 1.0/6.0)
skin        = 0.4
timestep    = 0.005

# set temperature to None for NVE-simulations
temperature = 1.0

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
print espressopp.Version().info()
print 'Setting up simulation ...'
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.convert.lammps.read('polymer_melt.lammps')
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=4, ydim=4, zdim=4)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)
system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
props = ['id', 'type', 'mass', 'pos']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, espressopp.Real3D(x[i], y[i], z[i])]
  new_particles.append(part)
  if i % 1000 == 0:
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
    new_particles = []
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

# Lennard-Jones with Verlet list
vl      = espressopp.VerletList(system, cutoff = rc)
potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# FENE bonds
fpl = espressopp.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
system.addInteraction(interCosine)

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'density             = ', density
print 'rc                  = ', rc
print 'dt                  = ', integrator.dt
print 'skin                = ', system.skin
print 'temperature         = ', temperature
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# espressopp.tools.decomp.tuneSkin(system, integrator)
#filena = "polymer_melt_%0i.bin" % integrator.step

filena = "polymer_melt_%0i.bin" % integrator.step
timing = True
filename = "polymer_melt_%0i.xyz" % integrator.step
filehdf5 = "polymer_melt_%0i.h5" % integrator.step
filehdf52 = "polymer_melt_%0i_test2.h5" % integrator.step
filehdf5parallel = "polymer_melt_%0i_paralleltest_another.h5" % integrator.step


espressopp.tools.analyse.info(system, integrator)
start_time = time.clock()
#jack = espressopp.io.DumpHDF5nton(system, integrator, filename="hdf5_NtoN_test.h5", unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
for k in range(nsteps):
  integrator.run(isteps)
  espressopp.tools.analyse.info(system, integrator)
  #if k % 2 == 0:
   # jack = espressopp.io.DumpHDF5nton(system, integrator, filename="hdf5_NtoN_test.h5", unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
    #jack = espressopp.io.DumpHDF5parallelckpt(system, integrator, filename=filehdf5parallel, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
    #jack.dump()

#end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)

#os.system("free -m && sync && echo 3 > /proc/sys/vm/drop_caches && free -m")
#os.system("free -m && sync && echo 3 > /proc/sys/vm/drop_caches && free -m")
#os.system("free -m && sync && echo 3 > /proc/sys/vm/drop_caches && free -m")

hdf5_file = True
#hdf5_file = False
txt_file = False
#txt_file = True
iomode_choice = 1 #1 Nto1
#iomode_choice = 2 # 2 NtoN
base_path_files = '/gpfs/fs2/project/zdvresearch/padua/polymer_melt_tests/'

if hdf5_file is True:
  if iomode_choice == 1:
    base_filename = os.path.join(base_path_files, 'N_to_1.h5')
    start_io_time_all = time.clock()
    jack = espressopp.io.HDF5File(system, integrator, filename=base_filename, iomode=iomode_choice, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
    start_io_time_write = time.clock()
    jack.write()
    end_io_time = time.clock()
  elif iomode_choice == 2:
    base_filename = os.path.join(base_path_files, 'N_to_N.h5')
    start_io_time_all = time.clock()
    jack = espressopp.io.HDF5File(system, integrator, filename=base_filename, iomode=iomode_choice, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
    start_io_time_write = time.clock()
    jack.write()
    end_io_time = time.clock()
elif hdf5_file is False and txt_file is True:
  base_filename = os.path.join(base_path_files, 'finiscila.xyz')
  start_io_time_all = time.clock()
  jack = espressopp.io.DumpXYZ(system, integrator, filename=base_filename, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
  start_io_time_write = time.clock()
  jack.dump()
  end_io_time = time.clock()


constr_write_time = "io time all (constr+write): %f\n" % (end_io_time - start_io_time_all)
write_time = "io time (only write): %f\n" % (end_io_time - start_io_time_write)

#print constr_write_time
#print write_time


end_time = time.clock()

espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)
total_time = "total time: %f\n" % (end_time - start_time)

if hdf5_file is True and iomode_choice == 1:
  mode = "hdf5_N_to_1"
elif hdf5_file is True and iomode_choice == 2:
  mode = "hdf5_N_to_N"
elif hdf5_file is False:
  mode = "XYZ"

cores = espressopp.MPI.COMM_WORLD.size

result_file = "result_cores_%i_mode_%s_4.dat" % (cores, mode)
with open(os.path.join(base_path_files, result_file), 'w') as fout:
  fout.write(constr_write_time)
  fout.write(write_time)
  fout.write(total_time)

os.system("sync")
cmd = "rm -f /gpfs/fs2/project/zdvresearch/padua/polymer_melt_tests/*.h5"
os.system(cmd)
#print constr_write_time
#print write_time
#print total_time

