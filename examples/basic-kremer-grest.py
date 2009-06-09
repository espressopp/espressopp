#!/usr/bin/env python
#
#  Python example script that does the same as forceloop but with
#  epxorted routines from the C++ classes
#
#  Authors: Jon + Thomas
#
#  Done:   
#     - bindings for particle/All and particel/Set
#     - addtional bindings for integrator/VerlocityVerlet 
#     - ParticleWriter as PythonComputer 
#
#  Furter ToDo's 
#     - Extend all C++ classes in xxx/__init__.py as Local + PMI classes
#     - writing UnitTest in Python

import espresso
import math
import logging

logging.getLogger("Langevin").setLevel(logging.INFO)
logging.getLogger("ForceComputer").setLevel(logging.INFO)
logging.getLogger("Integrator").setLevel(logging.INFO)

from _espresso import particles_Storage as Storage
from _espresso import Real3DProperty
from _espresso import Real3D
from _espresso import bc_PBC as PBC
from _espresso import particles_PythonComputer 
from _espresso import pairs_PythonComputer 
from _espresso import particles_All
from _espresso import pairs_All as All2
from _espresso import pairs_List as List
from _espresso import interaction_LennardJones as LennardJones
from _espresso import interaction_FENE as FENE
from _espresso import integrator_VelocityVerlet as VelocityVerlet
from _espresso import thermostat_Langevin as Langevin
from _espresso import force_ForceComputer as ForceComputer
 
##  Simple class to write position of particles

##  Global Parameters of the simulation

SIZE    = 16.0  # size of the simulation box
NCHAINS = 20     # number of polymers
NBEADS  = 5     # number of particles in one polymer

class ParticleWriter(particles_PythonComputer):

    def __init__(self, _property, storage):
        particles_PythonComputer.__init__(self, storage)
        self.property = _property

    def each(self, pid):
        print("%d %s" % (pid, self.property[pid]))

########################################################################
#
#  Class that allows to write XYZ coordinate files
#
########################################################################

class XYZWriter(particles_PythonComputer):

    """This class writes a XYZ coordinate file for all particles"""

    def __init__(self, storage, size, position, filename):

        particles_PythonComputer.__init__(self, storage)

        self.outfile = open(filename, "w")
        self.position = position
        self.size     = size

    def start(self):

        self.outfile.write("%d\nAtoms\n" % (self.size,))

    def each(self, pid):

        pos = self.position[pid]
        self.outfile.write("1 %f %f %f\n" % (pos[0], pos[1], pos[2]))

    def close(self):

        self.outfile.close()

########################################################################
#
#  Class that write PSF file                         
#
########################################################################

def writePSF(particleStorage, size, bondList, numbonds, filename):

   # Helper class to write particles
   
   class ParticleWriter(particles_PythonComputer):

      def __init__(self, storage, outfile):

          particles_PythonComputer.__init__(self, storage)
          self.outfile = outfile

      def each(self, pid):

          self.outfile.write("%8d" % (pid,))
          self.outfile.write("      1")
          self.outfile.write("%11d" % (pid,))
          self.outfile.write("   ")
          type   = 1
          self.outfile.write("%7ld" % (type,))
          self.outfile.write("0.000000       12.0000           0")
          self.outfile.write("\n")

   # Helper class to write particle pairs (here bonds)

   class BondWriter(pairs_PythonComputer):

      def __init__(self, pairs, outfile):

         pairs_PythonComputer.__init__(self, pairs)
         self.outfile = outfile
         # counter for bonds needed to add CR
         self.nbonds  = 0       

      def each(self, pid1, pid2):

          self.outfile.write("%8d%8d" % (pid1, pid2))
          self.nbonds = self.nbonds + 1
          if self.nbonds == 4:
             self.outfile.write("\n")
             self.nbonds = 0

   #  Open PSF file and write header

   outfile = open(filename, "w")

   outfile.write("PSF\n")
   outfile.write("\n")
   outfile.write("       1 !NTITLE\n")
   outfile.write(" REMARKS ESPResSo++ generated structure x-plor psf file\n")
   outfile.write("\n")
   outfile.write("%8d !NATOM\n" % (size,))

   #  Write the info about the particles

   writeParticles = ParticleWriter(particleStorage, outfile)
   particleStorage.foreach(writeParticles)

   # Write the header for the bonds

   outfile.write("\n")
   outfile.write("%8d !NBONDS\n" % (numbonds,))

   # write the bonds into the file via PairComputer

   writeBonds = BondWriter(bondList, outfile)
   bondList.foreach(writeBonds)

   outfile.close()

import random

def randomWalk(step):

    rsq = 2.0
    while rsq > 1.0:
       dx = 2.0*random.random() - 1.0
       dy = 2.0*random.random() - 1.0
       dz = 2.0*random.random() - 1.0
       rsq = dx*dx + dy*dy + dz*dz

    r  = math.sqrt(rsq)

    return Real3D(step / r * dx, step / r * dy, step / r * dz)

particleStorage = Storage()
position = Real3DProperty(particleStorage)
velocity = Real3DProperty(particleStorage)
force    = Real3DProperty(particleStorage)

pbc = PBC()
pbc.set(SIZE)

# Bond List for bonded interactions

bondList = List(pbc, particleStorage, position)

for chainid in range(NCHAINS):
    pid1 = particleStorage.addParticle()
    pos1 = Real3D(random.random() * SIZE, random.random() * SIZE, random.random() * SIZE)
    position[pid1] = pos1
    velocity[pid1] = Real3D(0.0)
    force[pid1]    = Real3D(0.0)

    for beadid in range(NBEADS-1):

        # create a new position by random walk
        pos2 = pos1 + randomWalk(step=1.0)

        pid2 = particleStorage.addParticle()

        position[pid2] = pos2
        velocity[pid2] = Real3D(0.0)
        force[pid2]    = Real3D(0.0)

        bondList.addPair(pid1, pid2)

        # save pos2 and pid2 for the next step of this loop

        pos1 = pos2
        pid1 = pid2

#  Write all the particles of the storage

#  We could print the particles already here
#  particleStorage.foreach(ParticleWriter(position, particleStorage))

allSet = particles_All(particleStorage)
allSet.foreach(ParticleWriter(position, particleStorage))

# Pair List for nonbonded interactions

allPairs = All2(pbc, allSet, position)

ljint = LennardJones()

fene = FENE()
fene.set(1.0, 0.5, 0.1)
# fene = FENE(r0 = 0.5, K=1.0, rMax = 0.1)

# allPairs.foreach(forcecompute)

#  Setting up the integrator

integrator = VelocityVerlet(allSet, position, velocity, force)

integrator.setTimeStep(0.005)

# Pair(ljint, allPairs).connect(integrator)
# Pair(fene, bondList).connect(integrator)

f1 = ForceComputer(ljint, allPairs)
f2 = ForceComputer(fene, bondList)
f1.connect(f1, integrator);
f2.connect(f2, integrator);

#  Adding a thermostat

thermostat = Langevin(1.0, 0.5)

thermostat.connect(thermostat, integrator)

writePSF(particleStorage, NCHAINS*NBEADS, 
         bondList, NCHAINS*(NBEADS-1), "dump.psf")

xyz = XYZWriter(particleStorage, NCHAINS*NBEADS, position, "dump.xyz")
xyz.start()
allSet.foreach(xyz)

# equilibration

N = 100
for k in range(N):
   # print 'pre-equilibration', k, 'of', N
   sigma = (1.0 * (k+1)) / N
   ljint.set(1.0, sigma, 2.5)
   integrator.run(20)
   xyz.start()
   allSet.foreach(xyz)

xyz.close()

print 'Ready, call now: vmd dump.psf dump.xyz'

# allSet.foreach(ParticleWriter(force, particleStorage))

