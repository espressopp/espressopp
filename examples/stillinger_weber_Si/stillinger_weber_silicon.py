#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

#################################
#                               #
#  ESPResSo++ Python script     #
#                               #
#################################

'''
  This is an example of using Stillinger-Weber nonbonded 
  3-body potential (Si melt).
  The initial configuration is a diamond structure.
  A result is radial distribution function 'rdf.dat' (units - [nm])
'''

import sys
import time
import math
import espresso
import mpi4py.MPI as MPI
from espresso import Real3D
from espresso.tools import decomp

from time import gmtime, strftime, localtime

# initial parameters for Stillinger-Weber potential
# 2 body
_A       = 7.049556277
_B       = 0.6022245584
_p       = 4
_q       = 0

# 3 body
_gamma   = 1.2
_theta   = 109.47
_lambda  = 21

_eps=1.0
_sig=1.0
rc   = 1.8

# general parameters
skin = 0.3
timestep = 0.001

# desired temperature in reduced units
desiredT = 0.15 # it corresponds to 3770 K

# units
sigma_real    = 0.20951 # nm
epsilon_real  = 209     # kJ/mol
mass_real     = 27.9861 # amu

time_real = sigma_real * math.sqrt(mass_real / epsilon_real) # ps

# reading initial file
file = open('ini_conf_diamond.xyz')
lines = file.readlines()
coord = []
vel = []
count = 0
for line in lines:
  if count==0:
    num_particles = int(line.split()[0])
    count=1
  elif count==1:
    Lx, Ly, Lz = map(float, line.split()[0:3])
    count=2
  elif(count==2):
    id1, type1, x, y, z, vx, vy, vz = map(float, line.split()[0:8])
    coord.append( Real3D(x,y,z) )
    vel.append( Real3D(vx,vy,vz) )

density = num_particles / (Lx * Ly * Lz)
box     = (Lx, Ly, Lz)

######################################################################
sys.stdout.write('Setting up simulation ...\n')
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(box, nodeGrid, rc, skin)

system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
######################################################################

# add particles to the system and then decompose
props = ['id', 'type', 'pos', 'v', 'mass']
new_particles = []
type = 0
mass = 1.0
for pid in range(num_particles):
  part = [pid, type, coord[pid], vel[pid], mass]
  new_particles.append(part)
  
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

print ''
print 'number of particles =', num_particles
print 'timestep: %f / %10.9f [ps]' % (timestep, time_real * timestep)
print 'box: ', '(%6.3f, %6.3f, %6.3f) / (%6.3f, %6.3f, %6.3f) [nm]' %  \
        (Lx, Ly, Lz, sigma_real*Lx, sigma_real*Ly, sigma_real*Lz)
print 'density = %.4f' % (density)
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# Stillinger Weber pair potential (similar to other short range pair nonbonded interactions)
vl = espresso.VerletList(system, cutoff=rc)
potSWpair = espresso.interaction.StillingerWeberPairTerm(A=_A, B=_B, p=_p, q=_q, epsilon=_eps, sigma=_sig, cutoff=rc)
interSWpair = espresso.interaction.VerletListStillingerWeberPairTerm(vl)
interSWpair.setPotential(type1=0, type2=0, potential=potSWpair)
system.addInteraction(interSWpair)

# Stillinger Weber triple potential (3-body nonbonded interaction)
vl3 = espresso.VerletListTriple(system, cutoff=rc)
potSWtriple = espresso.interaction.StillingerWeberTripleTerm(gamma=_gamma, theta0=_theta, lmbd=_lambda, epsilon=_eps, sigma=_sig, cutoff=rc)
interSWtriple = espresso.interaction.VerletListStillingerWeberTripleTerm(system, vl3)
interSWtriple.setPotential(type1=0, type2=0, type3=0, potential=potSWtriple)
system.addInteraction(interSWtriple)

#######################################################################
#                               IMPORTANT                             #
#  Currently in order to use many particle types one have to          #
#  introduce whole set of possible triples with appropriate           #
#  parameters. For example, we have two types, then:                  #
# >> interSWtriple.setPotential(type1=0, \                            #
#                               type2=0, \                            #
#                               type3=0, \                            #
#                               potential=potSWtriple000)             #
# >> interSWtriple.setPotential(type1=1, \                            #
#                               type2=0, \                            #
#                               type3=0, \                            #
#                               potential=potSWtriple100)             #
# >> interSWtriple.setPotential(type1=1, \                            #
#                               type2=1, \                            #
#                               type3=0, \                            #
#                               potential=potSWtriple110)             #
# >> interSWtriple.setPotential(type1=1, \                            #
#                               type2=1, \                            #
#                               type3=1, \                            #
#                               potential=potSWtriple111)             #
# >> interSWtriple.setPotential(type1=0, \                            #
#                               type2=1, \                            #
#                               type3=0, \                            #
#                               potential=potSWtriple010)             #
# >> interSWtriple.setPotential(type1=0, \                            #
#                               type2=1, \                            #
#                               type3=1, \                            #
#                               potential=potSWtriple011)             #
# >> interSWtriple.setPotential(type1=0, \                            #
#                               type2=0, \                            #
#                               type3=1, \                            #
#                               potential=potSWtriple001)             #
# >> interSWtriple.setPotential(type1=1, \                            #
#                               type2=0, \                            #
#                               type3=1, \                            #
#                               potential=potSWtriple101)             #
#######################################################################

# setup integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = timestep

# Langevin thermostat
lT             = espresso.integrator.LangevinThermostat(system)
lT.gamma       = 1.0
lT.temperature = desiredT
integrator.addExtension(lT)

print 'Equilibration.'
in_time = time.time()
espresso.tools.analyse.info(system, integrator)
for i in range (10):
  integrator.run(1000)
  espresso.tools.info(system, integrator)

##############################################
# calculate the radial distribution function #
##############################################

nprints = 5
ncycles = 2
nruns = 100

print '\n Measurements!'

rdf_array_total = []

for i in range (nprints):
  for j in range(ncycles):
    integrator.run(nruns)
    
    # calculate radial distribution function
    distance_step = 800
    rdf = espresso.analysis.RadialDistrF(system)
    rdf_array = rdf.compute(distance_step);
    
    for i in range( len(rdf_array) ):
      if(i>=len(rdf_array_total)):
        rdf_array_total.append(rdf_array[i])
      else:
        rdf_array_total[i] += rdf_array[i]
    
  espresso.tools.info(system, integrator)

fin_time = time.time()
print '\ngeneral running time: ', (fin_time-in_time)

# printing radial distribution function
nameFile = 'rdf.dat'
resFile = open (nameFile, 'w')
fmt = ' %12.6f %12.6f\n'

dr = Lx / 2. / float(distance_step) * sigma_real;

for i in range( len(rdf_array_total) ):
  resFile.write(fmt % ( (i+0.5)*dr, rdf_array_total[i] / float(ncycles*nprints) ))
resFile.close()
