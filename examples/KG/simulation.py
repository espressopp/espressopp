#! /usr/bin/env python

import espresso
import oldespresso
import MPI

#  This script reads an old ESPResSo blockfile and sets up a system
#  with the same data and the same parameters (as far as possible)


def translateParticle(particle, properties):
    """
    Translate partice from blockfile format into new format
    where each entry corresponds to one property
    """

    index = 0    #  used to index old particle

    newparticle = []

    for prop in properties:

        if prop == "id":
 
           id = int(particle[index])
           newparticle.append(id)
           index += 1

        elif prop == "type":
 
           type = int(particle[index])
           newparticle.append(type)
           index += 1

        elif prop == "q":
 
           q = float(particle[index])
           newparticle.append(q)
           index += 1

        elif prop == "pos":

           x = float(particle[index])
           y = float(particle[index+1])
           z = float(particle[index+2])
           newparticle.append(Real3D(x, y, z))
           index += 3

        elif prop == "v":

           vx = float(particle[index])
           vy = float(particle[index+1])
           vz = float(particle[index+2])
           newparticle.append(Real3D(vx, vy, vz))
           index += 3

        elif prop == "f":

           fx = float(particle[index])
           fy = float(particle[index+1])
           fz = float(particle[index+2])
           newparticle.append(Real3D(fx, fy, fz))
           index += 3

        else:

           raise "property %s unknown"%prop
    
    return newparticle 
           
reader = oldespresso.Reader("input.dat")

# Be careful: when asking for values, strings or list of strings are returned

size = reader["variable", "box_l"]
size = (float(size[0]), float(size[1]), float(size[2]))

N    = int(reader["variable", "n_part"][0])

system = espresso.System()

system.rng  = espresso.esutil.RNG()
system.bc   = espresso.bc.OrthorhombicBC(system.rng, size)

skin = reader["variable", "skin"]
system.skin = float(skin[0])

from espresso import Real3D, Int3D

nodeGrid = reader["variable", "node_grid"]
nodeGrid = Int3D(int(nodeGrid[0]), int(nodeGrid[1]), int(nodeGrid[2]))

cellGrid = reader["variable", "cell_grid"]
cellGrid = Int3D(int(cellGrid[0]), int(cellGrid[1]), int(cellGrid[2]))

system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# read in all particles

particles  = reader["particles"]
newparticles = []
properties = particles[0]

for particle in particles[1:]:

    newparticle = translateParticle(particle, properties)

    newparticles.append(newparticle)

system.storage.addParticles(newparticles, *properties)
print newparticles
print properties

# read in all bonds, type is ignored

bonds = reader["bonds"]
bondList = []

for bond in bonds:

    pid1 = int(bond[0])
    
    for partner in bond[1]:
   
       type = int(partner[0])
       pid2 = int(partner[1])

       bondList.append((pid1, pid2))

print "bondList = %s"%bondList

fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bondList)

system.storage.decompose()

# define integrator + timestep

integrator = espresso.integrator.VelocityVerlet(system)

integrator.dt = float(reader["variable", "time_step"][0])

# set all interactions

input_interactions = reader["interactions"]

for iia in input_interactions:

    if iia[0] == "ljforcecap":

       print "ljforcecap ignored"

    elif iia[1] == "FENE":

       print "set FENE : %s"%iia

       type = int(iia[0])
       K    = float(iia[2])
       r0   = float(iia[3])

       # TODO: similiar to LennardJones, but with FixedPairList

    elif iia[2] == "lennard-jones":

       print "set lennard-jones : %s"%iia

       type1  = int(iia[0])
       type2  = int(iia[1])
       sigma  = float(iia[3])
       eps    = float(iia[4])
       cutoff = float(iia[5])
       shift  = float(iia[6])

       potLJ = espresso.interaction.LennardJones(sigma, eps, cutoff = cutoff, shift = shift)

       vl = espresso.VerletList(system, cutoff = cutoff + system.skin)

       interLJ = espresso.interaction.VerletListLennardJones(vl)
       interLJ.setPotential(type1 = type1, type2 = type2, potential = potLJ)
       system.addInteraction(interLJ)

    else:

       print "interaction %s ignored"%iia

temp = espresso.analysis.Temperature(system)
press = espresso.analysis.Pressure(system)

temperature = temp.compute()
kineticEnergy = 0.5 * temperature * N
potentialEnergy = interLJ.computeEnergy()
print 'Start: tot energy = %10.6f pot = %10.6f kin = %10.f temp = %10.6f' \
        %(kineticEnergy + potentialEnergy,
          potentialEnergy, kineticEnergy, temperature)

nsteps = 10

integrator.run(nsteps)
temperature = temp.compute()
p = press.compute()
kineticEnergy = 0.5 * temperature * N
potentialEnergy = interLJ.computeEnergy()
print 'End: tot energy = %10.6f pot = %10.6f kin = %10.6f temp = %f p = %f' \
       %(kineticEnergy + potentialEnergy, 
         potentialEnergy, kineticEnergy, temperature, p)
