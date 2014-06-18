"""lennard_jones_93_wall.py

This program prepares a system of identical Lennard-Jones particles in a slab
geometry, in the x direction. The slab is imposed by Lennard-Jones 9-3 walls.

The system is warmed up with a velocity rescaling potential and capped
potentials. Then, a NVE run is performed during which the density along the x
direction is sampled.

The density is saved in the file 'lj_93wall_density.txt'

"""

import espresso
import mpi4py.MPI as MPI
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('num', type=int, help='Number of LJ units')
parser.add_argument('--dt', type=float, help='Timestep.', default=0.01)
parser.add_argument('--warmup-steps', type=int, help='Number of steps in warmup loop.', default=100)
parser.add_argument('--warmup-loops', type=int, help='Number of iterations in warmup loop.', default=20)
parser.add_argument('--steps', type=int, help='Number of steps in main loop.', default=100)
parser.add_argument('--loops', type=int, help='Number of iterations in main loop.', default=40)
parser.add_argument('--density', type=float, help='Number density.', default=0.8)
args = parser.parse_args()

def get_velocity(system, n):
    """Obtain total velocity of a espresso system."""
    total_v = espresso.Real3D(0.)
    total_m = 0.
    for i in range(n):
        p = system.storage.getParticle(i)
        total_v += p.v*p.mass
        total_m += p.mass
    return total_v/total_m

def reset_velocity(system, n):
    """Reset the total velocity of a espresso system."""
    excess_v = get_velocity(system, n)
    for i in range(n):
        v = system.storage.getParticle(i).v
        system.storage.modifyParticle(i, 'v', v-excess_v)

# LJ settins
sigma = 1.0
epsilon=1.0
caprad_LJ=0.85
rc = pow(2., 1./6.)

# General settings
skin = 0.3

# Wall settings
wall_sigma = 4.
wall_cutoff = wall_sigma*3.**(1./6.)
wall_r0  = 2.

# Box definition
L = pow(args.num/args.density, 1./3.)
box = (L+2*wall_r0, L, L)

# Initialize the espresso system
system         = espresso.System()
system.rng     = espresso.esutil.RNG()
system.bc      = espresso.bc.SlabBC(system.rng, box)
system.skin    = skin
nodeGrid       = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# Initialize particles
particles_list = []
for i in range(args.num):
    pos = system.bc.getRandomPos()
    # Restrict box size for particles in dimension 0
    pos[0] = (wall_sigma+wall_r0)+pos[0]/L*(L-2.*(wall_sigma+wall_r0))
    particles_list.append([i, pos])
system.storage.addParticles(particles_list, 'id', 'pos')

# Define capped LJ potential
verletList = espresso.VerletList(system, cutoff=rc)
LJCapped    = espresso.interaction.VerletListLennardJonesCapped(verletList)
LJCapped.setPotential(type1=0, type2=0, potential=espresso.interaction.LennardJonesCapped(epsilon=epsilon, sigma=sigma, cutoff=rc, caprad=caprad_LJ))
system.addInteraction(LJCapped)

# Define integrator and StochasticVelocityRescaling thermostat
integrator     = espresso.integrator.VelocityVerlet(system)
thermostat = espresso.integrator.StochasticVelocityRescaling(system)
thermostat.temperature = 1.0
integrator.addExtension(thermostat)

system.storage.decompose()

integrator.dt = args.dt

LJ93 = espresso.interaction.LennardJones93Wall()
LJ93.setParams(0, 6., 1., wall_cutoff, wall_r0)
SPLJ93 = espresso.interaction.SingleParticleLennardJones93Wall(system, LJ93)
system.addInteraction(SPLJ93)

espresso.tools.analyse.info(system, integrator, per_atom=True)

# initial value of epsilon for capped warmup
epsilon_start = 0.1
# Run system with capped potentials, thermostat and increasing LJ epsilon
for k in range(args.warmup_loops):
    LJCapped.setPotential(0,0,espresso.interaction.LennardJonesCapped(epsilon_start + (epsilon-epsilon_start)*k*1.0/(args.warmup_loops-1), sigma, rc, caprad=caprad_LJ))
    integrator.run(args.warmup_steps)
    espresso.tools.analyse.info(system, integrator, per_atom=True)

# Remove LJ93 interaction (to preserve later the ordering of energies)
system.removeInteraction(1)
# Remove LJ Capped potential
system.removeInteraction(0)

# Add non-capped LJ potential
LJ = espresso.interaction.VerletListLennardJones(verletList)
LJ.setPotential(0, 0, espresso.interaction.LennardJones(epsilon, sigma, rc))
system.addInteraction(LJ)

# Re-introduce LJ93
system.addInteraction(SPLJ93)

# Run system with non-capped potentials, thermostat and fixed LJ epsilon
for k in range(args.warmup_loops):
    integrator.run(args.warmup_steps)
    espresso.tools.analyse.info(system, integrator, per_atom=True)

# Disconnect the thermostat
thermostat.disconnect()
# Reset the total velocity of the system
reset_velocity(system, args.num)

xd = espresso.analysis.XDensity(system)
xd_data = []
# Run system with non-capped potentials, no thermostat, fixed LJ epsilon
for k in range(args.loops):
    integrator.run(args.steps)
    espresso.tools.analyse.info(system, integrator, per_atom=True)
    xd_data.append(xd.compute(100))

xd_data = np.array(xd_data)
np.savetxt('lj_93wall_density.txt', xd_data)
