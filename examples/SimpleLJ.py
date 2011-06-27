from espresso                import System, Real3D, VerletList
from espresso.bc             import OrthorhombicBC
from espresso.tools          import timers, replicate
from espresso.esutil         import RNG
from espresso.storage        import DomainDecomposition
from espresso.analysis       import Temperature
from espresso.integrator     import VelocityVerlet, Langevin
from espresso.interaction    import LennardJones, VerletListLennardJones
from espresso.tools.decomp   import nodeGrid, cellGrid
from espresso.tools.init_cfg import lattice, velocities     

import MPI
import logging
import random
import time

random.seed(12345)
rand1                  = random.random()

n_steps                = 10000
n_analyze              = 100
rho                    = 0.8442
num_particles          = 20**3 
skin                   = 0.3   
epsilon                = 1.0
sigma                  = 1.0
shift                  = False
dt                     = 0.01
temperature            = 0.6
gamma                  = 1.0

x, y, z, Lx, Ly, Lz    = lattice.create(num_particles, rho, perfect=False, RNG=rand1)
vx, vy, vz             = velocities.gaussian(T=temperature, N=num_particles)
box                    = (Lx, Ly, Lz)
rc                     = 2.5
s                      = System()
s.rng                  = RNG()
s.bc                   = OrthorhombicBC(s.rng,box)
s.skin                 = skin
nG                     = nodeGrid(MPI.COMM_WORLD.size)
cG                     = cellGrid(box, nG, rc, skin)
s.storage              = DomainDecomposition(s,nG,cG)

vl                     = VerletList(s,cutoff=rc+skin)
potLJ                  = LennardJones(epsilon, sigma, rc, shift)
interLJ                = VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
s.addInteraction(interLJ)

integrator             = VelocityVerlet(s)
integrator.dt          = dt

thermostat             = Langevin(s)
thermostat.gamma       = gamma
thermostat.temperature = temperature
integrator.langevin    = thermostat

props = ['id', 'type', 'mass', 'pos', 'v']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, Real3D(x[i], y[i], z[i]), Real3D(vx[i], vy[i], vz[i])]
  new_particles.append(part)
s.storage.addParticles(new_particles, *props)
s.storage.decompose()

t=0.0
print "t= %5.2f T= %5.2f E_pot= %8.5f" % (t, Temperature(s).compute(), interLJ.computeEnergy()/num_particles)
for k in range(0,n_steps/n_analyze):
  integrator.run(n_analyze)
  t+=dt*n_analyze
  print "t= %5.2f T= %5.2f E_pot= %8.5f" % (t, Temperature(s).compute(), interLJ.computeEnergy()/num_particles)

