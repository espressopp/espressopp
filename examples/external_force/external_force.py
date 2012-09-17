import espresso
import math

# specify number of particles in cubic box
N = 4
num_particles = pow(N, 3)
# create particle positions in 3d lattice
x, y, z, Lx, Ly, Lz = espresso.tools.lattice.create(N=num_particles, rho=0.5)

# setup default system (has no interactions or particles yet)
# system, integrator = espresso.standard_system.Default(box=(4*Lx, Ly, Lz))

# setup LennardJones system (with 0 particles)
system, integrator = espresso.standard_system.LennardJones(0, box=(4*Lx, Ly, Lz))

#add particles to the system
particles     = []
property_list = ['id', 'pos']
for k in range(num_particles):
  particles.append([k, espresso.Real3D(2*Lx+x[k], y[k], z[k])])
  # add particles in chunks of 1000 (this is faster than doing it one by one)
  if k % 1000 == 0:
    system.storage.addParticles(particles, *property_list)
    system.storage.decompose()
    particles = []
system.storage.addParticles(particles, *property_list)
system.storage.decompose()

# specify particles that will be affected by external field
ef_particleGroup  = espresso.ParticleGroup(system.storage)
for k in range(num_particles):
  particles.append([k, espresso.Real3D(x[k], y[k], z[k])])
  if k%N>0 and k%N<N-1 :
    ef_particleGroup.add(k)

# define external force
ext_force = espresso.integrator.ExtForce(system, espresso.Real3D(1,0,0), ef_particleGroup)
# add external force to the integrator
integrator.addExtension(ext_force)

sock = espresso.tools.vmd.connect(system)
for i in range(1000):
  # make 10 Velocity-Verlet integration steps
  integrator.run(1)
  # update postions in VMD
  espresso.tools.vmd.imd_positions(system, sock)
  # modify external force field
  ef = 5.0 * math.cos(0.01*i)
  ext_force.setExtForce(espresso.Real3D(ef, 0,0))

