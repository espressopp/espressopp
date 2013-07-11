import espresso
#from espresso import Real3D

# create default Lennard Jones (WCA) system with 0 particles and cubic box (L=10)
#system, integrator = espresso.standard_system.LennardJones(0, (10, 10, 10))
system, integrator = espresso.standard_system.LennardJones(0, box=(20, 20, 20), temperature=1.0)

# define a LB grid
lb = espresso.integrator.LatticeBoltzmann(system, 4, 4, 4)
#lb = espresso.integrator.LatticeBoltzmann(system, 4, 4, 4, 1.0, 1.0, 1.0, espresso.Real3D(0., 0., 0.))
#lb = espresso.integrator.LatticeBoltzmann(system, 4, 4, 4, 1.0, 1.0, 1.0, espresso.Real3D(0., 0., 0.), 3, 19)

# create a lattice grid for LB
integrator.addExtension(lb)

integrator.run(10)
