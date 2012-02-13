import espresso
from espresso.tools.convert import lammps

system, integrator = espresso.standard_system.PolymerMelt(num_chains         = 10,            \
                                                          monomers_per_chain = 50,            \
                                                          box                = (50, 50, 50),  \
                                                          bondlen            = 1.0,           \
                                                          rc                 = 2**(1.0/6.0),  \
                                                          dt                 = 0.0001,        \
                                                          epsilon            = 1.0,           \
                                                          sigma              = 1.0,           \
                                                          shift              = 'auto',        \
                                                          temperature        = 1.0)

lammps.write('lammps_data.start', system, writeVelocities=True)
