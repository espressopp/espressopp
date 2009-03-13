# - setup simulation system
#    - set system geometry
#    - set integrator
#    - set particles
#    - set interactions
#      - set bonded interactions
#      - set nonbonded interactions
#    - set decomposition type
# - setup topology
# - integrate
# - analysis

#####
system = ES.NVTSystem(T=1.0, langevin_gamma=0.5, timestep=0.001)

## SET GEOMETRY
system.geometry = ES.PBC(length=10)

## SET CHAINS
# default: first_pos=random
system.generateChains(number=100, length=64, 
                      interaction=ES.FENEInteraction(...))

# default: tuples=(suggested by tupleInteraction?) 
system.addInteraction(tupleInteraction=ES.LJPairInteraction(...))

## INTEGRATE
for sweeps=1..1000:
    system.integrate(steps=100)
    # - write system to disc
    # - VMD coupling
    # - analysis

