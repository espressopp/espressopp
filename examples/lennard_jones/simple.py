import espressopp

box = (10, 10, 10)
warmup_cutoff = pow(2.0, 1.0/6.0)
skin = 0.3

system  = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.rng.seed(90)
system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin = skin
NCPUs = espressopp.MPI.COMM_WORLD.size
nodeGrid = espressopp.tools.decomp.nodeGrid(NCPUs,box,warmup_cutoff, skin)
cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, warmup_cutoff, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)


integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = 0.001
thermostat = espressopp.integrator.LangevinThermostat(system)
thermostat.gamma = 1.0
thermostat.temperatur = 0.5
integrator.addExtension(thermostat)


pos = system.bc.getRandomPos()
print(system.storage)
system.storage.addParticle(1, pos)
