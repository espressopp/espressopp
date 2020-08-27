import espressopp
import mpi4py.MPI as MPI

#
NCPUs = espressopp.MPI.COMM_WORLD.size
box = [10.0, 10.0, 10.0]

warmup_cutoff = 1.0
skin = 0.3

system = espressopp.System()
system.skin = skin
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
system.comm = MPI.COMM_WORLD

# nodeGrid = espressopp.tools.decomp.nodeGrid(NCPUs, box, warmup_cutoff, skin)
# cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, warmup_cutoff, skin)


system.storage = espressopp.storage.DomainDecomposition(system)
system.storage.addParticle(pid=1, pos=espressopp.Real3D(0, 0, 0))


thermostat = espressopp.integrator.LangevinThermostat(system)
integrator = espressopp.integrator.VelocityVerlet(system)

integrator.addExtension(thermostat)

integrator.run(100)

espressopp.tools.fastwritexyz('out.xyz', system)
