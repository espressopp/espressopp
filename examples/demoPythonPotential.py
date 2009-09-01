import espresso.potential

if __name__ == 'espresso.pmi':
    import math

    class HarmonicLocal(espresso.potential.PythonCentralPotentialLocal):
        def __init__(self, D=1.0, r0=0.0):
            espresso.potential.PythonCentralPotentialLocal.__init__(self)
            self.D = D
            self.r0 = r0

        def computeEnergySqr(self, distSqr):
            return 0.5 * self.D * (math.sqrt(distSqr) - self.r0)**2

        def computeForce(self, *args):
            dist = espresso.toReal3D(*args)
            return self.D * (abs(dist) - self.r0) * dist

else:
    class Harmonic(espresso.potential.Potential):
        __metaclass__ = espresso.pmi.Proxy
        pmiproxydefs = espresso.potential.Potential.pmiproxydefs
        pmiproxydefs['pmiproperty'] = [ 'D', 'r0' ]
        pmiproxydefs['cls'] = 'HarmonicLocal'

        def __init__(self, D=1.0, r0=0.0):
            self.pmiinit('HarmonicLocal', D, r0)

    from espresso import pmi
    pmi.execfile_(__file__)

    from espresso import Real3D, Real3DProperty, RealProperty
    from espresso import boostmpi as mpi
    import espresso.storage
    import espresso.bc
    import espresso.pairs
    import random
    
    storage = espresso.storage.SingleNode(espresso.bc.PeriodicBC(), mpi.size-1)
    pos = Real3DProperty(storage)
    force = Real3DProperty(storage)
    energy = RealProperty(storage)
    bc = espresso.bc.PeriodicBC(10.0)

    for pid in range(10):
        storage.addParticle(pid)
        pos[pid] = bc.getRandomPos()
        force[pid] = (0.0, 0.0, 0.0)
        energy[pid] = 0.0

    potential = Harmonic()

    allpairs = espresso.pairs.All(set=storage, posProperty=pos, bc=bc)
    harmint = espresso.potential.Interaction(pairs=allpairs, potential=potential)

    harmint.addForces(force)
    harmint.addEnergies(energy)

    for pid in range(10):
        print('Particle %d: x=%s f=%s e=%f' % (pid, pos[pid], force[pid], energy[pid]))
