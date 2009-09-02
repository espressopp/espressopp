#######################################################
#
# Testing the decomposer
#

from espresso import pmi
from espresso import boostmpi as mpi

if __name__ == 'espresso.pmi':
    from espresso.particles import PythonComputerLocal
    # tag particles in a small sphere around (0,0,0) as demonstration
    class ParticleTesterLocal(PythonComputerLocal):
        def __init__(self, _position, _tag):
            PythonComputerLocal.__init__(self)
            self.position = _position
            self.tag = _tag

        def prepare(self, storage):
            self.position.checkFitsTo(storage)
            self.tag.checkFitsTo(storage)

        def apply(self, pid):
            pos = self.position[pid]
            self.tag[pid] = pos*pos < 0.5
            return True

    # write out tagged particles
    class ParticleWriterLocal(PythonComputerLocal):
        def __init__(self, _property, _tag):
            PythonComputerLocal.__init__(self)
            self.property = _property
            self.tag = _tag

        def prepare(self, storage):
            self.property.checkFitsTo(storage)
            self.tag.checkFitsTo(storage)
            self.total = 0
            self.sphere = 0

        def apply(self, pid) :
            self.total += 1
            if self.tag[pid]:
                print("%d %d %s" % (mpi.rank, pid, self.property[pid]))
                self.sphere +=1
            return True

        def _getCount(self):
            return self.sphere, self.total

        def getCount(self):
            return pmi.reduce("lambda x,y: (x[0]+y[0], x[1]+y[1])", 
                              self._getCount)

else:

    pmi.execfile_(__file__)

    from espresso import Real3D, Real3DProperty, IntegerProperty
    from espresso import bc
    from espresso.storage import SingleNode
    import random

    storage = SingleNode(bc.PeriodicBC(length=1), mpi.size-1)
    pos = storage.getPositionProperty()

    for count in range(0,100):
        storage.addParticle(count)
        pos[count] = Real3D(random.random(), random.random(), random.random())

    # add property a posteriori
    tag = IntegerProperty(storage)

    #tag particles
    storage.foreach(pmi.create("ParticleTesterLocal", pos.pmiobject, tag.pmiobject))

    # and print tagged ones
    writer=pmi.create("ParticleWriterLocal", pos.pmiobject, tag.pmiobject)
    storage.foreach(writer)
    count=writer.getCount()
    
    print("printed %d out of %d particles" % count)
