#######################################################
#
# Testing the decomposer
#

from espresso import pmi
from espresso import boostmpi as mpi

if __name__ == 'espresso.pmi':
    from espresso import pmi
    from espresso.particles import PythonComputerLocal
    import sys

    # tag particles in a small sphere around (0,0,0) as demonstration
    class ParticleTesterLocal(PythonComputerLocal):
        def __init__(self, _position, _tag):
            PythonComputerLocal.__init__(self, _position.decomposer)
            self.position = _position
            self.tag = _tag

        def prepare(self, storage):
            storage.checkProperty(self.position.cxxobject)
            storage.checkProperty(self.tag.cxxobject)

        def apply(self, pid):
            pos = self.position[pid]
            self.tag[pid] = pos*pos < 0.5

    # write out tagged particles
    class ParticleWriterLocal(PythonComputerLocal):
        def __init__(self, _property, _tag):
            PythonComputerLocal.__init__(self, _property.decomposer)
            self.property = _property
            self.tag = _tag

        def prepare(self, storage):
            storage.checkProperty(self.property.cxxobject)
            storage.checkProperty(self.tag.cxxobject)
            self.total = 0
            self.sphere = 0

        def apply(self, pid) :
            self.total += 1
            if self.tag[pid]:
                print("%d %d %s" % (mpi.rank, pid, self.property[pid]))
                self.sphere +=1

        def _getCount(self):
            return self.sphere, self.total

        def getCount(self):
            return pmi.reduce("lambda x,y: (x[0]+y[0], x[1]+y[1])", 
                              self._getCount)

else:

    pmi.execfile_(__file__)

    from espresso import Real3D, Real3DProperty
    from espresso.decomposition import SingleNode
    import random

    decomposer = SingleNode(mpi.size-1)
    pos = decomposer.createProperty("Real3D")

    for count in range(0,100):
        decomposer.addParticle(count)
        pos[count] = Real3D(random.random(), random.random(), random.random())

    # add property a posteriori
    tag = decomposer.createProperty("Integer")

    #tag particles
    decomposer.foreach(pmi.create("ParticleTesterLocal", pos.pmiobject, tag.pmiobject))

    # and print tagged ones
    writer=pmi.create("ParticleWriterLocal", pos.pmiobject, tag.pmiobject)
    decomposer.foreach(writer)
    count=writer.getCount()
    
    print("printed %d out of %d particles" % count)
