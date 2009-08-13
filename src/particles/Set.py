from espresso import pmi
import espresso.particles.Computer

class SetLocal(object):
    """Any derived object must define cxxobject, which must be a C++
    Set object. 
    """
    def foreach(self, computer):
        """
        apply the computer to each particle in the set (i.e. call
        \"apply\" with the particle identity as
        parameter). 
        
        The return value of foreach is the return value of
        \"finalize\" on the master node.
        
        Example:

            >>> class MyPythonComputer(espresso.particles.PythonComputer) :
            >>>    def __init__(self) :
            >>>        self.count = 0
            >>>    def apply(self, id) :
            >>>        self.count += 1
            >>>    def collect(self) :
            >>>        return mpi.reduce(self.count, x,y : return x+y, pmi.CONTROLLER)
            >>>
            >>> decomposer.foreach(pmi.create(\"MyPythonComputer\"))
            """
        self.cxxclass.foreach(self, computer)

    def __locateParticleLocal(self, particle):
        # to be used by Set.contains
        # don't call locally
        return collectives.locateItem(particle in self)

    def __contains__(self, particle):
        return self.cxxclass.__contains__(self, particle)

# TODO: implement:
#   * __iter__ (via foreach?)

if pmi.IS_CONTROLLER:
    class Set(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(pmicall = ['foreach'])

        def locateParticle(self, particle):
            return pmi.call('espresso.particles.SetLocal.__locateParticleLocal', self, particle)

        def __contains__(self, particle):
            try:
                self.locateParticle(particle)
            except KeyError:
                return false
            return true
            
