from espresso import pmi
from _espresso import particles_PythonComputer

# This is the abstract base class for a computer implemented in Python
class PythonComputerLocal(particles_PythonComputer):
    """
    a per-particle-computer using Python code. You mainly need to
    provide a function "apply", so that a
    PythonComputerLocal-object can be applied to a particle
    (identity). "apply" is mandatory to override. If any of two
    other functions is provided, namely "prepare" and/or "finalize",
    they are called right before calling "apply" and right after,
    respectively. The "prepare" function receives as input parameter a
    DecomposerLocal-object, the one that the subsequent
    "apply"-calls will operate on. Several calls to "apply"
    might occur before the call to "finalize", but the
    Decomposer-object will not change inbetween.

    For example, the Decomposer.foreach-method will first call
    "prepare" on each node, handing over its local instance. Then on
    each node, it will call "apply" for each particle, and after
    that, call "finalize" once on each node. In the case of
    Decomposer.foreach, the return value of "finalize" is passed back.
    """
    def __init__(self):
        if not hasattr(self, 'cxxinit'):
            particles_PythonComputer.__init__(self)
            self.cxxinit = True
        
    def prepare(self, storage):
        pass

    def apply(self, id): 
        """
        override this function to implement a per-particle computer. This routine will be
        called for each particle, with id being the identity of the particle.
        """
        pass

    def finalize(self): pass

if pmi.IS_CONTROLLER:
    class Computer(object):
        pass
