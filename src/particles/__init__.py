from espresso import pmi
from _espresso import particles_PythonComputer

__all__ = [ "PythonComputerLocal"]

class PythonComputerLocal(particles_PythonComputer) :
    """
    a per-particle-computer using Python code. You mainly need to
    provide two functions: "each", which is called for each particle
    and is mandatory to override, and "collect", which you only need
    to provide if you want something to happen after the function
    "each" has been called for each particle. The result of "collect"
    on the master node will be available for further processing, see
    the documentation of decomposition.Decomposer.foreach for further
    details.
    """
    def __init__(self) :
        particles_PythonComputer.__init__(self)

    def each(self, id) :
        """
        override this function to implement a per-particle computer. This routine will be
        called for each particle, with "id" being the identity of the particle.
        """
        raise RuntimeError("PythonComputer invoked without \"each\" being implemented.")
