from espresso import pmi
from _espresso import particles_PythonComputer as _PythonComputer

__all__ = [ "ComputerLocalBase"]

class ComputerLocalBase(_PythonComputer) :
    """
    a per-particle-computer using Python code. You mainly need to
    provide a function "__apply__", so that a
    ComputerLocalBase-object can be applied to a particle
    (identity). "__apply__" is mandatory to override. If any of two
    other functions is provided, namely "prepare" and/or "finalize",
    they are called right before calling "__apply__" and right after,
    respectively. The "prepare" function receives as input parameter a
    DecomposerLocal-object, the one that the subsequent
    "__apply__"-calls will operate on. Several calls to "__apply__"
    might occur before the call to "finalize", but the
    Decomposer-object will not change inbetween.

    For example, the Decomposer.foreach-method will first call
    "prepare" on each node, handing over its local instance. Then on
    each node, it will call "__apply__" for each particle, and after
    that, call "finalize" once on each node. In the case of
    Decomposer.foreach, the return value of "finalize" is passed back.
    """
    def __init__(self):
        _PythonComputer.__init__(self)

    def __apply__(self, id):
        """
        override this function to implement a per-particle computer. This routine will be
        called for each particle, with id being the identity of the particle.
        """
        raise RuntimeError("PythonComputer invoked without __apply__ being implemented.")
