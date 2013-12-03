"""
***********************************
**espresso.standard_system.KGMelt**
***********************************

"""
import espresso

class KGMelt:
  def __init__(self, num_chains, chain_len):
    self._num_chains    = num_chains
    self._chain_len     = chain_len
    self._num_particles = num_chains * chain_len
    self._density       = 0.8449
    self._L             = pow(self._num_particles / self._density, 1.0/3.0)
    self._box           = (self._L, self._L, self._L)
    self._system        = espresso.System()

  
