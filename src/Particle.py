#  Copyright (C) 2016,2017
#      Max Planck Institute for Polymer Research & JGU Mainz
#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""
*******************
espressopp.Particle
*******************


.. function:: espressopp.Particle(pid, storage)

		:param pid:
		:param storage:
		:type pid:
		:type storage:
"""
import _espressopp
import esutil
import pmi
from espressopp import toReal3DFromVector, toInt3DFromVector
import mpi4py.MPI as MPI
from espressopp.Exceptions import ParticleDoesNotExistHere


# Controller Particle:
# * requests are directly forwarded

# Parallel Particle:
# * throw exception if particle does not exist locally
# * otherwise do it

# Parallel Ghost Particle:
# * throw exception if particle does not exist locally
# * can not be written
# * should throw exception if data is not available

# The class _TmpParticle wraps the C++ internal pointer to a particle
# _TmpParticle should not be used, as it might die easily and will
# cause a SegFault when used after it has died.

class ParticleLocal(object):
    """The local particle.

    Throws an exception:
    * when the particle does not exists locally

    TODO: Should throw an exception:
    * when a ghost particle is to be written
    * when data is to be read from a ghost that is not available
    """
    def __init__(self, pid, storage):
      self.pid = pid
      self.storage = storage

    def __getTmp(self):
      return self.storage.lookupRealParticle(self.pid)

        #if tmp is None:
            # TODO: Exception
            # raise ParticleDoesNotExistHere('pid='+str(self.pid)+' rank='+str(pmi.rank) )
        #else:
        #  return tmp

    # Defining __getattr__ will make sure that you can use any
    # property defined in _TmpParticle
    def __getattr__(self, key):
      return getattr(self.__getTmp(), key)

#     def __setattr__(self, key, value):
#         return setattr(self.__getTmp(), key, value)

    # The following properties are modified between Python and C++
    @property
    def f(self): return self.__getTmp().f
    @f.setter
    def f(self, val): self.__getTmp().f = toReal3DFromVector(val)

    @property
    def v(self): return self.__getTmp().v
    @v.setter
    def v(self, val): self.__getTmp().v = toReal3DFromVector(val)

    @property
    def pos(self): return self.__getTmp().pos
    @pos.setter
    def pos(self, val): self.__getTmp().pos = toReal3DFromVector(val)

    @property
    def type(self): return self.__getTmp().type
    @type.setter
    def type(self, val): self.__getTmp().type = val

    @property
    def mass(self): return self.__getTmp().mass
    @mass.setter
    def mass(self, val): self.__getTmp().mass = val

    @property
    def q(self): return self.__getTmp().q
    @q.setter
    def q(self, val): self.__getTmp().q = val

    @property
    def radius(self): return self.__getTmp().radius
    @radius.setter
    def radius(self, val): self.__getTmp().radius = val

    @property
    def fradius(self): return self.__getTmp().fradius
    @radius.setter
    def fradius(self, val): self.__getTmp().fradius = val

    @property
    def vradius(self): return self.__getTmp().vradius
    @radius.setter
    def vradius(self, val): self.__getTmp().vradius = val

    @property
    def imageBox(self): return self.__getTmp().imageBox
    @imageBox.setter
    def imageBox(self, val): self.__getTmp().imageBox = toInt3DFromVector(val)

    @property
    def isGhost(self): return self.__getTmp().isGhost
    @isGhost.setter
    def isGhost(self, val): self.__getTmp().isGhost = val
    
    @property
    def isFixed(self): return self.__getTmp().isFixed
    @isFixed.setter
    def isFixed(self, val): self.__getTmp().isFixed = val

    @property
    def lambda_adr(self): return self.__getTmp().lambda_adr
    @isGhost.setter
    def lambda_adr(self, val): self.__getTmp().lambda_adr = val

    @property
    def drift_f(self): return self.__getTmp().drift_f
    @isGhost.setter
    def drift_f(self, val): self.__getTmp().drift_f = val

    @property
    def lambda_adrd(self): return self.__getTmp().lambda_adrd
    @isGhost.setter
    def lambda_adrd(self, val): self.__getTmp().lambda_adrd = val

    @property
    def extVar(self): return self.__getTmp().extVar
    @radius.setter
    def extVar(self, val): self.__getTmp().extVar = val

    @property
    def state(self): return self.__getTmp().state
    @state.setter
    def state(self, val): self.__getTmp().state = val

    def getLocalData(self, key):
        tmp = self.storage.lookupRealParticle(self.pid)
        if tmp is not None:
            return getattr(tmp, key)
        else:
            return None

    def locateParticle(self):
        tmp = self.storage.lookupRealParticle(self.pid)
        return (tmp is not None)

if pmi.isController:
    class Particle(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.ParticleLocal',
            pmiproperty = [ "id", "storage" ]
            )

        @property
        def node(self):
            value, node = pmi.reduce(pmi.MAXLOC, self, 'locateParticle')
            return node

        def __getattr__(self, key):
            value = pmi.reduce(pmi.MAX, self, 'getLocalData', key)
            return value
