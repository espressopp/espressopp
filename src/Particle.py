#  Copyright (C) 2017,2018
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2016,2018
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2016
#      Max Planck Institute for Polymer Research & JGU Mainz
#  Copyright (C) 2012-2015
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008-2011
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
**************************************
espressopp.Particle
**************************************

The Particle class. Particles are used to model atoms, coarse-grained beads, etc. and are the central part of all simulations. They are stored in the storage and their time evolution can be modeled using an integrator. Importantly, particles have various properties which the user and other modules can make use of and which can be accessed. They are listed below.

.. class:: espressopp.Particle(pid, storage)

        The particle class.

        :param pid: the particle id
        :param storage: the storage object
        :type pid: int
        :type storage: storage

.. py:data:: Real3D espressopp.Particle.pos

        position

.. py:data:: Real3D espressopp.Particle.v

        velocity

.. py:data:: Real3D espressopp.Particle.f

        force

.. py:data:: Real3D espressopp.Particle.modepos

        normal mode coordinate (position in normal mode space)

.. py:data:: Real3D espressopp.Particle.modemom

        normal mode momentum (momentum in normal mode space)

.. py:data:: Real3D espressopp.Particle.fm

        normal mode force (force in normal mode space)

.. py:data:: int espressopp.Particle.type

        particle type

.. py:data:: int espressopp.Particle.res_id

        molecule id (eg. chain id)

.. py:data:: int espressopp.Particle.pib

        path integral bead number (Trotter number)

.. py:data:: real espressopp.Particle.q

        charge

.. py:data:: real espressopp.Particle.mass

        mass

.. py:data:: real espressopp.Particle.varmass

        variable mass (for path integral-based adaptive resolution simulations)

.. py:data:: real espressopp.Particle.radius

        particle radius

.. py:data:: real espressopp.Particle.vradius

        radial velocity

.. py:data:: real espressopp.Particle.fradius

        radial force

.. py:data:: real espressopp.Particle.lambda_adr

        particle's resolution parameter (used in adaptive resolution simulations)

.. py:data:: real espressopp.Particle.lambda_adrd

        particle's gradient of the resolution function in the direction of resolution change (used in adaptive resolution simulations)

.. py:data:: bool espressopp.Particle.isGhost

        boolean flag to indicate whether particle is ghost particle or not

.. py:data:: Int3D espressopp.Particle.imageBox

        particle's image box

.. py:data:: real espressopp.Particle.extVar

        auxiliary variable associated with the particle (used in generalized Langevin friction)

.. py:data:: real espressopp.Particle.drift_f

        particle's drift force in the direction of resolution change (used in adaptive resolution simulations)

.. py:data:: real espressopp.Particle.state

        particle state (used in AssociationReaction)

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
    def fm(self): return self.__getTmp().fm
    @fm.setter
    def fm(self, val): self.__getTmp().fm = toReal3DFromVector(val)

    @property
    def v(self): return self.__getTmp().v
    @v.setter
    def v(self, val): self.__getTmp().v = toReal3DFromVector(val)

    @property
    def pos(self): return self.__getTmp().pos
    @pos.setter
    def pos(self, val): self.__getTmp().pos = toReal3DFromVector(val)

    @property
    def modepos(self): return self.__getTmp().modepos
    @modepos.setter
    def modepos(self, val): self.__getTmp().modepos = toReal3DFromVector(val)

    @property
    def modemom(self): return self.__getTmp().modemom
    @modemom.setter
    def modemom(self, val): self.__getTmp().modemom = toReal3DFromVector(val)

    @property
    def type(self): return self.__getTmp().type
    @type.setter
    def type(self, val): self.__getTmp().type = val

    @property
    def pib(self): return self.__getTmp().pib
    @pib.setter
    def pib(self, val): self.__getTmp().pib = val

    @property
    def mass(self): return self.__getTmp().mass
    @mass.setter
    def mass(self, val): self.__getTmp().mass = val

    @property
    def varmass(self): return self.__getTmp().varmass
    @varmass.setter
    def varmass(self, val): self.__getTmp().varmass = val

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
    @fradius.setter
    def fradius(self, val): self.__getTmp().fradius = val

    @property
    def vradius(self): return self.__getTmp().vradius
    @vradius.setter
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
    def lambda_adr(self): return self.__getTmp().lambda_adr
    @lambda_adr.setter
    def lambda_adr(self, val): self.__getTmp().lambda_adr = val

    @property
    def drift_f(self): return self.__getTmp().drift_f
    @drift_f.setter
    def drift_f(self, val): self.__getTmp().drift_f = val

    @property
    def lambda_adrd(self): return self.__getTmp().lambda_adrd
    @lambda_adrd.setter
    def lambda_adrd(self, val): self.__getTmp().lambda_adrd = val

    @property
    def extVar(self): return self.__getTmp().extVar
    @extVar.setter
    def extVar(self, val): self.__getTmp().extVar = val

    @property
    def state(self): return self.__getTmp().state
    @state.setter
    def state(self, val): self.__getTmp().state = val

    @property
    def res_id(self): return self.__getTmp().res_id
    @res_id.setter
    def res_id(self, val): self.__getTmp().res_id = val

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
