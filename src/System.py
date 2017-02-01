#  Copyright (C) 2015
#      Jakub Krajniak (jkrajniak at gmail.com)
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
*****************
espressopp.System
*****************

The main purpose of this class is to store pointers to some
important other classes and thus make them available to C++.
In a way the System class can be viewed as a container for
system wide global variables.
If you need to run more than one system at the same time you
can combine several systems with the help of the Multisystem
class.

**In detail the System class holds pointers to:**

* the `storage` (e.g. DomainDecomposition)
* the boundary conditions `bc` for the system (e.g. OrthorhombicBC)
* a random number generator `rng` which is for example used by a thermostat
* the `skin` which is needed for the Verlet lists and the cell grid
* a list of short range interactions that apply to the system these
  interactions are added with the `addInteraction()` method of the System

Example (not complete):

>>> LJSystem      = espressopp.System()
>>> LJSystem.bc   = espressopp.bc.OrthorhombicBC(rng, boxsize)
>>> LJSystem.rng
>>> LJSystem.skin = 0.4
>>> LJSystem.addInteraction(interLJ)


.. function:: espressopp.System()


.. function:: espressopp.System.addInteraction(interaction, name)

		:param interaction: 
		:type interaction: 
		:param name: The optional name of the interaction.
		:type name: string
		:rtype: bool

.. function:: espressopp.System.getInteraction(number)

		:param number: 
		:type number: 
		:rtype: 

.. function:: espressopp.System.getNumberOfInteractions()

		:rtype: 

.. function:: espressopp.System.removeInteraction(number)

		:param number: 
		:type number: 
		:rtype: 

.. function:: espressopp.System.removeInteractionByName(self, name)

		:param name: The name of the interaction to remove.
		:type name: str

.. function:: espressopp.System.getAllInteractions()

		:rtype: The dictionary with name as a key and Interaction object.

.. function:: espressopp.System.scaleVolume(\*args)

		:param \*args: 
		:type \*args: 
		:rtype: 

.. function:: espressopp.System.setTrace(switch)

		:param switch: 
		:type switch: 
"""

from espressopp import pmi, Real3D, toReal3DFromVector
from espressopp.esutil import cxxinit
from espressopp.Exceptions import Error

import _espressopp
import mpi4py.MPI as MPI


class SystemLocal(_espressopp.System):
    def __init__(self):

        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, _espressopp.System, pmi._PMIComm.getMPIsubcomm())
            else :
                pass
        else :
            cxxinit(self, _espressopp.System, pmi._MPIcomm)

        self._integrator = None
        self._interaction2id = {}
        self._interaction_pid = 0

    @property
    def integrator(self):
        return self._integrator

    @integrator.setter
    def integrator(self, _integrator):
        self._integrator = _integrator

    def addInteraction(self, interaction, name=None):

        if pmi.workerIsActive():
            ret_val = self.cxxclass.addInteraction(self, interaction)
            if name is not None:
                if name in self._interaction2id:
                    raise RuntimeError('Interaction with name {} already defined.'.format(name))
                self._interaction2id[name] = self._interaction_pid
            self._interaction_pid += 1
            return ret_val

    def removeInteraction(self, number):

        if pmi.workerIsActive():
            self.cxxclass.removeInteraction(self, number)

    def removeInteractionByName(self, name):
        if pmi.workerIsActive():
            if name not in self._interaction2id:
                raise RuntimeError('Interaction {} not found'.format(name))
            interaction_id = self._interaction2id[name]
            self.cxxclass.removeInteraction(self, interaction_id)
            self._interaction2id = {
                k: v if v < interaction_id else v - 1
                for k, v in self._interaction2id.iteritems()
                }
            self._interaction_pid = max(self._interaction2id.values()) + 1

    def getAllInteractions(self):
        if pmi.workerIsActive():
            return {k: self.getInteraction(v) for k, v in self._interaction2id.items()}

    def getNumberOfInteractions(self):

        if pmi.workerIsActive():
            return self.cxxclass.getNumberOfInteractions(self)

    def getInteraction(self, number):

        if pmi.workerIsActive():
            ni = self.getNumberOfInteractions()
            if ni > 0:
                if number >=0 and number < ni: 
                    return self.cxxclass.getInteraction(self, number)
                else:
                    raise Error("Interaction number %i does not exist" % number)
            else:
                raise Error("interaction list of system is empty")

    def getInteractionByName(self, name):
        if pmi.workerIsActive():
            return self.getInteraction(self._interaction2id[name])
            
    def scaleVolume(self, *args):

        if pmi.workerIsActive():
          if len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, Real3D):
              #print arg0," is a Real3D object"
              self.cxxclass.scaleVolume( arg0 )
            elif hasattr(arg0, '__iter__'):
              if len(arg0) == 3:
                #print args, " has iterator and length 3"
                self.cxxclass.scaleVolume(self, toReal3DFromVector(arg0) )
              elif len(arg0) == 1:
                #print args, " has iterator and length 1"
                self.cxxclass.scaleVolume(self, toReal3DFromVector(arg0[0], arg0[0], arg0[0]) )
              else:
                print args, " is invalid"
            else:
              #print args, " is scalar"
              self.cxxclass.scaleVolume(self, toReal3DFromVector( [arg0, arg0, arg0] ) )
          elif len(args) == 3:          
            #print args, " is 3 numbers"
            self.cxxclass.scaleVolume(self, toReal3DFromVector(*args) )
          else:
            print args, " is invalid"
          
    def setTrace(self, switch):

        if pmi.workerIsActive():
            self.cxxclass.setTrace(self, switch)

if pmi.isController:
  class System(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.SystemLocal',
      pmiproperty = ['storage', 'bc', 'rng', 'skin', 'maxCutoff', 'integrator'],
      pmicall = ['addInteraction','removeInteraction', 'removeInteractionByName',
            'getInteraction', 'getNumberOfInteractions','scaleVolume', 'setTrace',
            'getAllInteractions', 'getInteractionByName']
    )
