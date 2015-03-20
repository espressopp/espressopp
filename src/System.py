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


"""
************************************
**System** - Object
************************************

The main purpose of this class is to store pointers to some
important other classes and thus make them available to C++.
In a way the System class can be viewed as a container for
system wide global variables.
If you need to run more than one system at the same time you
can combine several systems with the help of the Multisystem
class.

In detail the System class holds pointers to:
---------------------------------------------

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
>>> LJSystem.addInteraction(interLJ, 'lj')

"""

from espressopp import pmi, Real3D, toReal3DFromVector
from espressopp.esutil import cxxinit
from espressopp.Exceptions import Error

import _espressopp
import mpi4py.MPI as MPI


class SystemLocal(_espressopp.System):
    def __init__(self):
        """Local construction of a System"""
        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, _espressopp.System, pmi._PMIComm.getMPIsubcomm())
            else :
                pass
        else :
            cxxinit(self, _espressopp.System, pmi._MPIcomm)

        self._interaction_id = 0
        self._interaction_label_id = {}

    def addInteraction(self, interaction, label=None):
        """Adds a short range list interaction"""
        if pmi.workerIsActive():
            ret_val = self.cxxclass.addInteraction(self, interaction)
            if label is None:
                lbl = 'e%04d' % self._interaction_id
            else:
                lbl = label
            if lbl in self._interaction_label_id:
                raise ValueError('Interaction with label {} exists'.format(lbl))
            self._interaction_label_id[lbl] = self._interaction_id
            self._interaction_id += 1
            return ret_val

    def removeInteraction(self, number):
        """Removes an interaction.

        Args:
            number: The interaction id
        """
        if pmi.workerIsActive():
            self.cxxclass.removeInteraction(self, number)
            # Find the key based on number
            if self._interaction_label_id:
                interaction_key = [
                    k for k, v in self._interaction_label_id.iteritems() if v == number
                ][0]
                # Renumerate
                del self._interaction_label_id[interaction_key]
                self._interaction_label_id = {
                    k: v if v < number else v-1 for k, v in self._interaction_label_id.iteritems()
                }
                self._interaction_id -= 1

    def removeInteractionByLabel(self, label):
        """Removes an interaction.

        Args:
            label: The label that define interaction.
        """
        if pmi.workerIsActive():
            interaction_id = self._interaction_label_id[label]
            self.removeInteraction(interaction_id)

    def getNumberOfInteractions(self):
        'get number of interactions of the system'
        if pmi.workerIsActive():
            return self.cxxclass.getNumberOfInteractions(self)

    def getInteraction(self, number):
        'get python object of the one single interaction number i'
        if pmi.workerIsActive():
            ni = self.getNumberOfInteractions()
            if ni > 0:
                if number >=0 and number < ni:
                    return self.cxxclass.getInteraction(self, number)
                else:
                    raise Error("Interaction number %i does not exist" % number)
            else:
                raise Error("interaction list of system is empty")

    def getInteractionByLabel(self, label):
        """Gets interaction by identifying it via label

        Args:
            label: The label that defines interaction.
        """
        if pmi.workerIsActive():
            interaction_id = self._interaction_label_id[label]
            return self.getInteraction(interaction_id)

    def getInteractionMap(self):
        """Returns the dictionary with interaction label and corresponding id."""
        return self._interaction_label_id

    def printInteractionLabels(self):
        """Prints the list of interactions registerd in the system."""
        for k, v in self._interaction_label_id.iteritems():
            print('e{} -> {}'.format(v, k))

    def scaleVolume(self, *args):
        'scale the Volume of the system, which means in detail: scale all particle coordinates, scale box length, scale cellgrid (if it exists)'
        if pmi.workerIsActive():
          if len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, Real3D):
              #print arg0," is a Real3D object"
              self.cxxclass.scaleVolume( arg0 )
            elif hasattr(arg0, '__iter__'):
              if len(arg0) == 3:
                #print args, " has iterator and length 2"
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
        'switch on or off VampirTrace'
        if pmi.workerIsActive():
            self.cxxclass.setTrace(self, switch)

if pmi.isController:
  class System(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.SystemLocal',
      pmiproperty = ['storage', 'bc', 'rng', 'skin', 'maxCutoff'],
      pmicall = [
          'addInteraction',
          'removeInteraction',
          'removeInteractioByLabel',
          'getInteraction',
          'getInteractionMap',
          'getNumberOfInteractions',
          'scaleVolume',
          'setTrace',
          'printInteractionLabels'
      ]
    )

