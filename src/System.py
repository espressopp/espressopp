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

>>> LJSystem      = espresso.System()
>>> LJSystem.bc   = espresso.bc.OrthorhombicBC(rng, boxsize)
>>> LJSystem.rng
>>> LJSystem.skin = 0.4
>>> LJSystem.addInteraction(interLJ)

"""

from espresso import pmi
from espresso.esutil import cxxinit
from espresso.Exceptions import Error

import _espresso
import MPI


class SystemLocal(_espresso.System):
    def __init__(self):
        'Local construction of a System'
        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, _espresso.System, pmi._PMIComm.getMPIsubcomm())
            else :
                pass
        else :
            cxxinit(self, _espresso.System, pmi._MPIcomm)

    def addInteraction(self, interaction):
        'add a short range list interaction'
        if pmi.workerIsActive():
            return self.cxxclass.addInteraction(self, interaction)

    def removeInteraction(self, number):
        'remove a short range interaction, number is the number of the interaction in the short range interaction list'
        if pmi.workerIsActive():
            self.cxxclass.removeInteraction(self, number)

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
            
    def scaleVolume(self, factor):
        'scale the Volume of the system, which means in detail: scale all particle coordinates, scale box length, scale cellgrid (if it exists)'
        if pmi.workerIsActive():
            self.cxxclass.scaleVolume(self, factor)

if pmi.isController:
    class System(object):
        'System object.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SystemLocal',
            pmiproperty = ['storage', 'bc', 'rng', 'skin'],
            pmicall = ['addInteraction','removeInteraction','getInteraction','getNumberOfInteractions','scaleVolume']
            )

