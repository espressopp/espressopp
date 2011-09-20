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

if pmi.isController:
    class System(object):
        'System object.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SystemLocal',
            pmiproperty = ['storage', 'bc', 'rng', 'skin'],
            pmicall = ['addInteraction','getInteraction','getNumberOfInteractions']
            )

