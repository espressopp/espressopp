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
**Version** - Object
************************************

Return version information of espresso module 

Example:

>>> version     = espresso.Version()
>>> print "Name                   = ", version.name
>>> print "Major version number   = ", version.major
>>> print "Minor version number   = ", version.minor
>>> print "Mercurial(hg) revision = ", version.hgrevision
>>> print "boost version          = ", version.boostversion
>>> print "Patchlevel             = ", version.patchlevel
>>> print "Compilation date       = ", version.date
>>> print "Compilation time       = ", version.time

to print a full version info string:

>>> print version.info()

"""

from espresso import pmi
from espresso.esutil import cxxinit

import _espresso
import mpi4py.MPI as MPI


class VersionLocal(_espresso.Version):
    def __init__(self):
        'Local Version object'
        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, _espresso.Version)
            else :
                pass
        else :
            cxxinit(self, _espresso.Version)

if pmi.isController:
    class Version(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.VersionLocal',
            pmiproperty = ['major', 'minor', 'hgrevision', 'boostversion', 'patchlevel', 'date', 'time', 'name'],
            pmicall = ['info']
            )

