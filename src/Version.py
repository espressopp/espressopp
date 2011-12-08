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
>>> print "Patchlevel             = ", version.patchlevel
>>> print "Compilation date       = ", version.date
>>> print "Compilation time       = ", version.time

to print a full version info string:

>>> print version.info()

"""

from espresso import pmi
from espresso.esutil import cxxinit

import _espresso
import MPI


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
            pmiproperty = ['major', 'minor', 'hgrevision', 'patchlevel', 'date', 'time', 'name'],
            pmicall = ['info']
            )

