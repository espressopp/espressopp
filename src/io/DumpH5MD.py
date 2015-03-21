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
*******************************************
**DumpH5MD** - IO Object
*******************************************
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.ParticleAccess import *
from _espresso import io_DumpH5MD

class DumpH5MDLocal(ParticleAccessLocal, io_DumpH5MD):
  'The (local) storage of configurations.'
  def __init__(self, system, integrator, filename='out.gro',
               h5md_group, unfolded=False, length_unit='LJ', author='XXX'):
    cxxinit(self, io_DumpH5MD, system, integrator, filename, h5md_group,
            unfolded, length_unit, author)

  def dump(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.dump(self)


if pmi.isController :
    class DumpH5MD(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.io.DumpH5MDLocal',
            pmicall = ['dump']
            )
