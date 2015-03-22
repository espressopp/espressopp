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


"""
*******************************************
**DumpH5MD** - IO Object
*******************************************

Methods:
* `dump()`
    writes configuration to h5md trajectory file. Coordinates are folded by default.

* `close()`
    Closes the H5MD file.

Properties:
* `hdf_file_id`
    The id of the H5MD file that could be used in h5py module to process file in Python
    during the invocation of scripts.

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *  #NOQA
from _espressopp import io_DumpH5MD


class DumpH5MDLocal(ParticleAccessLocal, io_DumpH5MD):
    'The (local) storage of configurations.'
    def __init__(self, system, integrator, filename='out.gro',
                 h5md_group='atoms', unfolded=False, author='XXX', email='xxx'):
        if pmi.workerIsActive():
            cxxinit(self, io_DumpH5MD, system, integrator, filename, h5md_group,
                    unfolded, author, email)

    def dump(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.dump(self)


if pmi.isController:
    class DumpH5MD(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.io.DumpH5MDLocal',
            pmicall=['dump', 'close'],
            pmiproperty=['hdf_file_id']
            )
