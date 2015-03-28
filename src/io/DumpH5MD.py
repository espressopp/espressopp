#  Copyright (C) 2015
#      Jakub Krajniak (jkrajniak at gmail.com)
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
H5MD file format based on HDF5 format. The goal is to have portable efficient storage for
molecular simulation. H5MD specifies the structure that can be used in many simulations
like molecular dynamics, monte carlo or computational fluid dynamics.
More information: http://h5md.nongnu.org/

Arguments for `DumpH5MD` object:
* `system` - The system object.
* `integrator` - The integrator object.
* `filename` - The name of the file. If the file exists, it will be renamed.
* `h5md_group` - The name of particles group.
* `unfolded` - If set to true then the coordinates will be saved as unfolded (default: false).
* `author` - The name of the author of h5md file.
* `email` - The e-mail of the author.

* `save_force` - If set to true then the forces acting on particles will be saved.
* `save_velocity` - If set to true then the velocities of particles will be saved.

Methods:
* `dump()`
    writes configuration to h5md trajectory file. Coordinates are folded by default.

* `close()`
    Closes the H5MD file.

* `flush()`
    Flush buffers to the storage.

Properties:
* `file_id`
    The id of the H5MD file that could be used in h5py module to process file in Python
    during the invocation of scripts.

>>> dump_h5md = espressopp.io.DumpH5MD(system, integrator, filename='out.h5', h5md_group='atoms',
>>>                                    author="John Lenon", email="john@lenon")
>>> for i in range(200):
>>>   integrator.run(10)
>>>   dump_h5md.dump()

The trajectory can also be written with ExtAnalyze extension
>>> dump_h5md = espressopp.io.DumpH5MD(system, integrator, filename='out.h5', h5md_group='atoms',
>>>                                    author="John Lenon", email="john@lenon")
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_h5md, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)
>>> dump_h5md.close()

Beware, you have to close the file at the end of the script. The other option is to
use context manager.

>>> with espressopp.io.DumpH5MD(system, integrator, filename='out.h5', h5md_group='atoms',
>>>                             author='John Lenon', email='john@lenon') as h5md:
>>>     ext_analyze = espressopp.integrator.ExtAnalyze(dump_h5md, 10)
>>>     integrator.run(100)
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *  #NOQA
from _espressopp import io_DumpH5MD


class DumpH5MDLocal(ParticleAccessLocal, io_DumpH5MD):
    'The (local) storage of configurations.'
    def __init__(self, system, integrator, filename='out.h5',
                 h5md_group='atoms', unfolded=False, author='XXX', email='xxx',
                 save_force=False, save_velocity=False):
        if pmi.workerIsActive():
            cxxinit(self, io_DumpH5MD, system, integrator, filename, h5md_group,
                    unfolded, author, email, save_force, save_velocity)

    def dump(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.dump(self)

    def close(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.close(self)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

if pmi.isController:
    class DumpH5MD(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.io.DumpH5MDLocal',
            pmicall=['dump', 'close', 'flush', '__enter__', '__exit__'],
            pmiproperty=['file_id']
            )
