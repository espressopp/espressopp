#Copyright (C) 2016
#      Max Planck Institute for Polymer Research & Johannes Gutenberg-Universitaet Mainz
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
*********************************************
**HDF5File** - IO Object
*********************************************

* `write()`
  write configuration to trajectory HDF5 file. By default filename is "out.xyz",
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.xyz"

* `iomode`
  Iomode: 0 serial, 1 Nto1, 2 NtoN

* `unfolded`
  False if coordinates are folded, True if unfolded. By default - False

* `append`
  True if new trajectory data is appended to existing trajectory file. By default - True

* `length_factor`
  If length dimension in current system is nm, and unit is 0.23 nm, for example, then
  length_factor should be 0.23

* `length_unit`
  It is length unit. Can be LJ, nm or A. By default - LJ

usage:

writing down trajectory

>>> dump_conf_hdf5 = espressopp.io.HDF5File(system, integrator, filename='trajectory.xyz', iomode=1)
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_hdf5.write()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_hdf5 = espressopp.io.HDF5File(system, integrator, filename='trajectory.xyz', iomode=1)
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_conf_xyz, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both exapmles will give the same result: 200 configurations in trajectory .h5 file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]`

>>> dump_conf_hdf5 = espressopp.io.HDF5File(system, integrator, filename='trj.h5', iomode=1, unfolded=False, length_factor=0.34, length_unit='nm', append=True)

will produce trj.h5 with  in nanometers // Federico P. comment: what in nanometers? It's clear coordinate but please don't leave sentence hanging!

.. function:: espressopp.io.HDF5File(system, integrator, filename, iomode, unfolded, length_factor, length_unit, append)

        :param system:
        :param integrator:
        :param filename: (default: 'out.h5')
        :param iomode: (default: 1)
        :param unfolded: (default: False)
        :param length_factor: (default: 1.0)
        :param length_unit: (default: 'LJ')
        :param append: (default: True)
        :type system:
        :type integrator:
        :type filename:
        :type iomode:
        :type unfolded:
        :type length_factor: real
        :type length_unit:
        :type append:

.. function:: espressopp.io.HDF5File.write()

        :rtype:
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *
from _espressopp import io_HDF5File

class HDF5FileLocal(ParticleAccessLocal, io_HDF5File):

  def __init__(self, system, integrator, filename='out.h5', iomode=1, unfolded=False, length_factor=1.0, length_unit='LJ', append=True):
    cxxinit(self, io_HDF5File, system, integrator, filename, iomode, unfolded, length_factor, length_unit, append)

  def write(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.write(self)


if pmi.isController :
  class HDF5File(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.io.HDF5FileLocal',
      pmicall = [ 'write' ],
      pmiproperty = ['filename', 'iomode', 'unfolded', 'length_factor', 'length_unit', 'append']
    )
