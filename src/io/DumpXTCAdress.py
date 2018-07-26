#  Copyright (C) 2017
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


r"""
*********************************************
**DumpXTCAdress** - IO Object
*********************************************

* `dump()`
  write configuration to trajectory file in XTC format (Gromacs/Gromos compressed coordinate trajectory). By default filename is "out.xtc",
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.xtc"

* `unfolded`
  False if coordinates are folded, True if unfolded. By default - False

* `append`
  True if new trajectory data is appended to existing trajectory file. By default - True

* `length_factor`
  If length dimension in current system is nm, and unit is 0.23 nm, for example, then
  length_factor should be 0.23

usage:

writing down trajectory

>>> dump_conf_xtc = espressopp.io.DumpXTCAdress(system, ftpl, integrator, filename='trajectory.xtc')
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_xtc.dump()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_xtc = espressopp.io.DumpXTCAdress(system, ftpl, integrator, filename='trajectory.xtc')
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_conf_xtc, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both exapmles will give the same result: 200 configurations in trajectory .xtc file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]`

>>> dump_conf_xtc = espressopp.io.DumpXTCAdress(system, integrator, filename='trj.xtc', unfolded=False, length_factor=0.34, append=True)

will produce trj.xtc with coordinates in nm

.. function:: espressopp.io.DumpXTCAdress(system, integrator, ftpl, filename, unfolded, length_factor, append)

        :param system: The Espressopp system object.
        :param ftpl: The FixedTupleListAdress object.
        :param integrator: The MD Integrator.
        :param filename: The output file name (default: 'out.xtc').
        :param unfolded: Store unfolded positions (default: False).
        :param length_factor: The length factor (default: 1.0).
        :param append: If True then new frames will be append to existing file (default: True).
        :type system: espressopp.System
        :type integrator: espressopp.integrator.MDIntegrator
        :type filename: str
        :type unfolded: bool
        :type length_factor: real
        :type append: bool

.. function:: espressopp.io.DumpXTCAdress.dump()

        :rtype:
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *
from _espressopp import io_DumpXTCAdress

class DumpXTCAdressLocal(ParticleAccessLocal, io_DumpXTCAdress):
    def __init__(self, system, ftpl, integrator, filename='out.xtc', unfolded=False, length_factor=1.0, append=True):
        cxxinit(self, io_DumpXTCAdress, system, ftpl, integrator, filename, unfolded, length_factor, append)

    def dump(self):
        if pmi.workerIsActive():
            self.cxxclass.dump(self)


if pmi.isController :
    class DumpXTCAdress(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.io.DumpXTCAdressLocal',
            pmicall = ('dump',),
            pmiproperty = ('filename', 'unfolded', 'length_factor', 'append')
        )
