#  Copyright (C) 2017,2018
#      Max Planck Institute for Polymer Research
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
********************************
espressopp.analysis.RadGyrXProfilePI
********************************

Class to compute the radius of gyration profile in adaptive path integral-based simulations along slabs in the x-direction of the system for specified particle type.

Examples:

>>> gyrationprofile_instance = espressopp.analysis.RadGyrXProfilePI(system=system)
>>> # creates instance of the class for calculating the radius of gyration profile

>>> gyrationprofile_list = gyrationprofile_instance.compute(bins=100, ntrotter=32, ptype=2)
>>> # computes the radius of gyration profile for particles of type 2 using 100 bins. The system uses 32 Trotter beads.

.. function:: espressopp.analysis.RadGyrXProfilePI(system)

                Constructs the RadGyrXProfilePI object.

                :param system: system object
                :type system: shared_ptr<System>

.. function:: espressopp.analysis.RadGyrXProfilePI.compute(bins, ntrotter, ptype):

                Calculates the radius of gyration profile in x-direction.

                :param bins: number of bins
                :param ntrotter: number of Trotter beads
                :param ptype: particle type
                :type bins: int
                :type ntrotter: int
                :type ptype: int
                :rtype: list of reals
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_RadGyrXProfilePI

class RadGyrXProfilePILocal(ObservableLocal, analysis_RadGyrXProfilePI):
  'The (local) class for computing the radius of gyration profile in x-direction of path integral ring polymers.'
  def __init__(self, system):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, analysis_RadGyrXProfilePI, system)

  def compute(self, bins, ntrotter, ptype):
    return self.cxxclass.compute(self, bins, ntrotter, ptype)

if pmi.isController :
  class RadGyrXProfilePI(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espressopp.analysis.RadGyrXProfilePILocal'
    )
