#  Copyright (C) 2016
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
*********************************
espressopp.analysis.AdressDensity
*********************************

Class to compute radial density profiles in adaptive resolution simulations based on distance to closest AdResS center. Works also for multiple overlapping AdResS regions.

Examples:

>>> densityprofile = espressopp.analysis.AdressDensity(system, verletlist)
>>> # creates the class

>>> densityprofile.addExclusions([1,2,3])
>>> # defines particle to be excluded from the calculation based on list of particle ids

>>> densityprofile.compute(100)
>>> # computes the densityprofile using 100 bins

.. function:: espressopp.analysis.AdressDensity(system, verletlist)

        :param system: system object
        :type system: shared_ptr<System>
        :param verletlist: verletlist object
        :type verletlist: shared_ptr<VerletListAdress>

.. function:: espressopp.analysis.AdressDensity.compute(bins)

        :param bins: number of bins
        :type bins: int
        :rtype: list of reals

.. function:: espressopp.analysis.AdressDensity.addExclusions(pidlist)

        :param pidlist: list of ids of particles to be excluded from the calculation
        :type pidlist: list of ints
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_AdressDensity

class AdressDensityLocal(ObservableLocal, analysis_AdressDensity):

  def __init__(self, system, verletlist):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, analysis_AdressDensity, system, verletlist)

  def addExclusions(self, pidlist):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      for pid in pidlist:
        self.cxxclass.addExclpid(self, pid)

  def compute(self, bins):
    return self.cxxclass.compute(self, bins)

if pmi.isController :
  class AdressDensity(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ 'addExclusions', 'compute' ],
      cls = 'espressopp.analysis.AdressDensityLocal'
    )
