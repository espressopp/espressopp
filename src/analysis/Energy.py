#  Copyright (C) 2012,2013,2016
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


r"""
**************************
espressopp.analysis.Energy
**************************


.. function:: espressopp.analysis.EnergyPot(system, per_atom)

		:param system:
		:param per_atom: (default: False)
		:type system:
		:type per_atom:

.. function:: espressopp.analysis.EnergyPot.compute()

		:rtype:

.. function:: espressopp.analysis.EnergyKin(system, per_atom)

		:param system:
		:param per_atom: (default: False)
		:type system:
		:type per_atom:

.. function:: espressopp.analysis.EnergyKin.compute()

		:rtype:

.. function:: espressopp.analysis.EnergyTot(system, per_atom)

		:param system:
		:param per_atom: (default: False)
		:type system:
		:type per_atom:

.. function:: espressopp.analysis.EnergyTot.compute()

		:rtype:
"""
import espressopp

class EnergyPot():

    def __init__(self, system, per_atom=False):
        self.system   = system
        self.per_atom = per_atom

    def compute(self):
        EPot = 0.0
        for k in xrange(self.system.getNumberOfInteractions()):
          EPot += self.system.getInteraction(k).computeEnergy()
        if self.per_atom:
          NPart  = espressopp.analysis.NPart(self.system).compute()
          return EPot / NPart
        else:
          return EPot

class EnergyKin():

    def __init__(self, system, per_atom=False):
        self.system   = system
        self.per_atom = per_atom

    def compute(self):
      NPart  = espressopp.analysis.NPart(self.system).compute()
      T      = espressopp.analysis.Temperature(self.system).compute()
      EKin   = (3.0/2.0) * NPart * T
      if self.per_atom:
          return EKin / NPart
      else:
          return EKin

class EnergyTot():

    def __init__(self, system, per_atom=False):
        self.system   = system
        self.per_atom = per_atom

    def compute(self):
      NPart  = espressopp.analysis.NPart(self.system).compute()
      T      = espressopp.analysis.Temperature(self.system).compute()
      EKin   = (3.0/2.0) * NPart * T
      EPot   = 0.0
      for k in xrange(self.system.getNumberOfInteractions()):
        EPot += self.system.getInteraction(k).computeEnergy()
      if self.per_atom:
        return (EPot + EKin) / NPart
      else:
        return (EPot + EKin)
