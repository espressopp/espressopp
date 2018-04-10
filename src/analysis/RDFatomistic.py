#  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
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
********************************
espressopp.analysis.RDFatomistic
********************************

Class to compute radial distribution functions in adaptive resolution simulations in subregions of the box. Can be used for regular atomistic/coarse-grained (AT/CG) adaptive resolution simulations as well as path integral-based adaptive resolution simulations. The two functions (compute, computePathIntegral) exhibit different behavior.

The regular compute function is used for regular AT/CG simulations and there are two options:

Option 1 (spanbased = True): the RDF can be calculated in a cuboid region in the center of the box (periodic in y,z, limited in x). In this case, particle pairs are considered for which at least one of them is in the defined cuboid region. This can be useful when the high resolution region has a slab geometry. No further normalization should be required.

Option 2 (spanbased = False): the routine can also calculate unnormalized RDFs using particle pairs with both particles being in the high resolution region (based on the resolution value lambda, the span parameter is not used then). This can be useful when atomistic region has complicated or spherical geometries.

In any case, only pairs of atomistic particles belonging to two different coarse-grained particles are considered. Furthermore, note that the routine uses L_y / half (L_y is the box length in y-direction) as the maximum distance for the RDF calculation, which is then binned according to rdfN during the computation. Hence, L_y should be the shortest box side (or, equally short as L_x and/or L_z).


The computePathIntegral function is used for path integral-based adaptive resolution functions. It calculates the radial distribution functions over pairs of particles between different atoms or coarse-grained beads. Note, however, that in these types of quantum/classical adaptive resolution simulations, regular coarse-grained espressopp particles are associated with each atom and the additional "AdResS" atomistic particles correspond to the different Trotter beads. This means that the routine will, for molecules consisting of multiple atoms, calculate intramolecular rdfs, averaging over the Trotter bead pairs of the ring polymers, which represent the atoms. In doing so, it considers only particles pair with matching Trotter number and with the correct atomistic types. The results are averaged over all Trotter beads. Also in this case L_y / half (L_y is the box length in y-direction) is used as the maximum distance for the RDF calculation, which is then binned according to rdfN during the computation. Furthermore, the calculation is always "spanbased" in x direction (the function ignores the spanbased flag), but in such a fashion that BOTH particles need to be in the defined cuboid region. Normalization is performed as derived in R. Potestio et al., Phys. Rev. Lett. 111, 060601 (2013), Supp. Info.

Examples:

>>> rdf_0_1 = espressopp.analysis.RDFatomistic(system = system, type1 = 0, type2 = 1, spanbased = True, span = 1.5)
>>> # creates the class for calculating the RDF between atomistic particles of type 1 and 0 between different molecules,
>>> # At least one of these particles has to be within plus/minus 1.5 from the center of the box in x-direction

>>> rdf_0_1.compute(100)
>>> # computes the rdf using 100 bins over a distance corresponding to L_y / 2.0

.. function:: espressopp.analysis.RDFatomistic(system, type1, type2, spanbased, span)

                Constructs the RDFatomistic object.

                :param system: system object
                :param type1: type of atom 1
                :param type2: type of atom 2
                :param spanbased: (default: True) If True, calculates RDFs in a cuboid region of radius span from the center (limited in x, periodic in y,z). If False, calculates RDFs with both particles being in the high resolution region (using lambda resolution values and ignoring span parameter).
                :param span: (default: 1.0) +/- distance from centre of box in x-direction of the cuboid region used for RDF calculation if spanbased == True. If spanbased == False, this parameter is not used.
                :type system: shared_ptr<System>
                :type type1: int
                :type type2: int
                :type spanbased: bool
                :type span: real

.. function:: espressopp.analysis.RDFatomistic.compute(rdfN)

                Calculates the atomistic RDF assuming a regular atomistic/coarse-grained adaptive resolution setup.

                :param rdfN: number of bins
                :type rdfN: int
                :rtype: list of reals

.. function:: espressopp.analysis.RDFatomistic.computePathIntegral(rdfN)

                Calculates the path integral-based RDF averaging over all Trotter bead pairs with the same Trotter bead number between different ring polymers assuming a path integral-based quantum/classical adaptive resolution setup.

                :param rdfN: number of bins
                :type rdfN: int
                :rtype: list of reals
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_RDFatomistic

class RDFatomisticLocal(ObservableLocal, analysis_RDFatomistic):

  def __init__(self, system, type1, type2, span = 1.0, spanbased = True):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, analysis_RDFatomistic, system, type1, type2, span, spanbased)

  def compute(self, rdfN):
    return self.cxxclass.compute(self, rdfN)

  def computePathIntegral(self, rdfN):
    return self.cxxclass.computePathIntegral(self, rdfN)

if pmi.isController :
  class RDFatomistic(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute", "computePathIntegral" ],
      cls = 'espressopp.analysis.RDFatomisticLocal'
    )
