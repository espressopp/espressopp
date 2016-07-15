#  Copyright (C) 2012,2013,2014,2015,2016
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
************************************
**espressopp.analysis.RDFatomistic**
************************************

Class to compute radial distribution functions in adaptive resolution simulations in subregions of the box.

On the one hand, the RDF can be calculated in a cuboid region in the center of the box (periodic in y,z, limited in x). In this case, particle pairs are considered for which at least one of them is in the defined cuboid region (spanbase = True). This can be useful when the high resolution region has a slab geometry. No further normalization should be required. On the other hand, the routine can also calculate unnormalized RDFs using particle pairs with both particles being in the high resolution region. This can be useful when atomistic region has complicated or spherical geometries. In any case, only atoms belonging to different molecules (i.e. different coarse-grained entities) are considered.

Examples:

>>> rdf_0_1 = espressopp.analysis.RDFatomistic(system = system, type1 = 0, type2 = 1, span = 1.5, spanbased = True)
>>> # creates the class for calculating the RDF between atomistic particles of type 1 and 0 between different molecules,
>>> # within plus/minus 1.5 from the center of the box in x-direction

>>> rdf_0_1.compute(100)
>>> # computes the rdf using 100 bins over a distance corresponding to L_y / 2.0

.. function:: espressopp.analysis.RDFatomistic(system, type1, type2, span, spanbased)

                :param system: system object
                :param type1: type of atom 1
                :param type2: type of atom 2
                :param span: (default: 1.0) radius of the cuboid region used for rdf calculation (if spanbased == True)
                :param spanbased: (default: True) calculates rdfs in a cuboid region of radius span from the center (limited in x, periodic in y,z)
                :type system: shared_ptr<System>
                :type type1: int
                :type type2: int
                :type span: real
                :type spanbase: bool

.. function:: espressopp.analysis.RDFatomistic.compute(rdfN)

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

if pmi.isController :
  class RDFatomistic(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espressopp.analysis.RDFatomisticLocal'
    )
