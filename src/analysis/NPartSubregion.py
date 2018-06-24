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
espressopp.analysis.NPartSubregion
********************************

Class to compute the number of (coarse-grained) particles in a subregion of the simulation box.

Examples:

>>> subregionparticles_instance = espressopp.analysis.NPartSubregion(system, parttype=1, span=0.75, geometry=1, center=[Lx/2, Ly/2, Lz/2])
>>> # creates instance of the class for calculating number of particles of type 1 in a subregion centered in the simulation box and bounded within +-0.75 in x-direction from the center

>>> number_of_particles = subregionparticles_instance.compute()
>>> # computes the number of particles in subregion of the simulation box

.. function:: espressopp.analysis.NPartSubregion(system, parttype, span, geometry, center)

                Constructs the NPartSubregion object.

                :param system: system object
                :param parttype: particle type to be considered for particle number calculation
                :param span: radius of the subregion to be considered
                :param geometry: geometry of the subregion. Can only be in ['spherical', 'bounded-x', 'bounded-y', 'bounded-z']
                :param center: center of the subregion
                :type system: shared_ptr<System>
                :type parttype: int
                :type span: real
                :type geometry: str in ['spherical', 'bounded-x', 'bounded-y', 'bounded-z']
                :type center: list of 3 reals (x,y,z coordinates of center)

.. function:: espressopp.analysis.NPartSubregion.compute():

                Calculates the number of particles in defined subregion.

                :rtype: int
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_NPartSubregion

class NPartSubregionLocal(ObservableLocal, analysis_NPartSubregion):
  'The (local) class for computing the number of particles in a subregion of the system.'
  def __init__(self, system, parttype, span, geometry, center):
    if geometry not in ['spherical', 'bounded-x', 'bounded-y', 'bounded-z']:
      raise ValueError('Error: Geometry must be in ["spherical", "bounded-x", "bounded-y", "bounded-z"]. Your input: {}'.format(geometry))
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      geometrydict = {'spherical': 0, 'bounded-x': 1, 'bounded-y': 2, 'bounded-z': 3}
      cxxinit(self, analysis_NPartSubregion, system, parttype, span, geometrydict[geometry])
      self.cxxclass.setCenter(self, center[0], center[1], center[2])

if pmi.isController :
  class NPartSubregion(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.analysis.NPartSubregionLocal'
    )
