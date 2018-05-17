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
espressopp.analysis.SubregionTracking
********************************

Class to compute the number of (coarse-grained) particles that belong to a specified particle list and that reside in a specified subregion of the simulation box (when specifying a list of particles that reside in a certain subregion at the beginning of the simulation, the routine can be used, for example, to track how many of these particles still stay in the same region after some simulation time).

Examples:

>>> subregiontracking_instance = espressopp.analysis.SubregionTracking(system, span=0.75, geometry=1, pidlist=tracklist, center=[Lx/2, Ly/2, Lz/2])
>>> # creates instance of the class for calculating number of particles that belong to particle id list tracklist and reside in a subregion which is centered in the simulation box and bounded within +-0.75 in x-direction from the center

>>> number_of_particles = subregiontracking_instance.compute()
>>> # computes the number of particles belonging to specified particle id list in specified subregion of the simulation box

.. function:: espressopp.analysis.SubregionTracking(self, system, span, geometry, center, pidlist)

                Constructs the SubregionTracking object.

                :param system: system object
                :param span: radius of the subregion to be considered
                :param geometry: geometry of the subregion. Can only be in ['spherical', 'bounded-x', 'bounded-y', 'bounded-z']
                :param center: center of the subregion
                :param pidlist: list of particle ids of coarse-grained particles that are counted in the specified subregion
                :type system: shared_ptr<System>
                :type span: real
                :type geometry: str in ['spherical', 'bounded-x', 'bounded-y', 'bounded-z']
                :type center: list of 3 reals (x,y,z coordinates of center)
                :type pidlist: list of ints

.. function:: espressopp.analysis.SubregionTracking.compute():

                Calculates the number of particles that are present in specified subregion and that belong to specified particle id list.

                :rtype: real
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_SubregionTracking

class SubregionTrackingLocal(ObservableLocal, analysis_SubregionTracking):
  'The (local) class for computing the number of particles that are present in a specified subregion of the system and that belong to a specified group of particles.'
  def __init__(self, system, span, geometry, center, pidlist):
    if geometry not in ['spherical', 'bounded-x', 'bounded-y', 'bounded-z']:
      raise ValueError('Error: Geometry must be in ["spherical", "bounded-x", "bounded-y", "bounded-z"]. Your input: {}'.format(geometry))
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      geometrydict = {'spherical': 0, 'bounded-x': 1, 'bounded-y': 2, 'bounded-z': 3}
      cxxinit(self, analysis_SubregionTracking, system, span, geometrydict[geometry])
      self.cxxclass.setCenter(self, center[0], center[1], center[2])
      for pid in pidlist:
        self.cxxclass.addPID(self, pid)

if pmi.isController :
  class SubregionTracking(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.analysis.SubregionTrackingLocal'
    )
