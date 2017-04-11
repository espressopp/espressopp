#  Copyright (C) 2014 Pierre de Buyl
#  Copyright (C) 2012,2013,2017
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
*********************************
espressopp.analysis.CMVelocity
*********************************

Compute and reset (set to zero) the center-of-mass (CM) velocity of the system. 

.. class:: espressopp.analysis.CMVelocity(system)

    :param system: system object
    :type system: espressopp.System
    
    .. note:: 
            
        :class:`CMVelocity <espressopp.analysis.CMVelocity>` can be attached to 
        the *integrator*. In this case the :meth:`reset` method is called, so 
        you will reset the CM-velocity every `n-th` steps.

    **Methods**

    .. function:: espressopp.analysis.CMVelocity.compute()

        Compute the CM-velocity of the system
        
        :rtype: Real3D

    .. function:: espressopp.analysis.CMVelocity.reset()

        Reset (set to zero) the CM-velocity of the system. Done by computing the
        CM-velocity of the system and subtracting it then from every particle.
        
        :rtype: void

    Example of resetting velocity

    >>> total_velocity = espressopp.analysis.CMVelocity(system)
    >>> total_velocity.reset()

    Example of attaching to integrator

    >>> # This extension can be attached to integrator
    >>> # and run `reset()` every `n-th` steps.
    >>> total_velocity = espressopp.analysis.CMVelocity(system)
    >>> ext_remove_com = espressopp.integrator.ExtAnalyze(total_velocity, 10)
    >>> integrator.addExtension(ext_remove_com)

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.AnalysisBase import *
from _espressopp import analysis_CMVelocity

class CMVelocityLocal(AnalysisBaseLocal, analysis_CMVelocity):

    def __init__(self, system):
	if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, analysis_CMVelocity, system)

if pmi.isController :
    class CMVelocity(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.CMVelocityLocal',
            pmiproperty = ["v"]
            )
