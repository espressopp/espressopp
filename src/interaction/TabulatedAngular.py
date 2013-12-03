#  Copyright (C) 2012,2013
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


"""
*****************************************
**espresso.interaction.TabulatedAngular**
*****************************************

"""
# -*- coding: iso-8859-1 -*-
from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_TabulatedAngular, \
                      interaction_FixedTripleListTabulatedAngular


class TabulatedAngularLocal(AngularPotentialLocal, interaction_TabulatedAngular):
    'The (local) tabulated angular potential.'
    def __init__(self, itype, filename):
        """Initialize the local TabulatedAngularLocal object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedAngular, itype, filename)

class FixedTripleListTabulatedAngularLocal(InteractionLocal, interaction_FixedTripleListTabulatedAngular):
    'The (local) tanulated angular interaction using FixedTriple lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListTabulatedAngular, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class TabulatedAngular(AngularPotential):
        'The TabulatedAngular potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.TabulatedAngularLocal',
            pmiproperty = ['itype', 'filename']
            )

    class FixedTripleListTabulatedAngular(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedTripleListTabulatedAngularLocal',
            pmicall = ['setPotential', 'getFixedTripleList']
            )
