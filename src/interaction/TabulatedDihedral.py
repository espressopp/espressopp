#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
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


r"""
******************************************************
**espressopp.interaction.TabulatedDihedral**
******************************************************







.. function:: espressopp.interaction.TabulatedDihedral(itype, filename)

		:param itype: 
		:param filename: 
		:type itype: 
		:type filename: 

.. function:: espressopp.interaction.FixedQuadrupleListTabulatedDihedral(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedQuadrupleListTabulatedDihedral.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 
"""
# -*- coding: iso-8859-1 -*-
# -*- coding: iso-8859-1 -*-
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.DihedralPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_TabulatedDihedral, \
                        interaction_FixedQuadrupleListTabulatedDihedral, \
                        interaction_FixedQuadrupleListTypesTabulatedDihedral


class TabulatedDihedralLocal(DihedralPotentialLocal, interaction_TabulatedDihedral):

    def __init__(self, itype, filename):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedDihedral, itype, filename)

class FixedQuadrupleListTabulatedDihedralLocal(InteractionLocal, interaction_FixedQuadrupleListTabulatedDihedral):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedQuadrupleListTabulatedDihedral, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedQuadrupleListTypesTabulatedDihedralLocal(InteractionLocal, interaction_FixedQuadrupleListTypesTabulatedDihedral):
    def __init__(self, system, vl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedQuadrupleListTypesTabulatedDihedral, system, vl)

    def setPotential(self, type1, type2, type3, type4, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, type1, type2, type3, type4, potential)

    def getPotential(self, type1, type2, type3):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self, type1, type2, type3)

    def setFixedQuadrupleList(self, fixedlist):
        if pmi.workerIsActive():
            self.cxxclass.setFixedQuadrupleList(self, fixedlist)

    def getFixedQuadrupleList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedQuadrupleList(self)

if pmi.isController:
    class TabulatedDihedral(DihedralPotential):
        'The TabulatedDihedral potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedDihedralLocal',
            pmiproperty = ['itype', 'filename']
            )

    class FixedQuadrupleListTabulatedDihedral(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedQuadrupleListTabulatedDihedralLocal',
            pmicall = ['setPotential', 'getFixedQuadrupleList']
            )

    class FixedQuadrupleListTypesTabulatedDihedral(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedralLocal',
            pmicall = ['setPotential','getPotential','setFixedQuadrupleList','getFixedQuadrupleList']
        )
