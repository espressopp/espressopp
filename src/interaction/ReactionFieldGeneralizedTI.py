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
*************************************************
espressopp.interaction.ReactionFieldGeneralizedTI
*************************************************

This module is for performing simulations (e.g. as part of Thermodynamic Integration) where some interactions are a linear function of a parameter :math:`\lambda`.

.. math::

  U(\lambda) = (1-\lambda)U_C^A

where :math:`U_C^A` is the standard Reaction Field interaction. This allows one to perform TI where the charges in TI state A (:math:`\lambda=0`) are the particle charges contained in the particle property ``charge`` and the charges in TI state B (:math:`\lambda=1`) are zero. 

The user specifies a list of particles, pidlist. For all pairs of particles with particletypes interacting via this potential, the RF interaction between two particles i and j is calculated as follows:

if (i not in pidlist) and (j not in pidlist):
  :math:`U_{RF}` (full RF interaction)
if (i in pidlist) and (j in pidlist):
  if annihilate==True:
    :math:`(1-\lambda)U_{RF}` (RF interaction scaled by 1-lambda)
  if annihilate==False:
    :math:`U_{RF}` (full RF interaction)
if (i in pidlist) xor (j in pidlist):
  :math:`(1-\lambda)U_{RF}` (RF interaction scaled by 1-lambda)

The default is annihilation (completely turning off charges of particles in pidlist in state B, so that interactions within pidlist are turned off and also cross-interactions between particles in pidlist and particles in the rest of the system). The alternative is decoupling (only cross-interactions between particles in pidlist and particles in the rest of the system are turned off. Interactions within pidlist are not affected.) If annihilation==False, then decoupling is performed. See: http://www.alchemistry.org/wiki/Decoupling_and_annihilation

Exclusions apply as normal, i.e. interactions are only calculated for pairs of particles not already excluded.

So far only VerletListAdressReactionFieldGeneralizedTI is implemented, however VerletListReactionFieldGeneralizedTI, VerletListHadressReactionFieldGeneralizedTI, etc. can also be easily implemented.

The :math:`\lambda` (``lambdaTI``) parameter used here should not be confused with the :math:`\lambda` (``lambda_adr``) particle property used in AdResS simulations.

See also the Thermodynamic Integration tutorial.

Example python script:

>>> #value of lambda
>>> lambdaTI = 0.3
>>> #construct RF potential with parameters prefactor,kappa,epsilon1,epsilon2,cutoff as in standard RF interaction
>>> pot = espressopp.interaction.ReactionFieldGeneralizedTI(prefactor=prefactor, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilon2, cutoff=rc, lambdaTI=lambdaTI, annihilate=False)
>>> #add list of indices of particles whose charge is 0 in TI state B
>>> pidlist = [1,2,3,4]
>>> pot.addPids(pidlist)
>>> #create interaction using VerletListAdress object and FixedTupleListAdress object
>>> qq_adres_interaction=espressopp.interaction.VerletListAdressReactionFieldGeneralizedTI(verletlist, ftpl)
>>> #loop over list of all types for particles interacting with this atomistic potential
>>> for i in types:
>>>   for k in types:
>>>     qq_adres_interaction.setPotentialAT(type1=i, type2=k, potential=pot)
>>> system.addInteraction(qq_adres_interaction)

During the MD run, one can then calculate the derivative of the RF energy wrt lambda

>>> #calculate dU/dlambda
>>> dUdl = qq_adres_interaction.computeEnergyDeriv()

.. function:: espressopppp.interaction.ReactionFieldGeneralizedTI(prefactor, kappa, epsilon1, epsilon2, cutoff, lambdaTI, annihilate)

		:param prefactor: (default: 1.0) RF parameter
		:param kappa: (default: 0.0) RF parameter
		:param epsilon1: (default: 1.0) RF parameter
		:param epsilon2: (default: 80.0) RF parameter
		:param cutoff: (default: infinity) interaction cutoff
		:param lambdaTI: (default: 0.0) TI lambda parameter
		:param annihilate: (default: True) switch between annihilation and decoupling
		:type prefactor: real
		:type kappa: real
		:type epsilon1: real
		:type epsilon2: real
		:type cutoff: real
		:type lambdaTI: real
		:type annihilate: bool

.. function:: espressopppp.interaction.ReactionFieldGeneralizedTI.addPids(pidlist)

		:param pidlist: list of particle ids of particles whose charge is zero in state B
		:type pidlist: python list

.. function:: espressopppp.interaction.VerletListAdressReactionFieldGeneralized(vl, fixedtupleList)

		:param vl: Verlet list
		:param fixedtupleList: list of tuples describing mapping between CG and AT particles 
		:type vl: VerletListAdress object
		:type fixedtupleList: FixedTupleListAdress object

.. function:: espressopppp.interaction.VerletListAdressReactionFieldGeneralized.setPotentialAT(type1, type2, potential)

		:param type1: atomtype
		:param type2: atomtype
		:param potential: espressopppp potential
		:type type1: int
		:type type2: int
		:type potential: Potential

.. function:: espressopppp.interaction.VerletListAdressReactionFieldGeneralized.setPotentialCG(type1, type2, potential)

		:param type1: atomtype
		:param type2: atomtype
		:param potential: espressopppp potential
		:type type1: int
		:type type2: int
		:type potential: Potential

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_ReactionFieldGeneralizedTI, \
                      interaction_VerletListAdressReactionFieldGeneralizedTI
                      #interaction_VerletListReactionFieldGeneralizedTI, \
                      #interaction_VerletListHadressReactionFieldGeneralizedTI

#NOTE: to use ReactionFieldGeneralizedTI with VerletList or VerletListHadress, uncomment and check the relevant code in this file and ReactionFieldGeneralizedTI.cpp, and implement computeEnergyDeriv in the relevant interaction template

class ReactionFieldGeneralizedTILocal(PotentialLocal, interaction_ReactionFieldGeneralizedTI):
    def __init__(self, prefactor=1.0, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff=infinity, lambdaTI=0.0, annihilate=True):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_ReactionFieldGeneralizedTI, prefactor, kappa, epsilon1, epsilon2, cutoff, lambdaTI, annihilate)

    def addPids(self, pidlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          for pid in pidlist:
            self.cxxclass.addPid(self, pid)

#class VerletListReactionFieldGeneralizedTILocal(InteractionLocal, interaction_VerletListReactionFieldGeneralizedTI):
#    def __init__(self, vl):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            cxxinit(self, interaction_VerletListReactionFieldGeneralizedTI, vl)
#
#    def setPotential(self, type1, type2, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotential(self, type1, type2, potential)
#
#    def getPotential(self, type1, type2):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            return self.cxxclass.getPotential(self, type1, type2)
    
class VerletListAdressReactionFieldGeneralizedTILocal(InteractionLocal, interaction_VerletListAdressReactionFieldGeneralizedTI):
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressReactionFieldGeneralizedTI, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
#class VerletListHadressReactionFieldGeneralizedTILocal(InteractionLocal, interaction_VerletListHadressReactionFieldGeneralizedTI):
#    def __init__(self, vl, fixedtupleList):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            cxxinit(self, interaction_VerletListHadressReactionFieldGeneralizedTI, vl, fixedtupleList)
#
#    def setPotentialAT(self, type1, type2, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotentialAT(self, type1, type2, potential)
#
#    def setPotentialCG(self, type1, type2, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
if pmi.isController:
    class ReactionFieldGeneralizedTI(Potential):
        'The ReactionFieldGeneralizedTI potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.ReactionFieldGeneralizedTILocal',
            pmiproperty = ['prefactor'],
            pmicall = [ 'addPids' ]
            )
        
    #class VerletListReactionFieldGeneralizedTI(Interaction):
    #    __metaclass__ = pmi.Proxy
    #    pmiproxydefs = dict(
    #        cls =  'espressopp.interaction.VerletListReactionFieldGeneralizedTILocal',
    #        pmicall = ['setPotential','getPotential']
    #        )
        
    class VerletListAdressReactionFieldGeneralizedTI(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressReactionFieldGeneralizedTILocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    #class VerletListHadressReactionFieldGeneralizedTI(Interaction):
    #    __metaclass__ = pmi.Proxy
    #    pmiproxydefs = dict(
    #        cls =  'espressopp.interaction.VerletListHadressReactionFieldGeneralizedTILocal',
    #        pmicall = ['setPotentialAT', 'setPotentialCG']
    #        )     
