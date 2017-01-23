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
*********************************************
espressopp.interaction.LennardJonesSoftcoreTI
*********************************************

This module is for performing simulations (e.g. as part of Thermodynamic Integration) where some Lennard-Jones interactions are a function of a parameter :math:`\lambda`, used to construct a pathway between states A and B.

For those interactions which are a function of :math:`\lambda`, the potential is softcore Lennard Jones with the following form:

.. math::

  U_S(r_{ij},\lambda) = (1-\lambda)U_H^A(r_A) + \lambda U_H^B(r_B)

  r_A=(\alpha\sigma^6_A\lambda^p+r_{ij}^6)^{1/6}

  r_B=(\alpha\sigma^6_B(1-\lambda)^p+r_{ij}^6)^{1/6}

where :math:`\epsilon_A`, :math:`\epsilon_B`, :math:`\sigma_A` and :math:`\sigma_B` are the parameters of states A and B, and :math:`\alpha` and :math:`p` are adjustable parameters of the softcore potential. The potentials :math:`U_H^A(r_A)` and :math:`U_H^B(r_B)` are the normal Lennard-Jones 12-6 hardcore potentials:

.. math::

  U_H^A(r_A) = 4.0\epsilon_A(\frac{\sigma_A}{r_A}^{12} - \frac{\sigma_A}{r_A}^6)

The user specifies a list of particles, pidlist. For all pairs of particles with particletypes interacting via this potential, the LJ interaction between two particles i and j is calculated as follows:

if (i not in pidlist) and (j not in pidlist):
  :math:`U_{H}^A` (full state A hardcore LJ interaction)
if (i in pidlist) and (j in pidlist):
  if annihilate==True:
    :math:`U_{S}` (softcore LJ interaction, function of lambda)
  if annihilate==False:
    :math:`U_{H}^A` (full state A hardcore LJ interaction)
if (i in pidlist) xor (j in pidlist):
  :math:`U_{S}` (softcore LJ interaction, function of lambda)

The default is annihilation (interactions within pidlist are coupled to lambda, and cross-interactions between particles in pidlist and particles in the rest of the system are also coupled to lambda). The alternative is decoupling (only cross-interactions between particles in pidlist and particles in the rest of the system are coupled to lambda. Interactions within pidlist are not affected by the value of lambda.) If annihilation==False, then decoupling is performed. See: http://www.alchemistry.org/wiki/Decoupling_and_annihilation

Exclusions apply as normal, i.e. interactions are only calculated for pairs of particles not already excluded.

This class does not do any automatic shifting of the potential.

So far only VerletListAdressLennardJonesSoftcoreTI is implemented, however VerletListLennardJonesSoftcoreTI, VerletListHadressLennardJonesSoftcoreTI, etc. can also be easily implemented.

The :math:`\lambda` (``lambdaTI``) parameter used here should not be confused with the :math:`\lambda` (``lambda_adr``) particle property used in AdResS simulations.

See also the Thermodynamic Integration tutorial.

Example python script:

>>> #value of lambda
>>> lambdaTI = 0.3
>>> #softcore parameters
>>> alphaSC = 0.5
>>> powerSC = 1.0
>>> #make list of indices of particles whose LJ parameters are different in TI states A and B
>>> pidlist = [1,2,3,4]
>>> #create interaction using VerletListAdress object and FixedTupleListAdress object
>>> lj_adres_interaction=espressopp.interaction.VerletListAdressLennardJonesSoftcoreTI(verletlist, ftpl)
>>> #loop over list of all types for particles interacting with this atomistic potential
>>> for i in types:
>>>   for k in types:
>>>     ljpot = espressopp.interaction.LennardJonesSoftcoreTI(epsilonA=epsA[i][k], sigmaA=sigA[i][k], epsilonB=epsB[i][k], sigmaB=sigB[i][k], alpha=alphaSC, power=powerSC, cutoff=cutoff, lambdaTI=lambdaTI, annihilate=False)
>>>     ljpot.addPids(pidlist)
>>>     lj_adres_interaction.setPotentialAT(type1=i, type2=k, potential=ljpot)
>>> system.addInteraction(lj_adres_interaction)

During the MD run, one can then calculate the derivative of the RF energy wrt lambda

>>> #calculate dU/dlambda
>>> dUdl = lj_adres_interaction.computeEnergyDeriv()

.. function:: espressopppp.interaction.LennardJonesSoftcoreTI(epsilonA, sigmaA, epsilonB, sigmaB, alpha, power, cutoff, lambdaTI, annihilate)

		:param epsilonA: (default: 1.0) LJ interaction parameter
		:param sigmaA: (default: 1.0) LJ interaction parameter
		:param epsilonB: (default: 0.0) LJ interaction parameter
		:param sigmaB: (default: 1.0) LJ interaction parameter
		:param alpha: (default: 1.0) softcore parameter
		:param power: (default: 1.0) softcore parameter
		:param cutoff: (default: infinity) interaction cutoff
		:param lambdaTI: (default: 0.0) TI lambda parameter
		:param annihilate: (default: True) switch between annihilation and decoupling
		:type epsilonA: real
		:type sigmaA: real
		:type epsilonB: real
		:type sigmaB: real
		:type alpha: real
		:type power: real
		:type cutoff: real
		:type lambdaTI: real
		:type annihilate: bool

.. function:: espressopppp.interaction.LennardJonesSoftcoreTI.addPids(pidlist)

		:param pidlist: list of particle ids of particles whose interaction parameters differ in state A and B
		:type pidlist: python list

.. function:: espressopppp.interaction.VerletListAdressLennardJones(vl, fixedtupleList)

		:param vl: Verlet list
		:param fixedtupleList: list of tuples describing mapping between CG and AT particles 
		:type vl: VerletListAdress object
		:type fixedtupleList: FixedTupleListAdress object

.. function:: espressopppp.interaction.VerletListAdressLennardJones.setPotentialAT(type1, type2, potential)

		:param type1: atomtype
		:param type2: atomtype
		:param potential: espressopppp potential
		:type type1: int
		:type type2: int
		:type potential: Potential

.. function:: espressopppp.interaction.VerletListAdressLennardJones.setPotentialCG(type1, type2, potential)

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
from _espressopp import interaction_LennardJonesSoftcoreTI, \
                      interaction_VerletListAdressLennardJonesSoftcoreTI
                      #interaction_VerletListLennardJonesSoftcoreTI, \
                      #interaction_VerletListHadressLennardJonesSoftcoreTI

#NOTE: to use LennardJonesSoftcoreTI with VerletList or VerletListHadress, uncomment and check the relevant code in this file and LennardJonesSoftcoreTI.cpp, and implement computeEnergyDeriv in the relevant interaction template

class LennardJonesSoftcoreTILocal(PotentialLocal, interaction_LennardJonesSoftcoreTI):
    def __init__(self, epsilonA=1.0, sigmaA=1.0, epsilonB=0.0, sigmaB=1.0, alpha=1.0, power=1.0, 
                 cutoff=infinity, lambdaTI=0.0, annihilate=True):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
           if sigmaB == 0.0:
             print 'Error in LennardJonesSoftcoreTI!\n sigmaB should never be 0.0 even when epsilonB is 0.0'
             return
           cxxinit(self, interaction_LennardJonesSoftcoreTI, 
                   epsilonA, sigmaA, epsilonB, sigmaB, alpha, power, cutoff, lambdaTI, annihilate)

    def addPids(self, pidlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          for pid in pidlist:
            self.cxxclass.addPid(self, pid)

#class VerletListLennardJonesSoftcoreTILocal(InteractionLocal, interaction_VerletListLennardJonesSoftcoreTI):
#    def __init__(self, vl):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            cxxinit(self, interaction_VerletListLennardJonesSoftcoreTI, vl)
#
#    def setPotential(self, type1, type2, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotential(self, type1, type2, potential)
#
#    def getPotential(self, type1, type2):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            return self.cxxclass.getPotential(self, type1, type2)
#
#    def getVerletListLocal(self):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            return self.cxxclass.getVerletList(self)

class VerletListAdressLennardJonesSoftcoreTILocal(InteractionLocal, interaction_VerletListAdressLennardJonesSoftcoreTI):
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesSoftcoreTI, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
#class VerletListHadressLennardJonesSoftcoreTILocal(InteractionLocal, interaction_VerletListHadressLennardJonesSoftcoreTI):
#    def __init__(self, vl, fixedtupleList):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            cxxinit(self, interaction_VerletListHadressLennardJonesSoftcoreTI, vl, fixedtupleList)
#
#    def setPotentialAT(self, type1, type2, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotentialAT(self, type1, type2, potential)
#
#    def setPotentialCG(self, type1, type2, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
if pmi.isController:
    class LennardJonesSoftcoreTI(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesSoftcoreTILocal',
            pmicall = [ 'addPids' ]
            )

    #class VerletListLennardJonesSoftcoreTI(Interaction):
    #    __metaclass__ = pmi.Proxy
    #    pmiproxydefs = dict(
    #        cls =  'espressopp.interaction.VerletListLennardJonesSoftcoreTILocal',
    #        pmicall = ['setPotential', 'getPotential', 'getVerletList']
    #        )

    class VerletListAdressLennardJonesSoftcoreTI(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesSoftcoreTILocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    #class VerletListHadressLennardJonesSoftcoreTI(Interaction):
    #    __metaclass__ = pmi.Proxy
    #    pmiproxydefs = dict(
    #        cls =  'espressopp.interaction.VerletListHadressLennardJonesSoftcoreTILocal',
    #        pmicall = ['setPotentialAT', 'setPotentialCG']
    #        )
            
