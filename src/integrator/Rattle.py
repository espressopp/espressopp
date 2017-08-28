#  Copyright (C) 2016
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
******************************
espressopp.integrator.Rattle
******************************

RATTLE algorithm for satisfying bond constraints and making the corresponding velocity corrections. 

Refs: 

Andersen, H. C. Rattle: A velocity version of the Shake algorithm for molecular dynamics calculations, J. Comp. Physics, 52, 24-34 (1983)

Allen & Tildesley, Computer Simulation of Liquids, OUP, 1987

RATTLE is implemented as an integrator extension, and takes as input a list of lists detailing, for each bond to be constrained: the indices of the two particles involved, the constraint distance, and the particle masses. 

This implementation is intended for use with hydrogen-heavy atom bonds, which form isolated groups of constrained bonds, e.g NH2 or CH3 groups. The particle which participates in only one constrained bond (i.e. the hydrogen) should be listed first. The particle listed second (the heavy atom) may participate in more than one constrained bond. This implementation will not work if both particles participate in more than one constrained bond.

Note: At the moment, the RATTLE implementation only works if all atoms in an isolated group of rigid bonds are on the same CPU. This can be achieved by grouping all the particles using DomainDecompositionAdress and FixedTupleListAdress. The groups of rigid bonds can be identified using the dictionary constrainedBondsDict (see example below).

Note: The constraints are not taken into account in other parts of the code, such as temperature or pressure calculation.

Python example script for one methanol molecule where atoms are indexed in the order C H1 H2 H3 OH HO:

>>> # list for each constrained bond which lists: heavy atom index, light atom index, bond length, heavy atom mass, light atom mass
>>> constrainedBondsList = [[1, 2, 0.109, 12.011, 1.008], [1, 3, 0.109, 12.011, 1.008], [1, 4, 0.109, 12.011, 1.008], [5, 6, 0.096, 15.9994, 1.008]]
>>> rattle = espressopp.integrator.Rattle(system, maxit = 1000, tol = 1e-6, rptol = 1e-6)
>>> rattle.addConstrainedBonds(constrainedBondsList)
>>> integrator.addExtension(rattle)

This list of lists of constrained bonds can be conveniently built using the espressopppp tool `findConstrainedBonds`.

>>> # Automatically identify hydrogen-containing bonds among the particles whose indices are in the list pidlist
>>> # pidlist - list of indices of particles in which to search for hydrogens (list of int)
>>> # masses - list of masses of all particles (list of real)
>>> # massCutoff - atoms with mass < massCutoff are identified as hydrogens (real)
>>> # bondtypes - dictionary (e.g. obtained using espressopppp.gromacs.read()), key: bondtype (int), value: list of tuples of the indices of the two particles in each bond of that bondtype (list of 2-tuples of integers)
>>> # bondtypeparams - dictionary (e.g. obtained using espressopppp.gromacs.read()), key: bondtype (int), value: espressopppp interaction potential instance
>>> hydrogenIDs, constrainedBondsDict, constrainedBondsList = espressopp.tools.findConstrainedBonds(pidlist, bondtypes, bondtypeparams, masses, massCutoff = 1.1)
>>> # hydrogenIDs - list of indices of hydrogen atoms
>>> # constrainedBondsDict - dictionary mapping from a heavy atom to all the light atoms it is bonded to, key: heavy atom index (int), value: list of light atom indices (list of int)
>>> # constrainedBondsList - list of lists, constrained bonds for use with Rattle.addConstrainedBonds()
>>> print "# found", len(hydrogenIDs)," hydrogens in the solute"
>>> print "# found", len(constrainedBondsDict)," heavy atoms involved in bonds to hydrogen"
>>> print "# will constrain", len(constrainedBondsList)," bonds using RATTLE"

.. function:: espressopppp.integrator.Rattle(system, maxit = 1000, tol = 1e-6, rptol = 1e-6)

                :param espressopp.System system: espressopp system
                :param int maxit: maximum number of iterations
                :param real tol: tolerance for deciding if constraint distance and current distance are similar enough
                :param real rptol: tolerance for deciding if the angle between the bond vector at end of previous timestep and current vector has become too large
  
.. function:: espressopppp.integrator.Rattle.addConstrainedBonds(bondDetailsLists)

                :param bondDetailsLists: list of lists, each list contains pid of heavy atom, pid of light atom, constraint distance, mass of heavy atom, mass of light atom
                :type bondDetailsLists: list of [int, int, real, real, real]
"""

from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_Rattle

class RattleLocal(ExtensionLocal, integrator_Rattle):

    def __init__(self, system, maxit = 1000, tol = 1e-6, rptol = 1e-6):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, integrator_Rattle, system, maxit, tol, rptol)

    def addConstrainedBonds(self, bondDetailsLists):
        """
        Each processor takes the broadcasted list.
        """
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          for blist in bondDetailsLists: #each list contains int pid1, int pid2, real constraintDist, real mass1, real mass2
            self.cxxclass.addBond(self, blist[0], blist[1], blist[2], blist[3], blist[4])

if pmi.isController:
    class Rattle(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.RattleLocal',
            pmicall = [ "addConstrainedBonds" ]
            )
