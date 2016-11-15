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

*****************************************
prepareComplexMolecules - set up proteins
*****************************************

various helper functions for setting up systems containing complex molecules such as proteins

.. function:: espressopp.tools.findConstrainedBonds(atomPids, bondtypes, bondtypeparams, masses, massCutoff = 1.1)

  Finds all heavyatom-hydrogen bonds in a given list of particle IDs, and outputs a list describing the bonds, in a format suitable for use with the RATTLE algorithm for constrained bonds
  
  :param atomPids: list of pids of atoms between which to search for bonds
  :type atomPids: list of int
  :param bondtypes: dictionary mapping from tuple of pids to bondtypeid, e.g. as returned by tools.convert.gromacs.read()
  :type bondtypes: dict, key: (int,int), value: int
  :param bondtypeparams: dictionary mapping from bondtypeid to class storing parameters of that bond type, e.g. as returned by tools.convert.gromacs.read()
  :type bondtypeparams: dict, key: int, value: espressopp bond type
  :param masses: list of masses, e.g. as returned by tools.convert.gromacs.read()
  :type masses: list of float
  :param massCutoff: for identifying light atoms (hydrogens), default 1.1 mass units, can also be increased e.g. for use with deuterated systems
  :type massCutoff: float
  :returns:
    | hydrogenIDs - list of pids (integer) of light atoms (hydrogens)
    | constrainedBondsDict - dict, keys: pid (integer) of heavy atom, values: list of pids of light atoms that are bonded to it
    | constrainedBondsList - list of lists, one entry for each constrained bond with format: [pid of heavy atom, pid of light atom, bond distance, mass of heavy atom ,mass of light atom]
  
  Can then be used with RATTLE, e.g.
  
  >>> rattle = espresso.integrator.Rattle(system, maxit = 1000, tol = 1e-6, rptol = 1e-6)
  >>> rattle.addConstrainedBonds(constrainedBondsList)
  >>> integrator.addExtension(rattle)


.. function:: espressopp.tools.getInternalNonbondedInteractions(atExclusions,pidlist)

  Gets the non-bonded pairs within a list of particle indices, excluding those which are in a supplied list of exclusions. Useful for example for getting the internal atomistic non-bonded interactions in a coarse-grained particle and adding them as a fixedpairlist
  
  :param atExclusions: list of excluded pairs
  :type atExclusions: list of 2-tuples of int
  :param pidlist: list of pids among which to create pairs
  :type pidlist: list of int
  :return: list of pairs which are not in atExclusions
  :rtype: list of 2-tuples of int
"""

def findConstrainedBonds(atomPids, bondtypes, bondtypeparams, masses, massCutoff = 1.1):

  hydrogenIDs = []
  constrainedBondsDict = {} 
  constraintDistances = {} #keys: (pid heavy, pid light), values: bond length

  # first find hydrogens
  for pid in atomPids: 
    if masses[pid-1] < massCutoff:
      hydrogenIDs.append(pid)
  # then find hydrogen-containing bonds
  for bid, bondlist in bondtypes.iteritems():
    for pair in bondlist:
      pidHyd = pidHea = 0
      if pair[0] in hydrogenIDs:
        pidHyd = pair[0]
        pidHea = pair[1]
      elif pair[1] in hydrogenIDs:
        pidHyd = pair[1]
        pidHea = pair[0]
      if (pidHyd > 0):
        if pidHea in constrainedBondsDict.keys():
          constrainedBondsDict[pidHea].append(pidHyd)
        else:
          constrainedBondsDict[pidHea] = [pidHyd]
        constraintDistances[(pidHea,pidHyd)] = bondtypeparams[bid].parameters['b0']
  
  constrainedBondsList = []
  for pidHea in constrainedBondsDict.keys():
    for pidHyd in constrainedBondsDict[pidHea]:
      constrainedBondsList.append([pidHea,pidHyd,constraintDistances[(pidHea,pidHyd)],masses[pidHea-1],masses[pidHyd-1]])

  return hydrogenIDs, constrainedBondsDict, constrainedBondsList

def getInternalNonbondedInteractions(atExclusions,pidlist):
  nonBondPairs = []
  for pid1 in pidlist:
    for pid2 in pidlist:
      if pid1==pid2: continue
      if (pid1,pid2) in atExclusions: continue
      if (pid2,pid1) in atExclusions: continue
      if (pid1,pid2) in nonBondPairs: continue #to avoid this we'd have to assume pids in pidlist were ordered and continuous
      if (pid2,pid1) in nonBondPairs: continue #to avoid this we'd have to assume pids in pidlist were ordered and continuous
      nonBondPairs.append((pid1,pid2))
  return nonBondPairs
