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
  :param bondtypes: dictionary mapping from tuple of pids to bondtypeid, e.g. as returned by tools.gromacs.read()
  :type bondtypes: dict, key: (int,int), value: int
  :param bondtypeparams: dictionary mapping from bondtypeid to class storing parameters of that bond type, e.g. as returned by tools.gromacs.read()
  :type bondtypeparams: dict, key: int, value: espressopp bond type
  :param masses: list of masses, e.g. as returned by tools.gromacs.read()
  :type masses: list of float
  :param massCutoff: for identifying light atoms (hydrogens), default 1.1 mass units, can also be increased e.g. for use with deuterated systems
  :type massCutoff: float
  :returns:
    | hydrogenIDs - list of pids (integer) of light atoms (hydrogens)
    | constrainedBondsDict - dict, keys: pid (integer) of heavy atom, values: list of pids of light atoms that are bonded to it
    | constrainedBondsList - list of lists, one entry for each constrained bond with format: [pid of heavy atom, pid of light atom, bond distance, mass of heavy atom ,mass of light atom]
  
  Can then be used with RATTLE, e.g.
  
  >>> rattle = espressopp.integrator.Rattle(system, maxit = 1000, tol = 1e-6, rptol = 1e-6)
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

.. function:: espressopp.tools.readSimpleSystem(filename,nparticles,header)

  Read in a column-formatted file containing information about the 
  particles in a simple system, for example a coarsegrained protein.

  This function expects the input file to have between 2 and 5 columns.
  The number of columns in the file is automatically detected.
  The function reads each column into a list and returns the lists.
  Column types are interpreted as follows:
  2 columns: float, float
  3 columns: float, float, int
  4 columns: float, float, int, str
  5 columns: float, float, int, str, str

  For example in the case of a coarsegrained protein model, these could be:
  mass, charge, corresponding atomistic index, beadname, beadtype

  :param filename: name of file to open and read
  :type filename: string
  :param nparticles: number of particles in file
  :type nparticles: int
  :param header: number of lines to skip at start of file (default 0)
  :type header: int

  Returns: between 2 and 5 lists
 
.. function:: espressopp.tools.applyBoreschRestraints(system,restraintAtoms,restraintK,restraintR0)

  Applies restraints between ligand and protein as defined in Boresch et al, JPCB 2003, 107, 9535-9551
  The restraints (one bond, two angles and three dihedrals) are applied between three ligand atoms A,B,C and three protein atoms a,b,c.

  In espressopp, the potential for harmonic bond and angles is k(x-x0)^2, but the harmonic dihedral potential is 0.5*k(x-x0)^2. This is taken care of in this function, i.e. the user should supply the force constants exactly as given in Boresch et al.

  :param system: espressopp system
  :type system: System object
  :param restraintAtoms: dictionary identifying the six atoms involved in the restraints, key = atom label (one of 'A','B','C','a','b','c'), value = atom index
  :type restraintAtoms: python dictionary
  :param restraintK: dictionary of force constants, key = restraint label ('aA' for the bond, 'baA' and 'aAB' for the angles, 'aABC', 'cbaA' and 'baAB' for the dihedrals), value = force constant (units as for the simulation forcefield, angle and dihedral constants should be in rad^-2)
  :type restraintK: python dictionary
  :param restraintR0: dictionary of equilibrium values (distance or angle), key = restraint label, value = distance (in distance units) or angle (in degrees)
  :type restraintR0: python dictionary

  Examples of the three dictionaries:

  >>> restraintAtoms = {'A':1981,'B':1966,'C':1993,'a':1588,'b':1581,'c':1567}
  >>> restraintK = {'aA':4184,'baA':41.84,'aAB':41.84,'aABC':41.84,'cbaA':41.84,'baAB':41.84} #kJ mol-1 nm-2,kJ mol-1 rad-2, all 10 kcal as in JPCB 2003
  >>> restraintR0 = {'aA':0.31,'baA':120.0,'aAB':90.0,'aABC':100.0,'cbaA':-170.0,'baAB':-105.0} #nm, degrees

"""

import espressopp
import math

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

def readSimpleSystem(filename,nparticles,header=0):
  f = open(filename,'r')
  mass=[]
  charge=[]
  index=[]
  name=[]
  ptype=[]
  for i in range(header):
    f.readline()
  for i in range(nparticles):
    line = f.readline()
    line = line.split()
    mass.append(float(line[0]))
    charge.append(float(line[1]))
    if len(line)>2:
      index.append(int(line[2]))
    if len(line)>3:
      name.append(line[3])
    if len(line)==5:
      ptype.append(line[4])
  f.close()
  if len(line)==2:
    return mass,charge
  elif len(line)==3:
    return mass,charge,index
  elif len(line)==4:
    return mass,charge,index,name
  elif len(line)==5:
    return mass,charge,index,name,ptype

def applyBoreschRestraints(system,restraintAtoms,restraintK,restraintR0):

  restraintinteraction = {}

  restraintBond = espressopp.FixedPairList(system.storage)
  restraintBond.addBonds([(restraintAtoms['a'],restraintAtoms['A'])])
  potint=espressopp.interaction.FixedPairListHarmonic(system, restraintBond, potential=espressopp.interaction.Harmonic(K=0.5*restraintK['aA'],r0=restraintR0['aA'], cutoff = 5.0, shift = 0.0))
  system.addInteraction(potint)
  restraintinteraction.update({0:potint})
  
  restraintAngle = espressopp.FixedTripleList(system.storage)
  restraintAngle.addTriples([(restraintAtoms['a'],restraintAtoms['A'],restraintAtoms['B'])])
  potint=espressopp.interaction.FixedTripleListAngularHarmonic(system, restraintAngle, potential=espressopp.interaction.AngularHarmonic(K=0.5*restraintK['aAB'],theta0=restraintR0['aAB']*math.pi/180.0))
  system.addInteraction(potint)
  restraintinteraction.update({1:potint})
  
  restraintAngle = espressopp.FixedTripleList(system.storage)
  restraintAngle.addTriples([(restraintAtoms['b'],restraintAtoms['a'],restraintAtoms['A'])])
  potint=espressopp.interaction.FixedTripleListAngularHarmonic(system, restraintAngle, potential=espressopp.interaction.AngularHarmonic(K=0.5*restraintK['baA'],theta0=restraintR0['baA']*math.pi/180.0))
  system.addInteraction(potint)
  restraintinteraction.update({2:potint})
  
  restraintDih = espressopp.FixedQuadrupleList(system.storage)
  restraintDih.addQuadruples([(restraintAtoms['c'],restraintAtoms['b'],restraintAtoms['a'],restraintAtoms['A'])])
  potint=espressopp.interaction.FixedQuadrupleListDihedralHarmonic(system, restraintDih, potential=espressopp.interaction.DihedralHarmonic(K=restraintK['cbaA'],phi0=restraintR0['cbaA']*math.pi/180.0))
  system.addInteraction(potint)
  restraintinteraction.update({3:potint})
  
  restraintDih = espressopp.FixedQuadrupleList(system.storage)
  restraintDih.addQuadruples([(restraintAtoms['b'],restraintAtoms['a'],restraintAtoms['A'],restraintAtoms['B'])])
  potint=espressopp.interaction.FixedQuadrupleListDihedralHarmonic(system, restraintDih, potential=espressopp.interaction.DihedralHarmonic(K=restraintK['baAB'],phi0=restraintR0['baAB']*math.pi/180.0))
  system.addInteraction(potint)
  restraintinteraction.update({4:potint})
  
  restraintDih = espressopp.FixedQuadrupleList(system.storage)
  restraintDih.addQuadruples([(restraintAtoms['a'],restraintAtoms['A'],restraintAtoms['B'],restraintAtoms['C'])])
  potint=espressopp.interaction.FixedQuadrupleListDihedralHarmonic(system, restraintDih, potential=espressopp.interaction.DihedralHarmonic(K=restraintK['aABC'],phi0=restraintR0['aABC']*math.pi/180.0))
  system.addInteraction(potint)
  restraintinteraction.update({5:potint})

  return restraintinteraction

