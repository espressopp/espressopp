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

*******************************
PSF - read and write psf format
*******************************

PSF file format given at http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node24.html

.. function:: espressopp.tools.psfwrite(filename, system, maxdist=None, molsize=4, typenames=None)

  Writes !NATOM and !NBOND sections of a psf format file

  :param filename: output file name
  :type filename: string
  :param system: espressopp system
  :type system: espressopp System object
  :param maxdist: if this is specified, only bonds in which the pair of particles are separated by a distance < maxdist are written to the !NBOND section
  :type maxdist: float
  :param molsize: if molsize>0, the molecule count is increased every molsize particles
  :type molsize: int
  :param typenames: dictionary used for mapping from espressopp's integer particle types to the particle type strings written in a psf file
  :type typenames: dict, key=int, value=string

.. function:: espressopp.tools.psfread(filename)

  Reads !NATOM section of a psf format file

  :param filename: input file name
  :type filename: string
  :returns: pid,segname,resindex,resname,atomname,atomtype,mass,charge (lists of type int,str,int,str,str,str,float,float)

"""

import espressopp
from math import sqrt

def psfwrite(filename, system, maxdist=None, molsize=4, typenames=None):
  #if not all molecules have the same number of atoms, set molsize=0
  file = open(filename,'w')
  maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
  nParticles    = int(espressopp.analysis.NPart(system).compute())
  file.write("PSF CMAP\n")
  file.write(" \n")
  file.write("       1 !NTITLE\n")
  file.write("  REMARK remark")
  file.write(" \n")
  st = "\n%8d !NATOM\n" % nParticles
  file.write(st)
  
  pid    = 0
  addToPid = 0 # if pid begins from 0, then addToPid should be +1
  if (molsize>0):
    mol    = 0
    molcnt = 0
  else: #not all molecules have same number of atoms
    mol    = 1
  name='FE' # default name, overwritten when typenames map is given
  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
      particle = system.storage.getParticle(pid)
      if(pid==0):
        addToPid = 1
      xpos   = particle.pos[0]
      ypos   = particle.pos[1]
      zpos   = particle.pos[2]
      type   = particle.type
      q      = particle.q
      mass   = particle.mass
      if typenames:
	  name=typenames[type]
      st = "%8d T%5d    UNX  %-4s %-4s  %9.6f%14.4f\n" % (pid+addToPid, mol, name, name, q, mass)
      file.write(st)
      pid    += 1
      if (molsize>0):
        molcnt += 1
        if molcnt == molsize:
          mol   += 1
          molcnt = 0
    else:
      pid += 1

  bond = []
  nInteractions = system.getNumberOfInteractions()
  for i in xrange(nInteractions):
      if system.getInteraction(i).bondType() == espressopp.interaction.Pair:
        try:
           
          FixedPairList = system.getInteraction(i).getFixedPairList().getBonds()
          j = 0
          while j < len(FixedPairList):
              fplb = FixedPairList[j]
              k = 0
              while k < len(fplb):
                if maxdist != None:
                  pid1 = fplb[k][0]
                  pid2 = fplb[k][1]
                  p1 = system.storage.getParticle(pid1)
                  p2 = system.storage.getParticle(pid2)
                  x1 = p1.pos[0]
                  y1 = p1.pos[1]
                  z1 = p1.pos[2]
                  x2 = p2.pos[0]
                  y2 = p2.pos[1]
                  z2 = p2.pos[2]
                  xx = (x1-x2) * (x1-x2)
                  yy = (y1-y2) * (y1-y2)
                  zz = (z1-z2) * (z1-z2)
                  d = sqrt( xx + yy + zz )
                  if (d <= maxdist):
                    bond.append(fplb[k])
                else:
                  bond.append(fplb[k])
                k += 1
                
              j += 1
              
        except:
          pass
              
  bond.sort()

  file.write("\n%8d !NBOND:\n" % (len(bond)))
  i = 0
  while i < len(bond):
    file.write("%8d%8d" % (bond[i][0]+addToPid, bond[i][1]+addToPid) ) #pid_count_translate[bond[i][1]]
    if ( ((i+1) % 4) == 0 and (i != 0) ) or i == len(bond)-1 :
      file.write("\n")
    i += 1

  file.write('END\n')
  file.close()
  
  
  
  """
  NOT finished yet - working on that see lammps.write how to do this!
   
  bonds     = []
  angles    = []
  dihedrals = [] 
  nInteractions = system.getNumberOfInteractions()
  for i in xrange(nInteractions):
      bT = system.getInteraction.bondType
      if   bT == espressopp.interaction.Pair:
             bl = system.getInteraction(i).getFixedPairList().getBonds
             for j in xrange(len(bl)):
               bonds.extend(bl[j])
      elif bT == espressopp.interaction.Angle:
             an = system.getInteraction(i).getFixedTripleList().getTriples
             for j in xrange(len(an)):
               angles.extend(an[j])
      elif bT == espressopp.interaction.Dihedral:
             di = system.getInteraction(i).getFixedQuadrupleList().getQuadruples
             for j in xrange(len(di)):
               dihedrals.extend(di[j])
  
  nbonds     = len(bonds)
  nangles    = len(angles)
  ndihedrals = len(dihedrals)
"""
  
def psfread(filename):

  pid=[]
  segname=[]
  resname=[]
  resindex=[]
  atomname=[]
  atomtype=[]
  mass=[]
  charge=[]

  file = open(filename,'r')
  line = file.readline()
  line = file.readline()
  line = file.readline()
  nskip = int(line.split()[0]) #get NREMARKS value
  for i in range(nskip):
    line = file.readline()
  line = file.readline()
  line = file.readline()
  natoms = int(line.split()[0]) #get NATOMS value
  for i in range(natoms):
    line = file.readline()
    pid.append(int(line[3:8]))
    segname.append(line[8:12])
    resindex.append(int(line[10:15]))
    resname.append(line[19:23])
    atomname.append(line[24:28])
    atomtype.append(line[29:33])
    charge.append(float(line[34:44]))
    mass.append(float(line[45:59]))
  file.close()
  return pid,segname,resindex,resname,atomname,atomtype,mass,charge
