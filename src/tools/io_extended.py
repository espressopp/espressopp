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


# -*- coding: utf-8 -*-

"""
**********************************************
io_extended - read/write configurational files
**********************************************

This Python module allows one to read and write configurational files.
One can choose folded or unfolded coordinates and write down velocities or not.
It is similar to lammps read and write, but it writes down only:
1) number of particles + types
2) number of bonds (number of pairs) + types
3) number of angles (number of triples) + types
4) number of dihedrals (number of quadruples) + types
5) system size (Lx,Ly,Lz)
6) p_id, p_type, p_positions
7) velocities (if true)
8) bonds (if exist)
9) angles (if exist)
10)dihedrals (if exist)

read returns:
Lx, Ly, Lz, p_ids, p_types, poss, vels, bonds, angles, dihedrals
if something does not exist then it will return the empty list
bonds, angles, dihedrals - will return list [type, (x,x,x,x)],
where type is the type of bond, angle or dihedral
(x,x,x,x) is (pid1,pid2) for bonds,

(pid1,pid2,pid3) for angles

(pid1,pid2,pid3,pid4) for dihedrals

"""   

import espressopp

def write(fileName, system, folded=True, writeVelocities=False):
    
  # first collect all the information that we need to write into the file
  numParticles  = int(espressopp.analysis.NPart(system).compute())
  box_x = system.bc.boxL[0]
  box_y = system.bc.boxL[1]
  box_z = system.bc.boxL[2]
  
  bonds          = []
  nbondtypes     = 0
  angles         = []
  nangletypes    = 0
  dihedrals      = []
  ndihedraltypes = 0
  
  nInteractions = system.getNumberOfInteractions()
  for i in xrange(nInteractions):
      bT = system.getInteraction(i).bondType()
      if   bT == espressopp.interaction.Pair:
             nbondtypes += 1
             bl  = system.getInteraction(i).getFixedPairList().getBonds()
             bln = []
             for j in xrange(len(bl)):
               bln.extend(bl[j])
             bonds.append(bln)
      elif bT == espressopp.interaction.Angular:
             nangletypes += 1
             an  = system.getInteraction(i).getFixedTripleList().getTriples()
             ann = []
             for j in xrange(len(an)):
               ann.extend(an[j]) 
             angles.append(ann)
      elif bT == espressopp.interaction.Dihedral:
             ndihedraltypes += 1
             di  = system.getInteraction(i).getFixedQuadrupleList().getQuadruples()
             din = []
             for j in xrange(len(di)):
               din.extend(di[j])   
             dihedrals.append(din)
  
  nbonds = 0
  for i in xrange(len(bonds)):
      nbonds += len(bonds[i])
  nangles = 0
  for i in xrange(len(angles)):
      nangles += len(angles[i])
  ndihedrals = 0
  for i in xrange(len(dihedrals)):
      ndihedrals += len(dihedrals[i])
      
  atomtypes = []
  maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
  pid   = 0
  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
        particle = system.storage.getParticle(pid)
        type   = particle.type
        if type in atomtypes:
          pid   += 1
        else:
          atomtypes.append(type)
          pid += 1
    else:
        pid   += 1
  natomtypes = len(atomtypes)
       
  # now we can write the file
  file = open(fileName,'w')
  file.write('io_extended\n\n')
  file.write('%5d atoms\n' % numParticles)
  file.write('%5d bonds\n' % nbonds)
  file.write('%5d angles\n' % nangles)
  file.write('%5d dihedrals\n' % ndihedrals)

  file.write('%5d atom types\n' % natomtypes)
  file.write('%5d bond types\n' % nbondtypes)
  file.write('%5d angle types\n' % nangletypes)
  file.write('%5d dihedral types\n' % ndihedraltypes)
  
  file.write('%.15f %.15f xlo xhi\n' % (0.0, box_x))
  file.write('%.15f %.15f ylo yhi\n' % (0.0, box_y))
  file.write('%.15f %.15f zlo zhi\n' % (0.0, box_z))
  
  file.write('\nAtoms\n\n');
  pid = 0
  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
      particle = system.storage.getParticle(pid)
      if folded:
        xpos   = particle.pos.x
        ypos   = particle.pos.y
        zpos   = particle.pos.z
      else:
        p = system.bc.getUnfoldedPosition(particle.pos, particle.imageBox)
        xpos   = p[0]
        ypos   = p[1]
        zpos   = p[2]
      type   = particle.type
      st = "%d %d %.15f %.15f %.15f\n"%(pid, type, xpos, ypos, zpos)
      file.write(st)
    pid   += 1

  # velocities are written in the same order as coordinates, thus it does not need ID.
  if writeVelocities:
    file.write('\nVelocities\n\n');
    pid   = 0
    while pid <= maxParticleID:
      if system.storage.particleExists(pid):
          particle = system.storage.getParticle(pid)
          xvel   = particle.v[0]
          yvel   = particle.v[1]
          zvel   = particle.v[2]
          st = "%.12f %.12f %.12f\n"%(xvel, yvel, zvel)
          file.write(st)
          pid   += 1
      else:
          pid   += 1

  if nbonds > 0:
    file.write('\nBonds\n\n')
    bn = 0
    for i in xrange(len(bonds)):
      for j in xrange(len(bonds[i])):
        file.write('%d %d %d %d\n' % (bn, i, bonds[i][j][0], bonds[i][j][1]))
        bn += 1

  if nangles > 0:
    file.write('\nAngles\n\n')
    an = 0
    for i in xrange(len(angles)):
      for j in xrange(len(angles[i])):
        file.write('%d %d %d %d %d\n' % (an, i, angles[i][j][1], angles[i][j][0], angles[i][j][2]))
        an += 1

  if ndihedrals > 0:
    file.write('\nDihedrals\n\n')
    dn = 0
    for i in xrange(len(dihedrals)):
      for j in xrange(len(dihedrals[i])):
        file.write('%d %d %d %d %d %d\n' % (dn, i, dihedrals[i][j][0], dihedrals[i][j][1], dihedrals[i][j][2], dihedrals[i][j][3]))
        dn += 1
  
  file.close()
  

def read(fileName, readVelocities=False):

  f = open(fileName)
  line = f.readline() # comment line
  while not 'atoms' in line: #skip possible blank line
    line = f.readline()
  num_particles = int(line.split()[0])
  num_bonds = int(f.readline().split()[0])
  num_angles = int(f.readline().split()[0])
  num_dihedrals = int(f.readline().split()[0])
 
  # find and store size of box
  line = ''
  while not 'xlo' in line:
    line = f.readline()
  xmin, xmax = map(float, line.split()[0:2])
  ymin, ymax = map(float, f.readline().split()[0:2])
  zmin, zmax = map(float, f.readline().split()[0:2])
  Lx = xmax - xmin
  Ly = ymax - ymin
  Lz = zmax - zmin

  # find and store coordinates
  line = ''
  while not 'Atoms' in line:
    line = f.readline()
  line = f.readline()

  p_ids = []
  p_types = []
  poss = []
  for i in xrange(num_particles):
    k, kk, rx, ry, rz = map(float, f.readline().split()[0:])
    p_ids.append(int(k))
    p_types.append(int(kk))
    poss.append(espressopp.Real3D(rx, ry, rz))

  vels = []
  if(readVelocities):
    # find and store velocities
    line = ''
    while not 'Velocities' in line:
      line = f.readline()
    line = f.readline() # blank line
    for i in xrange(num_particles):
      vx_, vy_, vz_ = map(float, f.readline().split()[0:])
      vels.append(espressopp.Real3D(vx_, vy_, vz_))

  bonds = []
  if(num_bonds != 0):
    # find and store bonds
    line = ''
    while not 'Bonds' in line:
      line = f.readline()
    line = f.readline()
    for i in xrange(num_bonds):
      bond_id, bond_type, pid1, pid2 = map(int, f.readline().split())
      bonds.append([bond_type, (pid1, pid2)])

  angles = []
  if(num_angles != 0):
    # find and store angles
    line = ''
    while not 'Angles' in line:
      line = f.readline()
    line = f.readline()
    for i in xrange(num_angles):
      angle_id, angle_type, pid1, pid2, pid3 = map(int, f.readline().split())
      angles.append([angle_type, (pid1, pid2, pid3)])


  dihedrals = []
  if(num_dihedrals != 0):
    # find and store angles
    line = ''
    while not 'Dihedrals' in line:
      line = f.readline()
    line = f.readline()
    for i in xrange(num_dihedrals):
      dihedral_id, dihedral_type, pid1, pid2, pid3, pid4 = map(int, f.readline().split())
      dihedrals.append([dihedral_type, (pid1, pid2, pid3, pid4)])

  f.close()

  return Lx, Ly, Lz, p_ids, p_types, poss, vels, bonds, angles, dihedrals
