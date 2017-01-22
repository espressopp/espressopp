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
**************************
lammps - read lammps files
**************************

This Python module allows one to use a LAMMPS data file as the
input to an ESPResSo++ simulation.

"""

import espressopp

def read(fin):

  f = open(fin)
  line = f.readline() # comment line
  while not 'atoms' in line: #skip possible blank line
    line = f.readline()
  num_particles = int(line.split()[0])
  num_bonds = int(f.readline().split()[0])
  num_angles = int(f.readline().split()[0])
  num_dihedrals = int(f.readline().split()[0])
  line = f.readline() # impropers
  line = f.readline() # blank line
  line = f.readline() # atom types and maybe the word "velocities"
  num_types = int(line.split()[0])
  velocities = True if 'velocities' in line else False #TODO fix this? why should there be the velocity keyword?
 
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

  if(num_types == 1):
    rstart = 3
    if(num_bonds == 0): rstart = 2
    
    x = []
    y = []
    z = []
    for i in xrange(num_particles):
      rx, ry, rz = map(float, f.readline().split()[rstart:])
      x.append(rx)
      y.append(ry)
      z.append(rz)
  else:
    p_type = []
    q = []
    x = []
    y = []
    z = []
    for i in xrange(num_particles):
      k, rq, rx, ry, rz = map(float, f.readline().split()[2:])
      p_type.append(int(k))
      q.append(rq)
      x.append(rx)
      y.append(ry)
      z.append(rz)

  if(num_bonds != 0):
    # find and store bonds
    line = ''
    while not 'Bonds' in line:
      line = f.readline()
    line = f.readline()
    bonds = []
    for i in xrange(num_bonds):
      bond_id, bond_type, pid1, pid2 = map(int, f.readline().split())
      bonds.append((pid1, pid2))

  if(num_angles != 0):
    # find and store angles
    line = ''
    while not 'Angles' in line:
      line = f.readline()
    line = f.readline()
    angles = []
    for i in xrange(num_angles):
      angle_id, angle_type, pid1, pid2, pid3 = map(int, f.readline().split())
      angles.append((pid1, pid2, pid3))


  if(num_dihedrals != 0):
    # find and store angles
    line = ''
    while not 'Dihedrals' in line:
      line = f.readline()
    line = f.readline()
    dihedrals = []
    for i in xrange(num_dihedrals):
      dihedral_id, dihedral_type, pid1, pid2, pid3, pid4 = map(int, f.readline().split())
      dihedrals.append((pid1, pid2, pid3, pid4))

  if(velocities):
    # find and store velocities
    line = ''
    while not 'Velocities' in line:
      line = f.readline()
    line = f.readline() # blank line
    vx = []
    vy = []
    vz = []
    for i in xrange(num_particles):
      vx_, vy_, vz_ = map(float, f.readline().split()[1:])
      vx.append(vx_)
      vy.append(vy_)
      vz.append(vz_)


  f.close()


  if(num_types != 1):
    return p_type, bonds, angles, q, x, y, z, Lx, Ly, Lz

  if(num_bonds == 0 and num_angles == 0 and num_dihedrals == 0 and not velocities):
    return x, y, z, Lx, Ly, Lz
  if(num_bonds == 0 and num_angles == 0 and num_dihedrals == 0 and velocities):
    return x, y, z, Lx, Ly, Lz, vx, vy, vz
  elif(num_bonds != 0 and num_angles == 0 and num_dihedrals == 0):
    return bonds, x, y, z, Lx, Ly, Lz
  elif(num_bonds != 0 and num_angles != 0 and num_dihedrals == 0):
    return bonds, angles, x, y, z, Lx, Ly, Lz
  else:
    return bonds, angles, dihedrals, x, y, z, Lx, Ly, Lz

def read_charmm(fin):

  f = open(fin)
  line = f.readline() # comment line
  while not 'atoms' in line: #skip possible blank line
    line = f.readline()
  num_particles = int(line.split()[0])
  num_bonds = int(f.readline().split()[0])
  num_angles = int(f.readline().split()[0])
  num_dihedrals = int(f.readline().split()[0])
  num_impropers = int(f.readline().split()[0])
  
  line = f.readline() # blank line
  
  line = f.readline() # atom types 
  num_types = int(line.split()[0])
  line = f.readline() # bond types 
  num_bond_types = int(line.split()[0])
  line = f.readline() # angle types 
  num_angle_types = int(line.split()[0])
  line = f.readline() # dihedral types 
  num_dihedral_types = int(line.split()[0])
  line = f.readline() # impropers types 
  num_improper_types = int(line.split()[0])

  velocities = True if 'velocities' in line else False #TODO fix this? why should there be the velocity keyword?
 
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

  #find and store masses
  line = ''
  while not 'Masses' in line:
    line = f.readline()
  line = f.readline()
  masses = []
  dumm = 0.0
  masses.append(dumm)
  for i in xrange(num_types):
      rmass = float(f.readline().split()[1])
      masses.append(rmass) 

#find and store LJ param
  line = ''
  while not 'Pair Coeffs' in line:
    line = f.readline()
  line = f.readline()
  epsilon = []
  sigma = []
  dumm = 0.0
  epsilon.append(dumm)
  sigma.append(dumm)
  for i in xrange(num_types):
      repsilon, rsigma = map(float, f.readline().split()[1:3])
      epsilon.append(repsilon) 
      sigma.append(rsigma)
      

  # find and store coordinates
  line = ''
  while not 'Atoms' in line:
    line = f.readline()
  line = f.readline()

  if(num_types == 1):
    rstart = 3
    if(num_bonds == 0): rstart = 2
    
    p_type = []
    q = []
    x = []
    y = []
    z = []
    for i in xrange(num_particles):
      rx, ry, rz = map(float,f.readline().split()[rstart:])
      x.append(rx)
      y.append(ry)
      z.append(rz)
  else:
    p_type = []
    q = []
    x = []
    y = []
    z = []
    for i in xrange(num_particles):
      k, rq, rx, ry, rz = map(float, f.readline().split()[2:])
      p_type.append(int(k))
      q.append(rq)
      x.append(rx)
      y.append(ry)
      z.append(rz)

  if(num_bonds != 0):
    # find and store bond coeff
    line = ''
    while not 'Bond Coeffs' in line:
      line = f.readline()
    line = f.readline()
    K = []
    r0 = []
    dumm = 0.0
    K.append(dumm)
    r0.append(dumm)
    for i in xrange(num_bond_types):
        rK, rr0 = map(float, f.readline().split()[1:3])
        K.append(rK) 
        r0.append(rr0)



  if(num_bonds != 0):
    # find and store bonds
    line = ''
    while not 'Bonds' in line:
      line = f.readline()
    line = f.readline()
    bonds = []
    bonds_type_arr = []
    for i in xrange(num_bonds):
      bond_id, bond_type, pid1, pid2 = map(int, f.readline().split())
      bonds.append((pid1, pid2))
      bonds_type_arr.append(bond_type)

  if(num_angles != 0):
    # find and store angle coeff
    line = ''
    while not 'Angle Coeffs' in line:
      line = f.readline()
    line = f.readline()
    Kt = []
    t0 = []
    dumm = 0.0
    Kt.append(dumm)
    t0.append(dumm)
    for i in xrange(num_bond_types):
        rKt, rt0 = map(float, f.readline().split()[1:3])
        Kt.append(rKt) 
        t0.append(rt0)

  if(num_angles != 0):
    # find and store angles
    line = ''
    while not 'Angles' in line:
      line = f.readline()
    line = f.readline()
    angles = []
    angles_type_arr = []
    for i in xrange(num_angles):
      angle_id, angle_type, pid1, pid2, pid3 = map(int, f.readline().split())
      angles.append((pid1, pid2, pid3))
      angles_type_arr.append(angle_type)

  if(num_dihedrals != 0):
  # find and store dihedrals coeff
    line = ''
    while not 'Dihedral Coeffs' in line:
      line = f.readline()
    line = f.readline()
    Kdh = []
    ndh = []
    ph0 = []
    dumm = 0.0
    Kdh.append(dumm)
    ndh.append(dumm)
    ph0.append(dumm)
    for i in xrange(num_bond_types):
        rKdh, rndh, rph0 = map(float, f.readline().split()[1:4])
        Kdh.append(rKdh) 
        ndh.append(rndh)
        ph0.append(rph0)

  if(num_dihedrals != 0):
    # find and store dihedrals
    line = ''
    while not 'Dihedrals' in line:
      line = f.readline()
    line = f.readline()
    dihedrals = []
    dihedrals_type_arr = []
    for i in xrange(num_dihedrals):
      dihedral_id, dihedral_type, pid1, pid2, pid3, pid4 = map(int, f.readline().split())
      dihedrals.append((pid1, pid2, pid3, pid4))
      dihedrals_type_arr.append(dihedral_type)

  if(velocities):
    # find and store velocities
    line = ''
    while not 'Velocities' in line:
      line = f.readline()
    line = f.readline() # blank line
    vx = []
    vy = []
    vz = []
    for i in xrange(num_particles):
      vx_, vy_, vz_ = map(float, f.readline().split()[1:])
      vx.append(vx_)
      vy.append(vy_)
      vz.append(vz_)


  f.close()




  if(num_bonds == 0 and num_angles == 0 and num_dihedrals == 0 and not velocities):
    return p_type, masses, epsilon, sigma, q, x, y, z, Lx, Ly, Lz
  if(num_bonds == 0 and num_angles == 0 and num_dihedrals == 0 and velocities):
    return p_type, masses, epsilon, sigma, q, x, y, z, Lx, Ly, Lz, vx, vy, vz
  if(num_bonds != 0 and num_angles == 0 and num_dihedrals == 0 and not velocities):
    return p_type, masses, epsilon, sigma, K, r0, bonds, bonds_type_arr, q, x, y, z, Lx, Ly, Lz
  if(num_bonds != 0 and num_angles == 0 and num_dihedrals == 0 and velocities):
    return p_type, masses, epsilon, sigma, K, r0, bonds, bonds_type_arr, q, x, y, z, Lx, Ly, Lz, vx, vy, vz
  if(num_bonds != 0 and num_angles != 0 and num_dihedrals == 0 and not velocities):
    return p_type, masses, epsilon, sigma, K, r0, bonds, bonds_type_arr, Kt, t0, angles, angles_type_arr, q, x, y, z, Lx, Ly, Lz
  if(num_bonds != 0 and num_angles != 0 and num_dihedrals == 0 and velocities):
    return p_type, masses, epsilon, sigma, K, r0, bonds, bonds_type_arr, Kt, t0, angles, angles_type_arr, q, x, y, z, Lx, Ly, Lz, vx, vy, vz
  if(num_bonds != 0 and num_angles != 0 and num_dihedrals != 0 and not velocities):
    return p_type, masses, epsilon, sigma, K, r0, bonds, bonds_type_arr, Kt, t0, angles, angles_type_arr, Kdh, ndh, ph0, dihedrals, dihedrals_type_arr, q, x, y, z, Lx, Ly, Lz
  if(num_bonds != 0 and num_angles != 0 and num_dihedrals != 0 and velocities):
    return p_type, masses, epsilon, sigma, K, r0, bonds, bonds_type_arr, Kt, t0, angles, angles_type_arr, Kdh, ndh, ph0, dihedrals, dihedrals_type_arr, q, x, y, z, Lx, Ly, Lz, vx, vy, vz

def write(fout, system, writeVelocities=False):  
    
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
  file = open(fout,'w')
  file.write('LAMMPS\n\n')
  file.write('%5d atoms\n' % numParticles)
  file.write('%5d bonds\n' % nbonds)
  file.write('%5d angles\n' % nangles)
  file.write('%5d dihedrals\n' % ndihedrals)
  #impropers are not supported yet
  file.write('%5d impropers\n\n' % 0)

  file.write('%5d atom types\n' % natomtypes)
  file.write('%5d bond types\n' % nbondtypes)
  file.write('%5d angle types\n' % nangletypes)
  file.write('%5d dihedral types\n' % ndihedraltypes)
  file.write('%5d improper types\n\n' % 0)
  
  file.write('%.4f %.4f xlo xhi\n' % (0.0, box_x))
  file.write('%.4f %.4f ylo yhi\n' % (0.0, box_y))
  file.write('%.4f %.4f zlo zhi\n' % (0.0, box_z))
  
  file.write('\nAtoms\n\n');
  pid   = 0
  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
        particle = system.storage.getParticle(pid)
        p = system.bc.getUnfoldedPosition(particle.pos, particle.imageBox)
        xpos   = p[0]
        ypos   = p[1]
        zpos   = p[2]
        type   = particle.type
        # we don't support molecule tags yet
        molecule_tag = (pid-1) / 200 + 1
        st = "%d %d %d %.3f %.3f %.3f\n"%(pid, molecule_tag, type+1, xpos, ypos, zpos)
        file.write(st)
        pid   += 1
    else:
        pid   += 1

  if writeVelocities:
    file.write('\nVelocities\n\n');
    pid   = 0
    while pid <= maxParticleID:
      if system.storage.particleExists(pid):
          particle = system.storage.getParticle(pid)
          xvel   = particle.v[0]
          yvel   = particle.v[1]
          zvel   = particle.v[2]
          st = "%d %15.10f %15.10f %15.10f\n"%(pid, xvel, yvel, zvel)
          file.write(st)
          pid   += 1
      else:
          pid   += 1

  if nbonds > 0:      
    file.write('\nBonds\n\n')
    bn = 1
    for i in xrange(len(bonds)):
        for j in xrange(len(bonds[i])):
            file.write('%d %d %d %d\n' % (bn, i+1, bonds[i][j][0], bonds[i][j][1]))
            bn += 1

  if nangles > 0:
    file.write('\nAngles\n\n')
    an = 1
    for i in xrange(len(angles)):
        for j in xrange(len(angles[i])):
            file.write('%d %d %d %d %d\n' % (an, i+1, angles[i][j][1], angles[i][j][0], angles[i][j][2]))
            an += 1

  if ndihedrals > 0:
    file.write('\nDihedrals\n\n')
    dn = 1
    for i in xrange(len(dihedrals)):
        for j in xrange(len(dihedrals[i])):
            file.write('%d %d %d %d %d %d\n' % (dn, i+1, dihedrals[i][j][0], dihedrals[i][j][1], dihedrals[i][j][2], dihedrals[i][j][3]))
            dn += 1
  
  file.close()
 
