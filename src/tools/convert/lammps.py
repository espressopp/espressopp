# -*- coding: utf-8 -*-

import espresso

"""This Python module allows one to use a LAMMPS data file as the
   input to an ESPResSo++ simulation."""

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
    for i in range(num_particles):
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
    for i in range(num_particles):
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
    for i in range(num_bonds):
      bond_id, bond_type, pid1, pid2 = map(int, f.readline().split())
      bonds.append((pid1, pid2))

  if(num_angles != 0):
    # find and store angles
    line = ''
    while not 'Angles' in line:
      line = f.readline()
    line = f.readline()
    angles = []
    for i in range(num_angles):
      angle_id, angle_type, pid1, pid2, pid3 = map(int, f.readline().split())
      angles.append((pid1, pid2, pid3))


  if(num_dihedrals != 0):
    # find and store angles
    line = ''
    while not 'Dihedrals' in line:
      line = f.readline()
    line = f.readline()
    dihedrals = []
    for i in range(num_dihedrals):
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
    for i in range(num_particles):
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

def write(fout, system, writeVelocities=False):  
    
  # first collect all the information that we need to write into the file
  numParticles  = int(espresso.analysis.NPart(system).compute())
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
  for i in range(nInteractions):
      bT = system.getInteraction(i).bondType()
      if   bT == espresso.interaction.Pair:
             nbondtypes += 1
             bl  = system.getInteraction(i).getFixedPairList().getBonds()
             bln = []
             for j in range(len(bl)):
               bln.extend(bl[j])
             bonds.append(bln)
      elif bT == espresso.interaction.Angular:
             nangletypes += 1
             an  = system.getInteraction(i).getFixedTripleList().getTriples()
             ann = []
             for j in range(len(an)):
               ann.extend(an[j]) 
             angles.append(ann)
      elif bT == espresso.interaction.Dihedral:
             ndihedraltypes += 1
             di  = system.getInteraction(i).getFixedQuadrupleList().getQuadruples()
             din = []
             for j in range(len(di)):
               din.extend(di[j])   
             dihedrals.append(din)
  
  nbonds = 0
  for i in range(len(bonds)):
      nbonds += len(bonds[i])
  nangles = 0
  for i in range(len(angles)):
      nangles += len(angles[i])
  ndihedrals = 0
  for i in range(len(dihedrals)):
      ndihedrals += len(dihedrals[i])
      
  atomtypes = []
  maxParticleID = int(espresso.analysis.MaxPID(system).compute())
  pid   = 0
  while pid <= maxParticleID:
    particle = system.storage.getParticle(pid)
    if particle.pos:
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
    particle = system.storage.getParticle(pid)
    if particle.pos:
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
      particle = system.storage.getParticle(pid)
      if particle.pos:
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
    for i in range(len(bonds)):
        for j in range(len(bonds[i])):
            file.write('%d %d %d %d\n' % (bn, i+1, bonds[i][j][0], bonds[i][j][1]))
            bn += 1

  if nangles > 0:
    file.write('\nAngles\n\n')
    an = 1
    for i in range(len(angles)):
        for j in range(len(angles[i])):
            file.write('%d %d %d %d %d\n' % (an, i+1, angles[i][j][1], angles[i][j][0], angles[i][j][2]))
            an += 1

  if ndihedrals > 0:
    file.write('\nDihedrals\n\n')
    dn = 1
    for i in range(len(dihedrals)):
        for j in range(len(dihedrals[i])):
            file.write('%d %d %d %d %d %d\n' % (dn, i+1, dihedrals[i][j][0], dihedrals[i][j][1], dihedrals[i][j][2], dihedrals[i][j][3]))
            dn += 1
  
  file.close()
 
