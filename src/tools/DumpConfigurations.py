import espresso

def writexyz(filename, system, velocities = True, unfolded = False, append = False):

  if append:
    file = open(filename,'a')
  else:
    file = open(filename,'w')
  numParticles  = int(espresso.analysis.NPart(system).compute())
  box_x = system.bc.boxL[0]
  box_y = system.bc.boxL[1]
  box_z = system.bc.boxL[2]
  st = "%d\n%15.10f %15.10f %15.10f\n" % (numParticles, box_x, box_y, box_z)
  file.write(st)
  maxParticleID = int(espresso.analysis.MaxPID(system).compute())
  pid   = 0
  while pid <= maxParticleID:
    if system.storage.particleExists(pid):
        particle = system.storage.getParticle(pid)
        if unfolded == False:
          xpos   = particle.pos[0]
          ypos   = particle.pos[1]
          zpos   = particle.pos[2]
        else:
          unfoldedpos = system.bc.getUnfoldedPosition(particle.pos, particle.imageBox)
          xpos   = unfoldedpos[0]
          ypos   = unfoldedpos[1]
          zpos   = unfoldedpos[2]
        xvel   = particle.v[0]
        yvel   = particle.v[1]
        zvel   = particle.v[2]
        type   = particle.type
        if velocities:
          st = "%d %d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n"%(pid, type, xpos, ypos, zpos, xvel, yvel, zvel)
        else:
          st = "%d %d %15.10f %15.10f %15.10f\n"%(pid, type, xpos, ypos, zpos)
        file.write(st)
        pid   += 1
    else:
        pid   += 1
  
  file.close()

def readxyz(filename):
  file = open(filename)
  line = file.readline()
  num_particles = int(line.split()[0])
  line = file.readline()
  Lx = float(line.split()[0])
  Ly = float(line.split()[1])
  Lz = float(line.split()[2])
  pid  = []
  type = []
  xpos = []
  ypos = []
  zpos = []
  xvel = []
  yvel = []
  zvel = []
  for i in range(num_particles):
    line = file.readline().split()
    pid.append(int(line[0]))
    type.append(int(line[1]))
    xpos.append(float(line[2]))
    ypos.append(float(line[3]))
    zpos.append(float(line[4]))
    if len(line) > 5:
      xvel.append(float(line[5]))
      yvel.append(float(line[6]))
      zvel.append(float(line[7]))
    else:
      xvel.append(0.0)
      yvel.append(0.0)
      zvel.append(0.0)
  return pid, type, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz
  file.close()
