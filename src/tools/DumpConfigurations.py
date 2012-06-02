import espresso

def writexyz(filename, system, velocities = True, unfolded = False):
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
    particle = system.storage.getParticle(pid)
    if particle.pos:
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
