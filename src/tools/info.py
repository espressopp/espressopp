import espresso

def getAllParticles(system, *properties):
  'returns a list of all particle properties of all particles of the system (currently no atomistic AdResS particles are included)'
  allParticles = []
  maxParticleID = int(espresso.analysis.MaxPID(system).compute())
  pid   = 0
  while pid <= maxParticleID:
    particle = system.storage.getParticle(pid)
    part = []
    if particle.pos:    
      for val in properties:
        if   val.lower() == "id"    : part.append(particle.id)
        elif val.lower() == "pos"   : part.append(particle.pos)
        elif val.lower() == "type"  : part.append(particle.type)
        elif val.lower() == "mass"  : part.append(particle.mass)
        elif val.lower() == "v"     : part.append(particle.v)
        elif val.lower() == "f"     : part.append(particle.f)
        elif val.lower() == "q"     : part.append(particle.q)
        elif val.lower() == "adrat" : part.append(particle.adrat)
        else: raise "unknown particle property: %s"%val
      allParticles.append(part)         
      pid   += 1
    else:
      pid   += 1    
  return allParticles
 
  
def getAllBonds(system):
  'return all bonds of the system (currently only FixedPairLists are supported)'
  bonds = []
  nInteractions = system.getNumberOfInteractions()
  for i in range(nInteractions):
      if system.getInteraction(i).isBonded():
          try:
              FixedPairList = system.getInteraction(i).getFixedPairList().getBonds()
              j = 0
              while j < len(FixedPairList):
                  fplb = FixedPairList[j]
                  k = 0
                  while k < len(fplb):
                      bonds.append(fplb[k])
                      k += 1
                  j += 1
          except:
              pass
  return bonds