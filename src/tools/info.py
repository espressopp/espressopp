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

import espressopp

def getAllParticles(system, *properties):
  """
  returns a list of all particle properties of all particles of the system (currently no atomistic AdResS particles are included)
  """
  allParticles = []
  maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
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
  """
  return all bonds of the system (currently only FixedPairLists are supported)
  """
  bonds = []
  nInteractions = system.getNumberOfInteractions()
  for i in xrange(nInteractions):
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
