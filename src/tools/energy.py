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
***************
** energy.py **
***************

helper functions for calculating certain types of energy on the python level

.. function:: espressopp.tools.getSelfExclEnergyReactionField(system,exclusions,prefactor,epsilon1,epsilon2,rc,pidlist=[],nParticles=0)

  When calculating the electrostatic energy using the reaction field approach, a correction should also be applied to all excluded atom pairs, including self pairs. This function is for calculating the additional contribution to the electrostatic energy from those 'self' and 'excluded' interactions.

 If a list of pids is supplied, the self and exclusion energy is calculated only for particles in the list. Otherwise, if the variable nParticles is supplied, the calculation is for all particles with pid 1 to nParticles inclusive.

  Note: this is not generalised reaction field (i.e. the generalised reaction field parameter kappa is set to zero).

  :param system: espressopp system
  :type system: System object
  :param exclusions: list of all pairs of particles with non-bonded exclusions in the system
  :type exclusions: list of 2-tuples
  :param prefactor: electrostatic prefactor in the desired units, e.g. 138.935485 kJ mol^-1 e^-2
  :type prefactor: real
  :param epsilon1: epsilon1 of reaction field potential, e.g. 1
  :type epsilon1: real
  :param epsilon2: epsilon2 of reaction field potential, e.g. approximately 80 for water
  :type epsilon2: real
  :param rc: interaction cutoff
  :type rc: real
  :param pidlist: (optional) list of pids for which the self and exclusion energy is calculated
  :type pidlist: list of integers
  :param nParticles: (optional) if pidlist is not supplied, the self and exclusion energy is calculated for all particles with pids from 1 to nParticles inclusive
  :type nParticles: integer
"""

import math
import sys

def getSelfExclEnergyReactionField(system,exclusions,prefactor,epsilon1,epsilon2,rc,pidlist=[],nParticles=0):

  #initialisation
  krc = 0 #to convert this function for generalized reac field, one would also need to pass a parameter kappa, krc=kappa*rc
  tmp1 = (epsilon1 - 4.0*epsilon2) * (1.0 + krc) - 2.0*epsilon2 * krc*krc;
  tmp2 = (epsilon1 + 2.0*epsilon2) * (1.0 + krc) + epsilon2 * krc*krc;
  rc3 = math.pow(rc,3);
  rc2 = math.pow(rc,2);
  B1 = tmp1/tmp2;
  B1 = (1.0+B1) / rc3;
  B1_half = B1/2.0;
  crf = 3*epsilon2/(rc*(2*epsilon2+1.0));
  prefactor/=epsilon1;

  # energy correction for self pairs
  selfEnergy = 0.0
  if pidlist == []: #if no list of pids supplied, calculate energy for all particles
    if nParticles == 0:
      print 'Error! In getSelfExclEnergyReactionField, you supplied neither a list of particles indices (pidlist) nor a maximum particle id (nParticles), so no self or excluded energy is calculated'
      quit()
    pidlist = xrange(1,nParticles+1)
  for pid in pidlist:
    part = system.storage.getParticle(pid)
    selfEnergy += part.q*part.q*prefactor*-1*crf

  # energy correction for exclusion pairs
  exclEnergy = 0.0
  for (pid1,pid2) in exclusions:
    if (pid1 in pidlist) and (pid2 in pidlist):
      part1 = system.storage.getParticle(pid1)
      part2 = system.storage.getParticle(pid2)
      dist = part1.pos - part2.pos
      distSqr = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]
      if distSqr>rc2:
        print 'Error, PBC problem in getSelfEnergyReactionField with particles ',pid1,pid2
        quit()
      exclEnergy += part1.q*part2.q*prefactor*-1.0*(B1_half*distSqr+crf)

  return exclEnergy + 0.5*selfEnergy
