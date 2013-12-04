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


def gaussian(T, N, zero_momentum=True, seed=7654321):

  """This Python module generates initial particle velocities
  with temperature T according to a Maxwell-Boltzmann distribution."""
  # TODO: account for mass

  import random

  sqrtT = T**0.5
  random.seed(seed)
  vx = []
  vy = []
  vz = []
  for i in range(N):
    vx.append(sqrtT * random.gauss(0.0, 1.0))
    vy.append(sqrtT * random.gauss(0.0, 1.0))
    vz.append(sqrtT * random.gauss(0.0, 1.0))

  if(zero_momentum):
    # remove net momentum
    sumvx = 0.0
    sumvy = 0.0
    sumvz = 0.0
    for vx_, vy_, vz_ in zip(vx, vy, vz):
      sumvx += vx_
      sumvy += vy_
      sumvz += vz_
    sumvx = sumvx / N
    sumvy = sumvy / N
    sumvz = sumvz / N
    for i in range(N):
      vx[i] = vx[i] - sumvx
      vy[i] = vy[i] - sumvy
      vz[i] = vz[i] - sumvz

  return vx, vy, vz
