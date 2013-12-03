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
******************************************
**decomp.py** - Auxiliary python functions
******************************************


*  `nodeGrid(n)`:

    It determines how the processors are distributed and how the cells are arranged.
    `n` - number of processes 

*  `cellGrid(box_size, node_grid, rc, skin)`:

    It returns an appropriate grid of cells.
    
*  `tuneSkin(system, integrator, minSkin=0.01, maxSkin=1.2, precision=0.001)`:

    It tunes the skin size for the current system
    
*  `printTimeVsSkin(system, integrator, minSkin=0.01, maxSkin=1.5, skinStep = 0.01)`:
    
    It prints time of running versus skin size in the range [minSkin, maxSkin] with
    the step skinStep
"""


import sys
import espresso

from espresso import Int3D
from espresso.Exceptions import Error

import math
import time

# need to factorize n to find optimum nodeGrid
def nodeGrid(n):
	ijkmax = 3*n*n + 1
	d1 = 1
	d2 = 1
	d3 = 1
	for i in range(1,n+1):
		for j in range(i,n+1):
			for k in range(j,n+1):
				if (i*j*k == n) and (i*i + j*j + k*k < ijkmax):
					d1 = i
					d2 = j
					d3 = k
					ijkmax = i*i + j*j + k*k
	return Int3D(d1,d2,d3)

def cellGrid(box_size, node_grid, rc, skin):
  rc_skin = rc + skin
  if rc_skin == 0:
    raise Error("interaction range (cutoff + skin) must be larger than 0")
  if (node_grid[0]<=0 or node_grid[1]<=0 or node_grid[2]<=0):
    raise Error("invalid node grid %s" % str(node_grid))
  ix = box_size[0] / (rc_skin * node_grid[0])
  if ix < 1:
    raise Error("local box size in direction 0 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime." % (ix, rc_skin))
  iy = box_size[1] / (rc_skin * node_grid[1])
  if iy < 1:
    raise Error("local box size in direction 1 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime." % (iy, rc_skin))
  iz = box_size[2] / (rc_skin * node_grid[2])
  if iz < 1:
    raise Error("local box size in direction 2 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime." % (iz, rc_skin))
  
  return Int3D(ix, iy, iz)

def tuneSkin(system, integrator, minSkin=0.01, maxSkin=1.5, precision=0.001, printInfo=True):
  if printInfo:
    print 'The tuning is started. It can take some time depending on your system.'
  
  fi = (1.0+math.sqrt(5.0))/2.0 # golden ratio
  
  npart = espresso.analysis.NPart(system).compute()
  
  # this is an empirical formula in order to get the appropriate number of steps
  nsteps = int( espresso.MPI.COMM_WORLD.size * 1000000.0 / float(npart) )
  
  if printInfo:
    print 'CellGrid before tuning: ', system.storage.getCellGrid()
    sys.stdout.write('\nSteps     = %d\n' % nsteps)
    sys.stdout.write('Precision = %g\n' % precision)
    sys.stdout.write('It runs till deltaSkin<precision\n')
  
  if printInfo:
    prnt_format1 = '\n%9s %10s %10s %10s %14s\n'
    sys.stdout.write(prnt_format1 % ('time1: ',' time2: ',' skin1: ', ' skin2: ', ' deltaSkin: '))
  
  while (maxSkin-minSkin>=precision):
    skin1 = maxSkin - (maxSkin-minSkin)/fi
    skin2 = minSkin + (maxSkin-minSkin)/fi

    system.skin = skin1
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(nsteps)
    end_time = time.time()
    time1 = end_time - start_time

    system.skin = skin2
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(nsteps)
    end_time = time.time()
    time2 = end_time - start_time

    if(time1>time2):
      minSkin = skin1
    else:
      maxSkin = skin2
      
    if printInfo:
      prnt_format2 = '%7.3f %10.3f %11.4f %10.4f %12.6f\n'
      sys.stdout.write(prnt_format2 % (time1, time2, minSkin, maxSkin, (maxSkin-minSkin)) )
    
      sys.stdout.write('\nNew skin: %g\n' % system.skin)
      sys.stdout.write('\nNew cell grid: %s\n' % system.storage.getCellGrid())
  
  system.skin = (maxSkin+minSkin)/2.0
  system.storage.cellAdjust()
  
  return (maxSkin+minSkin)/2.0

def printTimeVsSkin(system, integrator, minSkin=0.01, maxSkin=1.5, skinStep = 0.005):
  npart = espresso.analysis.NPart(system).compute()
  # this is an empirical formula in order to get the appropriate number of steps
  nsteps = int( espresso.MPI.COMM_WORLD.size * 20000000.0 / float(npart) )
  
  print '      Calculations is started. It will print out the dependece of time of \n\
      running of %d steps on the skin size into the file \'timeVSskin.dat\'.\n\
      The range of skin sizes is [%g, %g], skin step is %g. It can take some \n\
      time depending on your system.' % (nsteps, minSkin, maxSkin, skinStep)
  
  curSkin = minSkin
  
  fmt2 = ' %8.4f %8.4f\n' # format for writing the data
  nameFile = 'timeVSskin.dat'
  resFile = open (nameFile, 'w')

  count = 0
  while (curSkin < maxSkin):
    system.skin = curSkin
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(nsteps)
    end_time = time.time()
    time1 = end_time - start_time
    
    resFile.write(fmt2 % ( system.skin, time1 ))

    count = count +1
    if (count == 20):
      print 'skin: ', system.skin
      count = 0
    
    curSkin = curSkin + skinStep
  
  resFile.close()
  
  return
  
