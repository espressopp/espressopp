"""Python functions to determine how the processors are distributed
   and how the cells are arranged."""

from sys import exit
from espresso import Int3D
from espresso.Exceptions import Error

from espresso.tools import timers

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

def tuneSkin(system, integrator, precision):
  print 'The tuning is started. It can take some time. It depends on your system.'
  nodeGrid = system.storage.getNodeGrid()
  box      = system.bc.boxL
  minNodeLength = min( box[0]/nodeGrid[0], min( box[1]/nodeGrid[1], box[2]/nodeGrid[2] )  )
  maximumCut = system.maxCutoff # biggest value
  
  fi = (1.0+math.sqrt(5.0))/2.0 # golden ratio
  
  # initial data
  minSkin = 0.01 # lowest value
  maxSkin = 1.0 #minNodeLength - maximumCut
  
  while (maxSkin-minSkin>=precision):
    sk1 = maxSkin - (maxSkin-minSkin)/fi
    sk2 = minSkin + (maxSkin-minSkin)/fi

    system.skin = sk1
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(2000)
    end_time = time.time()
    time1 = end_time - start_time

    system.skin = sk2
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(2000)
    end_time = time.time()
    time2 = end_time - start_time

    if(time1>time2):
      minSkin = sk1
    else:
      maxSkin = sk2
      
    print 'time1: ', time1, '  time2: ', time2, '  deltaSkin: ', (maxSkin-minSkin)
    
  system.skin = (maxSkin+minSkin)/2.0
  system.storage.cellAdjust()
  
  return (maxSkin+minSkin)/2.0