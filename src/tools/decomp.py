"""Python functions to determine how the processors are distributed
   and how the cells are arranged."""

from sys import exit
from espresso import Int3D
from espresso.Exceptions import Error

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
    raise Error("local box size in direction 0 (=%6f) is smaller than interaction range (cutoff + skin = %6f)" % (ix, rc_skin))
  iy = box_size[1] / (rc_skin * node_grid[1])
  if iy < 1:
    raise Error("local box size in direction 1 (=%6f) is smaller than interaction range (cutoff + skin = %6f)" % (iy, rc_skin))
  iz = box_size[2] / (rc_skin * node_grid[2])
  if iz < 1:
    raise Error("local box size in direction 2 (=%6f) is smaller than interaction range (cutoff + skin = %6f)" % (iz, rc_skin))
  
  return Int3D(ix, iy, iz)
