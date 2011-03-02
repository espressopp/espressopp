"""Python functions to determine how the processors are distributed
   and how the cells are arranged."""

import espresso
from espresso import Int3D

# need to factorize n to find optimum nodeGrid
def nodeGrid(n):
  ijkmax = 3*n*n + 1
  d1 = 1
  d2 = 2
  d3 = n
  for i in range(1,n):
    for j in range(i,n):
      for k in range(j,n):
        if (i*j*k == n) and (i*i + j*j + k*k < ijkmax):
          d1 = k
          d2 = j
          d3 = i
          ijkmax = i*i + j*j + k*k
  return Int3D(d1,d2,d3)

def cellGrid(box_size, node_grid, rc, skin):
  rc_skin = rc + skin
  ix = max(1, box_size[0] / (rc_skin * node_grid[0]))
  iy = max(1, box_size[1] / (rc_skin * node_grid[1]))
  iz = max(1, box_size[2] / (rc_skin * node_grid[2]))
  return Int3D(ix, iy, iz)
