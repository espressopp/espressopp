"""Python functions to determine how the processors are distributed
   and how the cells are arranged."""

import espresso
from espresso import Int3D

# need to factorize p to find optimum nodeGrid
def nodeGrid(p):
  return Int3D(1, 1, p)

def cellGrid(box_size, node_grid, rc, skin):
  rc_skin = rc + skin
  ix = max(1, box_size[0] / (rc_skin * node_grid[0]))
  iy = max(1, box_size[1] / (rc_skin * node_grid[1]))
  iz = max(1, box_size[2] / (rc_skin * node_grid[2]))
  return Int3D(ix, iy, iz)
