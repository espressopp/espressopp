"""Python functions to determine how the processors are distributed
   and how the cells are arranged."""

import espresso
from espresso import Int3D

def nodeGrid(p):
  return Int3D(1, 1, p)

def cellGrid(box_size, node_grid, rc, skin):
  ix = max(1, box_size[0] * node_grid[0] / (rc + skin))
  iy = max(1, box_size[1] * node_grid[1] / (rc + skin))
  iz = max(1, box_size[2] * node_grid[2] / (rc + skin))
  return Int3D(ix, iy, iz)
