def replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=1, ydim=1, zdim=1):
  """
  Presently this routine works only for semiflexible polymers. A general
  class should be written to deal with files containing coordinates
  and topology data.

  This method takes the current configuration and replicates
  it by some number of boxes in each dimension. This may be used to
  increase the size of an equilibrated melt by a factor of 8 or more.

  xdim = ydim = zdim = 1 returns the original system not replicated.
  xdim = ydim = zdim = 2 returns the original system replicated to 8x.
  xdim = ydim = zdim = 3 returns the original system replicated to 27x.
  xdim = ydim = 1, zdim = 2 returns the original system replicated in the z-direction.
  """

  # replicate the particles
  x_replicated = x[:]
  y_replicated = y[:]
  z_replicated = z[:]
  for i in range(xdim):
    for j in range(ydim):
      for k in range(zdim):
        if(i + j + k != 0):
          for x_, y_, z_ in zip(x, y, z):
            x_replicated.append(x_ + i * Lx)
            y_replicated.append(y_ + j * Ly)
            z_replicated.append(z_ + k * Lz)

  # replicate the bonds and angles
  ct = 0
  num_particles_original = len(x)
  bonds_replicated = bonds[:]
  angles_replicated = angles[:]
  for i in range(xdim):
    for j in range(ydim):
      for k in range(zdim):
        if(i + j + k != 0):
          ct = ct + 1
          for p1, p2 in bonds:
            bonds_replicated.append((p1 + ct * num_particles_original, \
                                     p2 + ct * num_particles_original))
          for p1, p2, p3 in angles:
            angles_replicated.append((p1 + ct * num_particles_original, \
                                      p2 + ct * num_particles_original, \
                                      p3 + ct * num_particles_original))
 
  # modify the box size
  Lx = xdim * Lx
  Ly = ydim * Ly
  Lz = zdim * Lz

  return bonds_replicated, angles_replicated, x_replicated, y_replicated, z_replicated, Lx, Ly, Lz
