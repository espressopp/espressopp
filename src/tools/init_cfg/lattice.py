"""This Python module initializes particles on the sites
   of a simple cubic lattice. By setting perfect=False
   the particle positions will be given random displacements
   with a magnitude of one-tenth the lattice spacing."""

def create(N, density, perfect=True):

  import random

  cubes = []
  for i in range(100):
    cubes.append(i**3)
  if(cubes.count(N) != 1):
    print '\nWARNING: num_particles is not a perfect cube. Initial'
    print '         configuration may be inhomogeneous.\n'

  L = (N / density)**(1.0/3.0)
  a = int(N**(1.0/3.0))
  if(a**3 < N):
    a = a + 1

  lattice_spacing = L / a

  random.seed(10)
  def rnd(magn_):
    return magn_ * (2.0 * random.random() - 1.0)

  # magnitude of random displacements
  magn = lattice_spacing / 10.0
  if(perfect): magn = 0.0

  ct = 0
  x = []
  y = []
  z = []
  for i in range(a):
    for j in range(a):
      for k in range(a):
	if(ct < N):
	  x.append(0.5 * lattice_spacing + i * lattice_spacing + rnd(magn))
	  y.append(0.5 * lattice_spacing + j * lattice_spacing + rnd(magn))
	  z.append(0.5 * lattice_spacing + k * lattice_spacing + rnd(magn))
	  ct += 1

  return x, y, z, L, L, L
