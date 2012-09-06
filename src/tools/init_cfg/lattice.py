# -*- coding: iso-8859-1 -*-

"""This Python module initializes particles on the sites
   of a simple cubic lattice. By setting perfect=False
   the particle positions will be given random displacements
   with a magnitude of one-tenth the lattice spacing."""

def create(N, rho, perfect=True, RNG=None):

  if RNG == None:
    import random
  
  cubes = []
  for i in range(100):
    cubes.append(i**3)
  if(cubes.count(N) != 1):
    print '\nWARNING: num_particles is not a perfect cube. Initial'
    print '         configuration may be inhomogeneous.\n'

  L = (N / rho)**(1.0/3.0)
  a = int(N**(1.0/3.0))
  if(a**3 < N):
    a = a + 1

  lattice_spacing = L / a

  def rnd(magn_):
    if RNG == None:
      rand = random.random()
    else :
      rand = RNG()
    return magn_ * (2.0 * rand - 1.0)

  # magnitude of random displacements
  magn = 0.0 if perfect else lattice_spacing / 10.0

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

# TODO implement checking for a wrong number of particles, lightly nonideal lattice etc.
def createDiamond(N, rho, perfect=True, RNG=None):
  from espresso import Real3D
  
  #L = (N / 8.0 / rho)**(1.0/3.0)
  L = (N / rho)**(1.0/3.0)
  
  num_per_edge = int( (N/8.0)**(1.0/3.0) )
  
  if(8.0*num_per_edge**3 < N):
    num_per_edge = num_per_edge + 1

  #print 'num_per_site= ', num_per_edge
  
  a = L / num_per_edge
  #print 'a= ', a
  #print 'a1= ', (1.0 / rho)**(1.0/3.0)

  pos = []
  # in general structure is shifted relative to (0,0,0)
  R0 = Real3D(0.125 * a, 0.125 * a, 0.125 * a)
  
  R1 = Real3D(0.25 * a, 0.25 * a, 0.25 * a)
  a11 = a * Real3D(1,0,0)
  a22 = a * Real3D(0,1,0)
  a33 = a * Real3D(0,0,1)
  
  a1 = 0.5 * a * Real3D(0,1,1)
  a2 = 0.5 * a * Real3D(1,0,1)
  a3 = 0.5 * a * Real3D(1,1,0)
  for i in range(num_per_edge):
    for j in range(num_per_edge):
      for k in range(num_per_edge):
        Rijk = R0 + i*a11 + j*a22 + k*a33
        pos.append(Rijk)
        pos.append(Rijk+a1)
        pos.append(Rijk+a2)
        pos.append(Rijk+a3)

        pos.append(Rijk+R1)
        pos.append(Rijk+a1+R1)
        pos.append(Rijk+a2+R1)
        pos.append(Rijk+a3+R1)

  '''
  L1 = L-0.01
  pos.append( Real3D(0.01, 0.01, 0.01) )
  pos.append( Real3D(L1, 0.01, 0.01) )
  pos.append( Real3D(0.01, L1, 0.01) )
  pos.append( Real3D(0.01, 0.01, L1) )
  pos.append( Real3D(0.01, L1, L1) )
  pos.append( Real3D(L1, L1, 0.01) )
  pos.append( Real3D(L1, 0.01, L1) )
  pos.append( Real3D(L1, L1, L1) )
  '''
  
  return pos, L, L, L
