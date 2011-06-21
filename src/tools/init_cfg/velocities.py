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
