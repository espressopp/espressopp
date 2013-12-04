import espresso
from espresso import Real3D

d = 0.85
Nchains = 10

Mmonomers = 10
N = Nchains * Mmonomers
L = pow(N/d, 1.0/3)

system, integrator = espresso.standard_system.PolymerMelt(Nchains, Mmonomers,(10,10,10), dt = 0.005, temperature=1.0)


print "starting warmup"
org_dt = integrator.dt
pot = system.getInteraction(0).getPotential(0,0)
print pot
print "Nint = ", system.getNumberOfInteractions()
final_sigma = pot.sigma
final_epsilon = pot.epsilon
print "sigma=",pot.sigma, "epsilon=",pot.epsilon
maxParticleID = int(espresso.analysis.MaxPID(system).compute())
N = 1
number = 50

for k in range(number):
  if k < 10:
    continue
 
  force_capping = espresso.integrator.CapForce(system, 1000000.0/number*k)
  integrator.addExtension(force_capping)

  pot.sigma = final_sigma/number*k
  pot.epsilon = final_epsilon/number*k

  integrator.dt = 0.0001
  espresso.tools.analyse.info(system, integrator)
  integrator.run(N)
  espresso.tools.analyse.info(system, integrator)

integrator.dt = org_dt
pot.sigma = final_sigma
pot.epsilon = final_epsilon
force_capping.disconnect()
  
for k in range(10):
  integrator.run(70)
  espresso.tools.analyse.info(system, integrator)
integrator.step = 0

print "warmup finished"

for k in range(10):
  integrator.run(100)
  espresso.tools.analyse.info(system, integrator)

