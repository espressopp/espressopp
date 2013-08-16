"""
***************************************
**warmup** - method to warm up a system
***************************************

This method does a warm up for a system with a density of 0.85.

The method needs the following parameters:

* system, integrator
    ESPResSo system which schoul be warmed up and the correspondig integrator e.g.:
    >>>system, integrator = espresso.standard_system.LennardJones(100,(10,10,10))
*number
    number of steps of the warm up
  =80
    for a system with a desnsity of 0.85, if it explodes try a higher number
"""
import espresso
from espresso import Real3D

def warmup(system, integrator, number=80):

  print "starting warmup"

  org_dt = integrator.dt
  pot = system.getInteraction(0).clonePotential(0,0)
  final_sigma = pot.sigma
  final_epsilon = pot.epsilon
  N = 50

  integrator.dt = 0.0001
  force_capping = espresso.integrator.CapForce(system, 0.0)
  integrator.addExtension(force_capping)

  for k in range(11,number):
    force_capping.setAbsCapForce(1000000.0/number*k)
    pot.sigma = final_sigma/number*k
    pot.epsilon = final_epsilon/number*k
    system.getInteraction(0).setPotential(0,0,pot)
    espresso.tools.analyse.info(system, integrator)
    integrator.run(N)

  integrator.dt = org_dt
  pot.sigma = final_sigma
  pot.epsilon = final_epsilon
  force_capping.disconnect()
  
  for k in range(11):
    integrator.run(70)
    espresso.tools.analyse.info(system, integrator)

  integrator.step = 0

  print "warmup finished"
