import espresso

system, integrator = espresso.standard_system.PolymerMelt(10, 20, (10, 10, 10))
integrator.dt = 0.001
sock = espresso.tools.vmd.connect(system)
t = 0
while t < 1000:
  integrator.run(1)
  espresso.tools.vmd.imd_positions(system, sock)
  t += 1
