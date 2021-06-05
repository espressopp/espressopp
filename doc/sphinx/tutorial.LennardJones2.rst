Advanced Lennard Jones System
=============================

This tutorial needs the matplotlib.pyplot and numpy libraries and also VMD to be installed.

>>> import espressopp

After importing espressopp we import several other Python packages that we want to
use for graphical output of some system parameters (e.g. temperature and energy)

>>> import math 
>>> import time
>>> import matplotlib
>>> matplotlib.use('TkAgg')
>>> import matplotlib.pyplot as plt
>>> plt.ion()

We setup a standard Lennard-Jones system with 1000 particles and a density of 0.85
in a cubic siomulation box. ESPResSo++ provides a "shortcut" to setup such a system:

>>> N   = 1000
>>> rho = 0.85
>>> L   = pow(N/rho, 1.0/3)
>>> system, integrator = espressopp.standard_system.LennardJones(N,(L, L, L), dt=0.0001)

We also add a Langevin thermostat:

>>> langevin = espressopp.integrator.LangevinThermostat(system)
>>> langevin.gamma       = 1.0
>>> langevin.temperature = 1.0
>>> integrator.addExtension(langevin)

We do a very short warmup in the beginning to get rid of "extremely" high forces

>>> force_capping   = espressopp.integrator.CapForce(system, 1000000.0)
>>> integrator.addExtension(force_capping)
>>> espressopp.tools.analyse.info(system, integrator)
>>> for k in range(10):
>>>   integrator.run(100)
>>>   espressopp.tools.analyse.info(system, integrator)

Now let's initialize a graph. So that we can have a realtime-view on what is happening
in the simulation:

>>> plt.figure()

We want to observe temperature and energy of the system:

>>> T   = espressopp.analysis.Temperature(system)
>>> E   = espressopp.analysis.EnergyPot(system, per_atom=True)

x will be the x-axixs of the graph containg the time. yT and yE will be temperature
and energy as y-axes in 2 plots:

>>> x   = []
>>> yT  = []
>>> yE  = []
>>> yTmin = 0.0
>>> yEmin = 0.0
>>> x.append(integrator.dt * integrator.step)
>>> yT.append(T.compute())
>>> yE.append(E.compute())
>>> yTmax = max(yT)
>>> yEmax = max(yE)

Initialize the two graphs ('ro' means red circles, 'go' means green cirlces, see also pyplot documentation)

>>> plt.subplot(211)
>>> gT, = plt.plot(x, yT, 'ro')
>>> plt.subplot(212)
>>> gE, = plt.plot(x, yE, 'go')

We also want to observe the configuration with VMD. So we have to connect to vmd. This command
will automatically start vmd (vmd has to be found in your PATH environment for this to work)

>>> sock = espressopp.tools.vmd.connect(system)
>>> for k in range(200):
>>>   integrator.run(1000)
>>>   espressopp.tools.vmd.imd_positions(system, sock)

Update the x-, and y-axes:

>>>   x.append(integrator.dt * integrator.step)
>>>   yT.append(T.compute())
>>>   yE.append(E.compute())
>>>   yTmax = max(yT)
>>>   yEmax = max(yE)

Plot the temperature graph

>>>   plt.subplot(211)
>>>   plt.axis([x[0], x[-1], yTmin, yTmax*1.2 ])
>>>   gT.set_ydata(yT)
>>>   gT.set_xdata(x)
>>>   plt.draw()

Plot the energy graph

>>>   plt.subplot(212)
>>>   plt.axis([x[0], x[-1], yEmin, yEmax*1.2 ])
>>>   gE.set_ydata(yE)
>>>   gE.set_xdata(x)
>>>   plt.draw()

In the end save the equilibrated configurations as .eps and .pdf files

>>> plt.savefig('mypyplot.eps')
>>> plt.savefig('mypyplot.pdf')

