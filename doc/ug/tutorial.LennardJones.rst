Simple Lennard Jones System
===========================

Lets just copy and paste the beginning from the "System Setup" tutorial:

>>> import espressopp
>>> from espressopp import Real3D
>>> 
>>> system         = espressopp.System()
>>> box            = (10, 10, 10)
>>> rng            = espressopp.esutil.RNG()
>>> bc             = espressopp.bc.OrthorhombicBC(rng, box)
>>> system.bc      = bc
>>> system.rng     = rng
>>> maxcutoff      = pow(2.0, 1.0/6.0)
>>> skin           = 0.4
>>> system.skin    = skin
>>> nodeGrid       = (1,1,1)
>>> cellGrid       = (1,1,1)
>>> nodeGrid       = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,maxcutoff,skin)
>>> cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, maxcutoff, skin)
>>> ddstorage      = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
>>> system.storage = ddstorage
>>> 
>>> integrator     = espressopp.integrator.VelocityVerlet(system)
>>> dt             = 0.005
>>> integrator.dt  = dt

And lets add some random particles:

>>> num_particles = 20
>>> particle_list = []
>>> for k in range(num_particles):
>>>   pid  = k + 1
>>>   pos  = system.bc.getRandomPos()
>>>   v    = Real3D(0,0,0)
>>>   mass = system.rng()
>>>   type = 0
>>>   part = [pid, pos, type, v, mass]
>>>   particle_list.append(part)
>>> system.storage.addParticles(particle_list, 'id', 'pos', 'type', 'v', 'mass')
>>> system.storage.decompose()

All particles should interact via a Lennard Jones potential:

>>> LJPot = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=maxcutoff, shift='auto')

shift=True means that the potential will be shifted at the cutoff so that potLJ(cutoff)=0
Next we create a VerletList which will than be used in the interaction:
(the Verlet List object needs to know from which system to get its particles and which cutoff to use)

>>> verletlist = espressopp.VerletList(system, cutoff=maxcutoff)

Now create a non bonded interaction object and add the Lennard Jones potential to that:

>>> NonBondedInteraction = espressopp.interaction.VerletListLennardJones(verletlist)
>>> NonBondedInteraction.setPotential(type1=0, type2=0, potential=LJPot)

Tell the system about the newly created NonBondedInteraction object:

>>> system.addInteraction(NonBondedInteraction)

We should set the langevin thermostat in the integrator to cool down the random particle system:

>>> langevin             = espressopp.integrator.LangevinThermostat(system)
>>> langevin.gamma       = 1.0
>>> langevin.temperature = 1.0
>>> integrator.addExtension(langevin)

and finally let the system run and see how it relaxes or explodes:   

>>> espressopp.tools.analyse.info(system, integrator)
>>> for k in range(100):
>>>   integrator.run(10)
>>>   espressopp.tools.analyse.info(system, integrator)

Due to the random particle positions it may happen, that two or more particles are very close to
each other and the resulting repulsive force between them are so high that they 'shoot off' in
different directions with very high speed. Usually the numbers are then larger than the computer
can deal with. A typical error message you get could look like this:

.. note::
   ERROR: particle 5 has moved to outer space (one or more coordinates are nan)
   
In order to prevent this, systems that are setup in a random way and thus have strong overlaps between particels
have to be "warmed up" before they can be equilibrated. 

In ESPResSo++ there are several possible ways of warming up a system. As a first approach one could simply constrain
the forces in the integrator. For this purpose ESPResSo++ provides an integrator Extension
named CapForces. The two parameters of this Extension are the system and the maximum force that
a particle can get. The following python code shows how CapForces can be used. Add it to your
Lennard-Jones example just after adding the Langevin Extension:

>>> print "starting warmup with force capping ..."
>>> force_capping   = espressopp.integrator.CapForce(system, 1000000.0)
>>> integrator.addExtension(force_capping)
>>> # reduce the time step of the integrator to make the integration numerically more stable
>>> integrator.dt = 0.0001
>>> espressopp.tools.analyse.info(system, integrator)
>>> for k in range(10):
>>>   integrator.run(1000)
>>>   espressopp.tools.analyse.info(system, integrator)

After the warmup the time step of the integrator can be set to a larger value.
The CapForce extension can be disconnected after the warmup to get the original 
full Lennard-Jones potential back.

>>> integrator.dt   = 0.005
>>> integrator.step = 0
>>> force_capping.disconnect()
>>> print "warmup finished - force capping switched off."

Task 1: 
-------

write a python script that creates a random configuration of 1000 Lennard Jones
particles with a number density of 0.85 in a cubic simulation box.
Warm up and equilibrate this configuration.
Examine the output of the command

>>> espressopp.tools.analyse.info(system, integrator)

after each integration step. How fast is the energy of the system going down ?
How long do you have to warmup ? What are good parameters for
dt, force_capping and number of integration steps ?
