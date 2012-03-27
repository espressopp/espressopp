Simple Lennard Jones System
===========================

Lets just copy and paste the beginning from the "System Setup" tutorial:

>>> import espresso
>>> from espresso import Real3D
>>> 
>>> system         = espresso.System()
>>> box            = (10, 10, 10)
>>> rng            = espresso.esutil.RNG()
>>> bc             = espresso.bc.OrthorhombicBC(rng, box)
>>> system.bc      = bc
>>> system.rng     = rng
>>> maxcutoff      = pow(2.0, 1.0/6.0)
>>> skin           = 0.4
>>> system.skin    = skin
>>> nodeGrid       = (1,1,1)
>>> cellGrid       = (1,1,1)
>>> nodeGrid       = espresso.tools.decomp.nodeGrid(espresso.MPI.COMM_WORLD.size)
>>> cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, maxcutoff, skin)
>>> ddstorage      = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
>>> system.storage = ddstorage
>>> 
>>> integrator     = espresso.integrator.VelocityVerlet(system)
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

>>> LJPot = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=maxcutoff, shift=True)

shift=True means that the potential will be shifted at the cutoff so that potLJ(cutoff)=0
Next we create a VerletList which will than be used in the interaction:
(the Verlet List object needs to know from which system to get its particles and which cutoff to use)

>>> verletlist = espresso.VerletList(system, cutoff=maxcutoff+system.skin)

Now create a non bonded interaction object and add the Lennard Jones potential to that:

>>> NonBondedInteraction = espresso.interaction.VerletListLennardJones(verletlist)
>>> NonBondedInteraction.setPotential(type1=0, type2=0, potential=LJPot)

Tell the system about the newly created NonBondedInteraction object:

>>> system.addInteraction(NonBondedInteraction)

We should set the langevin thermostat in the integrator to cool down the random particles:

>>> langevin             = espresso.integrator.Langevin(system)
>>> langevin.gamma       = 1.0
>>> langevin.temperature = 1.0
>>> integrator.langevin  = langevin

and finally let the system run and see how it relaxes:

>>> espresso.tools.analyse.info(system, integrator)
>>> for k in range(100):
>>>   integrator.run(10)
>>>   espresso.tools.analyse.info(system, integrator)

