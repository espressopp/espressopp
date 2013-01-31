Polymer Melt
============




We first import espresso and then define all the parameters of the simulation:

>>> import espresso
>>> num_chains          = 10
>>> monomers_per_chain  = 10
>>> L                   = 10
>>> box                 = (L, L, L)
>>> bondlen             = 0.97
>>> rc                  = pow(2, 1.0/6.0)
>>> skin                = 0.3
>>> dt                  = 0.005
>>> epsilon             = 1.0
>>> sigma               = 1.0

Like in the simple Lennard Jones tutorial we setup the system and the integrator.
First the system with the excluded volume interaction (WCA, Lennard Jones type)

>>> system         = espresso.System()
>>> system.rng     = espresso.esutil.RNG()
>>> system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
>>> system.skin    = skin
>>> nodeGrid       = espresso.tools.decomp.nodeGrid(espresso.MPI.COMM_WORLD.size)
>>> cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
>>> system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
>>> interaction    = espresso.interaction.VerletListLennardJones(espresso.VerletList(system, cutoff=rc))
>>> potLJ          = espresso.interaction.LennardJones(epsilon, sigma, rc)
>>> interaction.setPotential(type1=0, type2=0, potential=potLJ)
>>> system.addInteraction(interaction)

Then the integrator with the Langevin extension

>>> integrator     = espresso.integrator.VelocityVerlet(system)  
>>> integrator.dt  = dt
>>> thermostat     = espresso.integrator.LangevinThermostat(system)
>>> thermostat.gamma  = 1.0
>>> thermostat.temperature = temperature
>>> integrator.addExtension(thermostat)

Know we add the particles. Keep in mind that we want to create a polymer melt. This means
that particles are "bonded" in chains. We setup each polymer chain as a random walk.

>>> props    = ['id', 'type', 'mass', 'pos', 'v']
>>> vel_zero = espresso.Real3D(0.0, 0.0, 0.0)

In providing bonding information for the particles we "setup" the bonded chains.
For this we use the FixedPairList object that needs to know where and in which storage
the particles can be found:

>>> bondlist = espresso.FixedPairList(system.storage)
>>> pid      = 1
>>> type     = 0
>>> mass     = 1.0  
>>> chain    = []

|espp| provides a function that will return position and bond information of a random walk.
You have to provide a start ID (particle id) and a starting position which we will get from the
random position generator of the boundary condition object:

>>> for i in range(num_chains):
>>>   startpos = system.bc.getRandomPos()
>>>   positions, bonds = espresso.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen)
>>>   for k in range(monomers_per_chain):  
>>>     part = [pid + k, type, mass, positions[k], vel_zero]
>>>     chain.append(part)
>>>   pid += monomers_per_chain
>>>   type += 1
>>>   system.storage.addParticles(chain, *props)
>>>   system.storage.decompose()
>>>   chain = []
>>>   bondlist.addBonds(bonds)

.. note::
   try out the command

   >>> espresso.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen)
  
   to see what it returns

Don't forget to distribute the particles and the bondlist to the CPUs in the end:

>>> system.storage.decompose()

Finally add the information about the bonding potential. In this example we are using
a FENE-potential between the bonded particles.

>>> potFENE   = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
>>> interFENE = espresso.interaction.FixedPairListFENE(system, bondlist, potFENE)
>>> system.addInteraction(interFENE)

Start the integrator and observe how the system explodes. Like in the random Lennard Jones
system, we have the same problem here: particles can strongly overlap and thus will
get very high forces accelerating them to infinite (for the computer) speed.

>>> espresso.tools.analyse.info(system, integrator)
>>> for k in range(nsteps):
>>>   integrator.run(isteps)
>>>   espresso.tools.analyse.info(system, integrator)
>>>   espresso.tools.analyse.info(system, integrator)

Task 2:
-------

Try to warmup and equilibrate a dense polymer melt (density=0.85) by using the warmup methods
that you have learned in the Lennard Jones tutorial.

Hint:
----
During warmup you can slowly switch on the excluded volume interaction by starting with a small
epsilon and increasing it during integration:
You can do this by continuously overwriting the interaction potential after some time interval.

>>> potLJ          = espresso.interaction.LennardJones(new_epsilon, sigma, rc)
>>> interaction.setPotential(type1=0, type2=0, potential=potLJ)


