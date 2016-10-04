Adaptive Resolution Simulations
===================================

.. |espp| replace:: ESPResSo++


Theory and Background
----------------------------------------------

ESPResSo++ provides functionality to run adaptive resolution simulations using the Adaptive Resolution Simulation Scheme (AdResS). In AdResS molecules in different regions in a simulation box are described by different non-bonded force fields, typically atomistic (AT) and coarse-grained (CG). These different subregions are interfaced and coupled via a hybrid region, where the interaction smoothly changes. Molecules can diffuse between the different regions and change their interaction on the fly.

There are two different AdResS approaches: The force-based scheme, in which forces are interpolated, as well as the energy-based scheme (Hamiltonian AdResS or H-AdResS) which interpolates on the level of potential energies. In force-based AdResS (see, for example, Praprotnik et al., J. Chem. Phys. 123, 224106 (2005) as well as Annu. Rev. Phys. Chem. 59, 545 (2008)), we have for the net force between the molecules :math:`\alpha` and :math:`\beta`

.. math::
  \mathbf{F}_{\alpha|\beta} = \lambda(\mathbf{R}_\alpha)\lambda(\mathbf{R}_\beta)\mathbf{F}_{\alpha|\beta}^{\text{AT}} + \left(1-\lambda(\mathbf{R}_\alpha)\lambda(\mathbf{R}_\beta)\right)\mathbf{F}_{\alpha|\beta}^{\text{CG}},
where :math:`\mathbf{F}_{\alpha|\beta}^{\text{AT}}` is an AT force field based on the individual atoms belonging to the molecules :math:`\alpha` and :math:`\beta` and :math:`\lambda` is a position dependent resolution function smoothly changing from 1 in the AT region to 0 in the CG region via the hybrid buffer region. It is evaluated based on the molecules' center of mass positions :math:`\mathbf{R_\alpha}`. Note that there can of course also be bonded interactions, but these are typically not interpolated, as they are computationally usually much cheaper to evaluate than the non-bonded forces. For the sake of clarity, we omit them here.

In H-AdResS (see Potestio et al., Phys. Rev. Lett. 110, 108301 (2013)), interpolation is performed directly on potential energies in the Hamiltonian as

.. math::
  H  = \sum_\alpha \sum_{i\in\alpha} \frac{\mathbf{p}_{\alpha i}^2}{2m_{\alpha i}}+\sum_{\alpha} \left\{\lambda({\mathbf{R}_\alpha}) V^{\text{AT}}_\alpha + \left(1 - \lambda({\mathbf{R}_\alpha})\right) V^{\text{CG}}_\alpha \right\},
where the first term corresponds to the kinetic energy and we again omitted intramolecular interactions. The forces obtained from this Hamiltonian are

.. math::
  \textbf{F}_{\alpha i} = & \sum_{\beta\neq\alpha}\sum_{j\in\beta}\left\{ \frac{\lambda_\alpha+\lambda_\beta}{2}\textbf{F}_{\alpha i|\beta j}^{\text{AT}} + \left(1-\frac{\lambda_\alpha+\lambda_\beta}{2}\right)\textbf{F}_{\alpha i|\beta}^{\text{CG}}  \right\} \\
  & - [V^{\text{AT}}_\alpha - V^{\text{CG}}_\alpha] \,\nabla_{\alpha i}\lambda_\alpha.
The last term, the so-called drift force, comes from applying the position gradient on the position-dependent resolution function :math:`\lambda`. It acts only in the hybrid region and unphysically pushes molecules from one region to the other. Therefore, it needs to be corrected. On the other hand, force-based AdResS, contrary to H-AdResS, does not allow a Hamiltonian formulation at all.

Usually, the force fields used in the different regions of the adaptive simulation setup have significantly different pressures given the same temperature and particle density. This pressure gradient leads, in addition to the drift force in H-AdResS, to particles being pushed across the hybrid region. Eventually, the system would evolve to an equilibrium state with a inhomogeneous density profile across the simulation box. Therefore, correction forces needs to be applied in the hybrid region to counter these effects. In H-AdResS one can use a so-called free energy correction (FEC), which on average cancels the drift force in the hybrid region (see Potestio et al., Phys. Rev. Lett. 110, 108301 (2013)). The FEC corresponds to the free energy difference between the subsystems and can, for instance, be derived from Kirkwood thermodynamic integration. An alternative approach which is typically used to cancel the pressure gradient in force-based AdResS is the so-called thermodynamic force (see Fritsch et al., Phys. Rev. Lett. 108, 170602 (2012)). It is derived by constructing the correction directly from the distorted density profile which is obtained without any correction and then refined iteratively.


|espp| code
----------------------

Several measures had to be taken to implement adaptive resolution simulations in ESPResSo++. On top of the normal particles, which serve as the CG particles in AdResS, another layer of extra AT particles is introduced such that one has access to both atomistic and CG particles throughout the whole system. A mapping between the two defines which atoms belong to which CG bead. The resolution function :math:`\lambda` is implemented as a particle property of the CG particles that is updated after each integration step based on the new positions. This happens in an extension to the Velocity Verlet integrator. The actual adaptive resolution scheme is then implemented via new interaction templates that define how forces and energies are computed in force-based and energy-based AdResS. These templates use for particle pairs in the atomistic region the actual atoms, this is the AT particles, for the force and energy computation while in the CG region they use the CG particles. In the hybrid region, both are used, as defined in the equations above. The drift term of H-AdResS is implemented similarly. Furthermore, the AdResS integrator extension makes sure that the atomistic particles in the CG region travel along with the CG particles and that similarly the CG particles in the AT region are properly updated according to the new atomistic positions after each integration step. The FEC as well as a module to apply the Thermodynamic Force are implemented as integrator extensions.

In the following, we explain the new features step by step (more details about parameters etc. can be found in the documention of the different classes).

Adress Domain Decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When setting up the storage we have to use an appropriate domain decomposition that accomodates storage and proper interprocessor communication of both AT and CG particles.

.. code-block:: python

  # (H-)AdResS domain decomposition
  system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

Atomistic and Coarse-Grained particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When adding particles to the storage, we have to define them as atomistic or coarse-grained. This has been implemented as the particle property "adrat". If it is 0, the particle is coarse-grained. If it is 1, it is an atomistic particle.

.. code-block:: python

  # add particles to system
  system.storage.addParticles(allParticles, "id", "pos", "v", "f", "type", "mass", "adrat")

When adding the particles as above, it is important that a set of atomistic particles belonging to one CG particle appears in the list of particles ``allParticles`` always after the corresponding CG particle.

Next, the FixedTupleListAdress defines which atomistic particles belong to which coarse-grained particles.

.. code-block:: python

  # create FixedTupleList object and add the tuples
  ftpl = espressopp.FixedTupleListAdress(system.storage)
  ftpl.addTuples(tuples)
  system.storage.setFixedTuplesAdress(ftpl)

In this example, ``tuples`` is a list of tuples, where each tuple itself is another short list in which the first element is the CG particle and the other elements are the AT particles belonging to it. Note that in ESPResSo++ the CG particle is positioned always in the center of mass of its atoms.

Having set up the FixedTupleList, we can also set up an AdResS fixed pair list that defines bonds between AT particles within individual molecules. This is done in the following way:

.. code-block:: python

  # add bonds between AT particles
  fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
  fpl.addBonds(bonds)

where ``bonds`` is a list of bonds between AT particles within CG molecules. Similarly, triple lists for angles, quadruple lists for dihedrals etc. are set up. Compared to conventional bonds, angles, etc. between different normal CG particles one just adds the suffix ``Adress`` to the appropriate list object and provides it also with the FixedTupleList (``ftpl`` in the example). Note that you can define several different such fixed pair lists and you can, for example, also in AdResS simulations still use the normal ``FixedPairList`` to define bonds between regular CG particles.

AdResS Verlet List
~~~~~~~~~~~~~~~~~~~~~~~~~
Next, we construct the AdResS Verlet list object for non-bonded interacting particle pairs:

.. code-block:: python

  # AdResS Verlet list
  vl = espressopp.VerletListAdress(system, cutoff=0.8, adrcut=1.4,
                                  dEx=1.5, dHy=1.0,
                                  adrCenter=[Lx/2, Ly/2, Lz/2], sphereAdr=False)

We have to provide the cutoffs of the list as well as the sizes of the atomistic and hybrid regions. The parameter ``cutoff`` corresponds to the cutoff used for CG particle pairs with both particles being in the CG region, while ``adrcut`` is the cutoff for all other particle pairs (at least one particle of the pair is in the AT or hybrid region). We want to stress that this  pair list is build based on the CG particles' positions. Hence, for the AT and hybrid region one needs in some situations to provide a Verlet list cutoff (adrcut) slightly larger than the actual maximum interaction range of the potential, in order to not lose interactions between some atom pairs. Let us clarify this with an example: Thinking of a pair of water molecules, both coarse-grained into single beads, these CG beads could be farther apart than the interaction cutoff. Two hydrogen atoms pointing towards each other, however, could in fact still be in interaction range. Therefore, an appropriate buffer needs to be provided.

The ``sphereAdr`` flag decides how to geometrically set up the change in resolution. If it's true, the AT region is a spherical region positioned at ``adrCenter`` with radius ``dEx``. If ``sphereAdr`` is false, the resolution changes along the x-axis of the system and ``dEx`` corresponds to half the width of the AT region. ``dHy`` always is the full width of the hybrid region. Instead of providing a 3D position for ``adrCenter`` as above, one can also provide a particle ID of a CG particle. In this case, the atomistic region will follow the movement of the particle. This should be only done, however, for force-based AdResS, since it would break the Hamiltonian character of H-AdResS, and also only when using a spherical adaptive geometry. Then, however, it is even possible to provide a list of particle IDs, in which case the AT region corresponds to the overlap of the spherical regions defined by the individual particles provided in the list. It will deform accordingly while these particle move.

Interactions
~~~~~~~~~~~~~~~~~~~~~~~~~
When adding interactions to the system we have to use the corresponding interaction templates. Here is how to set up a non-bonded interaction in a H-AdResS system:

.. code-block:: python

  # H-AdResS non-bonded interaction: WCA potential between AT particles
  # and tabulated potential between CG particles
  interNB = espressopp.interaction.VerletListHadressLennardJones(vl, ftpl)
  potWCA  = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift='auto',
                                    cutoff=rca)
  potCG = espressopp.interaction.Tabulated(itype=3, filename=tabCG, cutoff=rc) # CG
  interNB.setPotentialAT(type1=1, type2=1, potential=potWCA) # AT
  interNB.setPotentialCG(type1=0, type2=0, potential=potCG) # CG
  system.addInteraction(interNB)

First, we define the appropriate interaction type, in H-AdResS this is ``VerletListHadressLennardJones``. Next we define the actual potentials. Then we associate them with the H-AdResS interaction and add the interaction to the system. For force-based AdResS the only change required would be to use the ``VerletListAdressLennardJones`` interaction.

Note that the here used interaction, ``VerletListHadressLennardJones``, couples only Lennard-Jones-type potentials with tabulated ones. However, there exist more such interaction templates for other potentials and potential combinations.

AdResS Integrator Extension
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, we have to set up the AdResS integrator extension:

.. code-block:: python

  # AdResS integrator extension
  adress = espressopp.integrator.Adress(system, verletlist, ftpl, regionupdates = 1)
  integrator.addExtension(adress)

It takes as arguments the Verlet list and the fixed tuple list. Additionally, for the case of a moving and/or deforming AdResS region based on one or more particles, the parameter ``regionupdates`` specifies how regularly we want to update the shape of the AdResS region in number of steps. This is to avoid as much as possible of the additional communication required to inform different processors of the change of the AdResS region. The parameter defaults to 1 and is not used at all for static AdResS regions.

Having set up the AdResS extension, we can distribute all particles in the box and place the CG molecules in the centers of mass of the atoms which they belong to. This can be done conveniently via

.. code-block:: python

  # distribute atoms and CG molecules according to AdResS domain decomposition,
  # place CG molecules in the center of mass
  espressopp.tools.AdressDecomp(system, integrator)

Free Energy Compensation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When using H-AdResS, we probably want to also employ a FEC. This can be done as follows:

.. code-block:: python

  # set up FEC
  fec = espressopp.integrator.FreeEnergyCompensation(system, center=[Lx/2, Ly/2, Lz/2])
  fec.addForce(itype=3, filename="table_fec.tab", type=1)
  integrator.addExtension(fec)

The FEC takes as arguments the system object as well as the center of the AT region. Then we add the actual force, which needs to be provided in a table (first column: resolution :math:`\lambda`, second: energy, third: force). ``itype`` defines which type of interpolation should be used for values between the ones provided in the table. 1 corresponds to linear interpolation, 2 to akima splines, 3 to cubic splines. We suggest to use cubic splines. The FEC is applied on CG particles and distributed among the atoms belonging to the CG particle. ``type`` specifies the CG particle type for which this correction should be applied. One can, for example, use different FECs for different molecules types.

Thermodynamic Force
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When using force-based AdResS, or, alternatively, in addition to the FEC in H-AdResS, we can use the thermodynamic force. It can be set up in the following way, very similar to the FEC before:

.. code-block:: python

  # set up Thermodynamic Force
  thdforce = espressopp.integrator.TDforce(system, verletlist)
  thdforce.addForce(itype=3,filename="table_tf.tab",type=1)
  integrator.addExtension(thdforce)

It works largely as for the FEC with the following differences: The table should not provide resolution values in the first column but actual distance values, this is, the distance from the (closest) AT region center. This allows to extend the application of the thermodynamic force slightly beyond the borders of the hybrid region where the resolution is constant. Furthermore, the Thermodynamic Force needs the verletlist as argument.

It is also possible to define a thermodynamic force, which is suited for an adaptive resolution setup with an AT region that is constructed via the overlap of several spherical regions. In this case, the extension needs more information:

.. code-block:: python

  # set up Thermodynamic Force
  thdforce = espressopp.integrator.TDforce(system, verletlist, startdist = 0.9,
                                      enddist = 2.1, edgeweightmultiplier = 20)
  thdforce.addForce(itype=3,filename="table_tf.tab",type=1)
  integrator.addExtension(thdforce)

It gets three more parameters, ``startdist``, ``enddist`` and ``edgeweightmultiplier``. ``startdist`` explicitely says at which distance from the center of the closest AT region defining particle the thermodynamic force starts to act and ``enddist`` says where it ends. Hence, these value should correspond to what is actually written in the table. ``edgeweightmultiplier`` is a parameter that speficies how precisely the thermodynamic force should be applied in the overlap regions of different spheres. For most applications, however, 20 should provide reasonable results (for details, see Kreis et al., J. Chem. Theory Comput. 12, 4067 (2016)). The 3 additional parameters are of course also present with some default values in the basic case, but they are ignored unless we have an AT region that is constructed via the overlap of several spherical regions.

Examples
----------------------

We have provided several example scripts and setups that are available in the ESPResSo++ source code at ``examples/adress``. Most of them are based on published papers.

The reader is strongly encouraged to play around with them and test what happens when the setups are modified. Possible questions to ask are provided at the end of the following subsections, which explain the individual examples in more detail.


Force-AdResS: Tetrahedral Liquid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Subfolder: ``fadress_tetraliquid``. This example consists of the system that was used in the initial work introducing the force-based adaptive resolution method (see Praprotnik et al., J. Chem. Phys. 123, 224106 (2005) and Phys. Rev. E 73, 066701 (2006)). A liquid composed of artificial tetrahedral molecules, i.e. each molecule consists of 4 bonded atoms arranged in a tetrahedral geometry, is coupled to a CG model which describes the molecules as individual beads.

Questions: The geometry is set in such a way that the resolution changes along the x-axis of the box. Try changing the setup such that the AT region is of spherical shape. You can also try removing the thermostat. Does the system conserve energy? Also vary the size of the atomistic region and see what happens. Can you also make the system all-atomistic or all-CG? You can also try to compare computational times.

Force-AdResS: A Protein in Water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Subfolder: ``fadress_protein``. This system is an aqueous solution of the regulatory protein ubiquitin. The atomistic protein and the atomistic water around it is coupled to a coarse-grained water model, which maps water molecules farther away from the protein to single beads. The CG water interaction was parametrized with iterative Boltzmann inversion (IBI).  This system is similar to the setup which was used by Fogarty et al. (J. Chem. Phys. 142, 195101 (2015)) to study the structure and dynamics of a protein hydration shell.

Questions: The setup is significantly more complicated than the previous system. Try to understand the script. You can also have a look into the the actual source code and try to understand, for example, how the gromacs parser works. The example is set up as a fully atomistic simulation by setting the size of the atomistic region to a value larger than the simulation box. Try to change the script such that it is an actual adaptive setup. Do not forget the thermodynamic force! Furthermore, how is the high-resolution region positioned now?

Force-AdResS: Self-Adjusting Adaptive Resolution Simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Subfolder: ``fadress_selfadjusting``. This setup demonstrates how force-based adaptive resolution simulations with self-adjusting high-resolution regions can be set up (Kreis et al., J. Chem. Theory Comput. 12, 4067 (2016)). The system is a polyalanine-9 molecule in aqueous solution. A spherical AT region is associated with each atom of the peptide such that the overall AT region formed by the overlap of all these spheres elegantly envelops the peptide. The peptide starts in an extended configuration and as it folds, the AT region surrounding it adjusts itself accordingly. At the outside, we use again a coarse-grained IBI single-bead model for the water molecules.

Questions: Can you change the system such that fewer atoms are associated with AT region, for example, only the heavy atoms? Can you change the update frequency of the shape of the AT region?

H-AdResS: Tetrahedral Liquid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Subfolder: ``hadress_tetraliquid``. This is the system used by Potestio et al. in the paper that proposed the H-AdResS method (Phys. Rev. Lett. 110, 108301 (2013)). It is again a simple system composed of tetrahedral molecules that change their resolution and become individual beads in the CG region. The interpolation occurs along the x-axis. This example has three subfolders.

The first folder ``hadress_tetraliquid_plain`` runs a simple H-AdResS simulation without any free energy correction. Hence, the drift force strongly pushes molecules from one region to the other. The script contains analysis routines which measure both a density and a pressure profile along the direction of resolution change while the simulation is running. Gathering enough statistics takes a while, but we have also provided reference profiles which are obtained after a sufficiently long simulation. Have a look at them and try to interpret them.

The second folder ``hadress_tetraliquid_FEC`` contains the same setup but with a free energy correction. For this, two tables are provided, ``table_FEC_Helmholtz.dat`` and ``table_FEC_Gibbs.dat``. They were derived via Kirkwood thermodynamic integration. The first one is based on the Helmholtz free energy difference per particle between the two subsystem, and the second one corresponds to the Gibbs free energy difference per particle. Two density and pressure profiles obtained while applying these correction are also shown. Try to interpret them.

The third folder ``hadress_tetraliquid_KTI`` contains a simple implementation of Kirkwood thermodynamic integration (KTI) which could in principle, when run for long enough, be used to derive the FEC. This is not an adaptive resolution simulation. Instead, we tell the AdResS integrator extension that we want to run KTI. Then, the extension does not modify the resolution values associated with the different molecules and we can change them by hand during the simulation. In this way, we can set up a simulation in which we change the resolution of all molecules in the system every few steps and slowly proceed from a complete CG system to an all-atom one. Have a look and try to understand what is going on.

There are many more interesting things you can try out: Are the H-AdResS simulations energy conserving? Add the commented Langevin thermostat and compare. Also vary the timestep. Additionally, you can change the size of the hybrid region. What happens if it becomes smaller or larger? Furthermore, what happens if you change the system from H-AdResS to force-based AdResS?

H-AdResS: Water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Subfolder: ``hadress_water``. This is a slightly more advanced H-AdResS system in which an atomistic model is coupled to a coarse-grained one, mapping the three water atoms onto single beads.

Questions: Feel free to play around with the system. You could also try to figure out, how the gromacs parsers sets up the interactions and chooses the right H-AdResS interactions.
