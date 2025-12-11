---
title: 'Espresso++: A Fast and Extensible Molecular Simulation Package for Coarse-Grained Models'
tags:
  - molecular dynamics
  - coarse-grained simulations
  - soft matter
  - parallel computing
  - C++/Python
authors:
  - name: Sebastian Eibl
    orcid: 0000-0002-1069-2720
    affiliation: 2
  - name: Markus Rampp
    orcid: 0000-0001-8177-8698
    affiliation: 2
  - name: James Vance
    orcid: 0000-0001-7112-0382
    affiliation: 3
  - name: Nikita Tretyakov
    affiliation: 3
  - name: Tristan Bereau
    orcid: 0000-0001-9945-1271
    affiliation: 5
  - name: Horacio Vargas Guzman
    orcid: 0000-0003-2564-3005
    affiliation: 1
  - name: Bin Song
    affiliation: 1
  - name: Zhen-Hao Xu
    affiliation: 1
  - name: Pavel Kus
    affiliation: 2
  - name: Jakub Krajniak
    orcid: 0000-0001-9372-6975
    affiliation: 6
  - name: Torsten Stuehn
    orcid: 0009-0006-2144-2002
    affiliation: 1
  - name: Christoph Junghans
    orcid: 0000-0003-0925-1458
    affiliation: 3
affiliations:
  - name: Max Planck Institute for Polymer Research, Mainz, Germany
    index: 1
  - name: Max Planck Computing and Data Facility, Garchingen, Germany
    index: 2
  - name: Los Alamos National Laboratory, Los Alamos, USA
    index: 3
  - name: Johannes Gutenberg University of Mainz, Mainz, Germany
    index: 4
  - name: Heidelberg University, Heidelberg, Germany
    index: 5
  - name: Independent researcher, Poznań, Poland 
    index: 6


date: 2025-10-02
bibliography: paper.bib
---

# Summary

**Espresso++** is an open-source software package for **molecular dynamics (MD) simulations** with a particular emphasis on **coarse-grained models** of soft matter systems. Written in C++ with a flexible Python interface, it is designed for **high-performance computing (HPC)** environments and supports **massively parallel simulations** through MPI. The package enables simulations of polymers, membranes, colloids, complex fluids, and active matter with a wide range of interaction models and advanced algorithms.

Espresso++ builds upon the experience of its predecessor [ESPResSo](http://espressomd.org), but provides a cleaner, more modular codebase and enhanced extensibility. The software is actively developed by an international community of researchers in physics, chemistry, biology, and materials science.

# Statement of need

Molecular dynamics simulations are essential tools for exploring the behavior of soft matter systems at mesoscopic scales. Traditional all-atom MD codes (e.g., GROMACS, LAMMPS) are often not efficient or flexible enough for coarse-grained models that require custom interactions or specialized algorithms. Espresso++ addresses this gap by providing:

- A modular and extensible design, enabling researchers to easily implement new interaction potentials and integrators.
- Efficient parallelization for large-scale simulations of complex systems.
- A rich library of coarse-grained interaction models and algorithms tailored to soft matter.
- A Python-based scripting interface for ease of use, reproducibility, and coupling with external analysis tools.

Espresso++ is widely used in academic research for simulating phenomena such as polymer rheology, membrane dynamics, colloidal suspensions, and active particles.

# Functionality

Key features of Espresso++ include:

- **Inter-particle interactions**: Lennard-Jones, Coulomb, soft repulsive, bonded interactions, tabulated potentials, and more.
- **Algorithms**: Molecular dynamics, Langevin dynamics, dissipative particle dynamics (DPD), Brownian dynamics, Monte Carlo sampling.
- **Electrostatics**: Particle–particle particle–mesh (P3M), Ewald summation, and other long-range methods.
- **Parallelization**: Domain decomposition using MPI, optimized for massively parallel architectures.
- **Python interface**: Full simulation control and analysis scripting in Python.
- **Extensibility**: Modular design allows easy addition of new force fields, integrators, or analysis tools.

# New Feature since last release

Copy list from Markus' google doc

# Impact
Espresso++ has been applied in numerous scientific studies, including investigations of:

- Polymer rheology and entanglement effects
- Lipid membranes and vesicle dynamics
- Colloidal self-assembly
- Active matter and microswimmers

It is actively maintained and extended by a community of researchers, with contributions from multiple institutions. Its flexible architecture makes it a valuable tool for developing novel coarse-grained models and algorithms in soft matter research.

# Research Projects using ESPResSo++ in the last 5 years:

- Grommes, Dirk and Bruch, Olaf and Imhof, Wolfgang and Reith, Dirk, 
"Coarse-Grained Molecular Dynamics Study of the Melting Dynamics in Long Alkanes", 
POLYMERS (2025), 
doi:10.3390/polym17182500

- Gholami, Abbas and Kloth, Sebastian and Xu, Zhen-Hao and Kremer, Kurt and Vogel, Michael and Stuehn, Torsten and Rudzinski, Joseph F., 
"Structure and dynamics of ionic liquids under shear flow", 
JOURNAL OF CHEMICAL PHYSICS (2025), 
doi:10.1063/5.0279946

- Grommes, Dirk and Bruch, Olaf and Reith, Dirk, 
"Mimicking Polymer Processing Conditions on the Meso-Scale: Relaxation and Crystallization in Polyethylene Systems after Uni- and Biaxial Stretching", 
MOLECULES (2024), 
doi:10.3390/molecules29143391

- Hsu, Hsiao-Ping and Kremer, Kurt, 
"Entanglement-Stabilized Nanoporous Polymer Films Made by Mechanical Deformation", 
MACROMOLECULES (2024), 
doi:10.1021/acs.macromol.4c00187

- Hsu, Hsiao-Ping and Kremer, Kurt, 
"Glass transition temperature of (ultra-)thin polymer films", 
JOURNAL OF CHEMICAL PHYSICS (2023), 
doi:10.1063/5.0165902

- Vance, James and Xu, Zhen-Hao and Tretyakov, Nikita and Stuehn, Torsten and Rampp, Markus and Eibl, Sebastian and Junghans, Christoph and Brinkmann, Andre, 
"Code modernization strategies for short-range non-bonded molecular dynamics simulations", 
COMPUTER PHYSICS COMMUNICATIONS (2023),
doi:10.1016/j.cpc.2023.108760

- Ohkuma, Takahiro and Hagita, Katsumi and Murashima, Takahiro and Deguchi, Tetsuo, 
"Miscibility and exchange chemical potential of ring polymers in symmetric ring-ring blends", 
SOFT MATTER (2023), 
doi:10.1039/d3sm00108c

- Grommes, Dirk and Schenk, Martin R. and Bruch, Olaf and Reith, Dirk, 
"Initial Crystallization Effects in Coarse-Grained Polyethylene Systems after Uni- and Biaxial Stretching in Blow-Molding Cooling Scenarios", 
POLYMERS (2022), 
doi:10.3390/polym14235144

- Grommes, Dirk and Schenk, Martin R. and Bruch, Olaf and Reith, Dirk, 
"Investigation of Crystallization and Relaxation Effects in Coarse-Grained Polyethylene Systems after Uniaxial Stretching", 
POLYMERS (2021), 
doi:10.3390/polym13244466

- Brunk, Aaron and Duenweg, Burkhard and Egger, Herbert and Habrich, Oliver and Lukacova-Medvid'ova, Maria and Spiller, Dominic, 
"Analysis of a viscoelastic phase separation model", 
JOURNAL OF PHYSICS-CONDENSED MATTER (2021), 
doi:10.1088/1361-648X/abeb13

- Zhang, Zidan and Krajniak, Jakub and Ganesan, Venkat, 
"A Multiscale Simulation Study of Influence of Morphology on Ion Transport in Block Copolymeric Ionic Liquids", 
MACROMOLECULES (2021), 
doi:10.1021/acs.macromol.1c00025

- Tubiana, Luca and Kobayashi, Hideki and Potestio, Raffaello and Duenweg, Burkhard and Kremer, Kurt and Virnau, Peter and Daoulas, Kostas, 
"Comparing equilibration schemes of high-molecular-weight polymer melts with topological indicators", 
JOURNAL OF PHYSICS-CONDENSED MATTER (2021), 
doi:10.1088/1361-648X/abf20c

- Thaler, S. and Praprotnik, M. and Zavadlav, J., 
"Back-mapping augmented adaptive resolution simulation", 
JOURNAL OF CHEMICAL PHYSICS (2020), 
doi:10.1063/5.0025728

- Hsu, Hsiao-Ping and Kremer, Kurt, 
"Efficient equilibration of confined and free-standing films of highly entangled polymer melts", 
JOURNAL OF CHEMICAL PHYSICS (2020), 
doi:10.1063/5.0022781

- Fiorentini, Raffaele and Kremer, Kurt and Potestio, Raffaello, 
"Ligand-protein interactions in lysozyme investigated through a dual-resolution model", 
PROTEINS-STRUCTURE FUNCTION AND BIOINFORMATICS (2020), 
doi:10.1002/prot.25954

- Papez, Petra and Merzel, Franci and Praprotnik, Matej, 
"Rotational Dynamics of a Protein under Shear Flow Studied by the Eckart Frame Formalism", 
JOURNAL OF PHYSICAL CHEMISTRY B (2023), 
doi:10.1021/acs.jpcb.3c02324

- Smith, Spencer and Michalski, Peter and Carette, Jacques and Keshavarz-Motamed, Zahra, 
"State of the Practice for Lattice Boltzmann Method Software", 
ARCHIVES OF COMPUTATIONAL METHODS IN ENGINEERING (2024), 
doi:10.1007/s11831-023-09981-2

- Bause, Marius and Bereau, Tristan, 
"Reweighting non-equilibrium steady-state dynamics along collective variables", 
JOURNAL OF CHEMICAL PHYSICS (2021), 
doi:10.1063/5.0042972

- Rudzinski, Joseph F. and Bereau, Tristan, 
"Coarse-grained conformational surface hopping: Methodology and transferability", 
JOURNAL OF CHEMICAL PHYSICS (2020), 
doi:10.1063/5.0031249

- Singh, Manjesh K. and Hu, Minghan and Cang, Yu and Hsu, Hsiao-Ping and Therien-Aubin, Heloise and Koynov, Kaloian and Fytas, George and Landfester, Katharina and Kremer, Kurt, 
"Glass Transition of Disentangled and Entangled Polymer Melts: Single-Chain-Nanoparticles Approach", 
MACROMOLECULES (2020), 
doi:10.1021/acs.macromol.0c00550

- Zhao, Yani and Cortes-Huerto, Robinson and Kremer, Kurt and Rudzinski, Joseph F., 
"Investigating the Conformational Ensembles of Intrinsically Disordered Proteins with a Simple Physics-Based Model", 
JOURNAL OF PHYSICAL CHEMISTRY B (2020), 
doi:10.1021/acs.jpcb.0c01949

- Zhao, Yani and Singh, Manjesh K. and Kremer, Kurt and Cortes-Huerto, Robinson and Mukherji, Debashish, 
"Why Do Elastin-Like Polypeptides Possibly Have Different Solvation Behaviors in Water-Ethanol and Water-Urea Mixtures?", 
MACROMOLECULES (2020), 
doi:10.1021/acs.macromol.9b02123

- Lee, Eunsang and Paul, Wolfgang, 
"Additional Entanglement Effect Imposed by Small Sized Ring Aggregates in Supramolecular Polymer Melts: Molecular Dynamics Simulation Study", 
MACROMOLECULES (2020), 
doi:10.1021/acs.macromol.9b02209

- Grommes, Dirk and Reith, Dirk, 
"Determination of relevant mechanical properties for the production process of polyethylene by using mesoscale molecular simulation techniques", 
SOFT MATERIALS (2020), 
doi:10.1080/1539445X.2020.1722692


# Example usage

A minimal Python script to run a simple (repulsive only) Lennard-Jones type particle simulation in ESPResSo++ looks like:

```python
import espressopp

# simulation system parameters
num_particles = 10000      # total number of particles in the system
box           = (20,20,20) # size of the simulationbox (all length are in sigma)
rc            = 1.12246    # cut off for the short range non bonded potential
skin          = 0.3        # skin used for verlet neighbor list
dt            = 0.005      # time step for 1 md step
epsilon       = 1.0        # energy unit
sigma         = 1.0        # length unit
temperature   = 1.0        # temperature of the simulation, e.g. used in Langevin Thermostat

# system setup
system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin

# define underlying storage system for parallelisation
nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size,box,rc,skin)
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# interaction setup, here short range non-bonded Lennard Jones potential
interaction    = espressopp.interaction.VerletListLennardJones(espressopp.VerletList(system, cutoff=rc))
interaction.setPotential(type1=0, type2=0, potential=espressopp.interaction.LennardJones(epsilon, sigma, rc, shift='auto'))
system.addInteraction(interaction)

# integrator setup
integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = dt

# thermostat setup
thermostat             = espressopp.integrator.LangevinThermostat(system)
thermostat.gamma       = 1.0
thermostat.temperature = temperature
integrator.addExtension(thermostat)

# create random particle setup in the simulation box
props = ['id', 'type', 'mass', 'pos', 'v']
new_particles = []
pid = 1
while pid <= num_particles:
    type = 0
    mass = 1.0
    pos  = system.bc.getRandomPos()
    vel  = espressopp.Real3D(0.0, 0.0, 0.0)
    part = [pid, type, mass, pos, vel]
    new_particles.append(part)
    if pid % 1000 == 0:
        system.storage.addParticles(new_particles, *props)
        system.storage.decompose()
        new_particles = []
    pid += 1
system.storage.addParticles(new_particles, *props)
```

# Acknowledgements
We thank the Espresso++ developer community and all contributors listed in the AUTHORS file. We acknowledge funding from Los Alamos National Laboratory, Forschungszentrum Jülich, and other collaborating institutions.
