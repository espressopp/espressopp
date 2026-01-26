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
- Multi-scale simulation techniques such as AdResS and Lees-Edwards

Espresso++ is widely used in academic research for simulating phenomena such as polymer rheology, membrane dynamics, colloidal suspensions, and active particles. Recent research projects using ESPResSo++ include [@Grommes:2025]

[@Gholami2025]

[@Grommes2024]

[@Hsu2024]

[@Hsu2023]

[@Ohkuma2023]

[@Grommes2022]

[@Grommes2021]

[@Brunk2021]

[@Zhang2021]

[@Tubiana2021]

[@Thaler2020]

[@Hsu2020]

[@Fiorentini2020]

[@Pape2023]

[@Smith2023]

[@Bause2021]

[@Rudzinski2020]

[@Singh2020]

[@Zhao2020]

[@Zhao2020b]

[@Lee2020]

[@Grommes2020]

## Merge into above
Espresso++ has been applied in numerous scientific studies, including investigations of:

- Polymer rheology and entanglement effects
- Lipid membranes and vesicle dynamics
- Colloidal self-assembly
- Active matter and microswimmers

It is actively maintained and extended by a community of researchers, with contributions from multiple institutions. Its flexible architecture makes it a valuable tool for developing novel coarse-grained models and algorithms in soft matter research.


# Functionality

Key features of Espresso++ include:

- **Inter-particle interactions**: Lennard-Jones, Coulomb, soft repulsive, bonded interactions, tabulated potentials, and more.
- **Algorithms**: Molecular dynamics, Langevin dynamics, dissipative particle dynamics (DPD), Brownian dynamics, adpative resolution simulations (AdResS), Monte Carlo sampling.
- **Electrostatics**: Particle–particle particle–mesh (P3M), Ewald summation, and other long-range methods.
- **Parallelization**: Domain decomposition using MPI, optimized for massively parallel architectures.
- **Python interface**: Full simulation control and analysis scripting in Python.
- **Extensibility**: Modular design allows easy addition of new force fields, integrators, or analysis tools.

# New Feature since last release

Since the last major release of Espresso++ v2.0 in 2018 a number of new functionalities and features have been added, including:

- **SIMD vectorization and related optimizations**: enhance compute performance on modern CPUs [@Vance2023]
- **Cell decomposition**: allow sub-demcomposition into cells with a lenght of half or a third of the cutoff for direct force calaculations [@Yao]
- **HeSpaDDA**: heterogeneous spatial domain decomposition algorithm (HeSpaDDA) for \dots [@Guzman:2017]
- **new potentials and simulation methods**: AngularCosineSquared, TabulatedSubEnsAngular, surface hopping MD 
- **checkpoint the state of the random number generator (RNG)**: allow restaring from checkpointed state of RNG
- **I/O**: support for parallel writing and reading of H5MD checkpoints
- **Python 3 compatibility**

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

TODO add integrate() call
```

# Acknowledgements
We thank the Espresso++ developer community and all contributors listed in the AUTHORS file.
Espresso++ project is supported by the U.S. Department of Energy through Los Alamos National Laboratory (LANL). Los Alamos National Laboratory is operated by Triad National Security, LLC, for the National Nuclear Security Administration of the U.S. Department of Energy (contract no. 89233218CNA000001). This paper has been assigned a Los Alamos Unlimited Release number of LA-UR-26-XXXXX.
