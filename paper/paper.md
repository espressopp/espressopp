---
title: 'ESPResSo++: A Fast and Extensible Molecular Simulation Package for Coarse-Grained Models'
tags:
  - molecular dynamics
  - coarse-grained simulations
  - soft matter
  - parallel computing
  - C++/Python
authors:
  - name: Zhen-Hao Xu
    affiliation: 4
  - name: James Vance
    orcid: 0000-0001-7112-0382
    affiliation: 4
  - name: Nikita Tretyakov
    affiliation: 4
  - name: Sebastian Eibl
    orcid: 0000-0002-1069-2720
    affiliation: 2
  - name: Pavel Kus
    affiliation: 2
  - name: Jakub Krajniak
    orcid: 0000-0001-9372-6975
    affiliation: 6
  - name: Tristan Bereau
    orcid: 0000-0001-9945-1271
    affiliation: 5
  - name: Horacio Vargas Guzman
    orcid: 0000-0003-2564-3005
    affiliation: 1
  - name: Bin Song
    orcid: 0000-0003-4229-9242
    affiliation: 1
  - name: Markus Rampp
    orcid: 0000-0001-8177-8698
    affiliation: 2  
  - name: Torsten Stuehn
    orcid: 0009-0006-2144-2002
    affiliation: 1
  - name: Christoph Junghans
    orcid: 0000-0003-0925-1458
    affiliation: 3
affiliations:
  - name: Max Planck Institute for Polymer Research, Mainz, Germany
    index: 1
  - name: Max Planck Computing and Data Facility, Garching, Germany
    index: 2
  - name: Los Alamos National Laboratory, Los Alamos, USA
    index: 3
  - name: Johannes Gutenberg University of Mainz, Mainz, Germany
    index: 4
  - name: Heidelberg University, Heidelberg, Germany
    index: 5
  - name: Independent researcher, Poznań, Poland 
    index: 6


date: 2 October 2025
bibliography: paper.bib
---

# Summary

**ESPResSo++** is an open-source software package for molecular dynamics (MD) simulations with a particular emphasis on coarse-grained (CG) models of soft matter systems[@Praprotnik2008]. Written in C++ with a flexible Python interface, it is designed for high-performance computing (HPC) environments and supports massively parallel simulations through MPI. The package enables simulations of polymers, membranes, colloids and complex fluids with a wide range of interaction models and advanced algorithms.

ESPResSo++ builds upon the experience of its predecessor [ESPResSo](https://espressomd.org), but provides a cleaner, more modular codebase and enhanced extensibility. It is actively developed by an international community of researchers across multiple institutions and disciplines.

# Statement of need

Molecular dynamics simulations are essential tools for exploring the behavior of soft matter systems at mesoscopic scales. General-purpose MD codes such as GROMACS [@Abraham2015] and LAMMPS [@Thompson2022] are primarily optimized for all-atom simulations and may require significant customization for coarse-grained models that need non-standard interactions or specialized multiscale algorithms. ESPResSo++ addresses this gap by providing:

- A modular and extensible design, enabling researchers to easily implement new interaction potentials and integrators.
- Efficient parallelization for large-scale simulations of complex systems.
- A rich library of coarse-grained interaction models and algorithms tailored to soft matter.
- A Python-based scripting interface for ease of use, reproducibility, and coupling with external analysis tools.
- Multi-scale simulation techniques such as AdResS and Lees-Edwards.

# State of the field

Several established MD packages serve the soft matter and coarse-grained simulation communities. **LAMMPS** [@Thompson2022] is a general-purpose, highly scalable code with broad force field support and a flexible plugin architecture; however, its input scripting language is less suited to complex CG workflows, and adaptive resolution is not a core feature. **GROMACS** [@Abraham2015] provides exceptional performance for standard force fields and supports the Martini CG model ecosystem, but adding truly novel CG interactions requires modifying the C++ source, and its file-based workflow is less interactive than a Python API. **HOOMD-blue** [@Anderson2020] is the most directly comparable package, offering a Python-native API and GPU-accelerated CG simulations for soft matter; however, it lacks support for adaptive resolution methods. **OpenMM** [@Eastman2017] excels at GPU-accelerated simulations and allows custom force definitions via algebraic expressions, but targets primarily biomolecular systems and lacks MPI-based multi-node parallelization. The sibling project **ESPResSo** [@Weik2019] shares common roots but focuses on charged systems, hydrodynamic coupling (lattice Boltzmann), and electrokinetics rather than multiscale resolution bridging.

ESPResSo++ occupies a distinct niche as a package purpose-built for coarse-grained and multiscale soft matter simulations. Its key differentiator is first-class support for adaptive resolution simulation (AdResS), which allows seamless coupling of atomistic and coarse-grained representations within a single simulation. ESPResSo++ is currently the only actively maintained MD package that provides AdResS as a core feature; none of the other packages listed above offer this capability. It is worth noting that these capabilities were also contributed to existing general-purpose codes[@Fritsch2012; @Nagarajan2013; @Junghans2010], but got removed again later, hence a dedicated package allows tighter integration of multiscale algorithms with the domain decomposition, neighbor list construction, and force computation pipeline — design choices that would be difficult to retrofit into architectures optimized for all-atom throughput.

# Software design

ESPResSo++ uses a layered C++/Python architecture. The performance-critical computational kernel (force calculations, integration, domain decomposition, and communication) is implemented in C++17 and exposed to Python through Boost.Python bindings. Each C++ class provides a static `registerPython()` method, and a central registration dispatcher ensures that the full object hierarchy is available as the `espressopp` Python module. This design lets users assemble, configure, and control simulations entirely from Python while retaining the performance of compiled C++ for the inner loops.

**Modularity through templates and extensions.** Interactions are implemented using C++ class templates parameterized by a potential type (e.g.,\
`VerletListInteractionTemplate<_Potential>`), so that adding a new pairwise potential requires only implementing the `computeEnergy()` and `computeForce()` methods of a `Potential` subclass. The integrator follows a similar pattern: an `MDIntegrator` base class exposes an `addExtension()` mechanism through which thermostats, barostats, constraints, and adaptive resolution layers are attached via Boost.Signals2 callbacks. This signal-based coupling keeps extensions decoupled from the integration loop and from each other.

**Parallelization.** ESPResSo++ employs spatial domain decomposition with MPI. The simulation box is partitioned across MPI ranks using a `NodeGrid`, and each rank further subdivides its domain into linked cells via a `CellGrid`. Ghost particle exchange handles communication of boundary data. A non-blocking variant (`DomainDecompositionNonBlocking`) overlaps communication with computation, and the heterogeneous spatial domain decomposition algorithm (HeSpaDDA) [@Guzman:2017] provides load balancing for spatially inhomogeneous systems.

**Performance optimizations.** Since version 2.0 [@Guzman2019], ESPResSo++ has been modernized with SIMD vectorization through a structure-of-arrays (SOA) particle data layout (`ParticleArray`) with 64-byte alignment, yielding an overall three-times speedup for short-range non-bonded force calculations [@Vance:2023]. An improved cell decomposition scheme [@Yao:2004] allows sub-decomposition into cells with a length of half or a third of the cutoff radius, reducing the number of unnecessary distance calculations.

**Testing and documentation.** The code is tested through a combination of Boost.Test (C++) and Python `unittest` test suites, executed via CMake/CTest and continuous integration on GitHub Actions. User documentation is built with Sphinx; developer API documentation with Doxygen.

# Research impact statement

ESPResSo++ has enabled research across a broad range of soft matter topics over more than a decade of active development. It has been used in numerous peer-reviewed studies, including investigations of:

- Polymer rheology and entanglement effects [@Grommes:2025; @Grommes2024; @Hsu2023; @Hsu2024; @Ohkuma2023; @Grommes2022; @Grommes2021; @Tubiana2021; @Hsu2020; @Singh2020; @Zhao2020b; @Lee2020; @Grommes2020]
- Lipid membranes, protein and vesicle dynamics [@Pape2023; @Bause2021; @Zhao2020]
- Adaptive resolution simulations [@Thaler2020; @Fiorentini2020]
- Ionic liquids under shear flow [@Gholami2025; @Zhang2021]
- Coarse-grained conformational surface hopping [@Rudzinski2020]

The project is developed across multiple institutions — Max Planck Institute for Polymer Research, Max Planck Computing and Data Facility, Los Alamos National Laboratory, Johannes Gutenberg University of Mainz, and Heidelberg University — with contributions from both current and former developers documented in the AUTHORS file. The codebase has a public development history spanning more than ten years on GitHub, with continuous integration, code coverage tracking, and versioned releases.

# Functionality

Key features of ESPResSo++ include:

- **Inter-particle interactions**: Lennard-Jones, Coulomb, soft repulsive, bonded interactions, tabulated potentials, and more.
- **Algorithms**: Molecular dynamics, Langevin dynamics, dissipative particle dynamics (DPD), Brownian dynamics, adaptive resolution simulations (AdResS), Monte Carlo sampling.
- **Electrostatics**: Particle–particle particle–mesh (P3M), Ewald summation, and other long-range methods.
- **Parallelization**: Domain decomposition using MPI, optimized for massively parallel architectures.
- **Python interface**: Full simulation control and analysis scripting in Python.
- **Extensibility**: Modular design allows easy addition of new force fields, integrators, or analysis tools.

# New Features since last release

Since the last major release of ESPResSo++ v2.0 in 2018 [@Guzman2019] a number of new functionalities and features have been added, including:

- **SIMD vectorization and related optimizations**: enhance compute performance on modern CPUs [@Vance:2023]
- **Cell decomposition**: allow sub-decomposition into cells with a length of half or a third of the cutoff for direct force calculations [@Yao:2004]
- **HeSpaDDA**: heterogeneous spatial domain decomposition algorithm (HeSpaDDA) for larger scale simulations [@Guzman:2017]
- **new potentials and simulation methods**: AngularCosineSquared, TabulatedSubEnsAngular, surface hopping MD, Lees-Edwards boundary conditions 
- **Checkpoint the state of the random number generator (RNG)**: allow restarting from checkpointed state of RNG
- **I/O**: support for parallel writing and reading of H5MD checkpoints
- **Python 3 compatibility**

# Example usage

A minimal Python script to run a simple (repulsive only) Lennard-Jones type particle simulation in ESPResSo++ looks like:

```python
import espressopp
from mpi4py import MPI

# simulation system parameters
num_particles = 10000      # total number of particles in the system
box           = (20,20,20) # size of the simulationbox (all length are in sigma)
rc            = 1.12246    # cut off for the short range non bonded potential
skin          = 0.3        # skin used for verlet neighbor list
dt            = 0.005      # time step for 1 md step
epsilon       = 1.0        # energy unit
sigma         = 1.0        # length unit
temperature   = 1.0        # temperature of the simulation
LJcaprad      = 0.8        # inital capping radius for LJ interaction
                           # for random configurations

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
interaction    = espressopp.interaction.VerletListLennardJonesCapped(
                   espressopp.VerletList(system, cutoff=rc))
interaction.setPotential(type1=0, type2=0,
                         potential=espressopp.interaction.LennardJonesCapped(
                           epsilon, sigma, rc, shift='auto', caprad=LJcaprad))
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

integrator.dt = 0.0001 # very small timestep for initial warmup
for n in range(20):
  # warmup finished, switch to uncapped Lennard Jones potential and increase timestep
  if n == 10: 
    interaction = espressopp.interaction.VerletListLennardJones(
                    espressopp.VerletList(system, cutoff=rc))
    interaction.setPotential(type1=0, type2=0,
                             potential=espressopp.interaction.LennardJones(
                               epsilon, sigma, rc, shift='auto'))
    system.removeInteraction(0)
    system.addInteraction(interaction)
    integrator.dt = 0.005
  integrator.run(10000)
  Etot = system.getInteraction(0).computeEnergy()
  print("md time = {:4.1f},
        total energy of the system: {:10.2f}".format(integrator.dt*n*10000, Etot))

# write PDB file of (quasi) equilibrated LJ system.
# At this temperature it is more or less crystallized.
espressopp.tools.pdbwrite("simplelj.pdb", system, molsize=num_particles)

```

# AI usage disclosure

No generative AI tools were used in the creation of the ESPResSo++ software or its documentation. During the preparation of this manuscript, AI-assisted tools were used for copyediting and formatting. All content was reviewed and verified by the authors.

# Acknowledgements
We thank the ESPResSo++ developer community and all contributors listed in the AUTHORS file. ESPResSo++ has been supported by the Transregio TRR146 of the German Research Foundation. The ESPResSo++ project is supported by the U.S. Department of Energy through Los Alamos National Laboratory (LANL). Los Alamos National Laboratory is operated by Triad National Security, LLC, for the National Nuclear Security Administration of the U.S. Department of Energy (contract no. 89233218CNA000001). This paper has been assigned a Los Alamos Unlimited Release number of LA-UR-26-XXXXX.
