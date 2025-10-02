---
title: 'Espresso++: A Fast and Extensible Molecular Simulation Package for Coarse-Grained Models'
tags:
  - molecular dynamics
  - coarse-grained simulations
  - soft matter
  - parallel computing
  - C++/Python
authors:
  - name: Torsten Stuehn
    affiliation: 1
  - name: Sebastian Eibl
    affiliation: 2
  - name: Markus Rampp
    affiliation: 2
  - name: Nikita Tretyakov
    affiliation: 1
  - name: Tristan Bereau
    affiliation: 1
  - name: Horacio Vargas
    affiliation: 1
  - name: Bin Song
    affiliation: 1
  - name: James Vance
    affiliation: 1
  - name: Pavel Kus
    affiliation: 2
  - name: Others
    affiliation: 2
  - name: Christoph Junghans
    affiliation: 3
affiliations:
  - name: Max Planck Institute for Polymer Research, Mainz, Germany
    index: 1
  - name: Max Planck Computing and Data Facility, Garchingen, Germany
    index: 2
  - name: Los Alamos National Laboratory, Los Alamos, USA
    index: 3

# $ git log --since='Thu Jul 12 22:06:00 2018 +0200' | grep ^Author: | sort | uniq -c | sort -nr
# 201 Author: James Vance <vance@uni-mainz.de>
# 135 Author: Christoph Junghans <junghans@votca.org>
#  74 Author: Sebastian Eibl <XzzX@users.noreply.github.com>
#  49 Author: Sebastian Eibl <sebastian.eibl@mpcdf.mpg.de>
#  35 Author: niktre <niktre@gmail.com>
#  35 Author: Christoph Junghans <christoph.junghans@gmail.com>
#  34 Author: Jakub Krajniak <563684+jkrajniak@users.noreply.github.com>
#  15 Author: Tristan Bereau <bereau@mpip-mainz.mpg.de>
#  10 Author: Nikita Tretyakov <tretyakov@mpip-mainz.mpg.de>
#  10 Author: Jakub Krajniak <jkrajniak@gmail.com>
#   9 Author: hache <vargas@mpip-mainz.mpg.de>
#   8 Author: Bin Song <Bin.Song@mpip-mainz.mpg.de>
#   6 Author: Christoph Junghans <junghans@lanl.gov>
#   5 Author: Pavel Kus <pavel.kus@mpcdf.mpg.de>
#   5 Author: Bin Song <songbin6280@users.noreply.github.com>
#   4 Author: Torsten Stuehn <stuehn@mpip-mainz.mpg.de>
#   4 Author: Horacio Vargas <govarguz@users.noreply.github.com>
#   3 Author: Tristan Bereau <tristan.bereau@gmail.com>
#   2 Author: govarguz <horacio.v.g@gmail.com>
#   2 Author: espressopp-bot <espressopp-bot@users.noreply.github.com>
#   1 Author: Zhen-Hao Xu <zeh026@126.com>
#   1 Author: Pavel Kus <pavel.kus@gmail.com>
#   1 Author: Hache <hache@Haches-MacBook-Air.local>
#   1 Author: Bin Song <songbin@mpip-mainz.mpg.de>

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

# Example usage

A minimal Python script to run a coarse-grained polymer melt simulation in Espresso++ looks like:

```python
import espressopp

system = espressopp.System()
# configure box, particles, and interactions
# run integrator
Full examples and tutorials are available in the online documentation.

Impact
Espresso++ has been applied in numerous scientific studies, including investigations of:

- Polymer rheology and entanglement effects
- Lipid membranes and vesicle dynamics
- Colloidal self-assembly
- Active matter and microswimmers

It is actively maintained and extended by a community of researchers, with contributions from multiple institutions. Its flexible architecture makes it a valuable tool for developing novel coarse-grained models and algorithms in soft matter research.

# Acknowledgements
We thank the Espresso++ developer community and all contributors listed in the AUTHORS file. We acknowledge funding from Los Alamos National Laboratory, Forschungszentrum Jülich, and other collaborating institutions.
