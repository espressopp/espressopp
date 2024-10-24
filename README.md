# ESPResSo++

![Build Status](https://github.com/espressopp/espressopp/actions/workflows/validate.yml/badge.svg?branch=master)
[![Code Climate](https://codeclimate.com/github/espressopp/espressopp/badges/gpa.svg)](https://codeclimate.com/github/espressopp/espressopp)
[![codecov](https://codecov.io/gh/espressopp/espressopp/branch/master/graph/badge.svg?token=gx8YKTpfcR)](https://codecov.io/gh/espressopp/espressopp)

ESPResSo++ is an extensible, flexible, fast and parallel simulation software for
soft matter research. It is a highly versatile software package for the
scientific simulation and analysis of coarse-grained atomistic or bead-spring
models as they are used in soft matter research. ESPResSo and ESPResSo++ have
common roots and share parts of the developer/user community. However their
development is independent and they are different software packages. ESPResSo++
is free, open-source software published under the GNU General Public License
(GPL).

# Quick start:

To get a copy of the developer version (most recent version) of ESPResSo++, you can use git or docker. Using docker will give you a binary release (nothing to compile, but performance may not be optimal). If you use git clone or download a tarball, you will have to compile ESPResSo+ yourself, which might lead to better performance.

Using [docker](https://www.docker.com):
```sh
$ docker pull espressopp/espressopp
$ docker run -it espressopp/espressopp /bin/bash
```

Using git:
```sh
$ git clone https://github.com/espressopp/espressopp.git
```

Alternatively, you can download a tarball or zip file of [previous release versions](https://github.com/espressopp/espressopp/releases) of ESPResSo++.

# Dependencies
## C++ Dependencies
 - Boost ( >= 1.69.0)
 - MPI
 - FFTW3
 - GROMACS (required when `WITH_XTC` flag is enabled, GROMACS needs to be built with GMX_INSTALL_LEGACY_API)
 - [HPX 1.9.1](https://hpx-docs.stellar-group.org/tags/1.9.1/html/quickstart.html) (required when `HPX4ESPP_ENABLED=ON`)
 - HDF5

## Python Dependencies
ESPResSo++ requires Python 3.7 or newer. All required Python packages are listed in `requirements.txt`. You can install them via: `pip3 install -r requirements.txt`

# Quick install:

```sh
$ cd espressopp
$ cmake -B builddir -DCMAKE_INSTALL_PREFIX=/where/to/install/espressopp .
$ cmake --build builddir
$ cmake --install builddir
$ export PYTHONPATH=/where/to/install/espressopp/lib/python3*/site-packages:${PYTHONPATH}
```

After building go to the `examples` directory and have a look at the Python scripts.

You can also use [Pipenv](https://github.com/pypa/pipenv), simply after compilation call in the root directory
```sh
$ pipenv install
$ pipenv shell
```

then you can go to `examples` and have a look at the Python scripts.

## CMake options

You can customize the build process by applying following CMake flags

 - `WITH_XTC` - build E++ with support of dumping trajectory to GROMACS xtc files (default: OFF).
 - `CMAKE_INSTALL_PREFIX` - where the E++ should be installed.
 - `CMAKE_CXX_FLAGS` - put specific compilation flags.
 - `HPX4ESPP_ENABLED` - build E++ with hpx4espp module which integrates the HPX runtime (default: OFF).
 - `HPX_DIR` - path to `.../lib64/cmake/HPX` directory which is required when `HPX4ESPP_ENABLED=ON`.

Then, the flags can be used in `cmake`

```sh
$ cmake -B builddir -DWITH_XTC=ON -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_CXX_FLAGS=-O3 -DHPX4ESPP_ENABLED=ON -DHPX_DIR=/path/to/hpx/lib64/cmake/HPX .
$ make
```

## How to install E++ in some Linux distributions

### Ubuntu

```sh
$ apt-get -qq install -y build-essential openmpi-bin libfftw3-dev python3-dev libboost-all-dev git python3-mpi4py cmake wget python3-numpy ipython3 clang llvm ccache python3-pip doxygen sphinx-common python3-matplotlib graphviz texlive-latex-base texlive-latex-extra texlive-latex-recommended ghostscript libgromacs-dev clang-format curl latexmk libhdf5-dev python3-h5py sudo

$ cd espressopp
$ cmake -B builddir .
$ cmake --build builddir
```

### Fedora

```sh
$ dnf install -y make cmake wget git gcc-c++ doxygen python-devel openmpi-devel environment-modules python-pip clang llvm compiler-rt ccache findutils boost-devel boost-python3-devel python-sphinx fftw-devel python-matplotlib texlive-latex-bin graphviz boost-openmpi-devel ghostscript python3-mpi4py-openmpi texlive-hyphen-base texlive-cm texlive-cmap texlive-ucs texlive-ec gromacs-devel hwloc-devel lmfit-devel ocl-icd-devel hdf5-devel python-h5py atlas hdf5 liblzf python-six python-nose python-numpy
$ cd espressopp
$ cmake -B builddir .
$ cmake --build builddir
```

# Documentation

http://espressopp.github.io

# Reporting issues

Report bugs on the [GitHub issues site](https://github.com/espressopp/espressopp/issues)

# Using the hpx4espp module

## Installation

Download the latest release [HPX 1.9.1](https://github.com/STEllAR-GROUP/hpx/releases/tag/1.9.1).

Follow the installation procedure in [this documentation](https://hpx-docs.stellar-group.org/tags/1.9.1/html/quickstart.html).

Build espressopp:

```sh
$ cmake . -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_CXX_FLAGS=-O3 -DHPX4ESPP_ENABLED=ON -DHPX_DIR=/path/to/hpx/lib64/cmake/HPX
$ make
```

## Examples

Run the example script found in `examples/hpx4espp` using the following options:
```sh
[SLURM/MPI commands] python lennard_jones.py --hpx:threads=16 --hpx:print-bind [--no-hpx4espp] [--no-hpxStart]
```

For example, if you are running with SLURM on a node with 16 Skylake cores and want to compare across three versions
```
ORIG:
$ srun --nodes=1 -n 16 --cpus-per-task=2 python lennard_jones.py --no-hpx4espp --no-hpxStart

Vectorized kernels on hpx4espp w/o threading:
$ srun --nodes=1 -n 16 --cpus-per-task=2 python lennard_jones.py --hpx:threads=1 --no-hpxStart

HPX:
$ srun --nodes=1 -n 1 --cpus-per-task=32 python lennard_jones.py --hpx:threads=16 --hpx:print-bind
```
