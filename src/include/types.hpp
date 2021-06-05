/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <boost/mpi.hpp>
#include <exception>
#include <limits>
#include <memory>
#include "esconfig.hpp"

namespace espressopp
{
using std::const_pointer_cast;
using std::dynamic_pointer_cast;
using std::enable_shared_from_this;
using std::make_shared;
using std::shared_ptr;
using std::static_pointer_cast;
using std::weak_ptr;
namespace mpi
{
using namespace boost::mpi;
}

/* Forward declarations and typedefs. */
namespace esutil
{
class RNG;
}

class Real3D;
class Int3D;
class Tensor;

struct Particle;
struct ParticleList;
class ParticlePair;
struct PairList;

class ParticleTriple;
struct TripleList;
class ParticleQuadruple;
struct QuadrupleList;

struct Cell;
struct CellList;
struct NeighborCellList;
struct LocalCellList;

class VerletList;

class System;

namespace storage
{
class Storage;
class DomainDecomposition;
class DomainDecompositionNonBlocking;
}  // namespace storage

namespace bc
{
class BC;
}

namespace interaction
{
class Interaction;
struct InteractionList;
}  // namespace interaction

namespace vectorization
{
class ParticleArray;
class Vectorization;
class VerletList;
}  // namespace vectorization

class NoDefault : public std::exception
{
};

template <typename T>
int int_c(const T& val)
{
    return static_cast<T>(val);
}
}  // namespace espressopp

#endif
