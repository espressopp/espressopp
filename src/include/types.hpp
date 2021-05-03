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

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/mpi.hpp>
#include <exception>
#include <limits>
#include "esconfig.hpp"

namespace espressopp {
  using boost::shared_ptr;
  using boost::weak_ptr;
  using boost::make_shared;
  using boost::enable_shared_from_this;
  using boost::const_pointer_cast;
  using boost::static_pointer_cast;
  using boost::dynamic_pointer_cast;
  namespace mpi {
    using namespace boost::mpi;
  }

  /* Forward declarations and typedefs. */
  namespace esutil {
    class RNG;
  }

  class Real3D;
  class Int3D;
  class Tensor;

  class Particle;
  class ParticleList;
  class ParticlePair;
  class PairList;

  class ParticleTriple;
  class TripleList;
  class ParticleQuadruple;
  class QuadrupleList;

  class Cell;
  class CellList;
  class NeighborCellList;
  class LocalCellList;

  class VerletList;

  class System;

  namespace storage {
    class Storage;
    class DomainDecomposition;
    class DomainDecompositionNonBlocking;
  }

  namespace bc {
    class BC;
  }

  namespace interaction {
    class Interaction;
    class InteractionList;
  }

  class NoDefault: public std::exception {};
}

#endif
