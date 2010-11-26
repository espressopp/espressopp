#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/mpi.hpp>
#include <exception>
#include <limits>
#include "esconfig.hpp"

namespace espresso {
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
