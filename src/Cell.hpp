#ifndef _ESPRESSO_CELL_HPP
#define _ESPRESSO_CELL_HPP
#include <vector>
#include "Particle.hpp"
#include <esutil/ESPPIterator.hpp>

namespace espresso {
  /** A cell is a structure that manages a list of particles and a
      list of other cells that are active geometric neighbor cells.

      A word about the interacting neighbor cells:
      
      In a 3D lattice each cell has 27 neighbors (including
      itself!). Since we deal with pair forces, it is sufficient to
      calculate only half of the interactions (Newtons law: actio =
      reactio). For each cell 13+1=14 neighbors. This has only to be
      done for the inner cells.
      
      Caution: This implementation needs double sided ghost
      communication! For single sided ghost communication one would
      need some ghost-ghost cell interaction as well, which we do not
      need!
	
      It follows: inner cells: #neighbors = 14
      ghost cells:             #neighbors = 0
  */
  struct Cell;

  /** A NeighborCellInfo manages the information about a cell in a
      neighborhood cell list. It stores the cell, whether the cell is
      to be used when a loop over all pairs is done.
  */
  struct NeighborCellInfo {
    NeighborCellInfo(Cell *_cell, bool _useForAllPairs)
      : cell(_cell), useForAllPairs(_useForAllPairs) {}

    NeighborCellInfo(Cell &_cell, bool _useForAllPairs)
      : cell(&_cell), useForAllPairs(_useForAllPairs) {}

    Cell* cell;
    // stores whether the cell is to be used in an all pair loop
    bool useForAllPairs;
  };

  struct CellList : public std::vector< Cell* > {
    typedef esutil::ESPPIterator< std::vector< Cell* > > Iterator;
  };

  struct NeighborCellList : public std::vector< NeighborCellInfo > {
    typedef esutil::ESPPIterator< std::vector< NeighborCellInfo > > Iterator;
  };

  struct LocalCellList : public std::vector< Cell > {
    typedef esutil::ESPPIterator< std::vector< Cell > > Iterator;
  };

  struct Cell {
    ParticleList particles;
    NeighborCellList neighborCells;
  };
}
#endif
