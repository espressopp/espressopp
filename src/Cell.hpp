#ifndef _ESPRESSO_CELL_HPP
#define _ESPRESSO_CELL_HPP
#include <vector>
#include "Particle.hpp"

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

  struct CellList : public std::vector< Cell* > {
    typedef esutil::ESPPIterator< std::vector< Cell* > > Iterator;
  };

  struct LocalCellList : public std::vector< Cell > {
    typedef esutil::ESPPIterator< std::vector< Cell > > Iterator;
  };

  struct Cell {
    ParticleList        particles;
    CellList		neighborCells;
  };
  

}
#endif
