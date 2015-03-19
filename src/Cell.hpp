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

// ESPP_CLASS
#ifndef _ESPRESSO_CELL_HPP
#define _ESPRESSO_CELL_HPP
#include <vector>
#include "Particle.hpp"
#include <esutil/ESPPIterator.hpp>

namespace espressopp {
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

  struct CellList 
    : public esutil::ESPPContainer < std::vector< Cell* > > 
  {};

  struct NeighborCellList 
    : public esutil::ESPPContainer< std::vector< NeighborCellInfo > > 
  {};

  struct LocalCellList 
    : public esutil::ESPPContainer< std::vector< Cell > > 
  {};

  struct Cell {
    ParticleList particles;
    NeighborCellList neighborCells;
  };
}
#endif
