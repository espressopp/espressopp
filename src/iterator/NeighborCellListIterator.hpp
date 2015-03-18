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
#ifndef _ITERATOR_NEIGHBORCELLLISTITERATOR_HPP
#define _ITERATOR_NEIGHBORCELLLISTITERATOR_HPP

#include "Cell.hpp"

namespace espressopp {
  namespace iterator {
    /**
       Iterates all Particles in a list of cells.
    */
    class NeighborCellListIterator {
    public:
      NeighborCellListIterator();
      NeighborCellListIterator(NeighborCellList &cl, 
			       bool _useAllPairs);      

      NeighborCellListIterator &operator++();

      bool isValid() const;
      bool isDone() const;

      Particle &operator*() const;
      Particle *operator->() const;

      Cell &getCurrentCell();

    private:
      void findNonemptyCell();

      NeighborCellList::Iterator ncit;
      ParticleList::Iterator pit;
      bool useAllPairs;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    inline
    NeighborCellListIterator::
    NeighborCellListIterator() 
      : ncit(), pit(), useAllPairs(false)
    {}

    inline
    NeighborCellListIterator::
    NeighborCellListIterator(NeighborCellList &cl, bool _useAllPairs) 
      : ncit(cl), pit(), useAllPairs(_useAllPairs)  {
      // advance to the first valid neighbor list
      if (useAllPairs)
        while (ncit.isValid() && !ncit->useForAllPairs)
          ++ncit;
      
      if (ncit.isDone()) return;
      pit = ParticleList::Iterator(ncit->cell->particles);

      while (pit.isDone()) {
        ++ncit;
        if (useAllPairs)
          while (ncit.isValid() && !ncit->useForAllPairs)
            ++ncit;
        if (ncit.isDone()) break;
        pit = ParticleList::Iterator(ncit->cell->particles);
      }
    }

    inline NeighborCellListIterator &
    NeighborCellListIterator::operator++() {
      ++pit;
      while (pit.isDone()) {
        ++ncit;
        if (useAllPairs)
          while (ncit.isValid() && !(ncit->useForAllPairs))
            ++ncit;
        if (ncit.isDone()) break;
        pit = ParticleList::Iterator(ncit->cell->particles);
      }
      return *this;
    }

    inline bool
    NeighborCellListIterator::isValid() const 
    { return ncit.isValid(); }

    inline bool 
    NeighborCellListIterator::isDone() const 
    { return !isValid(); }

    inline Particle&
    NeighborCellListIterator::operator*() const 
    { return *pit; }

    inline Particle*
    NeighborCellListIterator::operator->() const 
    { return &**this; }

    inline Cell&
    NeighborCellListIterator::getCurrentCell()
    { return *ncit->cell; }
    
  }
}
#endif
