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
#ifndef _ITERATOR_CELLLISTALLPAIRSITERATOR_HPP
#define _ITERATOR_CELLLISTALLPAIRSITERATOR_HPP

#include "Cell.hpp"
#include "log4espp.hpp"
#include "iterator/CellListIterator.hpp"
//#include "iterator/NeighborCellListIterator.hpp"

namespace espressopp {
  namespace iterator {
    class CellListAllPairsIterator {
    public:
      CellListAllPairsIterator();
      CellListAllPairsIterator(CellList &cl); 
      
      CellListAllPairsIterator &operator++();
      
      bool isValid() const;
      bool isDone() const;
      
      const ParticlePair &operator*() const;
      const ParticlePair *operator->() const;
      
    private:
      static LOG4ESPP_DECL_LOGGER(theLogger);

      ParticlePair current;

      bool inSelfLoop;

      // current cell
      CellList::Iterator cit;
      // current particle
      ParticleList::Iterator pit;

      // current neighbor cell
      NeighborCellList::Iterator ncit;
      // current neighbor particle
      ParticleList::Iterator npit;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    inline 
    CellListAllPairsIterator::
    CellListAllPairsIterator()
    {}

    inline 
    CellListAllPairsIterator::
    CellListAllPairsIterator(CellList &cl) {
      cit = CellList::Iterator(cl);
      if (cit.isDone()) return;
      inSelfLoop = true;
      pit = ParticleList::Iterator((*cit)->particles);
      while (pit.isDone()) {
        ++cit;
        if (cit.isDone()) return;
        pit = ParticleList::Iterator((*cit)->particles);
      }
      
      npit = pit; 
      this->operator++();
    }

    inline CellListAllPairsIterator &
    CellListAllPairsIterator::
    operator++() {
      ++npit;
      while (npit.isDone()) {
        ++pit;
        while (pit.isDone()) {
          if (inSelfLoop) {
            LOG4ESPP_TRACE(theLogger, "pit.isDone(), inSelfLoop, starting neighbor loop");
            inSelfLoop = false;
            ncit = NeighborCellList::Iterator((*cit)->neighborCells);
          } else {
            LOG4ESPP_TRACE(theLogger, "pit.isDone(), !inSelfLoop, continuing neighbor loop");
            ++ncit;
          }

          while (ncit.isValid() && ncit->useForAllPairs)
            ++ncit;

          if (ncit.isDone()) {
            LOG4ESPP_TRACE(theLogger, "ncit.isDone(), go to next cell");
            ++cit;
            if (cit.isDone()) {
              LOG4ESPP_TRACE(theLogger, "cit.isDone(), LOOP FINISHED");
              return *this;
            }
            inSelfLoop = true;
          }
          //assert(inSelfLoop || ncit.isValid());
          pit = ParticleList::Iterator((*cit)->particles);
        }
        //assert(pit.isValid());

        if (inSelfLoop) {
          npit = pit;
          ++npit;
        } else {
          npit = ParticleList::Iterator(ncit->cell->particles);
        }
      }

      current.first = &*pit;
      current.second = &*npit;

      LOG4ESPP_TRACE(theLogger,
             "current pair: (" << current.first->id() <<
             ", " << current.second->id() << ")"
             );
      return *this;
    }

    inline bool 
    CellListAllPairsIterator::
    isValid() const 
    { return cit.isValid(); }

    inline bool 
    CellListAllPairsIterator::
    isDone() const 
    { return !isValid(); }
      
    inline const ParticlePair &
    CellListAllPairsIterator::
    operator*() const 
    { return current; }
    
    inline const ParticlePair *
    CellListAllPairsIterator::
    operator->() const 
    { return &(**this); }

  }
}

    // foreach(CellList &cl) {
    //   // loop over all cells
    //   for (cit = CellList::Iterator(cl); cit.isValid(); ++cit) {
    // 	// loop over pairs of particles where both particles are in
    // 	// the local cell
    // 	for (pit = ParticleList::Iterator(*cit); pit.isValid(); ++pit) {
    // 	  // loop over particles in the cell itself
    // 	  npit = pit; 
    // 	  ++npit;
    // 	  for (; npit.isValid(); ++npit)
    // 	    yield(*pit, *npit);
    // 	}

    // 	// now loop over pairs of particles in the cell that involve
    // 	// neighbor cells
    // 	for (ncit = NeighborCellList::Iterator(cit->neighborCells);
    // 	     ncit.isValid(); ++ncit) {
    // 	  if (ncit->useForAllPairs)
    // 	    for (pit = ParticleList::Iterator(*cit); pit.isValid(); ++pit)
    // 	      for (npit = ParticleList::Iterator(ncit->cell->particles); 
    // 		   npit.isValid(); ++npit)
    // 		yield(*pit, *npit);
    // 	}
    //   }
    // }

#endif
