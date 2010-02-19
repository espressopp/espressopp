#ifndef _ITERATOR_CELLLISTALLPAIRSITERATOR_HPP
#define _ITERATOR_CELLLISTALLPAIRSITERATOR_HPP

#include "Cell.hpp"
#include "log4espp.hpp"
#include "iterator/CellListIterator.hpp"
#include "iterator/NeighborCellListIterator.hpp"
#include <cassert>

namespace espresso {
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
    LOG4ESPP_LOGGER(CellListAllPairsIterator::theLogger, 
			   "CellListAllPairsIterator");

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
		     "current pair: (" << current.first->p.id <<
		     ", " << current.second->p.id << ")"
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
