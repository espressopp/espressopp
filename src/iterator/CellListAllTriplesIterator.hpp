// ESPP_CLASS
#ifndef _ITERATOR_CELLLISTALLTRIPLESITERATOR_HPP
#define _ITERATOR_CELLLISTALLTRIPLESITERATOR_HPP

#include "Cell.hpp"
#include "log4espp.hpp"
#include "iterator/CellListIterator.hpp"
#include "iterator/NeighborCellListIterator.hpp"

namespace espresso {
  namespace iterator {
    class CellListAllTriplesIterator {
    public:
      CellListAllTriplesIterator();
      CellListAllTriplesIterator(CellList &cl);
      
      CellListAllTriplesIterator &operator++();
      
      bool isValid() const;
      bool isDone() const;
      
      const ParticleTriple &operator*() const;
      const ParticleTriple *operator->() const;
      
    private:
      static LOG4ESPP_DECL_LOGGER(theLogger);

      ParticleTriple current;

      bool inSelfLoop1;
      bool inSelfLoop2;

      // current cell1
      CellList::Iterator cit1;
      // current particle1
      ParticleList::Iterator pit1;

      // current cell1
      NeighborCellList::Iterator cit2;
      // current particle2
      ParticleList::Iterator pit2;

      // current cell3
      NeighborCellList::Iterator cit3;
      // current particle3
      ParticleList::Iterator pit3;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    inline 
    CellListAllTriplesIterator::
    CellListAllTriplesIterator()
    {}

    inline 
    CellListAllTriplesIterator::
    CellListAllTriplesIterator(CellList &cl) {

      cit1 = CellList::Iterator(cl);
      if (cit1.isDone()) return;
      inSelfLoop1 = true;
      pit1 = ParticleList::Iterator((*cit1)->particles);
      while (pit1.isDone()) {
        ++cit1;
        if (cit1.isDone()) return;
        pit1 = ParticleList::Iterator((*cit1)->particles);
      }
      
      pit2 = pit1;
      this->operator++();
    }

    inline CellListAllTriplesIterator &
    CellListAllTriplesIterator::
    operator++() {
      ++pit2;
      while (pit2.isDone()) {
        ++pit1;
        while (pit1.isDone()) {
          if (inSelfLoop1) {
            LOG4ESPP_TRACE(theLogger, "pit1.isDone(), inSelfLoop1, starting neighbor loop");
            inSelfLoop1 = false;
            cit2 = NeighborCellList::Iterator((*cit1)->neighborCells);
          } else {
            LOG4ESPP_TRACE(theLogger, "pit1.isDone(), !inSelfLoop1, continuing neighbor loop");
            ++cit2;
          }

          while (cit2.isValid() && cit2->useForAllPairs)
            ++cit2;

          if (cit2.isDone()) {
            LOG4ESPP_TRACE(theLogger, "cit2.isDone(), go to next cell");
            ++cit1;
            if (cit1.isDone()) {
              LOG4ESPP_TRACE(theLogger, "cit1.isDone(), LOOP FINISHED");
              return *this;
            }
            inSelfLoop1 = true;
          }
          //assert(inSelfLoop1 || cit2.isValid());
          pit1 = ParticleList::Iterator((*cit1)->particles);
        }
        //assert(pit1.isValid());

        if (inSelfLoop1) {
          pit2 = pit1;
          ++pit2;
        } else {
          pit2 = ParticleList::Iterator(cit2->cell->particles);
        }
      }

      // WARNING : CellListAllTriplesIterator not finished yet DOES NOT WORK !
      // TODO: This still has to be implemented !
      pit3 = pit2;

      current.first  = &*pit1;
      current.second = &*pit2;
      current.third  = &*pit3;

      LOG4ESPP_TRACE(theLogger,
             "current triple: (" << current.first->p.id <<
             ", " << current.second->p.id << ", " << current.third->p.id << ")"
             );
      return *this;
    }

    inline bool 
    CellListAllTriplesIterator::
    isValid() const 
    { return cit1.isValid(); }

    inline bool 
    CellListAllTriplesIterator::
    isDone() const 
    { return !isValid(); }
      
    inline const ParticleTriple &
    CellListAllTriplesIterator::
    operator*() const 
    { return current; }
    
    inline const ParticleTriple *
    CellListAllTriplesIterator::
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
