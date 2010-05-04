// ESPP_CLASS
#ifndef _ITERATOR_NEIGHBORCELLLISTITERATOR_HPP
#define _ITERATOR_NEIGHBORCELLLISTITERATOR_HPP

#include "Cell.hpp"

namespace espresso {
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

    // INLINE IMPLEMENTATION
    NeighborCellListIterator::
    NeighborCellListIterator() 
      : ncit(), pit(), useAllPairs(false)
    {}

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
