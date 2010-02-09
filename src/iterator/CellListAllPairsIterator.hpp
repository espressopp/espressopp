#ifndef _ITERATOR_CELLLISTALLPAIRSITERATOR_HPP
#define _ITERATOR_CELLLISTALLPAIRSITERATOR_HPP

#include "Cell.hpp"
#include "iterator/CellListIterator.hpp"
#include "iterator/NeighborCellListIterator.hpp"

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
      ParticlePair current;
      CellListIterator clit;
      NeighborCellListIterator nclit;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    inline 
    CellListAllPairsIterator::
    CellListAllPairsIterator() 
    {}

    inline 
    CellListAllPairsIterator::
    CellListAllPairsIterator(CellList &cl)
      : clit(cl), nclit() {
      if (clit.isDone()) return;
      nclit = NeighborCellListIterator
	(clit.getCurrentCell().neighborCells(), true);
      while (nclit.isDone()) {
	++clit;
	if (clit.isDone()) return;
      }
      current.first = &*clit;
      current.seconds = &*nclit;
    }
      
    inline CellListAllPairsIterator &
    CellListAllPairsIterator::
    operator++() {
      ++nclit;
      while (nclit.isDone()) {
	++clit;
	if (clit.isDone()) return *this;
	nclit = NeighborCellListIterator
	  (clit.getCurrentCell().neighborCells(), true);
      }
      current.first = &*clit;
      current.seconds = &*nclit;
      return *this;
    }
      
    inline bool 
    CellListAllPairsIterator::
    isValid() const 
    { return clit.isValid(); }

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

// void loop(CellList &cl) {
//   // loop over the cell list
//   for (esutil::ESPPIterator< CellList > cit(cl);
//        !cit.isDone(); ++cit) {
//     // loop over the particles in the current cell
//     for (esutil::ESPPIterator< ParticleList > pit((*cit)->particles);
// 	 !pit.isDone(); ++pit) {
//       // loop over the neighbor cells
//       for (esutil::ESPPIterator< CellList > ncit((*cit)->neighborCells);
// 	   !ncit.isDone(); ++cit) {
// 	// compare cell ids
// 	if ((*ncit)->id < (*cit)->id) {
// 	  // do full loop over all particle pairs
// 	  for (esutil::ESPPIterator< ParticleList > npit((*ncit)->particles);
// 	       !npit.isDone(); ++npit) {
// 	    // use *pit, *npit
// 	  }
// 	} else if ((*ncit)->id == (*cit)->id) {
// 	  // this and neighbor cell are identical
// 	  // do half loop over all particles
// 	  esutil::ESPPIterator< ParticleList > npit(pit);
// 	  if (!npit.isDone()) {
// 	    ++npit;
// 	    for (; !npit.isDone(); ++npit) {
	      
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
// }

  }
}


#endif
