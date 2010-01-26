#ifndef _CELLLISTALLPAIRSITERATOR_HPP
#define _CELLLISTALLPAIRSITERATOR_HPP

#include <utility>
#include "Cell.hpp"
#include "esutil/ESPPIterator.hpp"

namespace espresso {
  typedef std::pair< Particle&, Particle& > ParticlePair;

  class CellListAllPairsIterator {
  public:
    CellListAllPairsIterator();
    CellListAllPairsIterator(CellListAllPairsIterator &clni);
    CellListAllPairsIterator(CellList &cl);

    CellListAllPairsIterator &operator++();

    bool isValid() const;
    bool isDone() const;

    ParticlePair &operator*() const;
    ParticlePair *operator->() const { return &(**this); }
   };

}

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



#endif
