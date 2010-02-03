#include "CellListAllPairsIterator.hpp"

namespace espresso {
  CellListAllPairsIterator &
  CellListAllPairsIterator::operator++() {
    // // advance the neighbor cell list iterator
    // ++nclit;
    // while (nclit.isDone()) {
    //   ++clit;
    //   if (clit.isDone()) return *this;
    //   current.first = *clit;
    //   nclit = CellListIterator(clit.getCurrentCell()->neighborCells);
    // }
    // current.second = *nclit;
    return *this;
  }
}
