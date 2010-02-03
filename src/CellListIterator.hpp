#ifndef _ESPRESSO_CELLLISTITERATOR_HPP
#define _ESPRESSO_CELLLISTITERATOR_HPP

#include <vector>
#include "Cell.hpp"
#include <iostream>

namespace espresso {
  /**
     Iterates all Particles in a list of cells. This is a Python-like,
     self-contained iterator: isValid() tells whether there are more
     particles to come.
  */
  class CellListIterator {
  public:
    CellListIterator() {}

    CellListIterator(CellList &cl) : cit(cl), pit() {
      if (cit.isDone()) return;
      pit = ParticleList::Iterator((*cit)->particles);
      findNonemptyCell();
    }

    CellListIterator &operator++() {
      ++pit;
      findNonemptyCell();
      return *this;
    }

    bool isValid() const { return cit.isValid(); }
    bool isDone() const { return !isValid(); }

    Particle &operator*() const { return *pit; }
    Particle *operator->() const { return &**this; }

    Cell *getCurrentCell() { return *cit; }

  private:
    void findNonemptyCell();

    CellList::Iterator cit;
    ParticleList::Iterator pit;
  };
}
#endif
