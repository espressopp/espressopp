#ifndef _ESPRESSO_CELLLISTITERATOR_HPP
#define _ESPRESSO_CELLLISTITERATOR_HPP

#include <vector>
#include "Cell.hpp"

namespace espresso {
  /**
     Iterates all Particles in a list of cells. This is a Python-like,
     self-contained iterator: isValid() tells whether there are more
     particles to come.
  */
  class CellListIterator {
  public:
    CellListIterator(std::vector<Cell *> &lst)
      : cCell(lst.begin()), endCell(lst.end()), part(0)
    {
      if (!isValid()) {
	end = 0;
	return;
      }
      end = (*cCell)->particles.size();
      if (part >= end) {
	findNonemptyCell();
      }
    }

    CellListIterator &operator++()
    {
      if (++part >= end) {
	findNonemptyCell();
      }
      return *this;
    }

    bool isValid() const { return cCell != endCell; }

    Particle &operator*() const { return (*cCell)->particles[part]; }
    Particle *operator->() const { return &((*cCell)->particles[part]); }

  private:
    void findNonemptyCell();

    std::vector<Cell *>::iterator cCell, endCell;
    int part, end;
  };
}
#endif
