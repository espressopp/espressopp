#include "CellListIterator.hpp"

using namespace boost;
using namespace espresso;

void CellListIterator::findNonemptyCell() {
  
  while (pit.isDone()) {
    ++cit;
    if (cit.isDone()) break;
    pit = ParticleList::Iterator((*cit)->particles);
  }
}

