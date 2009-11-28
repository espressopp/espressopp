#include "CellListIterator.hpp"

using namespace boost;
using namespace espresso;

void CellListIterator::findNonemptyCell()
{
  part = 0;
  while (++cCell != endCell) {
    end = (*cCell)->particles.size();
    if (end > 0) break;
  }      
}

