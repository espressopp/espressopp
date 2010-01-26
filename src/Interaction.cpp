
#include "Interaction.hpp"
#include "Storage.hpp"
#include "types.hpp"

using namespace espresso;
using namespace interaction;

LOG4ESPP_LOGGER(Interaction::theLogger, "Interaction");

double Interaction::computeStorageEnergy(shared_ptr<Storage> storage)
{
  std::vector<Cell>& localCells = storage->getLocalCells();

  double e = 0.0;

  for (size_t c = 0; c < localCells.size(); c++) {

    Cell* localCell = &localCells[c];

    ParticleList& list1 = localCell->particles;

    // Loop cell neighbors

    std::vector<NeighborCellInfo>& neighborCells = localCell->neighborCells;

    for (size_t n = 0; n < neighborCells.size(); n++) {

      Cell* neighborCell = neighborCells[n].cell;

      // avoid double cells

      if (localCell - neighborCell < 0) continue;

      // call now the efficient reimplemented method for derived interactions

      if (localCell == neighborCell) {

         e += computeCellEnergy(list1);

      } else {
      
         ParticleList& list2 = neighborCell->particles;
         e += computeCellEnergy(list1, list2);
      }
    }
  }

  return e;
}
