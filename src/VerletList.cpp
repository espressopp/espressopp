#include "Cell.hpp"

#include "System.hpp"
#include "VerletList.hpp"

using namespace espresso;

LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

VerletList::VerletList(System::SelfPtr system, double cut)

{
  cutsq = cut * cut;

  std::vector<Cell>& localCells = system->storage->getLocalCells();
 
  for (size_t c = 0; c < localCells.size(); c++) {

    Cell* localCell = &localCells[c];

    LOG4ESPP_INFO(theLogger, "build verlet list for local cell " << c);

    // Loop cell neighbors

    std::vector<NeighborCellInfo>& neighborCells = localCell->neighborCells;
 
    for (size_t n = 0; n < neighborCells.size(); n++) {

      Cell* neighborCell = neighborCells[n].cell;

      LOG4ESPP_DEBUG(theLogger, "loop cell pair " << c << " x " << n);
      LOG4ESPP_DEBUG(theLogger, c << " has " << localCell->particles.size() << " particles");
      LOG4ESPP_DEBUG(theLogger, n << " has " << neighborCell->particles.size() << " particles");

      // avoid double cells

      if (localCell - neighborCell < 0) continue;

      // now build all parir of localCell / neighborCell
 
      for (size_t p1 = 0; p1 < localCell->particles.size(); p1++) {

        Particle* pt1  = &localCell->particles[p1];
        double*   pos1 = pt1->r.p;

        size_t nparticles = neighborCell->particles.size();

        // Avoid double entries (i,j) and (j,i) in local cell

        if (localCell == neighborCell) {

          // avoid double particle pairs in a cell 

          nparticles = p1;

        }

        for (size_t p2 = 0; p2 < nparticles; p2++) {

          Particle* pt2  = &neighborCell->particles[p2];
          double*   pos2 = pt2->r.p;

          double d[3];

          d[0] = pos1[0] - pos2[0];
          d[1] = pos1[1] - pos2[1];
          d[2] = pos1[2] - pos2[2];

          double distsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];

          if (distsq >= cutsq) continue;

          myList.push_back(std::pair<Particle*, Particle*>(pt1, pt2));
        }
      }
    }
  }
}

VerletList::~VerletList()
{
}

