#include "VerletList.hpp"

#include "Cell.hpp"
#include "System.hpp"
#include "Storage.hpp"
#include "BC.hpp"

using namespace espresso;

LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/*-------------------------------------------------------------*/

VerletList::VerletList(shared_ptr< System > system, real cut) 
{
  LOG4ESPP_INFO(theLogger, "build VerletList, cut = " << cut);

  cutsq = cut * cut;

  std::vector< Cell* >& localCells = system->storage->getRealCells();
 
  bc = system->bc;

  for (size_t c = 0; c < localCells.size(); c++) {

    Cell* localCell = localCells[c];

    // Check for pairs in the cell itself

    size_t nparticles = localCell->particles.size();

    LOG4ESPP_INFO(theLogger, "find pairs of real cell " << c <<
                             " with " << nparticles << " particles");

    for (size_t p1 = 0; p1 < nparticles; p1++) {

      Particle* pt1  = &localCell->particles[p1];

      for (size_t p2 = p1 + 1; p2 < nparticles; p2++) {

        checkPair(pt1, &localCell->particles[p2]);
      }
    }

    // Loop cell neighbors

    std::vector<NeighborCellInfo>& neighborCells = localCell->neighborCells;
 
    LOG4ESPP_INFO(theLogger, "find pairs of real cell " << c << 
                " with " << neighborCells.size() << " neighbored cells");

    for (size_t n = 0; n < neighborCells.size(); n++) {

      Cell* neighborCell = neighborCells[n].cell;

      int nPNeighbor = neighborCell->particles.size();

      if (nPNeighbor == 0) continue;

      LOG4ESPP_DEBUG(theLogger, "Neighbor " << n << " has " <<
          neighborCell->particles.size() << " particles");

      // avoid double cells

      // if (localCell - neighborCell < 0) continue;

      // now build all parir of localCell / neighborCell
 
      for (size_t p1 = 0; p1 < nparticles; p1++) {

        Particle* pt1  = &localCell->particles[p1];

        for (size_t p2 = 0; p2 < nPNeighbor; p2++) {

          checkPair(pt1, &neighborCell->particles[p2]);
        }
      }
    }
  }
}

/*-------------------------------------------------------------*/

void VerletList::checkPair(Particle* pt1, Particle* pt2)
{
  real* pos1 = pt1->r.p;
  real* pos2 = pt2->r.p;

  real d[3];

  real distsq;

#define WORKAROUND

#ifdef WORKAROUND

  // This gives the right distance, but is much slower

  bc->getMinimumImageVector(d, distsq, pos1, pos2);

#else

  // This should also work, but doesn't yet
  //
  d[0] = pos1[0] - pos2[0];
  d[1] = pos1[1] - pos2[1];
  d[2] = pos1[2] - pos2[2];

  distsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
#endif

  LOG4ESPP_DEBUG(theLogger, "p1: " << pt1->p.id << " - p2: "
             << pt2->p.id << " -> distsq = " << distsq);

  if (pt1->p.id >= pt2->p.id) return;

  if (distsq > cutsq) return;

  myList.push_back(std::pair<Particle*, Particle*>(pt1, pt2));
}

/*-------------------------------------------------------------*/

VerletList::~VerletList()
{
}

