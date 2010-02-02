#include "VerletList.hpp"

#include "Cell.hpp"
#include "System.hpp"
#include "Storage.hpp"
#include "BC.hpp"

using namespace espresso;

LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/*-------------------------------------------------------------*/

VerletList::VerletList(shared_ptr< System > system, real cut) {
  cutsq = cut * cut;

  std::vector< Cell >& localCells = system->storage->getLocalCells();
 
  bc = system->bc;

  printf ("Build verlet list, here what local cells contain: \n");

  for (size_t c = 0; c < localCells.size(); c++) {

    Cell* localCell = &localCells[c];

    for (size_t p1 = 0; p1 < localCell->particles.size(); p1++) {

        Particle* pt1  = &localCell->particles[p1];
        printf ("Particle %lld is in celll %d\n", pt1->p.id, c);
    }
  }

  for (size_t c = 0; c < localCells.size(); c++) {

    Cell* localCell = &localCells[c];

    // Check for pairs in the local cell

    size_t nparticles = localCell->particles.size();

    LOG4ESPP_INFO(theLogger, "find pairs in local cell " << c <<
                             " with " << nparticles << " particles");

    for (size_t p1 = 0; p1 < nparticles; p1++) {

      Particle* pt1  = &localCell->particles[p1];

      for (size_t p2 = p1 + 1; p2 < nparticles; p2++) {

          checkPair(pt1, &localCell->particles[p2]);
      }
    }

    LOG4ESPP_INFO(theLogger, "find pairs in local cell " << c << " with neighbors");

    // Loop cell neighbors

    std::vector<NeighborCellInfo>& neighborCells = localCell->neighborCells;
 
    for (size_t n = 0; n < neighborCells.size(); n++) {

      Cell* neighborCell = neighborCells[n].cell;

      LOG4ESPP_DEBUG(theLogger, "loop cell pair " << c << " x " << n);
      LOG4ESPP_DEBUG(theLogger, c << " has " << 
          localCell->particles.size() << " particles");
      LOG4ESPP_DEBUG(theLogger, n << " has " <<
          neighborCell->particles.size() << " particles");

      // avoid double cells

      // if (localCell - neighborCell < 0) continue;

      // now build all parir of localCell / neighborCell
 
      for (size_t p1 = 0; p1 < localCell->particles.size(); p1++) {

        Particle* pt1  = &localCell->particles[p1];

        nparticles = neighborCell->particles.size();

        for (size_t p2 = 0; p2 < nparticles; p2++) {

          checkPair(pt1, &neighborCell->particles[p2]);
        }
      }
    }
  }
}

/*-------------------------------------------------------------*/

void VerletList::checkPair(Particle* pt1, Particle* pt2)
{
  double* pos1 = pt1->r.p;
  double* pos2 = pt2->r.p;

  double d[3];

  double distsq;

  bc->getMinimumImageVector(d, distsq, pos1, pos2);

  // d[0] = pos1[0] - pos2[0];
  // d[1] = pos1[1] - pos2[1];
  // d[2] = pos1[2] - pos2[2];

  // double distsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];

  printf("p1 %lld - p2 %lld, distsq = %f\n",
          pt1->p.id, pt2->p.id, distsq);

  if (pt1->p.id >= pt2->p.id) return;

  if (distsq > cutsq) return;

  myList.push_back(std::pair<Particle*, Particle*>(pt1, pt2));
}

VerletList::~VerletList()
{}

