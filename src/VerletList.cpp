#define LOG4ESPP_LEVEL_DEBUG

#include "VerletList.hpp"

#include "Real3D.hpp"
#include "Particle.hpp"
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

    LOG4ESPP_DEBUG(theLogger, "find pairs of real cell " << c
		   << " with local index " << localCell - system->storage->getFirstCell()
		   << " with " << nparticles << " particles");

    int c = myList.size();

    for (size_t p1 = 0; p1 < nparticles; p1++) {

      Particle& pt1  = localCell->particles[p1];

      for (size_t p2 = p1 + 1; p2 < nparticles; p2++) {

        checkPair(pt1, localCell->particles[p2]);
      }
      
    }

    LOG4ESPP_INFO(theLogger, "self size " << myList.size() - c);

    // Loop cell neighbors

    std::vector<NeighborCellInfo>& neighborCells = localCell->neighborCells;
 
    LOG4ESPP_DEBUG(theLogger, "find pairs of real cell " << c
		   << " with local index " << localCell - system->storage->getFirstCell()
		   << " with " << neighborCells.size() << " neighbored cells");

    for (size_t n = 0; n < neighborCells.size(); n++) {

      Cell* neighborCell = neighborCells[n].cell;

      int nPNeighbor = neighborCell->particles.size();

      LOG4ESPP_DEBUG(theLogger, "Neighbor " << n
		     << " with local index " << neighborCell - system->storage->getFirstCell()
		     << " has " << neighborCell->particles.size() << " particles");

      if (!neighborCells[n].useForAllPairs || nPNeighbor == 0) continue;

      // now build all parir of localCell / neighborCell

      int c = myList.size();
 
      for (size_t p1 = 0; p1 < nparticles; p1++) {

        Particle& pt1  = localCell->particles[p1];

        for (size_t p2 = 0; p2 < nPNeighbor; p2++) {

          checkPair(pt1, neighborCell->particles[p2]);
        }
      }

      LOG4ESPP_DEBUG(theLogger, "cross size " << myList.size() - c);
    }
  }
}

/*-------------------------------------------------------------*/

void VerletList::checkPair(Particle& pt1, Particle& pt2)
{
  Real3DRef pos1 = pt1.r.p;
  Real3DRef pos2 = pt2.r.p;

  Real3D d;
  real distsq;

  d = pos1; d -= pos2;
  distsq = d.sqr();

  LOG4ESPP_TRACE(theLogger, "p1: " << pt1.p.id << " @ " << pos1
		 << " - p2: " << pt2.p.id << " @ " << pos2
		 << " -> distsq = " << distsq);

  if (distsq > cutsq) return;

  myList.push_back(ParticlePair(&pt1, &pt2));
}

/*-------------------------------------------------------------*/

VerletList::~VerletList()
{
}

