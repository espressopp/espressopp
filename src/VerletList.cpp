#define LOG4ESPP_LEVEL_DEBUG

#include "VerletList.hpp"

#include "Real3D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "System.hpp"
#include "Storage.hpp"
#include "BC.hpp"
#include "iterator/CellListAllPairsIterator.hpp"

using namespace espresso;
using namespace espresso::iterator;

LOG4ESPP_LOGGER(VerletList::theLogger, "VerletList");

/*-------------------------------------------------------------*/

VerletList::VerletList(shared_ptr< System > system, real cut) 
{
  LOG4ESPP_INFO(theLogger, "build VerletList, cut = " << cut);

  cutsq = cut * cut;

  CellList cl = system->storage->getRealCells();

  for (CellListAllPairsIterator it(cl);it.isValid(); ++it) {

    checkPair(*it->first, *it->second);
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

