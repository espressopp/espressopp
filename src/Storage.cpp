#include <algorithm>

#define LOG4ESPP_SHORTNAMES
#define LOG4ESPP_LEVEL_DEBUG
#include "log4espp.hpp"

#include "Storage.hpp"
#include "CellListIterator.hpp"

using namespace boost;
using namespace espresso;

LOG4ESPP_LOGGER(Storage::logger, "Storage");


const int STORAGE_COMM_TAG = 0xaa;

Storage::Storage(System *_system,
                 const boost::mpi::communicator &_comm,
                 bool _useVList)
  : comm(_comm), system(_system), useVList(_useVList)
{
}

Storage::~Storage() {}

longint Storage::getNActiveParticles() const {
  longint cnt = 0;
  for (std::vector<Cell *>::const_iterator it = activeCells.begin(),
	 end = activeCells.end();
       it != end; ++it) {
    longint size = (*it)->getNParticles();
    if (size) {
      LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
    }
    cnt += size;
  }
  return cnt;
}

void Storage::updateLocalParticles(ParticleList &l) {
  for (ParticleList::iterator
	 it  = l.begin(),
	 end = l.end();
       it != end; ++it) {
    localParticles[it->p.identity] = &(*it);
  }
}

void Storage::addParticle(longint id, const real p[3])
{
  Cell *cell;

  Particle n;
  n.init();
  n.p.identity = id;
  for (int i = 0; i < 3 ; ++i) {
    n.r.p[i] = p[i];
    n.l.i[i] = 0;
  }
  system->foldPosition(n.r.p, n.l.i);
  cell = mapPositionToCellClipped(n.r.p);

  appendIndexedParticle(cell->particles, n);

  LOG4ESPP_TRACE(logger, "got particle id="
		 << id << " @ " << p[0] << " " << p[1] << " " << p[2] << " ; put it into cell " << cell - getFirstCell());
  LOG4ESPP_TRACE(logger, "cell size is now " << cell->size());
}

Particle *Storage::appendUnindexedParticle(ParticleList &l, Particle &part)
{
  l.push_back(part);
  return &l.back();
}

Particle *Storage::appendIndexedParticle(ParticleList &l, Particle &part)
{
  // see whether the array was resized; STL hack
  Particle *begin = &l.front();

  l.push_back(part);
  Particle *p = &l.back();

  if (begin != &l.front())
    updateLocalParticles(l);
  else
    localParticles[p->p.identity] = p;
  return p;
}

Particle *Storage::moveUnindexedParticle(ParticleList &dl, ParticleList &sl, int i)
{
  dl.push_back(sl[i]);
  int newSize = sl.size() - 1;
  if (i != newSize) {
    sl[i] = sl.back();
  }
  sl.resize(newSize);
  return &dl.back();
}

Particle *Storage::moveIndexedParticle(ParticleList &dl, ParticleList &sl, int i)
{
  // see whether the arrays were resized; STL hack
  Particle *dbegin = &dl.front();
  Particle *sbegin = &sl.front();

  dl.push_back(sl[i]);
  int newSize = sl.size() - 1;
  if (i != newSize) {
    sl[i] = sl.back();
  }
  sl.resize(newSize);

  Particle *dst = &dl.back();
  Particle *src = &(sl[i]);

  // fix up destination list
  if (dbegin !=  &dl.front()) {
    updateLocalParticles(dl);
  }
  else {
    localParticles[dst->p.identity] = dst;
  }
  // fix up resorted source list; due to moving, the last particle
  // might have been moved to the position of the actually moved one
  if (sbegin != &sl.front()) {
    updateLocalParticles(sl);
  }
  else if (i != newSize) {
    localParticles[src->p.identity] = src;
  }

  return dst;
}

void Storage::fetchParticles(Storage &old)
{
  LOG4ESPP_DEBUG(logger, "number of received cells = "
		 << old.getActiveCells().size());

  for (CellListIterator it(old.getActiveCells());
       it.isValid(); ++it) {
    Particle &part = *it;
    Cell *nc = mapPositionToCellClipped(part.r.p);
    appendUnindexedParticle(nc->particles, part);
  }

  // update localParticles
  for(std::vector<Cell *>::iterator
	it = activeCells.begin(),
	end = activeCells.end();
      it != end; ++it) {
    updateLocalParticles((*it)->particles);
  }
}

void Storage::sendParticles(ParticleList &l, longint node)
{
  LOG4ESPP_DEBUG(logger, "send " << l.size() << " particles to " << node);

  longint nPart = l.size();
  comm.send(node, STORAGE_COMM_TAG, nPart);
  if (nPart > 0) {
    comm.send(node, STORAGE_COMM_TAG,
	      static_cast<char *>(static_cast<void *>(&(l[0]))), nPart*sizeof(Particle));
  }
  l.clear();

  LOG4ESPP_DEBUG(logger, "done");
}

void Storage::recvParticles(ParticleList &l, longint node)
{
  LOG4ESPP_DEBUG(logger, "recv from " << node);

  longint nPart;
  comm.recv(node, STORAGE_COMM_TAG, nPart);

  longint curSize = l.size();
  LOG4ESPP_DEBUG(logger, "will get " << nPart << " particles, have " << curSize);

  if (nPart > 0) {
    l.resize(curSize + nPart);

    comm.recv(node, STORAGE_COMM_TAG, static_cast<char *>(static_cast<void *>(&(l[curSize]))),
	      nPart*sizeof(Particle));
  }
  LOG4ESPP_DEBUG(logger, "done");
}
