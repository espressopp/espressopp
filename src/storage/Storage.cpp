#include "python.hpp"

#include <algorithm>

#include "log4espp.hpp"

#include "System.hpp"
#include "Storage.hpp"
#include "bc/BC.hpp"
#include "iterator/CellListIterator.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "Particle.hpp"
#include "Buffer.hpp"

using namespace boost;
using namespace espresso::iterator;

namespace espresso {
 
  namespace storage {
    LOG4ESPP_LOGGER(Storage::logger, "Storage");

    const int STORAGE_COMM_TAG = 0xaa;

    const int Storage::dataOfUpdateGhosts = 0;
    const int Storage::dataOfExchangeGhosts = DATA_PROPERTIES;

    Storage::Storage(shared_ptr< System > system)
      : SystemAccess(system)
    {
      //logger.setLevel(log4espp::Logger::TRACE);
    }

    Storage::~Storage() {}

    longint Storage::getNRealParticles() const {
      longint cnt = 0;
      for (CellList::const_iterator 
	     it = realCells.begin(),
	     end = realCells.end();
	   it != end; ++it) {
	longint size = (*it)->particles.size();
	if (size) {
	  LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
	}
	cnt += size;
      }
      return cnt;
    }

    inline void Storage::removeFromLocalParticles(Particle *p, bool weak) {
      /* no pointer left, can happen for ghosts when the real particle
	 e has already been removed */
      if (localParticles.find(p->id()) == localParticles.end())
	return;

      if (!weak || localParticles[p->id()] == p) {
	LOG4ESPP_TRACE(logger, "removing local pointer for particle id="
		       << p->id() << " @ " << p);
	localParticles.erase(p->id());
      }
      else {
	LOG4ESPP_TRACE(logger, "NOT removing local pointer for particle id="
		       << p->id() << " @ " << p << " since pointer is @ "
		       << localParticles[p->id()]);
      }
    }

    inline void Storage::updateInLocalParticles(Particle *p, bool weak) {
      if (!weak || localParticles.find(p->id()) == localParticles.end()) {
	LOG4ESPP_TRACE(logger, "updating local pointer for particle id="
		       << p->id() << " @ " << p);
	localParticles[p->id()] = p;
      }
      else {
	LOG4ESPP_TRACE(logger, "NOT updating local pointer for particle id="
		       << p->id() << " @ " << p << " has already pointer @ "
		       << localParticles[p->id()]);
      }
    }

    void Storage::updateLocalParticles(ParticleList &list) {
      for (ParticleList::Iterator it(list); it.isValid(); ++it) {
	updateInLocalParticles(&(*it));
      }
    }

    void Storage::resizeCells(longint nCells) {
      cells.resize(nCells);
      localCells.reserve(nCells);
      for (LocalCellList::iterator
	     it = cells.begin(),
	     end = cells.end(); it != end; ++it) {
	localCells.push_back(&(*it));
      }
    }

    Particle* Storage::
    addParticle(longint id, const Real3D& p) {
      if (!checkIsRealParticle(id, p))
	return static_cast< Particle* >(0);

      Cell *cell;

      Particle n;
      n.init();
      n.id() = id;
      n.position()= p;
      n.image() = Int3D(0);
      getSystem()->bc->foldPosition(n.position(), n.image());
      cell = mapPositionToCellClipped(n.position());

      appendIndexedParticle(cell->particles, n);

      LOG4ESPP_TRACE(logger, "got particle id ="
		     << id << " @ " << p << " ; put it into cell " << cell - getFirstCell());
      LOG4ESPP_TRACE(logger, "folded it to "
		     << n.r.p[0] << " " << n.r.p[1] << " " << n.r.p[2] );
      LOG4ESPP_TRACE(logger, "cell size is now " << cell->particles.size());

      return &cell->particles.back();
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
	updateInLocalParticles(p);
      return p;
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
      if (dbegin != &dl.front()) {
	updateLocalParticles(dl);
      }
      else {
	updateInLocalParticles(dst);
      }
      // fix up resorted source list; due to moving, the last particle
      // might have been moved to the position of the actually moved one
      if (sbegin != &sl.front()) {
	updateLocalParticles(sl);
      }
      else if (i != newSize) {
	updateInLocalParticles(src);
      }

      return dst;
    }

    void Storage::fetchParticles(Storage &old)
    {
      LOG4ESPP_DEBUG(logger, "number of received cells = "
		     << old.getRealCells().size());

      for (CellListIterator it(old.getRealCells());
	   it.isValid(); ++it) {
	Particle &part = *it;
	Cell *nc = mapPositionToCellClipped(part.position());
	appendUnindexedParticle(nc->particles, part);
      }

      // update localParticles
      for(CellList::Iterator it(realCells); it.isValid(); ++it) {
	updateLocalParticles((*it)->particles);
      }
    }

    void Storage::sendParticles(ParticleList &list, longint node)
    {
      LOG4ESPP_DEBUG(logger, "send " << list.size() << " particles to " << node);

      // pack for transport
      OutBuffer data(*getSystem()->comm);
      int size = list.size();
      data.write(size);
      for (ParticleList::Iterator it(list); it.isValid(); ++it) {
	removeFromLocalParticles(&(*it));
	data.write(*it);
      }

      beforeSendParticles(list, data);

      list.clear();

      // ... and send
      data.send(node, STORAGE_COMM_TAG);

      LOG4ESPP_DEBUG(logger, "done");
    }

    void Storage::recvParticles(ParticleList &list, longint node)
    {
      LOG4ESPP_DEBUG(logger, "recv from " << node);

      // receive packed data

      InBuffer data(*getSystem()->comm);

      data.recv(node, STORAGE_COMM_TAG);

      // ... and unpack
      int size;
      data.read(size);
      int curSize = list.size();
      LOG4ESPP_DEBUG(logger, "got " << size << " particles, have " << curSize);

      if (size > 0) {
	list.resize(curSize + size);

	for (int i = 0; i < size; ++i) {
	  Particle *p = &list[curSize + i];
	  data.read(*p);
	  updateInLocalParticles(p);
	}
      }
      afterRecvParticles(list, data);

      LOG4ESPP_DEBUG(logger, "done");
    }

    void Storage::invalidateGhosts()
    {
      for(CellListIterator it(getGhostCells());
	  it.isValid(); ++it) {
	/* remove only ghosts from the hash if the localParticles hash
	   actually points to the ghost.  If there are local ghost cells
	   to implement pbc, the real particle will be the one accessible
	   via localParticles.
	*/
	removeFromLocalParticles(&(*it), true);
      }
    }

    void Storage::decompose() {
      invalidateGhosts();
      decomposeRealParticles();
      exchangeGhosts();
      onParticlesChanged();
    }

    void Storage::packPositionsEtc(OutBuffer &buf,
				   Cell &_reals, int extradata, const Real3D& shift)
    {
      ParticleList &reals  = _reals.particles;

      LOG4ESPP_DEBUG(logger, "pack data from reals in "
		     << (&_reals - getFirstCell()));
      LOG4ESPP_DEBUG(logger, "also packing "
		     << ((extradata & DATA_PROPERTIES) ? "properties " : "")
		     << ((extradata & DATA_MOMENTUM) ? "momentum " : "")
		     << ((extradata & DATA_LOCAL) ? "local " : ""));
      LOG4ESPP_DEBUG(logger, "positions are shifted by "
		     << shift[0] << "," << shift[1] << "," << shift[2]);

      for(ParticleList::iterator src = reals.begin(), end = reals.end(); src != end; ++src) {

        buf.write(*src, extradata, shift);
      }
    }

    void Storage::unpackPositionsEtc(Cell &_ghosts, InBuffer &buf, int extradata)
    {
      ParticleList &ghosts  = _ghosts.particles;

      LOG4ESPP_DEBUG(logger, "unpack data to ghosts in "
		     << (&_ghosts - getFirstCell()));
      LOG4ESPP_DEBUG(logger, "also unpacking "
		     << ((extradata & DATA_PROPERTIES) ? "properties " : "")
		     << ((extradata & DATA_MOMENTUM) ? "momentum " : "")
		     << ((extradata & DATA_LOCAL) ? "local " : ""));

      for(ParticleList::iterator dst = ghosts.begin(), end = ghosts.end(); dst != end; ++dst) {
        buf.read(*dst, extradata);
	if (extradata & DATA_PROPERTIES) {
	  updateInLocalParticles(&(*dst), true);
	}
	dst->ghost() = 1;
      }
    }

    void Storage::copyRealsToGhosts(Cell &_reals, Cell &_ghosts,
				    int extradata,
				    const Real3D& shift)
    {
      ParticleList &reals  = _reals.particles;
      ParticleList &ghosts = _ghosts.particles;

      LOG4ESPP_DEBUG(logger, "copy data from reals in "
		     << (&_reals - getFirstCell()) << " to ghosts in "
		     << (&_ghosts - getFirstCell()));
      LOG4ESPP_DEBUG(logger, "also copying "
		     << ((extradata & DATA_PROPERTIES) ? "properties " : "")
		     << ((extradata & DATA_MOMENTUM) ? "momentum " : "")
		     << ((extradata & DATA_LOCAL) ? "local " : ""));
      LOG4ESPP_DEBUG(logger, "positions are shifted by "
		     << shift[0] << "," << shift[1] << "," << shift[2]);

      ghosts.resize(reals.size());

      for(ParticleList::iterator src = reals.begin(), end = reals.end(), dst = ghosts.begin();
	  src != end; ++src, ++dst) {

        dst->copyAsGhost(*src, extradata, shift);
      }
    }

    void Storage::packForces(OutBuffer &buf, Cell &_ghosts)
    {
      LOG4ESPP_DEBUG(logger, "pack ghost forces to buffer from cell "
		     << (&_ghosts - getFirstCell()));

      ParticleList &ghosts = _ghosts.particles;
  
      for(ParticleList::iterator src = ghosts.begin(), end = ghosts.end(); src != end; ++src) {
        buf.write(src->particleForce());
	LOG4ESPP_TRACE(logger, "from particle " << src->id() 
                               << ": packing force " << src->force());
      }
    }

    void Storage::unpackForces(Cell &_reals, InBuffer &buf)
    {
      LOG4ESPP_DEBUG(logger, "add forces from buffer to cell "
		     << (&_reals - getFirstCell()));

      ParticleList &reals = _reals.particles;

      for(ParticleList::iterator dst = reals.begin(), end = reals.end(); dst != end; ++dst) {
	ParticleForce f;
        buf.read(f);
	LOG4ESPP_TRACE(logger, "for particle " << dst->id() 
                        << ": unpacking force " << f.force());
	dst->particleForce() = f;
      }
    }

    void Storage::unpackAndAddForces(Cell &_reals, InBuffer &buf)
    {
      LOG4ESPP_DEBUG(logger, "add forces from buffer to cell "
		     << (&_reals - getFirstCell()));

      ParticleList &reals = _reals.particles;

      for(ParticleList::iterator dst = reals.begin(), end = reals.end(); dst != end; ++dst) {
	ParticleForce f;
        buf.read(f);
	LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": unpacking force "
		       << f.f() << " and adding to " << dst->force());
	dst->particleForce() += f;
      }
    }

    void Storage::
    addGhostForcesToReals(Cell &_ghosts, Cell &_reals)
    {
      LOG4ESPP_DEBUG(logger, "add forces from ghosts in cell "
		     << (&_ghosts - getFirstCell()) << " to reals in cell "
		     << (&_reals - getFirstCell()));

      ParticleList &reals  = _reals.particles;
      ParticleList &ghosts = _ghosts.particles;

      for(ParticleList::iterator dst = reals.begin(), end = reals.end(), src = ghosts.begin();
	  dst != end; ++dst, ++src) {

	LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": adding force "
		       << src->force() << " to " << dst->force()); 
	dst->particleForce() += src->particleForce();
      }
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    Storage::registerPython() {
      using namespace espresso::python;
      class_< Storage, boost::noncopyable >("storage_Storage", no_init)
	.def("addParticle", &Storage::addParticle, 
	     return_value_policy< reference_existing_object >())
	.def("lookupLocalParticle", &Storage::lookupLocalParticle,
	     return_value_policy< reference_existing_object >())
	.def("lookupRealParticle", &Storage::lookupRealParticle,
	     return_value_policy< reference_existing_object >())
	.def("decompose", &Storage::decompose)
        .add_property("system", &Storage::getSystem)
	;

    }

  }
}
