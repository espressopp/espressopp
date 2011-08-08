#include "python.hpp"

//#include <algorithm>

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
      : SystemAccess(system),
        inBuffer(*system->comm),
        outBuffer(*system->comm)
    {
      //logger.setLevel(log4espp::Logger::TRACE);
      LOG4ESPP_INFO(logger, "Created new storage object for a system, has buffers");
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


    void Storage::removeAdrATParticle(longint id) {

    	if (localAdrATParticles.find(id) == localAdrATParticles.end()) {
    		std::cout << "not removing AT particle "<< id << ", since not found \n";
    		return;
    	}

    	//std::cout << " removing AT particle "<< id << "\n";


		// remove from ParticleList (vector)
        Particle *dbegin = &AdrATParticles.front(); // see whether the array was moved
    	Particle* p = lookupAdrATParticle(id);
    	int i = p - &AdrATParticles[0];
    	int newSize = AdrATParticles.size() - 1;
    	if (i != newSize) { // if we are not removing the particle in last place
    		AdrATParticles[i] = AdrATParticles.back();
    	}
    	AdrATParticles.resize(newSize);


    	// remove from particle map
    	localAdrATParticles.erase(id);


    	// update particle map with particle list
    	if (dbegin != &AdrATParticles.front()) {
    		updateLocalParticles(AdrATParticles, true);
    	}
    	else if (i != newSize) {
    		Particle *np = &(AdrATParticles[i]);
    		updateInLocalAdrATParticles(np);
    	}
	}


    inline void Storage::updateInLocalParticles(Particle *p, bool weak) {
      if (!weak || localParticles.find(p->id()) == localParticles.end()) {
          LOG4ESPP_TRACE(logger, "updating local pointer for particle id="
		       << p->id() << " @ " << p);
          localParticles[p->id()] = p;

          /*std::cout << "updating local pointer for VP particle id "
                                << p->id() << " @ " << p << "\n";*/
      }
      else {
          LOG4ESPP_TRACE(logger, "NOT updating local pointer for particle id="
		       << p->id() << " @ " << p << " has already pointer @ "
		       << localParticles[p->id()]);
      }
    }

    inline void Storage::updateInLocalAdrATParticles(Particle *p) {
          localAdrATParticles[p->id()] = p;
          /*std::cout << " updating local pointer for AT particle id "
                  << p->id() << " @ " << p << "\n";*/
    }

    void Storage::updateLocalParticles(ParticleList &list, bool adress) {
      if (adress) {
          for (ParticleList::Iterator it(list); it.isValid(); ++it) {
              updateInLocalAdrATParticles(&(*it));
          }
      }
      else {
          for (ParticleList::Iterator it(list); it.isValid(); ++it) {
              updateInLocalParticles(&(*it));
          }
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

    Particle* Storage::addParticle(longint id, const Real3D& p) {
      if (!checkIsRealParticle(id, p)) {
    	//std::cout << getSystem()->comm->rank() << ": " << "VP particle: " << id << " does not belong here!\n";
        return static_cast< Particle* >(0);
      }

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

      //std::cout << getSystem()->comm->rank() << ": " << "add VP particle: " << id << "\n";

      return &cell->particles.back();
    }


    Particle* Storage::addAdrATParticle(longint id, const Real3D& p, const Real3D& vpp) {

      if (!checkIsRealParticle(id, vpp)) {
    	//std::cout << getSystem()->comm->rank() << ": " << "AT particle: " << id << " does not belong here!\n";
    	return static_cast< Particle* >(0);
      }

      Particle n;
      n.init();
      n.id() = id;
      n.position()= p;

      // see whether the array was resized; STL hack
      Particle *begin = &AdrATParticles.front();

      AdrATParticles.push_back(n);
      Particle* local = &AdrATParticles.back();

      if (begin != &AdrATParticles.front()) {
          updateLocalParticles(AdrATParticles, true);
      }
      else {
          updateInLocalAdrATParticles(local);
      }

      //std::cout << getSystem()->comm->rank() << ": " << "add AT particle: " << id << "\n";

      return local;
    }

    // this is called from fixedtuplelist only!
    Particle* Storage::addAdrATParticleFTPL(Particle n) {

	  //Particle n;
	  //n.init();
	  //n.id() = id;

	  // see whether the array was resized; STL hack
	  Particle *begin = &AdrATParticles.front();

	  AdrATParticles.push_back(n);
	  Particle* local = &AdrATParticles.back();

	  if (begin != &AdrATParticles.front()) {
		  updateLocalParticles(AdrATParticles, true);
	  }
	  else {
		  updateInLocalAdrATParticles(local);
	  }

	  //std::cout << getSystem()->comm->rank() << ": " << "add AT particle: " << id << "\n";

	  return local;

	}


    /*Particle* Storage::addParticle(longint id, const Real3D& p, int type) {
        Particle* pt = addParticle(id, p);
        pt->setType(type);
        return pt;
    }*/

    Particle *Storage::appendUnindexedParticle(ParticleList &l, Particle &part)
    {
      l.push_back(part);
      return &l.back();
    }

    Particle *Storage::appendUnindexedAdrParticle(ParticleListAdr &l, Particle &part)
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

      if (begin != &l.front()) {
          updateLocalParticles(l);
      }
      else {
          updateInLocalParticles(p);
      }

      return p;
    }

    Particle *Storage::moveIndexedParticle(ParticleList &dl, ParticleList &sl, int i)
    {

      //std::cout << getSystem()->comm->rank() << ": " << "move indexed particles: " << sl.size() << "\n";

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



      // for AdResS -- not used anymore, since I avoid it with globaltuples
      // TODO reorganize later
      // fix tuples
      /*
      std::cout << " moving src particle " << src->getId() << "-" << src->ghost() <<
    		  " from " << src << " to dst " << dst << "\n";

      FixedTupleList::iterator its;
	  its = fixedtupleList->find(src);
	  if (its != fixedtupleList->end()) {
		  std::cout << " src particle found in tuples \n";
	  } else {
		  std::cout << " src particle NOT found in tuples \n";
	  }*/

	  /*
	  std::vector<Particle*> tmp; // = it->second;

	  for (std::vector<Particle*>::iterator ita = its->second.begin(); ita != its->second.end(); ++ita) {
		  tmp.push_back(*ita); // copy vector
	  }

	  fixedtupleList->insert(std::make_pair(dst, tmp)); // copy vector to new position
	  fixedtupleList->erase(its); // erase old position
	  tmp.clear();
	  */

	  /*
	  itd = fixedtupleList->find(dst);
	  if (itd != fixedtupleList->end()) {
		  std::cout << " dst particle found in tuples \n";
	  } else {
		  std::cout << " dst particle NOT found in tuples \n";
	  }*/


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

      //std::cout << getSystem()->comm->rank() << ": " << " send " << list.size() << " particles to node " << node <<"\n";

      // pack for transport

      OutBuffer& data = outBuffer;

      data.reset();
      int size = list.size();
      data.write(size);
      for (ParticleList::Iterator it(list); it.isValid(); ++it) {
          removeFromLocalParticles(&(*it));
          data.write(*it);
      }

      //std::cout << getSystem()->comm->rank() << ": beforeSendParticles " << list.size() << "\n";
      beforeSendParticles(list, data); // this also takes care of AdResS AT Particles

      list.clear();

      // ... and send
      data.send(node, STORAGE_COMM_TAG);

      LOG4ESPP_DEBUG(logger, "done");
    }

    void Storage::recvParticles(ParticleList &list, longint node)
    {
      LOG4ESPP_DEBUG(logger, "recv from " << node);

      InBuffer& data = inBuffer;  // reuse storage buffer

      data.recv(node, STORAGE_COMM_TAG);

      // ... and unpack
      int size;
      //std::cout << getSystem()->comm->rank() << ": data.read(size) (recvParticles) \n";
      data.read(size);
      int curSize = list.size();
      LOG4ESPP_DEBUG(logger, "got " << size << " particles, have " << curSize);

      //std::cout << getSystem()->comm->rank() << ": " << " recv " << size << " particles from node " << node <<"\n";

      if (size > 0) {
        list.resize(curSize + size);

        for (int i = 0; i < size; ++i) {
          Particle *p = &list[curSize + i];
          //std::cout << getSystem()->comm->rank() << ": data.read(particle) (recvParticles) \n";
          data.read(*p);
          updateInLocalParticles(p);
        }

        //std::cout << getSystem()->comm->rank() << ": ";
        afterRecvParticles(list, data); // this also takes care of AdResS AT Particles

      }

      //std::cout << getSystem()->comm->rank() << ": ";
      //afterRecvParticles(list, data); // this also takes care of AdResS AT Particles

      LOG4ESPP_DEBUG(logger, "done");
    }

    void Storage::invalidateGhosts() {
      for(CellListIterator it(getGhostCells());
	  it.isValid(); ++it) {
	/* remove only ghosts from the hash if the localParticles hash
	   actually points to the ghost.  If there are local ghost cells
	   to implement pbc, the real particle will be the one accessible
	   via localParticles.
	*/
          removeFromLocalParticles(&(*it), true);
      }

      // TODO reorganize later
      // for AdResS
      //std::cout << getSystem()->comm->rank() << ": ----- CLEAR GHOST ADRESS PARTICLES (AdrATParticlesG) ----\n";
      AdrATParticlesG.clear();
      // clear ghost tuples
      FixedTupleList::iterator it = fixedtupleList->begin();
      for (;it != fixedtupleList->end(); ++it) {
    	  Particle* vp = it->first;
    	  if (vp->ghost()) {
    		  //std::cout << "erasing ghost particle in tuple: " << vp->id() << "-" << vp->ghost() << "\n";
    		  fixedtupleList->erase(it);
    	  }
      }

    }

    void Storage::decompose() {
      invalidateGhosts();
      decomposeRealParticles();
      //std::cout << getSystem()->comm->rank() << ": (onTuplesChanged) ";
      onTuplesChanged(); // for AdResS, renamed to not confuse with bonds
      exchangeGhosts();
      //std::cout << getSystem()->comm->rank() << ": ";
      onParticlesChanged();
    }

    void Storage::packPositionsEtc(OutBuffer &buf,
				   Cell &_reals, int extradata, const Real3D& shift) {
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


        /*std::cout << getSystem()->comm->rank() << ": -> " << src->id() << "-" << src->ghost() <<
        		" (" << src->getPos() << ") type " << src->getType() << " @ " << &(*src) << " \n";*/


        // TODO reorganize later
        // for AdResS
		// write all AT particles belonging to this VP into the buffer
        FixedTupleList::iterator it;
		it = fixedtupleList->find(&(*src));
		if (it != fixedtupleList->end()) {
			std::vector<Particle*> atList;
			atList = it->second;

			for (std::vector<Particle*>::iterator itv = atList.begin();
				  itv != atList.end(); ++itv) {
				Particle &at = **itv;

				/*std::cout << getSystem()->comm->rank() << ":  ---> " << at.id() << "-" << at.ghost() <<
						" (" << at.position() << ") type " << at.type() << " \n";*/

				if (at.type() == 0) { // TODO just for testing
					std::cout << "SERIOUS ERROR, particle is of wrong type\n";
					exit(1);
					return;
				}

				buf.write(at, extradata, shift);
			}
		}
		else {
			std::cout << getSystem()->comm->rank() << ": packposetc " << "VP particle "<< src->id() << "-" << src->ghost() << " not found in tuples!\n";
			exit(1);
			return;
		}

      }
    }

    void Storage::unpackPositionsEtc(Cell &_ghosts, InBuffer &buf, int extradata) {
      ParticleList &ghosts  = _ghosts.particles;

      LOG4ESPP_DEBUG(logger, "unpack data to ghosts in "
		     << (&_ghosts - getFirstCell()));
      LOG4ESPP_DEBUG(logger, "also unpacking "
		     << ((extradata & DATA_PROPERTIES) ? "properties " : "")
		     << ((extradata & DATA_MOMENTUM) ? "momentum " : "")
		     << ((extradata & DATA_LOCAL) ? "local " : ""));

      for(ParticleList::iterator dst = ghosts.begin(), end = ghosts.end(); dst != end; ++dst) {

    	//std::cout << getSystem()->comm->rank() << ": buf.read(particle, extradata) (unpackPosEtc) \n";
        buf.read(*dst, extradata);

        if (extradata & DATA_PROPERTIES) {
        	updateInLocalParticles(&(*dst), true);
        }

        dst->ghost() = 1;


        // TODO reorganize later
        // for AdResS

        /*std::cout << getSystem()->comm->rank() << ": <- " << dst->id() << "-" << dst->ghost() <<
        		" (" << dst->getPos() << ") type " << dst->getType() << " \n";*/

        std::vector<Particle*> tmp;
        Particle* atg; // atom, ghost
        Particle tmpatg; // temporary particle, to be inserted into adr. at. ghost part.


        for (int i = 1; i <= 3; ++i) { // TODO change i to actual value
        	// read the AT partcles
        	//std::cout << getSystem()->comm->rank() << ": buf.read(AT particle, extradata) (unpackPosEtc), i: " << i << "\n";
			buf.read(tmpatg, extradata);

			atg = appendUnindexedAdrParticle(AdrATParticlesG, tmpatg);
			atg->ghost() = 1;
			atg->setType(1); // TODO it's fixed now, for testing purposes
			atg->id() = dst->id()+i;

			tmp.push_back(atg);

			/*std::cout << getSystem()->comm->rank() << ":  <--- " << atg->id() << "-" << atg->ghost() <<
       	                		" (" << atg->getPos() << ") type " << atg->getType() << " \n";*/

			if (atg->getType() == 0) { // TODO it's fixed now, for testing purposes
				std::cout << "SERIOUS ERROR, particle is of wrong type\n";
				exit(-1);
				return;
			}

        }

        fixedtupleList->insert(std::make_pair(&(*dst), tmp));
        tmp.clear();

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

      //std::cout << "Copy reals to ghosts ... \n";

      for(ParticleList::iterator src = reals.begin(), end = reals.end(), dst = ghosts.begin();
              src != end; ++src, ++dst) {
        dst->copyAsGhost(*src, extradata, shift);

        // for AdResS
        copyGhostTuples(*src, *dst, extradata, shift); //TODO reorganize later
      }
    }

    void Storage::copyGhostTuples(Particle& src, Particle& dst, int extradata, const Real3D& shift) {

        /*std::cout << getSystem()->comm->rank() << ": " << src.id() << "-" << src.ghost() << " (" << src.position() << ") -> "
                        << dst.id() << "-" << dst.ghost() << " (" << dst.position() << ") @ " << &dst << " \n";*/

        // create ghosts of particles in tuples
        //std::cout << "VP particle is " << src->id() << " @ " << lookupRealParticle(src->id()) << "\n";
        //std::cout << "VP particle is " << src.id() << " @ " << &src << "\n";
        //std::cout << "Creating ghost tuples...\n";
        FixedTupleList::iterator it;
        it = fixedtupleList->find(&src);
        if (it != fixedtupleList->end()) {
            std::vector<Particle*> atList;
            atList = it->second;

            //std::cout << " size of AT vector: " << atList.size() << "\n";

            Particle* atg; // atom, ghost
            std::vector<Particle*> tmp; // temporary vector

            for (std::vector<Particle*>::iterator itv = atList.begin();
                  itv != atList.end(); ++itv) {
                Particle &at = **itv;
                //std::cout << " AT id: " << at.id() << "\n";

                Particle n; // temporary particle, to be inserted into adr. at. ghost part.
                n.id() = at.getId();
                n.type() = at.getType();


                // see whether the array was resized; STL hack
                //Particle *begin = &AdrATParticlesG.front();

                atg = appendUnindexedAdrParticle(AdrATParticlesG, n);
                atg->copyAsGhost(at, extradata, shift);
                tmp.push_back(atg);

                /*if (begin != &AdrATParticlesG.front())
                    std::cout << "\n  -- AdrATParticlesG array resized!! --\n\n";*/

                /*std::cout << " " << at.id() << "-" << at.ghost() << " (" << at.position() << ") -> "
                        <<  atg->id() << "-" << atg->ghost() << " (" << atg->getPos() << ")\n";*/

            }
            fixedtupleList->insert(std::make_pair(&dst, tmp));
            tmp.clear();
            //std::cout << "\n";
        }
        else {
        	std::cout << "copyGhostTuples: VP particle "<< src.id() << "-" << src.ghost() << " (" << src.position() << ")" << " not found in tuples!\n";
        	exit(1);
        	return;
        }

    }

    void Storage::packForces(OutBuffer &buf, Cell &_ghosts)
    {
      LOG4ESPP_DEBUG(logger, "pack ghost forces to buffer from cell "
		     << (&_ghosts - getFirstCell()));

      ParticleList &ghosts = _ghosts.particles;
  
      for(ParticleList::iterator src = ghosts.begin(), end = ghosts.end(); src != end; ++src) {

        buf.write(src->particleForce());

        LOG4ESPP_TRACE(logger, "from particle " << src->id() << ": packing force " << src->force());


        /*std::cout << getSystem()->comm->rank() << ": " << " from particle " << src->id()
        		<< " packing force " << src->force() << "\n";*/



        // TODO reorganize later
        // for AdResS
		// write all AT particle forces belonging to this VP into the buffer
        FixedTupleList::iterator it;
		it = fixedtupleList->find(&(*src));
		if (it != fixedtupleList->end()) {
			std::vector<Particle*> atList;
			atList = it->second;

			for (std::vector<Particle*>::iterator itv = atList.begin();
				  itv != atList.end(); ++itv) {
				Particle &at = **itv;

				/*std::cout << getSystem()->comm->rank() << ":  ---> " << at.id() << "-" << at.ghost() <<
						" (" << at.position() << ") type " << at.type() << " \n";*/

				buf.write(at.particleForce());
			}
		}
		else {
			std::cout << getSystem()->comm->rank() << ": packforces " << "VP particle "<< src->id() << "-" << src->ghost() << " not found in tuples!\n";
			exit(1);
			return;
		}




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
    	  LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": unpacking force " << f.force());
    	  dst->particleForce() = f;


    	  // no AdResS code, since this is not used anywhere

      }
    }

    void Storage::unpackAndAddForces(Cell &_reals, InBuffer &buf)
    {
      LOG4ESPP_DEBUG(logger, "add forces from buffer to cell "
		     << (&_reals - getFirstCell()));

      ParticleList &reals = _reals.particles;

      for(ParticleList::iterator dst = reals.begin(), end = reals.end(); dst != end; ++dst) {
    	  ParticleForce f;
    	  //std::cout << getSystem()->comm->rank() << ": buf.read(force) (unpackAndAddForces) \n";
    	  buf.read(f);
    	  LOG4ESPP_TRACE(logger, "for particle " << dst->id() << ": unpacking force "
		       << f.f() << " and adding to " << dst->force());
    	  dst->particleForce() += f;



    	  /*std::cout << getSystem()->comm->rank() << ": " << " for particle " << dst->id() << " unpacking force "
    		  << " and adding to " << dst->force() << "\n";*/


    	  // TODO reorganize later
		  // for AdResS

    	  // iterate through atomistic particles in fixedtuplelist
		  FixedTupleList::iterator it;
		  it = fixedtupleList->find(&(*dst));

		  if (it != fixedtupleList->end()) {

			 std::vector<Particle*> atList1;
			 atList1 = it->second;

			 //std::cout << "AT forces ...\n";
			 for (std::vector<Particle*>::iterator itv = atList1.begin(); itv != atList1.end(); ++itv) {
				 Particle &p3 = **itv;
				 //std::cout << getSystem()->comm->rank() << ": buf.read(AT force) (unpackAndAddForces) \n";
				 buf.read(f);
				 p3.particleForce() += f;
			 }
		  }
		  else {
			 std::cout << " unpackForces: one of the VP particles not found in tuples: " << dst->id() << "-" << dst->ghost();
			 exit(1);
			 return;
		  }

      }
    }

    void Storage::addGhostForcesToReals(Cell &_ghosts, Cell &_reals)
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

          // for AdResS
          addAdrGhostForcesToReals(*src, *dst); // TODO reorganize later
      }
    }


    void Storage::addAdrGhostForcesToReals(Particle& src, Particle& dst) {
        /*std::cout << "adding force from " << src.id() << "-" << src.ghost() <<
              " to " << dst.id() << "-" << dst.ghost() << "\n";*/

        // iterate through atomistic particles in fixedtuplelist
        FixedTupleList::iterator it3;
        FixedTupleList::iterator it4;
        it3 = fixedtupleList->find(&src);
        it4 = fixedtupleList->find(&dst);

        //std::cout << "\nInteraction " << p1.id() << " - " << p2.id() << "\n";
        if (it3 != fixedtupleList->end() && it4 != fixedtupleList->end()) {

           std::vector<Particle*> atList1;
           std::vector<Particle*> atList2;
           atList1 = it3->second;
           atList2 = it4->second;

           //std::cout << "AT forces ...\n";
           for (std::vector<Particle*>::iterator itv = atList1.begin(),
                   itv2 = atList2.begin(); itv != atList1.end(); ++itv, ++itv2) {

               Particle &p3 = **itv;
               Particle &p4 = **itv2;

               /*std::cout << " from " << p3.id() << "-" << p3.ghost() << " to " <<
                       p4.id() << "-" << p4.ghost() << "\n";*/

               p4.particleForce() += p3.particleForce();
           }
        }
        else {
           std::cout << " one of the VP particles not found in tuples: " << src.id() << "-" <<
                   src.ghost() << ", " << dst.id() << "-" << dst.ghost();
           exit(1);
           return;
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

    .def("addAdrATParticle", &Storage::addAdrATParticle,
         return_value_policy< reference_existing_object >())

    .def("setFixedTuples", &Storage::setFixedTuples)


  //
  //.def("addParticle", &Storage::addParticle,
  //     return_value_policy< reference_existing_object >())

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
