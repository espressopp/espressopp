/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

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
#include "esutil/Error.hpp"

#include <iostream>
#include <boost/unordered/unordered_map.hpp>
using namespace std;

using namespace boost;
using namespace espressopp::iterator;

namespace espressopp {
 
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
      for (CellList::const_iterator it = realCells.begin(), end = realCells.end(); it != end; ++it) {
        longint size = (*it)->particles.size();
        if (size) {
          LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
        }
        cnt += size;
      }
      return cnt;
    }
    
    longint Storage::getNLocalParticles() const {
      longint cnt = 0;
      for (CellList::const_iterator it = localCells.begin(), end = localCells.end(); it != end; ++it) {
        longint size = (*it)->particles.size();
        if (size) {
          LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
        }
        cnt += size;
      }
      return cnt;
    }
    longint Storage::getNGhostParticles() const {
      longint cnt = 0;
      for (CellList::const_iterator it = ghostCells.begin(), end = ghostCells.end(); it != end; ++it) {
        longint size = (*it)->particles.size();
        if (size) {
          LOG4ESPP_TRACE(logger, "cell " << ((*it) - getFirstCell()) << " size " << size);
        }
        cnt += size;
      }
      return cnt;
    }
    longint Storage::getNAdressParticles() const {
          longint cnt = 0;
            return localAdrATParticles.size();
    }

    python::list Storage::getRealParticleIDs() {
      python::list pids;
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
    	pids.append(cit->getId());
      }
      return pids;
    }

    // TODO find out why python crashes if inlined
    //inline
    void Storage::removeFromLocalParticles(Particle *p, bool weak) {
      /* no pointer left, can happen for ghosts when the real particle
	 e has already been removed */
      if (localParticles.find(p->id()) == localParticles.end()){
        return;
      }

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

    /** Scale coordinates of all real particles by factor s */
    void Storage::scaleVolume(real s) {
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D pos = cit->getPos();
        pos *= s;
        cit->setPos(pos);
      }
    }
    /** Scale coordinates of all real particles by factor s. Anisotropic case */
    void Storage::scaleVolume(Real3D s) {
      for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
        Real3D pos = cit->getPos();
        pos[0] *= s[0];
        pos[1] *= s[1];
        pos[2] *= s[2];
        cit->setPos(pos);
      }
    }

    void Storage::removeAdrATParticle(longint id) {

    	if (localAdrATParticles.find(id) == localAdrATParticles.end()) {
    		std::cout << "not removing AT particle "<< id << ", since not found \n";
    		return;
    	}


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

    // TODO find out why python crashes if inlined
    //inline
    void Storage::updateInLocalParticles(Particle *p, bool weak) {
      if (!weak || localParticles.find(p->id()) == localParticles.end()) {
          LOG4ESPP_TRACE(logger, "updating local pointer for particle id="
		       << p->id() << " @ " << p);


          localParticles[p->id()] = p;

          /*
          // AdResS testing TODO
          Particle* oldp = localParticles[p->id()];
          Particle* newp = p;

          localParticles[p->id()] = p;

          //std::cout << "old *p " << oldp << "\n";
          //std::cout << "new *p " << newp << "\n";

          // TODO reorganize
          if (oldp && fixedtupleList) {
              FixedTupleList::iterator it;
              it = fixedtupleList->find(oldp);
              if (it != fixedtupleList->end()) {
                  std::vector<Particle*> tmp;
                  tmp = it->second;
                  fixedtupleList->insert(std::make_pair(newp, tmp));
                  fixedtupleList->erase(it);
              }
              else {
                  //std::cout << "updateInLocalParticles: Particle not found in tuples!\n";
              }
          }
          */
      }
      else {
          LOG4ESPP_TRACE(logger, "NOT updating local pointer for particle id="
		       << p->id() << " @ " << p << " has already pointer @ "
		       << localParticles[p->id()]);
      }
    }

    inline
    void Storage::updateInLocalAdrATParticles(Particle *p) {
          localAdrATParticles[p->id()] = p;
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
      for (LocalCellList::iterator it = cells.begin(), end = cells.end(); it != end; ++it) {
        localCells.push_back(&(*it));
      }
    }

    Particle* Storage::addParticle(longint id, const Real3D& p) {
      if (!checkIsRealParticle(id, p)) {
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

      //std::cout << "add particle: " << n.id() << " (" << n.position() << ")\n";

      appendIndexedParticle(cell->particles, n);

      LOG4ESPP_TRACE(logger, "got particle id ="
		     << id << " @ " << p << " ; put it into cell " << cell - getFirstCell());
      LOG4ESPP_TRACE(logger, "folded it to "
		     << n.r.p[0] << " " << n.r.p[1] << " " << n.r.p[2] );
      LOG4ESPP_TRACE(logger, "cell size is now " << cell->particles.size());

      return &cell->particles.back();
    }
    
    int Storage::removeParticle(longint id){
      Particle* p = lookupRealParticle(id);
      if(p){
        Cell *cell = mapPositionToCellChecked(p->position());
        
        removeFromLocalParticles( p );

        int rem_part_ind = p - &cell->particles[0];
    	int newSize = cell->particles.size() - 1;
    	if (rem_part_ind != newSize) {
          cell->particles[rem_part_ind] = cell->particles.back();
    	}
        cell->particles.resize(newSize);
        
        updateLocalParticles( cell->particles );

        onParticlesChanged();
        Particle* p1 = lookupRealParticle(id);
        if(p1){
          // it should not be printed out
          std::cout<< "Part still exst. pid: " << id << "  rpid: " << p1->id() << std::endl;
        }
        return 1;
      }
      else {
        return 0;
      }
      // TODO particle should be removed from different particle groups and lists too
      
    }
    
    void Storage::removeAllParticles(){
      localParticles.clear();
      for (CellList::iterator it = localCells.begin(), end = localCells.end(); it != end; ++it) {
        (*it)->particles.clear();
      }
      onParticlesChanged();
    }
    

    Particle* Storage::addAdrATParticle(longint id, const Real3D& p, const Real3D& _vpp) {

      if (!checkIsRealParticle(id, _vpp)) {
    	return static_cast< Particle* >(0);
      }

      Particle n;
      n.init();
      n.id() = id;
      n.position()= p;
      n.image() = Int3D(0);

      //std::cout << "add ATparticle: " << n.id() << " (" << n.position() << ")\n";

      // fold AT particles for same amount as VP
      Real3D vpp_old = _vpp;
      Real3D vpp_new = _vpp;

      getSystem()->bc->foldPosition(vpp_new);

      if (vpp_old != vpp_new) {
          //std::cout << "VP old pos (" << vpp_old << ") ";
          //std::cout << "new pos (" << vpp_new << ")\n";

          Real3D moved = vpp_old - vpp_new;
          n.position() = n.position() - moved;

          //std::cout << " Moved AT particle to " << n.position() << "\n";
      }


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

      return local;
    }

    // this is called from fixedtuplelist only!
    Particle* Storage::addAdrATParticleFTPL(Particle n) {

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

    //Particle *Storage::appendUnindexedAdrParticle(ParticleListAdr &l, Particle &part)
    Particle *Storage::appendUnindexedAdrParticle(ParticleList &l, Particle &part) {
      l.push_back(part);
      return &l.back();
    }

    void Storage::appendParticleListToGhosts(ParticleList &l) {
      AdrATParticlesG.push_back(l);
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

      for (CellListIterator it(old.getRealCells()); it.isValid(); ++it) {
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

      OutBuffer& data = outBuffer;

      data.reset();
      int size = list.size();
      data.write(size);
      for (ParticleList::Iterator it(list); it.isValid(); ++it) {
          removeFromLocalParticles(&(*it));
          data.write(*it);
      }

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

        afterRecvParticles(list, data); // this also takes care of AdResS AT Particles
      }

      LOG4ESPP_DEBUG(logger, "done");
    }

    void Storage::invalidateGhosts() {
      for(CellListIterator it(getGhostCells()); it.isValid(); ++it) {
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

        LOG4ESPP_TRACE(logger, "from particle " << src->id() << ": packing force " << src->force());
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
      }
    }
    

    void Storage::clearSavedPositions(){
      savedRealPositions.clear();
      savedImages.clear();
    }
    void Storage::savePosition(size_t id){
      if( getSystemRef().comm->size()==1 ){
        Particle* p = lookupRealParticle(id);
        if(p){
          savedRealPositions[id] = p->position();
          savedImages[id] = p->image();
        }
      }
      else{
        esutil::Error err(getSystem()->comm);
        stringstream msg;
        msg<<" At the moment it works only for one CPU. One can not store old positions"
                " for several CPUs";
        err.setException( msg.str() );
      }
    }
    
    void Storage::restorePositions(){
      int count = 0;
      int totCount = 0;
      if( !savedRealPositions.empty() ){
        for (map<size_t,Real3D>::iterator itr=savedRealPositions.begin(); itr != savedRealPositions.end(); ++itr) {
          size_t id = itr->first;
          Particle* p = lookupRealParticle(id);
          if(p){
            p->position() = itr->second;
          }
        }
        for (map<size_t,Int3D>::iterator itr=savedImages.begin(); itr != savedImages.end(); ++itr) {
          size_t id = itr->first;
          Particle* p = lookupRealParticle(id);
          if(p){
            p->image() = itr->second;
          }
        }
        count++;
      }
      
      mpi::all_reduce(*getSystem()->comm, count, totCount, std::plus<int>());
      
      if(totCount == 0){
        esutil::Error err(getSystem()->comm);
        stringstream msg;
        msg<<" There is nothing to restore. Check whether you saved positions";
        err.setException( msg.str() );
      }
    }
    

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    Storage::registerPython() {
      using namespace espressopp::python;
      class_< Storage, boost::noncopyable >("storage_Storage", no_init)
	    .def("clearSavedPositions", &Storage::clearSavedPositions)
	    .def("savePosition", &Storage::savePosition)
	    .def("restorePositions", &Storage::restorePositions)
	    .def("addParticle", &Storage::addParticle, return_value_policy< reference_existing_object >())
	    .def("removeParticle", &Storage::removeParticle)
	    .def("removeAllParticles", &Storage::removeAllParticles)
        .def("addAdrATParticle", &Storage::addAdrATParticle, return_value_policy< reference_existing_object >())
        .def("setFixedTuplesAdress", &Storage::setFixedTuplesAdress)
        //.def("addParticle", &Storage::addParticle, return_value_policy< reference_existing_object >())
	    .def("lookupLocalParticle", &Storage::lookupLocalParticle, return_value_policy< reference_existing_object >())
	    .def("lookupRealParticle", &Storage::lookupRealParticle, return_value_policy< reference_existing_object >())
	    .def("decompose", &Storage::decompose)
	    .def("getRealParticleIDs", &Storage::getRealParticleIDs)
        .add_property("system", &Storage::getSystem)
	    ;
    }
  }
}
