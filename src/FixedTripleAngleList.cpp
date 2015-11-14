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
#include <sstream>
#include "FixedTripleAngleList.hpp"
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"

#include <cmath>

namespace espressopp {

  LOG4ESPP_LOGGER(FixedTripleAngleList::theLogger, "FixedTripleAngleList");

  FixedTripleAngleList::FixedTripleAngleList(shared_ptr< storage::Storage > _storage)
  : storage(_storage), triplesAngles()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedTripleAngleList");
    
    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedTripleAngleList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedTripleAngleList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedTripleAngleList::onParticlesChanged, this));
  }

  FixedTripleAngleList::~FixedTripleAngleList() {

    LOG4ESPP_INFO(theLogger, "~FixedTripleAngleList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedTripleAngleList::
  add(longint pid1, longint pid2, longint pid3) {
    bool returnVal = true;
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);

    // ADD THE LOCAL TRIPLET
    Particle *p1 = storage->lookupLocalParticle(pid1);
    Particle *p2 = storage->lookupRealParticle(pid2);
    Particle *p3 = storage->lookupLocalParticle(pid3);

    // middle particle is the reference particle and must exist here
    if (!p2){
      // particle does not exists here (some other CPU must have it)
      returnVal = false;
    }
    else{
      if (!p1) {
        std::stringstream msg;
        msg << "triple particle p1 " << pid1 << " does not exists here and cannot be added";
        msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
        err.setException( msg.str() );
      }

      if (!p3) {
        std::stringstream msg;
        msg << "triple particle p3 " << pid3 << " does not exists here and cannot be added";
        msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
        err.setException( msg.str() );
      }
    }
    err.checkException();
    
    if(returnVal){
      // add the triple locally
      this->add(p1, p2, p3);
      //printf("me = %d: pid1 %d, pid2 %d, pid3 %d\n", mpiWorld->rank(), pid1, pid2, pid3);

      Real3D pos1 = p1->position();
      Real3D pos2 = p2->position();
      Real3D pos3 = p3->position();

      Real3D r12(0.0,0.0,0.0);
      system.bc->getMinimumImageVectorBox(r12, pos2, pos1);
      Real3D r32(0.0,0.0,0.0);
      system.bc->getMinimumImageVectorBox(r32, pos2, pos3);

      real cosVal = r12*r32 / (r12.abs() * r32.abs());

      if(cosVal>1) cosVal=1;
      if(cosVal<-1) cosVal=-1;

      triplesAngles.insert(std::make_pair(pid2, 
              std::make_pair(std::make_pair(pid1, pid3), acos(cosVal) )));
    }
    LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
    return returnVal;
  }

  python::list FixedTripleAngleList::getTriples(){
	python::tuple triple;
	python::list triples;
	for (TriplesAngles::const_iterator it=triplesAngles.begin(); it!=triplesAngles.end(); it++) {
      triple = python::make_tuple(it->first, it->second.first.first, 
              it->second.first.second);
      triples.append(triple);
    }

	return triples;
  }
  
  python::list FixedTripleAngleList::getTriplesAngles(){
	python::tuple triple;
	python::list triples;
	for (TriplesAngles::const_iterator it=triplesAngles.begin(); it!=triplesAngles.end(); it++) {
      triple = python::make_tuple(it->first, it->second.first.first,
              it->second.first.second, it->second.second);
      triples.append(triple);
    }

	return triples;
  }
  
  real FixedTripleAngleList::getAngle(int pid1, int pid2, int pid3){
    real returnVal = -3;
    
    TriplesAngles::iterator itr;
	TriplesAngles::iterator lastElement;
	
	// locate an iterator to the first pair object associated with key
	itr = triplesAngles.find(pid2);
	if (itr == triplesAngles.end())
		return returnVal; // no elements associated with key, so return immediately

	// get an iterator to the element that is one past the last element associated with key
	lastElement = triplesAngles.upper_bound(pid2);
    
    std::pair<longint, longint> neededPair = std::make_pair(pid1, pid3);

	for ( ; itr != lastElement; ++itr){
      if(neededPair==itr->second.first){
        returnVal = itr->second.second;
        break;
      }
    }
    
	return returnVal;
  }
  
  void FixedTripleAngleList::beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
    std::vector< longint > toSendInt;
    std::vector< real > toSendReal;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      // find all triples that involve this particle
      int n = triplesAngles.count(pid);

      if (n > 0) {
        std::pair<TriplesAngles::const_iterator,
        TriplesAngles::const_iterator> equalRange = triplesAngles.equal_range(pid);

        // first write the pid of this particle
        // then the number of partners (n)
        // and then the pids of the partners
        toSendInt.reserve(toSendInt.size()+2*n+1);
        toSendReal.reserve(toSendReal.size()+n);
        toSendInt.push_back(pid);
        toSendInt.push_back(n);
        for(TriplesAngles::const_iterator it = equalRange.first; it!=equalRange.second; ++it) {
          toSendInt.push_back(it->second.first.first);
          toSendInt.push_back(it->second.first.second);
          toSendReal.push_back(it->second.second);
        }

        // delete all of these triples from the global list
        triplesAngles.erase(pid);
      }
    }
    // send the list
    buf.write(toSendInt);
    buf.write(toSendReal);
    LOG4ESPP_INFO(theLogger, "prepared fixed triple list before send particles");
  }

  void FixedTripleAngleList::afterRecvParticles(ParticleList &pl, InBuffer& buf) {

    std::vector< longint > receivedInt;
    std::vector< real > receivedReal;
    int n;
    longint pid1, pid2, pid3;
    real angleVal;
    TriplesAngles::iterator it = triplesAngles.begin();
    // receive the triple list
    buf.read(receivedInt);
    buf.read(receivedReal);
    int size = receivedInt.size(); int i = 0;  int j = 0;
    while (i < size) {
      // unpack the list
      pid2 = receivedInt[i++];
      n = receivedInt[i++];
      for (; n > 0; --n) {
        pid1 = receivedInt[i++];
        pid3 = receivedInt[i++];
        angleVal = receivedReal[j++];
        
        // add the triple
        it = triplesAngles.insert(it, std::make_pair(pid2, 
                std::make_pair(std::make_pair(pid1, pid3), angleVal)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed triple list after receive particles");
  }

  void FixedTripleAngleList::onParticlesChanged() {
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    // (re-)generate the local triple list from the global list
    //printf("FixedTripleAngleList: rebuild local triple list from global\n");
    this->clear();
    longint lastpid2 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    for (TriplesAngles::const_iterator it = triplesAngles.begin(); it!=triplesAngles.end(); ++it) {
      //printf("lookup global triple %d %d %d\n", it->first, it->second.first, it->second.second);
      if (it->first != lastpid2) {
        p2 = storage->lookupRealParticle(it->first);
        if (p2 == NULL) {
          std::stringstream msg;
          msg << "triple particle p2 " << it->first << " does not exists here";
          err.setException( msg.str() );
        }
	    lastpid2 = it->first;
      }
      p1 = storage->lookupLocalParticle(it->second.first.first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "triple particle p1 " << it->second.first.first << " does not exists here";
        err.setException( msg.str() );
      }
      p3 = storage->lookupLocalParticle(it->second.first.second);
      if (p3 == NULL) {
        std::stringstream msg;
        msg<< "triple particle p3 " << it->second.first.second << " does not exists here";
        err.setException( msg.str() );
      }
      this->add(p1, p2, p3);
    }
    err.checkException();
    LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTripleAngleList::registerPython() {

    using namespace espressopp::python;

    bool (FixedTripleAngleList::*pyAdd)(longint pid1, longint pid2, longint pid3)
      = &FixedTripleAngleList::add;
    //bool (FixedTripleAngleList::*pyAdd)(pvec pids)
    //      = &FixedTripleAngleList::add;

    class_< FixedTripleAngleList, shared_ptr< FixedTripleAngleList > >
      ("FixedTripleAngleList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      .def("size", &FixedTripleAngleList::size)
      .def("getTriples",  &FixedTripleAngleList::getTriples)
      .def("getTriplesAngles",  &FixedTripleAngleList::getTriplesAngles)
      .def("getAngle",  &FixedTripleAngleList::getAngle)
     ;
  }
}
