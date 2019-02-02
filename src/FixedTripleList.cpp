/*
  Copyright (C) 2012,2013,2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2017
      Jakub Krajniak (jkrajniak at gmail.com)
  
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
#include "FixedTripleList.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "esutil/Error.hpp"


namespace espressopp {

  /*
  FixedTripleList::FixedTripleList(shared_ptr< storage::Storage > _storage)
  : FixedListComm (_storage){}
  */


  LOG4ESPP_LOGGER(FixedTripleList::theLogger, "FixedTripleList");

  FixedTripleList::FixedTripleList(shared_ptr< storage::Storage > _storage)
    : storage(_storage), globalTriples()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedTripleList");

    sigBeforeSend = storage->beforeSendParticles.connect
      (boost::bind(&FixedTripleList::beforeSendParticles, this, _1, _2));
    sigAfterRecv = storage->afterRecvParticles.connect
      (boost::bind(&FixedTripleList::afterRecvParticles, this, _1, _2));
    sigOnParticleChanged = storage->onParticlesChanged.connect
      (boost::bind(&FixedTripleList::onParticlesChanged, this));
  }

  FixedTripleList::~FixedTripleList() {

    LOG4ESPP_INFO(theLogger, "~FixedTripleList");

    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
    sigOnParticleChanged.disconnect();
  }

  /*
  bool FixedTripleList::add(longint pid1, longint pid2, longint pid3) {
      std::vector<longint> tmp;
      tmp.push_back(pid1);
      tmp.push_back(pid3);
      tmp.push_back(pid2); // this is used as key

      return FixedListComm::add(tmp);
  }*/


  bool FixedTripleList::
  add(longint pid1, longint pid2, longint pid3) {
    // three swaps needed for (1, 2, 3) == (1, 3, 2)
    // this is important for the line:
    // if (it->second == std::pair<longint, longint>(pid2, pid3))

    //if (pid3 < pid2)
    //  std::swap(pid2, pid3);
    //if (pid2 < pid1)
    //  std::swap(pid1, pid2);
    //if (pid3 < pid2)
    //  std::swap(pid2, pid3);
    
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
          msg << "adding error: triple particle p1 " << pid1 << " does not exists here and cannot be added";
          msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
          //std::runtime_error(err.str());
          err.setException( msg.str() );
      }
      if (!p3) {
          std::stringstream msg;
          msg << "adding error: triple particle p3 " << pid3 << " does not exists here and cannot be added";
          msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
          //std::runtime_error(err.str());
          err.setException( msg.str() );
      }
    }

    err.checkException();
    
    if(returnVal){
      // add the triple locally
      this->add(p1, p2, p3);
      //printf("me = %d: pid1 %d, pid2 %d, pid3 %d\n", mpiWorld->rank(), pid1, pid2, pid3);

      // ADD THE GLOBAL TRIPLET
      // see whether the particle already has triples
      std::pair<GlobalTriples::const_iterator,
                GlobalTriples::const_iterator> equalRange
        = globalTriples.equal_range(pid2);
      if (equalRange.first == globalTriples.end()) {
        // if it hasn't, insert the new triple
        globalTriples.insert(std::make_pair(pid2,
                             std::pair<longint, longint>(pid1, pid3)));
      }
      else {
        // otherwise test whether the triple already exists
        for (GlobalTriples::const_iterator it = equalRange.first; it != equalRange.second; ++it) {
          if (it->second == std::pair<longint, longint>(pid1, pid3)) {
            // TODO: Triple already exists, generate error!
  	      ;
          }
        }
        // if not, insert the new triple
        globalTriples.insert(equalRange.first, std::make_pair(pid2, std::pair<longint, longint>(pid1, pid3)));
      }
      LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
    }
    return returnVal;
  }

  python::list FixedTripleList::getTriples()
  {
	python::tuple triple;
	python::list triples;
	for (GlobalTriples::const_iterator it=globalTriples.begin(); it != globalTriples.end(); it++) {
      triple = python::make_tuple( it->second.first,it->first, it->second.second);
      triples.append(triple);
    }

	return triples;
  }

  std::vector<longint> FixedTripleList::getTripleList() {
    std::vector<longint> ret;
    for (GlobalTriples::const_iterator it=globalTriples.begin(); it != globalTriples.end(); it++) {
      ret.push_back(it->second.first);
      ret.push_back(it->first);
      ret.push_back(it->second.second);
    }
    return ret;
  }

  void FixedTripleList::
  beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
    std::vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      // LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find triples");
      //printf ("me = %d: send particle with pid %d find triples\n", mpiWorld->rank(), pid);

      // find all triples that involve this particle
      int n = globalTriples.count(pid);
      //printf ("me = %d: send particle with pid %d, has %d global triples\n", 
                //mpiWorld->rank(), pid, n);

      if (n > 0) {
          std::pair<GlobalTriples::const_iterator,
          GlobalTriples::const_iterator> equalRange
              = globalTriples.equal_range(pid);

          // first write the pid of this particle
          // then the number of partners (n)
          // and then the pids of the partners
          toSend.reserve(toSend.size()+2*n+1);
          toSend.push_back(pid);
          toSend.push_back(n);
          for (GlobalTriples::const_iterator it = equalRange.first;
                  it != equalRange.second; ++it) {
              toSend.push_back(it->second.first);
              toSend.push_back(it->second.second);
                  //printf ("send global triple: pid %d and partner %d\n", pid, it->second.first);
                  //printf ("send global triple: pid %d and partner %d\n", pid, it->second.second);
          }

          // delete all of these triples from the global list
          globalTriples.erase(equalRange.first, equalRange.second);
      }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed triple list before send particles");
  }

  void FixedTripleList::
  afterRecvParticles(ParticleList &pl, InBuffer& buf) {

    std::vector< longint > received;
    int n;
    longint pid1, pid2, pid3;
    GlobalTriples::iterator it = globalTriples.begin();
    // receive the triple list
    buf.read(received);
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid2 = received[i++];
      n = received[i++];
      //printf ("me = %d: recv particle with pid %d, has %d global triples\n",
                //mpiWorld->rank(), pid1, n);
      for (; n > 0; --n) {
	    pid1 = received[i++];
	    pid3 = received[i++];
	    // add the triple to the global list
        //printf("received triple %d %d %d, add triple to global list\n", pid1, pid2, pid3);
	    it = globalTriples.insert(it, std::make_pair(pid2,std::pair<longint, longint>(pid1, pid3)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed triple list after receive particles");
  }

  void FixedTripleList::onParticlesChanged() {
    
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    // (re-)generate the local triple list from the global list
    //printf("FixedTripleList: rebuild local triple list from global\n");
    this->clear();
    longint lastpid2 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    for (GlobalTriples::const_iterator it = globalTriples.begin();
	it != globalTriples.end(); ++it) {
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
      p1 = storage->lookupLocalParticle(it->second.first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "triple particle p1 " << it->second.first << " does not exists here";
        err.setException( msg.str() );
      }
      p3 = storage->lookupLocalParticle(it->second.second);
      if (p3 == NULL) {
        std::stringstream msg;
        msg << "triple particle p3 " << it->second.second << " does not exists here";
        err.setException( msg.str() );
      }
      this->add(p1, p2, p3);
    }
    err.checkException();
    
    LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
  }

  void FixedTripleList::remove() {
      this->clear();
      globalTriples.clear();
      sigBeforeSend.disconnect();
      sigAfterRecv.disconnect();
      sigOnParticleChanged.disconnect();
  }
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTripleList::registerPython() {

    using namespace espressopp::python;

    bool (FixedTripleList::*pyAdd)(longint pid1, longint pid2, longint pid3)
      = &FixedTripleList::add;
    //bool (FixedTripleList::*pyAdd)(pvec pids)
    //      = &FixedTripleList::add;

    class_< FixedTripleList, shared_ptr< FixedTripleList > >
      ("FixedTripleList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      .def("size", &FixedTripleList::size)
      .def("remove",  &FixedTripleList::remove)
      .def("getTriples",  &FixedTripleList::getTriples)
     ;
  }
}
