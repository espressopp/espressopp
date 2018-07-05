/*
  Copyright (C) 2016,2017
      Jakub Krajniak (jkrajniak at gmail.com)
  Copyright (C) 2012,2013,2016
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
#include "FixedQuadrupleList.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "esutil/Error.hpp"

//using namespace std;

namespace espressopp {

  /*
  FixedQuadrupleList::FixedQuadrupleList(shared_ptr< storage::Storage > _storage)
  : FixedListComm (_storage){}
  */


  LOG4ESPP_LOGGER(FixedQuadrupleList::theLogger, "FixedQuadrupleList");

  FixedQuadrupleList::FixedQuadrupleList(shared_ptr< storage::Storage > _storage) 
    : storage(_storage), globalQuadruples()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedQuadrupleList");

    sigBeforeSend = storage->beforeSendParticles.connect
      (boost::bind(&FixedQuadrupleList::beforeSendParticles, this, _1, _2));
    sigAfterRecv = storage->afterRecvParticles.connect
      (boost::bind(&FixedQuadrupleList::afterRecvParticles, this, _1, _2));
    sigOnParticlesChanged = storage->onParticlesChanged.connect
      (boost::bind(&FixedQuadrupleList::onParticlesChanged, this));
  }

  FixedQuadrupleList::~FixedQuadrupleList() {

    LOG4ESPP_INFO(theLogger, "~FixedQuadrupleList");

    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
    sigOnParticlesChanged.disconnect();
  }


  /*
  bool FixedQuadrupleList::add(longint pid1, longint pid2, longint pid3, longint pid4) {
      std::vector<longint> tmp;
      tmp.push_back(pid2);
      tmp.push_back(pid3);
      tmp.push_back(pid4);
      tmp.push_back(pid1); // this is used as key

      return FixedListComm::add(tmp);
  }*/


  bool FixedQuadrupleList::
  add(longint pid1, longint pid2, longint pid3, longint pid4) {
    // here we assume pid1 < pid2 < pid3 < pid4
    bool returnVal = true;
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);

    // ADD THE LOCAL QUADRUPLET
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    Particle *p3 = storage->lookupLocalParticle(pid3);
    Particle *p4 = storage->lookupLocalParticle(pid4);
    if (!p1){
      // Particle does not exist here, return false
      returnVal = false;
    }
    else{
      if (!p2) {
        std::stringstream msg;
        msg << "quadruple particle p2 " << pid2 << " does not exists here and cannot be added";
        err.setException( msg.str() );
      }
      if (!p3) {
        std::stringstream msg;
        msg << "quadruple particle p3 " << pid3 << " does not exists here and cannot be added";
        err.setException( msg.str() );
      }
      if (!p4) {
        std::stringstream msg;
        msg << "quadruple particle p4 " << pid4 << " does not exists here and cannot be added";
        err.setException( msg.str() );
      }
    }
    err.checkException();
    
    if(returnVal){
      // add the quadruple locally
      this->add(p1, p2, p3, p4);
      // ADD THE GLOBAL QUADRUPLET
      // see whether the particle already has quadruples
      std::pair<GlobalQuadruples::const_iterator,
                GlobalQuadruples::const_iterator> equalRange
        = globalQuadruples.equal_range(pid1);
      if (equalRange.first == globalQuadruples.end()) {
        // if it hasn't, insert the new quadruple
        globalQuadruples.insert(std::make_pair(pid1,
          Triple<longint, longint, longint>(pid2, pid3, pid4)));
      }
      else {
        // otherwise test whether the quadruple already exists
        for (GlobalQuadruples::const_iterator it = equalRange.first; it != equalRange.second; ++it)
  	      if (it->second == Triple<longint, longint, longint>(pid2, pid3, pid4))
  	        // TODO: Quadruple already exists, generate error!
  	    	;
        // if not, insert the new quadruple
        globalQuadruples.insert(equalRange.first,
          std::make_pair(pid1, Triple<longint, longint, longint>(pid2, pid3, pid4)));
      }
    }

    LOG4ESPP_INFO(theLogger, "added fixed quadruple to global quadruple list");
    return returnVal;
  }

  python::list FixedQuadrupleList::getQuadruples()
  {
	python::tuple quadruple;
	python::list quadruples;
	for (GlobalQuadruples::const_iterator it=globalQuadruples.begin(); it != globalQuadruples.end(); it++) {
      quadruple = python::make_tuple(it->first, it->second.first, it->second.second, it->second.third);
      quadruples.append(quadruple);
    }

	return quadruples;
  }

  std::vector<longint> FixedQuadrupleList::getQuadrupleList() {
    std::vector<longint> ret;
    for (GlobalQuadruples::const_iterator it=globalQuadruples.begin(); it != globalQuadruples.end(); it++) {
      ret.push_back(it->second.first);
      ret.push_back(it->first);
      ret.push_back(it->second.second);
      ret.push_back(it->second.third);
    }
    return ret;
  }

  void FixedQuadrupleList::
  beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
    
    std::vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      // LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find quadruples");
      //printf ("me = %d: send particle with pid %d find quadruples\n", mpiWorld->rank(), pid);

      // find all quadruples that involve this particle
      int n = globalQuadruples.count(pid);
      //printf ("me = %d: send particle with pid %d, has %d global quadruples\n", 
                //mpiWorld->rank(), pid, n);

      if (n > 0) {
	    std::pair<GlobalQuadruples::const_iterator, GlobalQuadruples::const_iterator> equalRange = globalQuadruples.equal_range(pid);
	    // first write the pid of this particle
	    // then the number of partners (n)
	    // and then the pids of the partners
	    toSend.reserve(toSend.size()+3*n+1);
	    toSend.push_back(pid);
	    toSend.push_back(n);
	    for (GlobalQuadruples::const_iterator it = equalRange.first; it != equalRange.second; ++it) {
	      toSend.push_back(it->second.first);
	      toSend.push_back(it->second.second);
	      toSend.push_back(it->second.third);
          //printf ("send global quadruple: pid %d and partner %d\n", pid, it->second.first);
          //printf ("send global quadruple: pid %d and partner %d\n", pid, it->second.second);
          //printf ("send global quadruple: pid %d and partner %d\n", pid, it->second.third);
        }
	    // delete all of these quadruples from the global list
	    globalQuadruples.erase(equalRange.first, equalRange.second);
      }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed quadruple list before send particles");
  }

  void FixedQuadrupleList::
  afterRecvParticles(ParticleList &pl, InBuffer& buf) {

    std::vector< longint > received;
    int n;
    longint pid1, pid2, pid3, pid4;
    GlobalQuadruples::iterator it = globalQuadruples.begin();
    // receive the quadruple list
    buf.read(received);
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid1 = received[i++];
      n = received[i++];
      //printf ("me = %d: recv particle with pid %d, has %d global quadruples\n",
                //mpiWorld->rank(), pid1, n);
      for (; n > 0; --n) {
	pid2 = received[i++];
	pid3 = received[i++];
	pid4 = received[i++];
	// add the quadruple to the global list
        //printf("received quadruple %d %d %d %d, add quadruple to global list\n", pid1, pid2, pid3, pid4);
	it = globalQuadruples.insert(it, std::make_pair(pid1,
          Triple<longint, longint, longint>(pid2, pid3, pid4)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed quadruple list after receive particles");
  }

  void FixedQuadrupleList::onParticlesChanged() {
    
    // (re-)generate the local quadruple list from the global list
    //printf("FixedQuadrupleList: rebuild local quadruple list from global\n");
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    Particle *p4;
    for (GlobalQuadruples::const_iterator it = globalQuadruples.begin(); it != globalQuadruples.end(); ++it) {
      //printf("lookup global quadruple %d %d %d %d\n",
      //it->first, it->second.first, it->second.second, it->second.third);
      if (it->first != lastpid1) {
	  p1 = storage->lookupRealParticle(it->first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p1 " << it->first << " does not exists here";
        err.setException( msg.str() );
      }
	  lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second.first);
      if (p2 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p2 " << it->second.first << " does not exists here";
        err.setException( msg.str() );
      }
      p3 = storage->lookupLocalParticle(it->second.second);
      if (p3 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p3 " << it->second.second << " does not exists here";
        err.setException( msg.str() );
      }
      p4 = storage->lookupLocalParticle(it->second.third);
      if (p4 == NULL) {
        std::stringstream msg;
        msg << "quadruple particle p4 " << it->second.third << " does not exists here";
        err.setException( msg.str() );
      }
      this->add(p1, p2, p3, p4);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed quadruple list from global list");
  }

  void FixedQuadrupleList::remove() {
      this->clear();
      globalQuadruples.clear();
      sigBeforeSend.disconnect();
      sigAfterRecv.disconnect();
  }
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedQuadrupleList::registerPython() {

    using namespace espressopp::python;

    bool (FixedQuadrupleList::*pyAdd)(longint pid1, longint pid2,
           longint pid3, longint pid4) = &FixedQuadrupleList::add;
    //bool (FixedQuadrupleList::*pyAdd)(pvec pids)
    //          = &FixedQuadrupleList::add;

    class_< FixedQuadrupleList, shared_ptr< FixedQuadrupleList > >
      ("FixedQuadrupleList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      .def("size", &FixedQuadrupleList::size)
      .def("remove",  &FixedQuadrupleList::remove)
      .def("getQuadruples",  &FixedQuadrupleList::getQuadruples)
     ;
  }
}
