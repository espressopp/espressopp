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
#include "FixedPairDistList.hpp"
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"
#include "bc/BC.hpp"
#include "esutil/Error.hpp"

namespace espressopp {

  LOG4ESPP_LOGGER(FixedPairDistList::theLogger, "FixedPairDistList");


  FixedPairDistList::FixedPairDistList(shared_ptr< storage::Storage > _storage) 
    : storage(_storage), pairsDist()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedPairDistList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairDistList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairDistList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedPairDistList::onParticlesChanged, this));
  }

  FixedPairDistList::~FixedPairDistList() {

    LOG4ESPP_INFO(theLogger, "~FixedPairDistList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedPairDistList::
  add(longint pid1, longint pid2) {
    bool returnVal = true;
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    if (pid1 > pid2)
      std::swap(pid1, pid2);

    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    if (!p1){
      // Particle does not exist here, return false
      returnVal = false;
    }
    else{
      if (!p2) {
        std::stringstream msg;
        msg << "bond particle p2 " << pid2 << " does not exists here and cannot be added";
        //std::runtime_error(err.str());
        err.setException( msg.str() );
      }
    }
    err.checkException();
    
    if(returnVal){
      // add the pair locally
      this->add(p1, p2);

      Real3D pos1 = p1->position();
      Real3D pos2 = p2->position();

      Real3D r21(0.0,0.0,0.0);
      system.bc->getMinimumImageVectorBox(r21, pos1, pos2);

      // add the particle pair to the globalPairs list
      pairsDist.insert(std::make_pair(pid1, std::make_pair(pid2, r21.abs()) ));
    }

    LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    return returnVal;
  }

  python::list FixedPairDistList::getPairs(){
	python::tuple pair;
	python::list pairs;
	for (PairsDist::const_iterator it=pairsDist.begin(); it!=pairsDist.end(); it++) {
      pair = python::make_tuple(it->first, it->second.first);
      pairs.append(pair);
    }

	return pairs;
  }
  
  python::list FixedPairDistList::getPairsDist(){
	python::tuple pair;
	python::list pairs;
	for (PairsDist::const_iterator it=pairsDist.begin(); it!=pairsDist.end(); it++) {
      pair = python::make_tuple(it->first, it->second.first, it->second.second);
      pairs.append(pair);
    }

	return pairs;
  }
  
  real FixedPairDistList::getDist(int pid1, int pid2){
    real returnVal = -3;
    
    PairsDist::iterator itr;
	PairsDist::iterator lastElement;
	
	// locate an iterator to the first pair object associated with key
	itr = pairsDist.find(pid1);
	if (itr == pairsDist.end())
      return returnVal; // no elements associated with key, so return immediately

	// get an iterator to the element that is one past the last element associated with key
	lastElement = pairsDist.upper_bound(pid1);
    
	for ( ; itr != lastElement; ++itr){
      if(pid2==itr->second.first){
        returnVal = itr->second.second;
        break;
      }
    }
    
	return returnVal;
  }

  void FixedPairDistList::
  beforeSendParticles(ParticleList& pl, 
		      OutBuffer& buf) {
    std::vector< longint > toSendInt;
    std::vector< real > toSendReal;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

      // find all pairs that involve this particle
      
      int n = pairsDist.count(pid);

      if (n > 0) {
        std::pair<PairsDist::const_iterator,
          PairsDist::const_iterator> equalRange
          = pairsDist.equal_range(pid);

        // first write the pid of the first particle
        // then the number of partners
        // and then the pids of the partners
        toSendInt.reserve(toSendInt.size()+n+1);
        toSendReal.reserve(toSendReal.size()+n);
        toSendInt.push_back(pid);
        toSendInt.push_back(n);
        for (PairsDist::const_iterator it = equalRange.first;
             it != equalRange.second; ++it) {
          toSendInt.push_back(it->second.first);
          toSendReal.push_back(it->second.second);
        }

        // delete all of these pairs from the global list
        //globalPairs.erase(equalRange.first->first, equalRange.second->first);
        pairsDist.erase(pid);
        // std::cout << "erasing particle " << pid << " from here" << std::endl;
      }
    }
    // send the list
    buf.write(toSendInt);
    buf.write(toSendReal);
    LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
  }

  void FixedPairDistList::
  afterRecvParticles(ParticleList &pl, 
		     InBuffer& buf) {
    std::vector< longint > receivedInt;
    std::vector< real > receivedReal;
    int n;
    longint pid1, pid2;
    real distVal;
    PairsDist::iterator it = pairsDist.begin();
    // receive the bond list
    buf.read(receivedInt);
    buf.read(receivedReal);
    int size = receivedInt.size(); int i = 0; int j = 0;
    while (i < size) {
      // unpack the list
      pid1 = receivedInt[i++];
      n = receivedInt[i++];
      for (; n > 0; --n) {
        pid2 = receivedInt[i++];
        distVal = receivedReal[j++];
        it = pairsDist.insert(it, std::make_pair(pid1, std::make_pair(pid2, distVal)));
      }
    }
    if (i != size) {
      LOG4ESPP_ERROR(theLogger, 
        "ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
  }

  void FixedPairDistList::
  onParticlesChanged() {
    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");
    
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    for (PairsDist::const_iterator it = pairsDist.begin();
	 it != pairsDist.end(); ++it) {
      if (it->first != lastpid1) {
	    p1 = storage->lookupRealParticle(it->first);
        if (p1 == NULL) {
          std::stringstream msg;
          msg << "bond particle p1 " << it->first << " does not exists here";
          //std::runtime_error(err.str());
          err.setException( msg.str() );
        }
	    lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second.first);
      if (p2 == NULL) {
          std::stringstream msg;
          msg << "bond particle p2 " << it->second.first << " does not exists here";
          //std::runtime_error(err.str());
          err.setException( msg.str() );
      }
      this->add(p1, p2);
    }
    err.checkException();
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairDistList::registerPython() {

    using namespace espressopp::python;

    bool (FixedPairDistList::*pyAdd)(longint pid1, longint pid2)
      = &FixedPairDistList::add;
    //bool (FixedPairDistList::*pyAdd)(pvec pids) = &FixedPairDistList::add;

    class_<FixedPairDistList, shared_ptr<FixedPairDistList> >
      ("FixedPairDistList", init <shared_ptr<storage::Storage> >())
      .def("add", pyAdd)
      .def("size", &FixedPairDistList::size)
      .def("getPairs",  &FixedPairDistList::getPairs)
      .def("getPairsDist",  &FixedPairDistList::getPairsDist)
      .def("getDist",  &FixedPairDistList::getDist)
      ;
  }
}
