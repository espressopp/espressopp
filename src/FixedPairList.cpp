/*
  Copyright (C) 2016
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
#include "FixedPairList.hpp"
//#include <utility>
//#include <algorithm>
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "boost/serialization/vector.hpp"
#include "Buffer.hpp"

#include "esutil/Error.hpp"

using namespace std;

namespace espressopp {

  /*
  FixedPairList::FixedPairList(shared_ptr< storage::Storage > _storage)
  : FixedListComm (_storage){}

  FixedPairList::~FixedPairList() {
    //std::cout << "~fixedpairlist" << std::endl;
    //FixedListComm::~FixedListComm();
  }*/


  LOG4ESPP_LOGGER(FixedPairList::theLogger, "FixedPairList");


  FixedPairList::FixedPairList(shared_ptr< storage::Storage > _storage)
    : storage(_storage), globalPairs()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedPairList");

    sigBeforeSend = storage->beforeSendParticles.connect
      (boost::bind(&FixedPairList::beforeSendParticles, this, _1, _2));
    sigAfterRecv = storage->afterRecvParticles.connect
      (boost::bind(&FixedPairList::afterRecvParticles, this, _1, _2));
    sigOnParticlesChanged = storage->onParticlesChanged.connect
      (boost::bind(&FixedPairList::onParticlesChanged, this));
  }

  FixedPairList::~FixedPairList() {

    LOG4ESPP_INFO(theLogger, "~FixedPairList");
    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
    sigOnParticlesChanged.disconnect();
  }


  /*
  bool FixedPairList::add(longint pid1, longint pid2) {
    std::vector<longint> tmp;
    tmp.push_back(pid2);
    tmp.push_back(pid1); // this is used as key

    return FixedListComm::add(tmp);
  }*/

  real FixedPairList::getLongtimeMaxBondSqr() {
	  return longtimeMaxBondSqr;
  }

  void FixedPairList::setLongtimeMaxBondSqr(real d) {
	  longtimeMaxBondSqr = d;
  }

  void FixedPairList::resetLongtimeMaxBondSqr() {
	  longtimeMaxBondSqr = 0.0;
  }

  bool FixedPairList::
  add(longint pid1, longint pid2) {
    bool returnVal = true;
    if (pid1 > pid2)
      std::swap(pid1, pid2);

    System& system = storage->getSystemRef();
    
    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    
    if (!p1){
      // Particle does not exist here, return false
      returnVal=false;
    }
    else{
      if (!p2) {
	LOG4ESPP_DEBUG(theLogger, "Particle p2 " << pid2 << " not found");
      }
    }
    
    if(returnVal){
      // add the pair locally
      this->add(p1, p2);
      // ADD THE GLOBAL PAIR
      // see whether the particle already has pairs
      std::pair<GlobalPairs::const_iterator,
        GlobalPairs::const_iterator> equalRange
        = globalPairs.equal_range(pid1);
      if (equalRange.first == globalPairs.end()) {
        // if it hasn't, insert the new pair
        globalPairs.insert(std::make_pair(pid1, pid2));
      }
      else {
        // otherwise test whether the pair already exists
        for (GlobalPairs::const_iterator it = equalRange.first; it != equalRange.second; ++it) {
  	    if (it->second == pid2) {
  	      // TODO: Pair already exists, generate error!
  	      ;
  	    }
        }
        // if not, insert the new pair
        globalPairs.insert(equalRange.first, std::make_pair(pid1, pid2));
      }
      LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    }
    LOG4ESPP_DEBUG(theLogger, "Leaving add with returnVal " << returnVal);
    return returnVal;
  }

  python::list FixedPairList::getBonds()
  {
	python::tuple bond;
	python::list bonds;
	for (GlobalPairs::const_iterator it=globalPairs.begin(); it != globalPairs.end(); it++) {
      bond = python::make_tuple(it->first, it->second);
      bonds.append(bond);
    }

	return bonds;
  }

  std::vector<longint> FixedPairList::getPairList() {
    std::vector<longint> ret;
    for (GlobalPairs::const_iterator it = globalPairs.begin(); it != globalPairs.end(); it++) {
      ret.push_back(it->first);
      ret.push_back(it->second);
    }
    return ret;
  }

  python::list FixedPairList::getAllBonds() {
    std::vector<longint> local_bonds;
    std::vector<std::vector<longint> > global_bonds;
    python::list bonds;

    for (GlobalPairs::const_iterator it = globalPairs.begin(); it != globalPairs.end(); it++) {
      local_bonds.push_back(it->first);
      local_bonds.push_back(it->second);
    }
    System& system = storage->getSystemRef();
    if (system.comm->rank() == 0) {
      mpi::gather(*system.comm, local_bonds, global_bonds, 0);
      python::tuple bond;

      for (std::vector<std::vector<longint> >::iterator it = global_bonds.begin();
           it != global_bonds.end(); ++it) {
        for (std::vector<longint>::iterator iit = it->begin(); iit != it->end();) {
          longint pid1 = *(iit++);
          longint pid2 = *(iit++);
          bonds.append(python::make_tuple(pid1, pid2));
        }
      }
    } else {
      mpi::gather(*system.comm, local_bonds, global_bonds, 0);
    }
    return bonds;
  }

  void FixedPairList::beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
    std::vector< longint > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      
      LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pid << ", find pairs");

      // find all pairs that involve this particle
      
      int n = globalPairs.count(pid);

      if (n > 0) {
        std::pair<GlobalPairs::const_iterator,
          GlobalPairs::const_iterator> equalRange
          = globalPairs.equal_range(pid);

        // first write the pid of the first particle
        // then the number of partners
        // and then the pids of the partners
        toSend.reserve(toSend.size()+n+1);
        toSend.push_back(pid);
        toSend.push_back(n);
        for (GlobalPairs::const_iterator it = equalRange.first;
             it != equalRange.second; ++it) {
          toSend.push_back(it->second);
          LOG4ESPP_DEBUG(theLogger, "send global bond: pid "
                       << pid << " and partner " << it->second);
        }

        // delete all of these pairs from the global list
        globalPairs.erase(equalRange.first, equalRange.second);
        // std::cout << "erasing particle " << pid << " from here" << std::endl;
      }
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
  }

  void FixedPairList::afterRecvParticles(ParticleList &pl, InBuffer& buf) {
    std::vector< longint > received;
    int n;
    longint pid1, pid2;
    GlobalPairs::iterator it = globalPairs.begin();
    // receive the bond list
    buf.read(received);
    int size = received.size(); int i = 0;
    while (i < size) {
      // unpack the list
      pid1 = received[i++];
      n = received[i++];
      LOG4ESPP_DEBUG(theLogger, "recv particle " << pid1 << 
                                ", has " << n << " global pairs");
      for (; n > 0; --n) {
	pid2 = received[i++];
	// add the bond to the global list
        LOG4ESPP_DEBUG(theLogger, "received pair " << pid1 << " , " << pid2);
	it = globalPairs.insert(it, std::make_pair(pid1, pid2));
      }
    }
    if (i != size) {
      LOG4ESPP_ERROR(theLogger, 
        "ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed pair list after receive particles");
  }

  void FixedPairList::onParticlesChanged() {
    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    for (GlobalPairs::const_iterator it = globalPairs.begin(); it != globalPairs.end(); ++it) {
      if (it->first != lastpid1) {
	    p1 = storage->lookupRealParticle(it->first);
        if (p1 == NULL) {
          std::stringstream msg;
          msg << "onParticlesChanged error. Fixed Pair List particle p1 " << it->first << " does not exists here";
          err.setException( msg.str() );
          //std::runtime_error(err.str());
        }
	    lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second);
      if (p2 == NULL) {
          std::stringstream msg;
          msg << "onParticlesChanged error. Fixed Pair List particle p2 " << it->second << " does not exists here";
          //std::runtime_error(err.str());
          err.setException( msg.str() );
      }
      this->add(p1, p2);
    }
    err.checkException();
    
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  void FixedPairList::remove() {
      this->clear();
      globalPairs.clear();
      sigBeforeSend.disconnect();
      sigAfterRecv.disconnect();
      sigOnParticlesChanged.disconnect();
  }

  int FixedPairList::totalSize() {
    int local_size = globalPairs.size();
    int global_size;
    System& system = storage->getSystemRef();
    mpi::all_reduce(*system.comm, local_size, global_size, std::plus<int>());
    return global_size;
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairList::registerPython() {

    using namespace espressopp::python;

    bool (FixedPairList::*pyAdd)(longint pid1, longint pid2)
      = &FixedPairList::add;
    //bool (FixedPairList::*pyAdd)(pvec pids) = &FixedPairList::add;

    class_<FixedPairList, shared_ptr<FixedPairList> >
      ("FixedPairList", init <shared_ptr<storage::Storage> >())
      .def("add", pyAdd)
      .def("size", &FixedPairList::size)
      .def("totalSize", &FixedPairList::totalSize)
      .def("getBonds",  &FixedPairList::getBonds)
      .def("remove",  &FixedPairList::remove)
      .def("getAllBonds", &FixedPairList::getAllBonds)
      .def("resetLongtimeMaxBondSqr", &FixedPairList::resetLongtimeMaxBondSqr)
      .def("getLongtimeMaxBondSqr", &FixedPairList::getLongtimeMaxBondSqr)
      ;
  }
}
