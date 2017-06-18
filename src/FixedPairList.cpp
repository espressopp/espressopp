/*
  Copyright (C) 2012,2013,2016
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  Copyright (C) 2017
      Jakub Krajniak
  
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
#include "FixedPairList.hpp"
#include "storage/Storage.hpp"

#include "esutil/Error.hpp"

using namespace std;

namespace espressopp {

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
    esutil::Error err(system.comm);

    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);
    
    if (!p1){
      // Particle does not exist here, return false
      returnVal=false;
    }
    else{
      if (!p2) {
        std::stringstream msg;
        msg << "bond particle p2 " << pid2 << " does not exists here and cannot be added";
        err.setException(msg.str());
      }
    }

    err.checkException();

    if (returnVal) {
      // ADD THE GLOBAL PAIR
      // see whether the particle already has pairs
      bool found = false;
      std::pair<GlobalPairs::const_iterator, GlobalPairs::const_iterator> equalRange =
          globalPairs.equal_range(pid1);
      if (equalRange.first != globalPairs.end()) {
        // otherwise test whether the pair already exists
        for (GlobalPairs::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
          if (it->second == pid2)
            found = true;
        }
      }
      returnVal = !found;
      if (!found) {
        // add the pair locally
        this->add(p1, p2);
        // Update list of integers.
        globalPairs.insert(equalRange.first, std::make_pair(pid1, pid2));
        // Throw signal onTupleAdded.
        onTupleAdded(pid1, pid2);
        LOG4ESPP_INFO(theLogger, "added fixed pair " << pid1 << "-" << pid2 << " to global pair list");
      }
    }
    return returnVal;
  }

  bool FixedPairList::iadd(longint pid1, longint pid2) {
    bool returnVal = true;
    if (pid1 > pid2)
      std::swap(pid1, pid2);

    System& system = storage->getSystemRef();

    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupRealParticle(pid1);
    Particle *p2 = storage->lookupLocalParticle(pid2);

    if (!p1) {
      // Particle does not exist here, return false
      returnVal = false;
    } else if (!p2) {
        returnVal = false;
    }

    if (returnVal) {
      // ADD THE GLOBAL PAIR
      // see whether the particle already has pairs
      bool found = false;
      std::pair<GlobalPairs::const_iterator, GlobalPairs::const_iterator> equalRange =
          globalPairs.equal_range(pid1);
      if (equalRange.first != globalPairs.end()) {
        // otherwise test whether the pair already exists
        for (GlobalPairs::const_iterator it = equalRange.first; it != equalRange.second && !found; ++it) {
          if (it->second == pid2)
            found = true;
        }
      }
      returnVal = !found;
      if (!found) {
        // add the pair locally
        this->add(p1, p2);
        // Update list of integers.
        globalPairs.insert(equalRange.first, std::make_pair(pid1, pid2));
        // Throw signal onTupleAdded.
        onTupleAdded(pid1, pid2);
        LOG4ESPP_INFO(theLogger, "added fixed pair " << pid1 << "-" << pid2 << " to global pair list");
      }
    }
    return returnVal;
  }

  bool FixedPairList::remove(longint pid1, longint pid2, bool no_signal) {
    LOG4ESPP_DEBUG(theLogger, "FPL remove " << pid1 << "-" << pid2);
    bool returnValue = false;
    std::pair<GlobalPairs::iterator, GlobalPairs::iterator> equalRange, equalRange_rev;
    equalRange = globalPairs.equal_range(pid1);
    if (equalRange.first != globalPairs.end()) {
      for (GlobalPairs::iterator it = equalRange.first; it != equalRange.second;) {
        if (it->second == pid2) {
          LOG4ESPP_DEBUG(theLogger, "FPL, found " << it->first << " - " << it->second);
          if (!no_signal)
            onTupleRemoved(pid1, pid2);
          it = globalPairs.erase(it);
          returnValue = true;
        } else {
          it++;
        }
      }
    }
    if (!returnValue) {
      equalRange_rev = globalPairs.equal_range(pid2);
      if (equalRange_rev.first != globalPairs.end()) {
        for (GlobalPairs::iterator it = equalRange_rev.first; it != equalRange_rev.second;) {
          if (it->second == pid1) {
            LOG4ESPP_DEBUG(theLogger, "FPL, found " << it->first << " - " << it->second);
            if (!no_signal)
              onTupleRemoved(pid1, pid2);
            it = globalPairs.erase(it);
            returnValue = true;
          } else {
            it++;
          }
        }
      }
    }
    return returnValue;
  }

  bool FixedPairList::removeByPid1(longint pid1, bool noSignal, bool removeAll, longint removeCounter) {
    bool returnValue = false;
    std::pair<GlobalPairs::iterator, GlobalPairs::iterator> equalRange = globalPairs.equal_range(pid1);
    if (equalRange.first == globalPairs.end())
      return returnValue;

    if (removeAll) {
      for(GlobalPairs::iterator it = equalRange.first; it != equalRange.second;) {
        if (!noSignal)
          onTupleRemoved(it->first, it->second);
        it = globalPairs.erase(it);
        returnValue = true;
      }
    } else {
      longint num_removed = 0;
      for(GlobalPairs::iterator it = equalRange.first; it != equalRange.second && num_removed < removeCounter;) {
        if (!noSignal)
          onTupleRemoved(it->first, it->second);
        it = globalPairs.erase(it);
        returnValue = true;
        num_removed++;
      }
    }
    return returnValue;
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
          msg << "onParticlesChanged error. Fixed Pair List particle p1 " << it->first << " does not exists here.";
          msg << " pair: " << it->first << "-" << it->second;
          err.setException( msg.str() );
          //std::runtime_error(err.str());
        }
	    lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second);
      if (p2 == NULL) {
          std::stringstream msg;
          msg << "onParticlesChanged error. Fixed Pair List particle p2 " << it->second << " does not exists here.";
          msg << " p1: " << *p1;
          msg << " pair: " << it->first << "-" << it->second;
          //std::runtime_error(err.str());
          err.setException( msg.str() );
      }
      this->add(p1, p2);
    }
    err.checkException();
    
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  void FixedPairList::updateParticlesStorage() {
    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

    System& system = storage->getSystemRef();

    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;
    for (GlobalPairs::const_iterator it = globalPairs.begin(); it != globalPairs.end(); ++it) {
      if (it->first != lastpid1) {
        p1 = storage->lookupRealParticle(it->first);
        if (p1 == NULL) {
          std::stringstream msg;
          msg << "updateParticlesStorage error. Fixed Pair List particle p1 " << it->first << " does not exists here.";
          msg << " p1: " << *p1;
          msg << " pair: " << it->first << "-" << it->second;
          throw std::runtime_error(msg.str());
        }
        lastpid1 = it->first;
      }
      p2 = storage->lookupLocalParticle(it->second);
      if (p2 == NULL) {
        std::stringstream msg;
        msg << "updateParticlesStorage error. Fixed Pair List particle p2 " << it->second << " does not exists here.";
        msg << " p1: " << *p1;
        msg << " pair: " << it->first << "-" << it->second;
        throw std::runtime_error(msg.str());
      }
      this->add(p1, p2);
    }
    LOG4ESPP_INFO(theLogger, "regenerated local fixed pair list from global list");
  }

  int FixedPairList::totalSize() {
    int local_size = globalPairs.size();
    int global_size;
    System& system = storage->getSystemRef();
    mpi::all_reduce(*system.comm, local_size, global_size, std::plus<int>());
    return global_size;
  }

  void FixedPairList::clearAndRemove() {
      this->clear();
      globalPairs.clear();
      sigBeforeSend.disconnect();
      sigAfterRecv.disconnect();
      sigOnParticlesChanged.disconnect();
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairList::registerPython() {

    using namespace espressopp::python;

    bool (FixedPairList::*pyAdd)(longint pid1, longint pid2)
      = &FixedPairList::add;
    //bool (FixedPairList::*pyAdd)(pvec pids) = &FixedPairList::add;

    class_<FixedPairList, shared_ptr<FixedPairList>, boost::noncopyable >
      ("FixedPairList", init <shared_ptr<storage::Storage> >())
      .def("add", pyAdd)
      .def("remove", &FixedPairList::remove)
      .def("removeByPid1", &FixedPairList::removeByPid1)
      .def("size", &FixedPairList::size)
      .def("totalSize", &FixedPairList::totalSize)
      .def("getBonds",  &FixedPairList::getBonds)
      .def("getAllBonds", &FixedPairList::getAllBonds)
      .def("clearAndRemove",  &FixedPairList::clearAndRemove)
      .def("resetLongtimeMaxBondSqr", &FixedPairList::resetLongtimeMaxBondSqr)
      .def("getLongtimeMaxBondSqr", &FixedPairList::getLongtimeMaxBondSqr)
      ;
  }
}
