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

  bool FixedTripleList::iadd(longint pid1, longint pid2, longint pid3) {
    bool returnVal = true;
    System& system = storage->getSystemRef();

    // ADD THE LOCAL TRIPLET
    Particle *p1 = storage->lookupLocalParticle(pid1);
    Particle *p2 = storage->lookupRealParticle(pid2);
    Particle *p3 = storage->lookupLocalParticle(pid3);

    // middle particle is the reference particle and must exist here
    if (!p2){
      // particle does not exists here (some other CPU must have it)
      returnVal = false;
    } else {
      std::stringstream msg;
      if (!p1) {
        msg << "adding error: triple particle p1 " << pid1 <<
            " does not exists here and cannot be added";
        msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
        throw std::runtime_error(msg.str());
      }
      if (!p3) {
        msg << "adding error: triple particle p3 " << pid3 <<
            " does not exists here and cannot be added";
        msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
        throw std::runtime_error(msg.str());
      }
    }

    if (returnVal) {
      // ADD THE GLOBAL TRIPLET
      // see whether the particle already has triples
      bool found = false;
      std::pair<GlobalTriples::const_iterator,
                GlobalTriples::const_iterator> equalRange = globalTriples.equal_range(pid2);
      if (equalRange.first != globalTriples.end()) {
        // otherwise test whether the triple already exists
        for (GlobalTriples::const_iterator it = equalRange.first;
             it != equalRange.second && !found; ++it) {
          if (it->second == std::pair<longint, longint>(pid1, pid3) ||
              it->second == std::pair<longint, longint>(pid3, pid1))
            found = true;
        }
      }
      returnVal = !found;
      if (!found) {
        // add the triple locally
        this->add(p1, p2, p3);
        globalTriples.insert(equalRange.first,
                             std::make_pair(pid2, std::pair<longint, longint>(pid1, pid3)));
        onTupleAdded(pid1, pid2, pid3);
      }
      LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
    }
    return returnVal;
  }

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
    
    if (returnVal) {
      // ADD THE GLOBAL TRIPLET
      // see whether the particle already has triples
      bool found = false;
      std::pair<GlobalTriples::const_iterator,
                GlobalTriples::const_iterator> equalRange = globalTriples.equal_range(pid2);
      if (equalRange.first != globalTriples.end()) {
        // otherwise test whether the triple already exists
        for (GlobalTriples::const_iterator it = equalRange.first;
             it != equalRange.second && !found; ++it) {
          if (it->second == std::pair<longint, longint>(pid1, pid3) ||
              it->second == std::pair<longint, longint>(pid3, pid1)) {
            found = true;
          }
        }
      }
      returnVal = !found;
      if (!found) {
        // add the triple locally
        this->add(p1, p2, p3);
        globalTriples.insert(equalRange.first,
            std::make_pair(pid2, std::pair<longint, longint>(pid1, pid3)));
        onTupleAdded(pid1, pid2, pid3);
      }
      LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
    }
    return returnVal;
  }

  bool FixedTripleList::remove(longint pid1, longint pid2, longint pid3, bool no_signal) {
    bool returnVal = false;
    // Remove entries.
    std::pair<GlobalTriples::iterator, GlobalTriples::iterator> equalRange =
        globalTriples.equal_range(pid2);
    if (equalRange.first != globalTriples.end()) {
      // otherwise test whether the triple already exists
      for (GlobalTriples::iterator it = equalRange.first; it != equalRange.second;) {
        if (it->second == std::pair<longint, longint>(pid1, pid3) ||
            it->second == std::pair<longint, longint>(pid3, pid1)) {
          LOG4ESPP_DEBUG(theLogger, "removed triple " << it->first << "-" << it->second.first
              << "-" << it->second.second
              << " bond: " << pid1 << "-" << pid2);
          if (!no_signal)
            onTupleRemoved(it->second.first, it->first, it->second.second);
          it = globalTriples.erase(it);
          returnVal = true;
        } else {
          ++it;
        }
      }
    }
    return returnVal;
  }

  bool FixedTripleList::removeByBond(longint pid1, longint pid2) {
    bool return_val = false;

    for (GlobalTriples::iterator it = globalTriples.begin(); it != globalTriples.end();) {
      longint a1 = it->second.first;
      longint a2 = it->first;
      longint a3 = it->second.second;
      if ((a1 == pid1 && a2 == pid2) || (a2 == pid1 && a3 == pid2) ||
          (a1 == pid2 && a2 == pid1) || (a2 == pid2 && a3 == pid1)) {
        it = globalTriples.erase(it);
        return_val = true;
      } else {
        ++it;
      }
    }
    return return_val;
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

  python::list FixedTripleList::getAllTriples() {
    std::vector<longint> local_triplets;
    std::vector<std::vector<longint> > global_triplets;
    python::list triplets;

    for (GlobalTriples::const_iterator it=globalTriples.begin(); it != globalTriples.end(); it++) {
      local_triplets.push_back(it->second.first);
      local_triplets.push_back(it->first);
      local_triplets.push_back(it->second.second);
    }

    System& system = storage->getSystemRef();
    if (system.comm->rank() == 0) {
      mpi::gather(*system.comm, local_triplets, global_triplets, 0);

      for (std::vector<std::vector<longint> >::iterator it = global_triplets.begin();
           it != global_triplets.end(); ++it) {
        for (std::vector<longint>::iterator iit = it->begin(); iit != it->end();) {
          longint pid1 = *(iit++);
          longint pid2 = *(iit++);
          longint pid3 = *(iit++);
          triplets.append(python::make_tuple(pid1, pid2, pid3));
        }
      }
    } else {
      mpi::gather(*system.comm, local_triplets, global_triplets, 0);
    }
    return triplets;
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

  void FixedTripleList::updateParticlesStorage() {
    System& system = storage->getSystemRef();

    // (re-)generate the local triple list from the global list
    this->clear();
    longint lastpid2 = -1;
    Particle *p1;
    Particle *p2;
    Particle *p3;
    for (GlobalTriples::const_iterator it = globalTriples.begin(); it != globalTriples.end(); ++it) {
      if (it->first != lastpid2) {
        p2 = storage->lookupRealParticle(it->first);
        if (p2 == NULL) {
          std::stringstream msg;
          msg << "triple particle p2 " << it->first << " does not exists here";
          throw std::runtime_error(msg.str());
        }
        lastpid2 = it->first;
      }
      p1 = storage->lookupLocalParticle(it->second.first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "triple particle p1 " << it->second.first << " does not exists here";
        throw std::runtime_error(msg.str());
      }
      p3 = storage->lookupLocalParticle(it->second.second);
      if (p3 == NULL) {
        std::stringstream msg;
        msg << "triple particle p3 " << it->second.second << " does not exists here";
        throw std::runtime_error(msg.str());
      }
      this->add(p1, p2, p3);
    }

    LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
  }

  int FixedTripleList::totalSize() {
    int local_size = globalTriples.size();
    int global_size;
    System& system = storage->getSystemRef();
    mpi::all_reduce(*system.comm, local_size, global_size, std::plus<int>());
    return global_size;
  }

  void FixedTripleList::clearAndRemove() {
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

    class_< FixedTripleList, shared_ptr< FixedTripleList >, boost::noncopyable  >
      ("FixedTripleList", init< shared_ptr< storage::Storage > >())
      .def("add", pyAdd)
      .def("size", &FixedTripleList::size)
      .def("totalSize", &FixedTripleList::totalSize)
      .def("getTriples",  &FixedTripleList::getTriples)
      .def("getAllTriples", &FixedTripleList::getAllTriples)
      .def("remove", &FixedTripleList::remove)
      .def("clearAndRemove", &FixedTripleList::clearAndRemove)
     ;
  }
}
