/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz
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

#include "vec/Vectorization.hpp"
#include "vec/FixedTripleList.hpp"
#include "vec/storage/StorageVec.hpp"

#include "python.hpp"
#include "storage/Storage.hpp"
#include "Buffer.hpp"
#include "esutil/Error.hpp"

#include <sstream>
#include <boost/bind.hpp>

namespace espressopp { namespace vec {

  LOG4ESPP_LOGGER(FixedTripleList::theLogger, "FixedTripleList");

  FixedTripleList::FixedTripleList(std::shared_ptr<espressopp::storage::Storage> storage)
    : globalTriples()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedTripleList");

    if(!storage->getSystem()->vectorization) {
      throw std::runtime_error("system has no vectorization");
    }
    vectorization = storage->getSystem()->vectorization;

    if(!(vectorization->storageVec))
      throw std::runtime_error("vectorization->storageVec cannot be null");
    auto& storageVec = vectorization->storageVec;
    storageVec->enableLocalParticles();

    sigBeforeSend = storage->beforeSendParticles.connect
      (boost::bind(&FixedTripleList::beforeSendParticles, this, _1, _2));
    sigAfterRecv = storage->afterRecvParticles.connect
      (boost::bind(&FixedTripleList::afterRecvParticles, this, _1, _2));
    sigOnParticlesChanged = storage->onParticlesChanged.connect(
      boost::bind(&FixedTripleList::onParticlesChanged, this));
  }

  FixedTripleList::~FixedTripleList()
  {
    LOG4ESPP_INFO(theLogger, "~FixedTripleList");

    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
    sigOnParticlesChanged.disconnect();
  }

  bool FixedTripleList::
  add(size_t pid1, size_t pid2, size_t pid3)
  {
    bool returnVal = true;
    auto& system  = vectorization->getSystemRef();
    esutil::Error err(system.comm);

    auto const& storageVec = vectorization->storageVec;
    size_t const p1 = storageVec->lookupLocalParticleVec(pid1);
    size_t const p2 = storageVec->lookupRealParticleVec(pid2);
    size_t const p3 = storageVec->lookupLocalParticleVec(pid3);

    // middle particle is the reference particle and must exist here
    if (p2==VEC_PARTICLE_NOT_FOUND){
      // particle does not exists here (some other CPU must have it)
      returnVal = false;
    } else {
      if (p1==VEC_PARTICLE_NOT_FOUND) {
        std::stringstream msg;
        msg << "adding error: triple particle p1 " << pid1 << " does not exists here and cannot be added";
        msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
        err.setException(msg.str());
      }
      if (p3==VEC_PARTICLE_NOT_FOUND) {
        std::stringstream msg;
        msg << "adding error: triple particle p3 " << pid3 << " does not exists here and cannot be added";
        msg << " triplet: " << pid1 << "-" << pid2 << "-" << pid3;
        err.setException(msg.str());
      }
    }
    err.checkException();

    if(returnVal){
      // add the triple locally
      this->push_back({p1, p2, p3});

      // ADD THE GLOBAL TRIPLET
      // see whether the particle already has triples
      const auto equalRange = globalTriples.equal_range(pid2);
      if (equalRange.first == globalTriples.end()) {
        // if it hasn't, insert the new triple
        globalTriples.insert(std::make_pair(pid2, std::pair<size_t, size_t>(pid1, pid3)));
      } else {
        // otherwise test whether the triple already exists
        for (auto it = equalRange.first; it != equalRange.second; ++it) {
          if (it->second == std::pair<size_t, size_t>(pid1, pid3)) {
            // TODO: Triple already exists, generate error!
            ;
          }
        }
        // if not, insert the new triple
        globalTriples.insert(equalRange.first, std::make_pair(pid2, std::pair<size_t, size_t>(pid1, pid3)));
      }
      LOG4ESPP_INFO(theLogger, "added fixed triple to global triple list");
    }
    return returnVal;
  }

  python::list FixedTripleList::getTriples()
  {
    python::tuple triple;
    python::list triples;
    for (auto it = globalTriples.cbegin(); it != globalTriples.cend(); it++) {
      triple = python::make_tuple( it->second.first,it->first, it->second.second);
      triples.append(triple);
    }
    return triples;
  }

  std::vector<size_t> FixedTripleList::getTripleList()
  {
    std::vector<size_t> ret;
    for (auto it = globalTriples.cbegin(); it != globalTriples.cend(); it++) {
      ret.push_back(it->second.first);
      ret.push_back(it->first);
      ret.push_back(it->second.second);
    }
    return ret;
  }

  void FixedTripleList::
  beforeSendParticles(ParticleList& pl, OutBuffer& buf)
  {
    std::vector< size_t > toSend;
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      const size_t pid = pit->id();

      // find all triples that involve this particle
      size_t n = globalTriples.count(pid);

      if (n > 0)
      {
        const auto equalRange = globalTriples.equal_range(pid);

        // first write the pid of this particle
        // then the number of partners (n)
        // and then the pids of the partners
        toSend.reserve(toSend.size()+2*n+1);
        toSend.push_back(pid);
        toSend.push_back(n);
        for (auto it = equalRange.first; it != equalRange.second; ++it)
        {
          toSend.push_back(it->second.first);
          toSend.push_back(it->second.second);
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
  afterRecvParticles(ParticleList &pl, InBuffer& buf)
  {
    std::vector< size_t > received;
    size_t n;
    size_t pid1, pid2, pid3;
    auto it = globalTriples.begin();
    // receive the triple list
    buf.read(received);
    size_t size = received.size(); size_t i = 0;
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
        it = globalTriples.insert(it, std::make_pair(pid2,std::pair<size_t, size_t>(pid1, pid3)));
      }
    }
    if (i != size) {
      printf("ATTETNTION:  recv particles might have read garbage\n");
    }
    LOG4ESPP_INFO(theLogger, "received fixed triple list after receive particles");
  }

  void FixedTripleList::onParticlesChanged()
  {
    auto& system  = vectorization->getSystemRef();
    esutil::Error err(system.comm);

    // (re-)generate the local triple list from the global list
    //printf("FixedTripleList: rebuild local triple list from global\n");
    this->clear();
    size_t lastpid2 = VEC_PARTICLE_NOT_FOUND;
    auto const& storageVec = vectorization->storageVec;
    size_t p1, p2, p3;
    for (auto it = globalTriples.cbegin(); it != globalTriples.cend(); ++it)
    {
      //printf("lookup global triple %d %d %d\n", it->first, it->second.first, it->second.second);
      if (it->first != lastpid2) {
        p2 = storageVec->lookupRealParticleVec(it->first);
        if (p2 == VEC_PARTICLE_NOT_FOUND) {
          std::stringstream msg;
          msg << "triple particle p2 " << it->first << " does not exists here";
          err.setException( msg.str() );
        }
        lastpid2 = it->first;
      }
      p1 = storageVec->lookupLocalParticleVec(it->second.first);
      if (p1 == VEC_PARTICLE_NOT_FOUND) {
        std::stringstream msg;
        msg << "triple particle p1 " << it->second.first << " does not exists here";
        err.setException( msg.str() );
      }
      p3 = storageVec->lookupLocalParticleVec(it->second.second);
      if (p3 == VEC_PARTICLE_NOT_FOUND) {
        std::stringstream msg;
        msg << "triple particle p3 " << it->second.second << " does not exists here";
        err.setException( msg.str() );
      }
      this->push_back({p1, p2, p3});
    }
    err.checkException();

    LOG4ESPP_INFO(theLogger, "regenerated local fixed triple list from global list");
  }

  void FixedTripleList::remove()
  {
    this->clear();
    globalTriples.clear();
    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
    sigOnParticlesChanged.disconnect();
  }
  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedTripleList::registerPython()
  {
    using namespace espressopp::python;

    bool (FixedTripleList::*pyAdd)(size_t pid1, size_t pid2, size_t pid3)
      = &FixedTripleList::add;

    class_< FixedTripleList, std::shared_ptr< FixedTripleList > >
      ("vec_FixedTripleList", init< std::shared_ptr<espressopp::storage::Storage> >())
      .def("add", pyAdd)
      .def("size", &FixedTripleList::size)
      .def("remove",  &FixedTripleList::remove)
      .def("getTriples",  &FixedTripleList::getTriples)
     ;
  }
}}
