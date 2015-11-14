/*
  Copyright (c) 2014
      Jakub Krajniak (jkrajniak at gmail.com)
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
#include "FixedQuadrupleListAdress.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "esutil/Error.hpp"


namespace espressopp {
LOG4ESPP_LOGGER(FixedQuadrupleListAdress::theLogger, "FixedQuadrupleListAdress");

FixedQuadrupleListAdress::FixedQuadrupleListAdress(
    shared_ptr< storage::Storage > _storage, shared_ptr<FixedTupleListAdress> _fixedtupleList)
    : FixedQuadrupleList(_storage), fixedtupleList(_fixedtupleList) {
  LOG4ESPP_INFO(theLogger, "construct FixedQuadrupleListAdress");

  sigBeforeSendAT = fixedtupleList->beforeSendATParticles.connect(
    boost::bind(&FixedQuadrupleListAdress::beforeSendATParticles, this, _1, _2));

  sigAfterRecvAT = fixedtupleList->afterRecvATParticles.connect
    (boost::bind(&FixedQuadrupleListAdress::afterRecvParticles, this, _1, _2));

  // We do not need those signals.
  sigBeforeSend.disconnect();
  sigAfterRecv.disconnect();
}

FixedQuadrupleListAdress::~FixedQuadrupleListAdress() {
  LOG4ESPP_INFO(theLogger, "~FixedQuadrupleListAdress");
  sigBeforeSendAT.disconnect();
  sigAfterRecvAT.disconnect();
}

bool FixedQuadrupleListAdress::add(longint pid1, longint pid2, longint pid3, longint pid4) {
  // here we assume pid1 < pid2 < pid3 < pid4
  bool returnVal = true;
  System& system = storage->getSystemRef();
  esutil::Error err(system.comm);

  // ADD THE LOCAL QUADRUPLET
  Particle *p1 = storage->lookupAdrATParticle(pid1);
  Particle *p2 = storage->lookupAdrATParticle(pid2);
  Particle *p3 = storage->lookupAdrATParticle(pid3);
  Particle *p4 = storage->lookupAdrATParticle(pid4);
  if (!p1){
    // Particle does not exist here, return false
    returnVal = false;
  } else {
    if (!p2) {
      std::stringstream msg;
      msg << "Quadruple particle p2 " << pid2 << " does not exists here and cannot be added.";
      err.setException( msg.str() );
    }
    if (!p3) {
      std::stringstream msg;
      msg << "Quadruple particle p3 " << pid3 << " does not exists here and cannot be added.";
      err.setException( msg.str() );
    }
    if (!p4) {
      std::stringstream msg;
      msg << "Quadruple particle p4 " << pid4 << " does not exists here and cannot be added.";
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
    } else {
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
  LOG4ESPP_INFO(theLogger, "Added fixed quadruple to local quadruple list.");
  return returnVal;
}

void FixedQuadrupleListAdress::beforeSendATParticles(std::vector<longint>& atpl, OutBuffer& buf) {
  std::vector< longint > toSend;

  // loop over the VP particle list
  for (std::vector<longint>::iterator it = atpl.begin(); it != atpl.end(); ++it) {
    longint pid = *it;

    LOG4ESPP_DEBUG(theLogger, "Send particle with pid " << pid << ", find pairs.");

    // find all triples that involve this particle
    int n = globalQuadruples.count(pid);

    if (n > 0) {
      std::pair<GlobalQuadruples::const_iterator,
        GlobalQuadruples::const_iterator> equalRange = globalQuadruples.equal_range(pid);

      // Data structure:
      // pid1 - N - n0 - n1 - ... - nN - pid2 - N ....
      toSend.reserve(toSend.size()+3*n+1);
      toSend.push_back(pid);
      toSend.push_back(n);
      for (GlobalQuadruples::const_iterator it = equalRange.first;
           it != equalRange.second; ++it) {
        toSend.push_back(it->second.first);
        toSend.push_back(it->second.second);
        toSend.push_back(it->second.third);
      }

      // Delete all of these triples from the global list.
      globalQuadruples.erase(equalRange.first, equalRange.second);
    }
  }

  // send the list
  buf.write(toSend);
  LOG4ESPP_INFO(theLogger, "Prepared fixed triple list before send particles.");
}

// Override parent function, this one should be empty
void FixedQuadrupleListAdress::beforeSendParticles(ParticleList& pl, OutBuffer& buf) {

}

void FixedQuadrupleListAdress::onParticlesChanged() {
  // (re-)generate the local quadruple list from the global list
  LOG4ESPP_INFO(theLogger, "Rebuild local bond list from global\n");

  System& system = storage->getSystemRef();
  esutil::Error err(system.comm);

  this->clear();
  longint lastpid1 = -1;
  Particle *p1;
  Particle *p2;
  Particle *p3;
  Particle *p4;
  for (GlobalQuadruples::const_iterator it = globalQuadruples.begin();
       it != globalQuadruples.end(); ++it) {
    if (it->first != lastpid1) {
      p1 = storage->lookupAdrATParticle(it->first);
      if (p1 == NULL) {
        std::stringstream msg;
        msg << "Quadruple particle p1 " << it->first << " does not exists here.";
        err.setException( msg.str() );
      }
      lastpid1 = it->first;
    }
    p2 = storage->lookupAdrATParticle(it->second.first);
    if (p2 == NULL) {
      std::stringstream msg;
      msg << "Quadruple particle p2 " << it->second.first << " does not exists here.";
      err.setException( msg.str() );
    }
    p3 = storage->lookupAdrATParticle(it->second.second);
    if (p3 == NULL) {
      std::stringstream msg;
      msg << "Quadruple particle p3 " << it->second.second << " does not exists here.";
      err.setException( msg.str() );
    }
    p4 = storage->lookupAdrATParticle(it->second.third);
    if (p4 == NULL) {
      std::stringstream msg;
      msg << "Quadruple particle p4 " << it->second.third << " does not exists here.";
      err.setException( msg.str() );
    }
    this->add(p1, p2, p3, p4);
  }
  err.checkException();
  LOG4ESPP_INFO(theLogger, "Regenerated local fixed quadruple list from global list.");
}

/****************************************************
** REGISTRATION WITH PYTHON
****************************************************/
void FixedQuadrupleListAdress::registerPython() {
  using namespace espressopp::python;

  bool (FixedQuadrupleListAdress::*pyAdd)(longint pid1, longint pid2,
         longint pid3, longint pid4) = &FixedQuadrupleListAdress::add;

  class_< FixedQuadrupleListAdress, shared_ptr< FixedQuadrupleListAdress > >
    ("FixedQuadrupleListAdress",
        init<shared_ptr< storage::Storage>, shared_ptr<FixedTupleListAdress> >())
    .def("add", pyAdd)
    .def("size", &FixedQuadrupleListAdress::size)
    .def("getQuadruples",  &FixedQuadrupleListAdress::getQuadruples)
   ;
}
}  // end namespace espressopp
