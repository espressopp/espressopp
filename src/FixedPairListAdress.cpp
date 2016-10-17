/*
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

#include "storage/Storage.hpp"
#include <boost/bind.hpp>
#include "FixedPairListAdress.hpp"
#include "Buffer.hpp"
#include "esutil/Error.hpp"

//using namespace std;

namespace espressopp {

  LOG4ESPP_LOGGER(FixedPairListAdress::theLogger, "FixedPairListAdress");


  FixedPairListAdress::
  FixedPairListAdress(shared_ptr< storage::Storage > _storage,
          shared_ptr<FixedTupleListAdress> _fixedtupleList)
          : FixedPairList(_storage), fixedtupleList(_fixedtupleList) {
    LOG4ESPP_INFO(theLogger, "construct FixedPairListAdress");
    
    sigBeforeSendAT = fixedtupleList->beforeSendATParticles.connect
          (boost::bind(&FixedPairListAdress::beforeSendATParticles, this, _1, _2));

    sigAfterRecvAT = fixedtupleList->afterRecvATParticles.connect
      (boost::bind(&FixedPairListAdress::afterRecvParticles, this, _1, _2));

    // We do not need those signals here.
    sigBeforeSend.disconnect();
    sigAfterRecv.disconnect();
  }

  FixedPairListAdress::~FixedPairListAdress() {

    LOG4ESPP_INFO(theLogger, "~FixedPairListAdress");

    sigBeforeSendAT.disconnect();
    sigAfterRecvAT.disconnect();
  }


  // override parent function (use lookupAdrATParticle)
  bool FixedPairListAdress::add(longint pid1, longint pid2) {

    if (pid1 > pid2)
      std::swap(pid1, pid2);

    bool returnVal = true;
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    // ADD THE LOCAL PAIR
    Particle *p1 = storage->lookupAdrATParticle(pid1);
    Particle *p2 = storage->lookupAdrATParticle(pid2);
    if (!p1){
      // Particle does not exist here, return false
      returnVal = false;
    }
    else{
      if (!p2) {
        std::stringstream msg;
        msg << "Atomistic bond particle p2 (id=" << pid2 << ") does not exists here "
            << "and cannot be added. "
            << " The pair " << pid1 << " - " << pid2 << " could not be created.";
        err.setException( msg.str() );
      }
    }
    err.checkException();
    
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
    return returnVal;
  }


  void FixedPairListAdress::beforeSendATParticles(std::vector<longint>& atpl,
          OutBuffer& buf) {

        //std::cout << "beforeSendATParticles() fixed pl (size " << atpl.size() << ")\n";

        std::vector< longint > toSend;

        // loop over the VP particle list
        for (std::vector<longint>::iterator it = atpl.begin();
                it != atpl.end(); ++it) {
          longint pid = *it;

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
            globalPairs.erase(pid);
          }
        }

        // send the list
        buf.write(toSend);
        LOG4ESPP_INFO(theLogger, "prepared fixed pair list before send particles");
  }


  // override parent function, this one should be empty
  void FixedPairListAdress::beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
  }


  // override parent function (use lookupAdrATParticle())
  void FixedPairListAdress::onParticlesChanged() {

    LOG4ESPP_INFO(theLogger, "rebuild local bond list from global\n");

    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    this->clear();
    longint lastpid1 = -1;
    Particle *p1;
    Particle *p2;

    for (GlobalPairs::const_iterator it = globalPairs.begin();
	 it != globalPairs.end(); ++it) {

        if (it->first != lastpid1) {
            p1 = storage->lookupAdrATParticle(it->first);
            if (p1 == NULL) {
                std::stringstream msg;
                msg << "FixedPairListAdress ";
                msg << "Atomistic bond particle p1 (id=" << it->first << ") does not exists here.";
                err.setException( msg.str() );
            }
            lastpid1 = it->first;
        }

        p2 = storage->lookupAdrATParticle(it->second);
        if (p2 == NULL) {
            std::stringstream msg;
            msg << "FixedPairListAdress ";
            msg << "Atomistic bond particle p2 (id=" << it->second << ") does not exists here.";
            err.setException( msg.str() );
        }

        //std::cout << " adding (" << p1->getId() << ", " << p2->getId() << ")\n";
        this->add(p1, p2);
    }
    err.checkException();
    
    LOG4ESPP_INFO(theLogger, "Regenerated local fixed pair list from global list");
  }

  python::list FixedPairListAdress::getBonds()
  {
	python::tuple bond;
	python::list bonds;
	for (GlobalPairs::const_iterator it=globalPairs.begin(); it != globalPairs.end(); it++) {
      bond = python::make_tuple(it->first, it->second);
      bonds.append(bond);
    }

	return bonds;
  }

  void FixedPairListAdress::remove(void) {
      this->clear();
      globalPairs.clear();
      sigBeforeSend.disconnect();
      sigAfterRecv.disconnect();
      sigBeforeSendAT.disconnect();
      sigAfterRecvAT.disconnect();
      sigOnParticlesChanged.disconnect();
  }


  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedPairListAdress::registerPython() {

    using namespace espressopp::python;

    bool (FixedPairListAdress::*pyAdd)(longint pid1, longint pid2)
      = &FixedPairListAdress::add;

    class_<FixedPairListAdress, shared_ptr<FixedPairListAdress> >
      ("FixedPairListAdress",
              init <shared_ptr<storage::Storage>,
                     shared_ptr<FixedTupleListAdress> >())
      .def("add", pyAdd)
      .def("remove",  &FixedPairListAdress::remove)
      .def("getBonds",  &FixedPairListAdress::getBonds)
      ;
  }

}
