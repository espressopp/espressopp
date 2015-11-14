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
#include "FixedSingleList.hpp"
#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "esutil/Error.hpp"

using namespace std;

namespace espressopp {

  LOG4ESPP_LOGGER(FixedSingleList::theLogger, "FixedSingleList");


  FixedSingleList::FixedSingleList(shared_ptr< storage::Storage > _storage)
    : storage(_storage), globalSingles()
  {
    LOG4ESPP_INFO(theLogger, "construct FixedSingleList");

    con1 = storage->beforeSendParticles.connect
      (boost::bind(&FixedSingleList::beforeSendParticles, this, _1, _2));
    con2 = storage->afterRecvParticles.connect
      (boost::bind(&FixedSingleList::afterRecvParticles, this, _1, _2));
    con3 = storage->onParticlesChanged.connect
      (boost::bind(&FixedSingleList::onParticlesChanged, this));
  }

  FixedSingleList::~FixedSingleList() {

    LOG4ESPP_INFO(theLogger, "~FixedSingleList");

    con1.disconnect();
    con2.disconnect();
    con3.disconnect();
  }

  bool FixedSingleList::
  add(longint pid1) {
    bool returnVal = true;
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    
    // ADD THE LOCAL SINGLE
    Particle *p1 = storage->lookupRealParticle(pid1);
    
    if (p1){
      // add the particle locally
      this->add(p1);
      // add the particle to the pid list if it doesn't exists yet
      GlobalSingles::const_iterator it = globalSingles.find(pid1);
      if (it == globalSingles.end()) {
        globalSingles.insert(pid1);
      } else {
  	      // TODO: Single already exists, generate error!
  	      ;
      }
      LOG4ESPP_INFO(theLogger, "added fixed pair to global pair list");
    } else {
    	returnVal = false;
    }
    
    return returnVal;
  }

  python::list FixedSingleList::getSingles()
  {
	python::list sl;
	longint pid;
	for (GlobalSingles::iterator it=globalSingles.begin(); it != globalSingles.end(); it++) {
	  pid = *it;
      sl.append(pid);
    }

	return sl;
  }

  void FixedSingleList::beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
	std::vector< longint > toSend(pl.size());
    // loop over the particle list
    for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
      longint pid = pit->id();
      toSend.push_back(pid);
      globalSingles.erase(pid);
      LOG4ESPP_DEBUG(theLogger, "erase and send particle with pid from FixedSingleList" << pid);
    }
    // send the list
    buf.write(toSend);
    LOG4ESPP_INFO(theLogger, "prepared fixed single list before send particles");
  }

  void FixedSingleList::afterRecvParticles(ParticleList &pl, InBuffer& buf) {
	longint size = pl.size();
    std::vector< longint > received(size);
    buf.read(received);
    for (std::vector< longint >::iterator it=received.begin(); it!=received.end(); it++) {
    	globalSingles.insert(*it);
    }
    LOG4ESPP_INFO(theLogger, "received fixed single list after receive particles");
  }

  void FixedSingleList::
  onParticlesChanged() {
    LOG4ESPP_INFO(theLogger, "rebuild local fixed single list from global\n");
    System& system = storage->getSystemRef();
    esutil::Error err(system.comm);
    Particle *p1;
    
    this->clear();
    for (GlobalSingles::const_iterator it=globalSingles.begin(); it != globalSingles.end(); ++it) {
    	p1 = storage->lookupRealParticle(*it);
        if (p1 == NULL) {
          std::stringstream msg;
          msg << "onParticlesChanged error. Fixed Single List particle p1 " << *it << " does not exists here";
          err.setException( msg.str() );
        }
        this->add(p1);
    }
    err.checkException();
    LOG4ESPP_INFO(theLogger, "regenerated local fixed single list from global list");
  }

  /****************************************************
  ** REGISTRATION WITH PYTHON
  ****************************************************/

  void FixedSingleList::registerPython() {

    using namespace espressopp::python;

    bool (FixedSingleList::*pyAdd)(longint pid1) = &FixedSingleList::add;

    class_<FixedSingleList, shared_ptr<FixedSingleList> >
      ("FixedSingleList", init <shared_ptr<storage::Storage> >())
      .def("add", pyAdd)
      .def("size", &FixedSingleList::size)
      .def("getSingles",  &FixedSingleList::getSingles)
      ;
  }
}
