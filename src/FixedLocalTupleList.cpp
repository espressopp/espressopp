/*
  Copyright (C) 2017
      Max Planck Institute for Polymer Research
  
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
#include "FixedLocalTupleList.hpp"

#include <boost/bind.hpp>
#include "storage/Storage.hpp"
#include "Buffer.hpp"

#include "esutil/Error.hpp"


namespace espressopp {
    
    LOG4ESPP_LOGGER(FixedLocalTupleList::theLogger, "FixedLocalTupleList");
    
    FixedLocalTupleList::FixedLocalTupleList(shared_ptr< storage::Storage > _storage)
	: storage(_storage), globalTuples()
    {
	LOG4ESPP_INFO(theLogger, "construct FixedLocalTupleList");
	
	con1 = storage->beforeSendParticles.connect
	    (boost::bind(&FixedLocalTupleList::beforeSendParticles, this, _1, _2));
	con2 = storage->afterRecvParticles.connect
	    (boost::bind(&FixedLocalTupleList::afterRecvParticles, this, _1, _2));
	con3 = storage->onParticlesChanged.connect
	    (boost::bind(&FixedLocalTupleList::onParticlesChanged, this));
    }

    FixedLocalTupleList::~FixedLocalTupleList() {
	
	LOG4ESPP_INFO(theLogger, "~FixedLocalTupleList");
	
	con1.disconnect();
	con2.disconnect();
	con3.disconnect();
    }
    
    bool FixedLocalTupleList::
    addTuple(boost::python::list& tuple) {
        bool returnVal = true;
        System& system = storage->getSystemRef();
        esutil::Error err(system.comm);
        
        Particle* vp, *at;
        longint pidK; // the pid used as key
        std::vector<Particle*> tmp; // temporary vector
        std::vector<longint> pids; //used to extract from tuple;
        std::vector<longint> pidstmp; // temporary vector
        
        for (int i = 0; i < boost::python::len(tuple); ++i) {
            pids.push_back(boost::python::extract<int>(tuple[i]));
        }
        
        tuple::iterator it = pids.begin();
        vp = storage->lookupRealParticle(*it);
        if (!vp) {// Particle does not exist here, return false
            returnVal = false;
        } else {
            pidK = *it; // first pid is key
            for (++it; it != pids.end(); ++it) {
                at = storage->lookupLocalParticle(*it);
                if (!at) { // Particle does not exist here, return false
                    std::stringstream msg;
                    msg << "ERROR: particle " << *it << " not found \n";
                    err.setException(msg.str());
                    returnVal = false;
                    break;
                }
                tmp.push_back(at);
                pidstmp.push_back(*it); // pidK is not in this vector
            }
        }
        err.checkException();
	
        if (returnVal) {
            this->add(vp, tmp); // add to TupleList
            globalTuples.insert(make_pair(pidK, pidstmp));
	    LOG4ESPP_INFO(theLogger, "added fixed tuple to global tuples");
        }
	
        tmp.clear();
        pids.clear();
        pidstmp.clear();
	

        return returnVal;
    }

    python::list FixedLocalTupleList::getTuples() {
        python::list tuple;
        python::list alltuples;
        for (GlobalTuples::iterator it = globalTuples.begin();
                it != globalTuples.end(); it++) {
            python::list tuple;
            tuple.append((*it).first); // key is also part of the tuple!
            for (tuple::iterator it2=(*it).second.begin(); it2 !=(*it).second.end(); it2++){
                tuple.append(*it2);
            }
            alltuples.append(tuple);
        }
	
        return alltuples;
    }
    
    void FixedLocalTupleList::
    beforeSendParticles(ParticleList& pl, OutBuffer& buf) {
        std::vector<longint> toSend;
        // loop over the particle list
        for (ParticleList::Iterator pit(pl); pit.isValid(); ++pit) {
            longint pidK = pit->id();
            LOG4ESPP_DEBUG(theLogger, "send particle with pid " << pidK << ", find tuples");
	    
            // find particle that involves this particle id
            GlobalTuples::const_iterator it = globalTuples.find(pidK);
            if (it != globalTuples.end()) {
		int s = it->second.size();
		toSend.reserve(toSend.size() + s + 2);
                // first write the pid of the first particle
		toSend.push_back(pidK);
                // write the size of the vector
		toSend.push_back(s);

                // iterate through vector and add pids
                for (tuple::const_iterator it2 = it->second.begin();
                        it2 != it->second.end(); ++it2) {
                    toSend.push_back(*it2);
                }
                // delete this pid from the global list
                globalTuples.erase(pidK);
		LOG4ESPP_DEBUG(theLogger, "Erase pid " << pidK);
            }
        }
	buf.write(toSend);
    }

    void FixedLocalTupleList::
    afterRecvParticles(ParticleList &pl, InBuffer& buf) {
	std::vector<longint> received;
        std::vector<longint> pids;
        int i, n;
        longint pidK;
        GlobalTuples::iterator it = globalTuples.begin();

	buf.read(received);

	int size = received.size();

	i = 0;
	while (i < size) {
	    pidK = received[i++];
	    if (pidK > 0) {
		n = received[i++];
		for (; n > 0; --n) {
		    LOG4ESPP_DEBUG(theLogger, "received vector for pid " << pidK);
		    int pidV = received[i++];
		    pids.push_back(pidV);
		}
		it = globalTuples.insert(it, std::make_pair(pidK, pids));
		pids.clear();
		LOG4ESPP_DEBUG(theLogger, "Insert pid " << pidK);
	    }
	}
    }


    void FixedLocalTupleList::onParticlesChanged() {
	// TODO errors should be thrown in a nice way
        LOG4ESPP_INFO(theLogger, "rebuild local particle list from global tuples\n");
        System& system = storage->getSystemRef();
        esutil::Error err(system.comm);
        
        this->clear();
        Particle* vp, * at;
        std::vector<Particle*> tmp;
	
        GlobalTuples::const_iterator it = globalTuples.begin();
	
        for (;it != globalTuples.end(); ++it) {
	    if (globalTuples.count(it->first) > 1) {
		std::stringstream msg;
                msg << " particle in tuple " << it->first << " is duplicated";
                err.setException( msg.str() );
	    }
            vp = storage->lookupRealParticle(it->first);
            if (vp == NULL) {
            	std::stringstream msg;
                msg << " particle in tuple " << it->first << " does not exists here";
                err.setException( msg.str() );
            }
	    
            // iterate through vector in map
            for (tuple::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                at = storage->lookupLocalParticle(*it2);
                if (at == NULL) {
		    std::stringstream msg;
		    msg << " particle in tuple " << *it2 << " does not exists here";
		    err.setException( msg.str() );
                }
                tmp.push_back(at);
            }
            this->add(vp, tmp);
            tmp.clear();
        }
        LOG4ESPP_INFO(theLogger, "regenerated local fixed list from global tuples");
    }
    
    /****************************************************
     ** REGISTRATION WITH PYTHON
     ****************************************************/
    
    void FixedLocalTupleList::registerPython() {
	
	using namespace espressopp::python;
	
	class_< FixedLocalTupleList, shared_ptr< FixedLocalTupleList > >
	    ("FixedLocalTupleList", init< shared_ptr< storage::Storage > >())
	    .def("addTuple", &FixedLocalTupleList::addTuple)
	    .def("getTuples", &FixedLocalTupleList::getTuples)
	    .def("size", &FixedLocalTupleList::size)
	    ;
    }
}
