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

// ESPP_CLASS
#ifndef _FIXEDLOCALTUPLELIST_HPP
#define _FIXEDLOCALTUPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espressopp {
    class FixedLocalTupleList : public TupleList {
    protected:
	boost::signals2::connection con1, con2, con3;
	shared_ptr<storage::Storage> storage;
	typedef std::vector<longint> tuple;
	typedef std::multimap <longint,tuple > GlobalTuples;
	GlobalTuples globalTuples;
	using TupleList::add;
	
    public:
	FixedLocalTupleList(shared_ptr<storage::Storage> _storage);
	virtual ~FixedLocalTupleList();
	virtual bool addTuple(boost::python::list& tuple);
	virtual void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
	void afterRecvParticles(ParticleList& pl, class InBuffer &buf);
	virtual void onParticlesChanged();
	
	python::list getTuples();

	/** Get the number of triples in the GlobalTriples list */
	int size() {
	    return globalTuples.size();
	}
	
	static void registerPython();
	
	
    private:
	static LOG4ESPP_DECL_LOGGER(theLogger);
    };
}

#endif
