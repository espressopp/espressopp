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

// ESPP_CLASS
#ifndef _FIXEDTUPLELIST_HPP
#define _FIXEDTUPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include "Real3D.hpp"


namespace espressopp {
	/* A list of particles, where the first particle is a dependent (virtual) particle.
	 *  All particles reside on the same node of the first particle.
	 */
  class FixedTupleList : public TupleList {
      protected:
		boost::signals2::connection con1, con2, con3;
		shared_ptr<storage::Storage> storage;
        typedef std::vector<longint> tuple;
		typedef std::multimap <longint,tuple > GlobalTuples;
		GlobalTuples globalTuples;
		using TupleList::add;

	  public:
		FixedTupleList(shared_ptr<storage::Storage> _storage);
		virtual ~FixedTupleList();
		virtual bool addTuple(boost::python::list& tuple);
		virtual void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
		void afterRecvParticles(ParticleList& pl, class InBuffer &buf);
		virtual void onParticlesChanged();
		void unwrapMinimumImage(int id);

		python::list getTuples();

		Real3D calcTupleCOM(int tupleid);

		std::vector<Particle *> getTupleByID(int id);


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
