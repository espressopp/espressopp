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

// ESPP_CLASS
#ifndef _FIXEDTRIPLELIST_HPP
#define _FIXEDTRIPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

namespace espressopp {
  class FixedTripleList : public TripleList {
      protected:
		boost::signals2::connection sigAfterRecv, sigOnParticleChanged, sigBeforeSend;
		shared_ptr<storage::Storage> storage;
		typedef boost::unordered_multimap <longint,std::pair <longint, longint> > GlobalTriples;
		GlobalTriples globalTriples;
		using TripleList::add;

      //FixedListComm<FixedTripleList, 3> _comm;

	  public:
		FixedTripleList(shared_ptr<storage::Storage> _storage);
		virtual ~FixedTripleList();
		//bool add(pvec pids) { _comm.add(pids); }
		/** Add the given particle triple to the list on this processor if the
		particle with the lower id belongs to this processor.  Note that
		this routine does not check whether the triple is inserted on
		another processor as well.
		\return whether the triple was inserted on this processor.
		*/
		virtual bool add(longint pid1, longint pid2, longint pid3);
		virtual void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
		void afterRecvParticles(ParticleList& pl, class InBuffer &buf);
		virtual void onParticlesChanged();
		virtual std::vector<longint> getTripleList();
		python::list getTriples();

	    /** Get the number of triples in the GlobalTriples list */
	    int size() {
	    	return globalTriples.size();
	    }
		
		void remove();
		static void registerPython();

		
	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);


  };
}

#endif
