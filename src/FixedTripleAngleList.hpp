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
#ifndef _FIXEDTRIPLEANGLELIST_HPP
#define _FIXEDTRIPLEANGLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
//#include <boost/unordered_map.hpp>
#include <map>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

/*
 * This list is temporary solution. In general it should be derived from FixedTripleList.
 * in order not to generate similar code.
 * 
 * It will store the initial angle for each triple
 */

using namespace std;

namespace espressopp {
  class FixedTripleAngleList: public TripleList{
      protected:
		boost::signals2::connection con1, con2, con3;
		typedef multimap <longint,pair<pair<longint, longint>, real> > TriplesAngles;
		shared_ptr <storage::Storage> storage;
		TriplesAngles triplesAngles;
		using TripleList::add;

	  public:
		FixedTripleAngleList( shared_ptr<storage::Storage> _storage );
		virtual ~FixedTripleAngleList();
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

		// get only triples list
        python::list getTriples();
		// get triples and corresponding angles
		python::list getTriplesAngles();
        
        // get angle value for current triple
        real getAngle(int, int, int);

	    /** Get the number of triples in the GlobalTriples list */
	    int size() { return triplesAngles.size(); }

		static void registerPython();
	

	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);
  };
}

#endif
