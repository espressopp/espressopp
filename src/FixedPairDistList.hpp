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
#ifndef _FIXEDPAIRDISTLIST_HPP
#define _FIXEDPAIRDISTLIST_HPP

#include "log4espp.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
//#include <boost/unordered_map.hpp>
#include <map>
#include <boost/signals2.hpp>
#include "types.hpp"

//#include "FixedListComm.hpp"

namespace espressopp {
	class FixedPairDistList : public PairList{
	  protected:
	    typedef std::multimap<longint, std::pair<longint, real> > PairsDist;
		boost::signals2::connection con1, con2, con3;
		shared_ptr <storage::Storage> storage;
		PairsDist pairsDist;
		using PairList::add;

	  public:
		FixedPairDistList(shared_ptr< storage::Storage > _storage);
		virtual ~FixedPairDistList();

		/** Add the given particle pair to the list on this processor if the
		particle with the lower id belongs to this processor.  Note that
		this routine does not check whether the pair is inserted on
		another processor as well.
		\return whether the particle was inserted on this processor.
		*/
		virtual bool add(longint pid1, longint pid2);
		virtual void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
		void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
		virtual void onParticlesChanged();

	    python::list getPairs();
	    //PairsDist* getGlobalPairs() {return &globalPairs;};
	    python::list getPairsDist();

        real getDist(int, int);
	    /** Get the number of bonds in the GlobalPairs list */
	    int size() {return pairsDist.size();}

	    static void registerPython();

	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);
	};
}

#endif

