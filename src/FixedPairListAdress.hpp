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

// ESPP_CLASS
#ifndef _FIXEDPAIRLISTADRESS_HPP
#define _FIXEDPAIRLISTADRESS_HPP

#include "log4espp.hpp"

//#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "FixedTupleListAdress.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include "types.hpp"

namespace espressopp {

    /**
     * This is a subclass of FixedPairList. It should be used for AdResS fixed
     * pairs. It overrides some parent functions, to use AT particles.
     *
     */
	class FixedPairListAdress : public FixedPairList {
	  public:
		FixedPairListAdress(shared_ptr<storage::Storage> _storage,
		        shared_ptr<FixedTupleListAdress> _fixedtupleList);
		virtual ~FixedPairListAdress();

		/** Add the given particle pair to the list on this processor if the
		particle with the lower id belongs to this processor.  Note that
		this routine does not check whether the pair is inserted on
		another processor as well.
		\return whether the particle was inserted on this processor.
		*/
		bool add(longint pid1, longint pid2);
		void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
		void beforeSendATParticles(std::vector<longint>& atpl, class OutBuffer& buf);
		void onParticlesChanged();
		void remove();
		python::list getBonds();
		static void registerPython();

	  protected:
		// fixedtuple list connects to this and triggers beforeSendATParticles()
		boost::signals2::connection sigBeforeSendAT, sigAfterRecvAT;

	  private:
		shared_ptr<FixedTupleListAdress> fixedtupleList;
		using PairList::add;
		static LOG4ESPP_DECL_LOGGER(theLogger);
	};
}

#endif

