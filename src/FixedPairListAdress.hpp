// ESPP_CLASS
#ifndef _FIXEDPAIRLISTADRESS_HPP
#define _FIXEDPAIRLISTADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"

//#include "Particle.hpp"
#include "FixedPairList.hpp"
#include "FixedTupleList.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>


namespace espresso {
	class FixedPairListAdress : public FixedPairList {
	  public:
		FixedPairListAdress(shared_ptr<storage::Storage> _storage,
		        shared_ptr<FixedTupleList> _fixedtupleList);
		~FixedPairListAdress();

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

		static void registerPython();

	  protected:
		// this connects to fixedtuple list and triggers beforeSendATParticles()
		boost::signals2::connection con4;

	  private:
		shared_ptr<FixedTupleList> fixedtupleList;
		using PairList::add;
		static LOG4ESPP_DECL_LOGGER(theLogger);
	};
}

#endif

