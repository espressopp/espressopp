// ESPP_CLASS
#ifndef _FIXEDTRIPLELISTADRESS_HPP
#define _FIXEDTRIPLELISTADRESS_HPP

#include "log4espp.hpp"
#include "types.hpp"

//#include "Particle.hpp"
#include "FixedTripleList.hpp"
#include "FixedTupleListAdress.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>


namespace espresso {

    /**
     * This is a subclass of FixedTripleList. It should be used for AdResS fixed
     * triples. It overrides some parent functions, to use AT particles.
     *
     */
	class FixedTripleListAdress : public FixedTripleList {
	  public:
		FixedTripleListAdress(shared_ptr<storage::Storage> _storage,
		        shared_ptr<FixedTupleListAdress> _fixedtupleList);
		~FixedTripleListAdress();

		/** Add the given particle triple to the list on this processor if the
		particle with the lower id belongs to this processor.  Note that
		this routine does not check whether the pair is inserted on
		another processor as well.
		\return whether the particle was inserted on this processor.
		*/
		bool add(longint pid1, longint pid2, longint pid3);
		void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
		void beforeSendATParticles(std::vector<longint>& atpl, class OutBuffer& buf);
		void onParticlesChanged();

		static void registerPython();

	  protected:
		// fixedtuple list connects to this and triggers beforeSendATParticles()
		boost::signals2::connection con;

	  private:
		shared_ptr<FixedTupleListAdress> fixedtupleList;
		using TripleList::add;
		static LOG4ESPP_DECL_LOGGER(theLogger);
	};
}

#endif

