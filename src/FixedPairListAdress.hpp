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

namespace espresso {

    /**
     * This is a subclass of FixedPairList. It should be used for AdResS fixed
     * pairs. It overrides some parent functions, to use AT particles.
     *
     */
	class FixedPairListAdress : public FixedPairList{
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
		python::list getBonds();

		static void registerPython();

	  protected:
		// fixedtuple list connects to this and triggers beforeSendATParticles()
		boost::signals2::connection con;

	  private:
		shared_ptr<FixedTupleListAdress> fixedtupleList;
		using PairList::add;
		static LOG4ESPP_DECL_LOGGER(theLogger);
	};
}

#endif

