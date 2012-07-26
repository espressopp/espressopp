// ESPP_CLASS
#ifndef _FIXEDPAIRLIST_HPP
#define _FIXEDPAIRLIST_HPP

#include "log4espp.hpp"
#include "types.hpp"
#include "python.hpp"
#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

//#include "FixedListComm.hpp"

namespace espresso {
	class FixedPairList : public PairList {
	  public:
	    typedef boost::unordered_multimap <longint, longint> GlobalPairs;

	  protected:
		boost::signals2::connection con1, con2, con3;
		shared_ptr <storage::Storage> storage;
		GlobalPairs globalPairs;
		using PairList::add;

	  public:
		FixedPairList(shared_ptr <storage::Storage> _storage);
		virtual ~FixedPairList();

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

	    python::list getBonds();
	    GlobalPairs* getGlobalPairs() {return &globalPairs;};

	    /** Get the number of bonds in the GlobalPairs list */
	    int size() {
	    	return globalPairs.size();
	    }

	    static void registerPython();

	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);
	};
}

#endif

