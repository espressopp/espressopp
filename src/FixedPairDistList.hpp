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

namespace espresso {
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

