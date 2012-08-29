// ESPP_CLASS
#ifndef _FIXEDTRIPLELIST_HPP
#define _FIXEDTRIPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
//#include <boost/unordered_map.hpp>
#include <map>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

namespace espresso {
  class FixedTripleList : public TripleList {
      protected:
		boost::signals2::connection con1, con2, con3;
		shared_ptr<storage::Storage> storage;
		typedef std::multimap <longint,std::pair <longint, longint> > GlobalTriples;
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

		python::list getTriples();

	    /** Get the number of triples in the GlobalTriples list */
	    int size() {
	    	return globalTriples.size();
	    }

		static void registerPython();
	

	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);


  };
}

#endif
