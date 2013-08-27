// ESPP_CLASS
#ifndef _FIXEDTUPLELIST_HPP
#define _FIXEDTUPLELIST_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
#include "Real3D.hpp"

namespace espresso {
  class FixedTupleList : public TupleList {
      protected:
		boost::signals2::connection con1, con2, con3;
		shared_ptr<storage::Storage> storage;
        typedef std::vector<longint> tuple;
		typedef std::multimap <longint,tuple > GlobalTuples;
		GlobalTuples globalTuples;
		using TupleList::add;

      //FixedListComm<FixedTupleList, 3> _comm;

	  public:
		FixedTupleList(shared_ptr<storage::Storage> _storage);
		virtual ~FixedTupleList();
		virtual bool addTuple(boost::python::list& tuple);
		virtual void beforeSendParticles(ParticleList& pl, class OutBuffer &buf);
		void afterRecvParticles(ParticleList& pl, class InBuffer &buf);
		virtual void onParticlesChanged();

		python::list getTuples();

		Real3D calcTupleCOM(int tupleid);

	    /** Get the number of triples in the GlobalTriples list */
	    int size() {
	    	return globalTuples.size();
	    }

		static void registerPython();
	

	  private:
		static LOG4ESPP_DECL_LOGGER(theLogger);


  };
}

#endif
