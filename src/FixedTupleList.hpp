// ESPP_CLASS
#ifndef _FIXEDTUPLELIST_HPP
#define _FIXEDTUPLELIST_HPP

//#include "log4espp.hpp"
//#include "types.hpp"

//#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
//#include <boost/unordered_map.hpp>
//#include <boost/signals2.hpp>
#include "FixedListComm.hpp"

namespace espresso {
    class FixedTupleList: public FixedListComm  {
        protected:
            //boost::signals2::connection con1, con2, con3;
            //shared_ptr< storage::Storage > storage;
            //typedef boost::unordered_multimap< longint, longint > GlobalPairs;
            //GlobalPairs globalPairs;
            //using PairList::add;

        public:
            FixedTupleList(shared_ptr< storage::Storage > _storage);
            //~FixedTupleList();

            //bool add(longint pid1, longint pid2);
            //void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
            //void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
            //void onParticlesChanged();

            static void registerPython();

        private:
            //static LOG4ESPP_DECL_LOGGER(theLogger);
    };
}

#endif

