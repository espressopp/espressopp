// ESPP_CLASS
#ifndef _FIXEDTUPLELIST_HPP
#define _FIXEDTUPLELIST_HPP

#include "log4espp.hpp"
//#include "types.hpp"

#include "Particle.hpp"
#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "FixedListComm.hpp"

namespace espresso {
    class FixedTupleList: public TupleList  {
        protected:
            boost::signals2::connection con1, con2, con3;
            shared_ptr<storage::Storage> storage;
            typedef std::vector<longint> tuple;
            typedef boost::unordered_multimap<longint, tuple> GlobalTuples;
            GlobalTuples globalTuples;
            using TupleList::add;

        public:
            FixedTupleList(shared_ptr<storage::Storage> _storage);
            ~FixedTupleList();

            void add(longint pid) { tmppids.push_back(pid); } // add particle id (called from python)
            void addTs() { addT(tmppids); tmppids.clear(); } // add tuple
            void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
            void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
            void onParticlesChanged();

            static void registerPython();

        private:
            tuple tmppids;
            bool addT(tuple pids); // add tuple
            static LOG4ESPP_DECL_LOGGER(theLogger);
    };
}

#endif

