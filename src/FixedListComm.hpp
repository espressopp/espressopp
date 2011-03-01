// ESPP_CLASS
#ifndef _FIXEDLISTCOMM_HPP
#define _FIXEDLISTCOMM_HPP

#include "log4espp.hpp"
#include "types.hpp"

#include "Particle.hpp"
//#include "esutil/ESPPIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>

namespace espresso {
    class FixedListComm: public PairList, public TripleList, public QuadrupleList {
        protected:
        boost::signals2::connection con1, con2, con3;
        shared_ptr<storage::Storage> storage;
        typedef std::vector<longint> pvec;
        typedef boost::unordered_multimap<longint, pvec> GlobalList;
        GlobalList globalLists;
        //using PairList::add;

        // if one wants to get rid of the virtual,
        // move this to a template
        //virtual void add(std::vector<Particle*> tmp)=0;

        public:
        FixedListComm(shared_ptr<storage::Storage> _storage);
        ~FixedListComm();
        bool add(pvec pids);
        void beforeSendParticles(ParticleList& pl, class OutBuffer& buf);
        void afterRecvParticles(ParticleList& pl, class InBuffer& buf);
        void onParticlesChanged();

        private:
        static LOG4ESPP_DECL_LOGGER(theLogger);
    }; // class
} //ns espresso

#endif
