// ESPP_CLASS
#ifndef _SETTLE_HPP
#define _SETTLE_HPP

#include "log4espp.hpp"

#include "iterator/CellListIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "integrator/VelocityVerletAdress.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "Triple.hpp"
#include "FixedTupleList.hpp"

namespace espresso {
    class Settle  {

        public:
            Settle(shared_ptr<storage::Storage> _storage,
            		shared_ptr<integrator::VelocityVerlet> _integrator,
            		real mO, real mH,
            		real distHH, real distOH);
            ~Settle();

            void add(longint pid) { molIDs.insert(pid); } // add molecule id (called from python)
            void saveOldPos();
            void applyConstraints();
            void settlep(longint molID);

            static void registerPython();

        private:
            std::set<longint> molIDs; // IDs of water molecules

            real mO, mH, distHH, distOH;
    		real mOrmT, mHrmT;
    		real rc, ra, rb;
    		real rra; // inverse ra

    		// positions in previous timestep
    		// key is molecule ID (VP particle ID), value is OHH coordinates
    		typedef boost::unordered_map<longint, Triple<Real3D, Real3D, Real3D> > OldPos;
    		OldPos oldPos;

            boost::signals2::connection con1, con2;
			shared_ptr<storage::Storage> storage;
			shared_ptr<integrator::VelocityVerlet> integrator; // this is needed for signal connection
			shared_ptr<FixedTupleList> fixedtupleList;

            static LOG4ESPP_DECL_LOGGER(theLogger);
    };
}

#endif

