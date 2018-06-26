/*
  Copyright (C) 2012-2018
      Max Planck Institute for Polymer Research
  Copyright (C) 2008-2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// ESPP_CLASS
#ifndef _SETTLE_HPP
#define _SETTLE_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Extension.hpp"
//#include "iterator/CellListIterator.hpp"
#include <boost/unordered_map.hpp>
#include <boost/signals2.hpp>
//#include "integrator/VelocityVerlet.hpp"
//#include "Triple.hpp"
#include "FixedTupleList.hpp"
#include "FixedTupleListAdress.hpp"
#include "boost/signals2.hpp"

namespace espressopp {
  namespace integrator {
    class Settle : public Extension {

        public:
            Settle(shared_ptr<System> _system, shared_ptr<FixedTupleListAdress> _fixedTupleList,
            		real mO, real mH, real distHH, real distOH);
            ~Settle();

            void add(longint pid) { molIDs.insert(pid); } // add molecule id (called from python)
            void saveOldPos();
            void applyConstraints();
            void correctVelocities();
            void settlep(longint molID);
            void settlev(longint molID);

            static void registerPython();

        private:
            boost::signals2::connection _befIntP, _aftIntP, _aftIntV, _aftIntSlow;
            std::set<longint> molIDs; // IDs of water molecules

            real mO, mH, distHH, distOH;
    	    real mOrmT, mHrmT;
    	    real rc, ra, rb;
    	    real rra; // inverse ra
            real mOmH, mOmH2;
            real twicemO,twicemH,mH2;

    	    // positions in previous timestep
    	    // key is molecule ID (first member of fixedtuplelist), value is OHH coordinates
    	    typedef boost::unordered_map<longint, Triple<Real3D, Real3D, Real3D> > OldPos;
    	    OldPos oldPos;

	    shared_ptr<FixedTupleListAdress> fixedTupleList;
	    void connect();
	    void disconnect();
            static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif

