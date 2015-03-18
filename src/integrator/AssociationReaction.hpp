/*
  Copyright (C) 2014
      Pierre de Buyl
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
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
#ifndef _INTEGRATOR_ASSOCIATIONREACTION_HPP
#define _INTEGRATOR_ASSOCIATIONREACTION_HPP

#include "types.hpp"
#include "logging.hpp"
#include "VerletList.hpp"
#include "FixedPairList.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"
#include "interaction/Potential.hpp"

#include "boost/signals2.hpp"


namespace espressopp {
  namespace integrator {

    const int AR_COMM_TAG = 0xac;
    using namespace storage;

    /** Reaction scheme for polymer growth and curing/crosslinking

	This extension enables the rate-controlled stochastic curing of polymer
	systems, either for chain growth of step growth, depending on the
	parameters.

	The variables typeA, typeB, stateAMin control the particles that enter
	the curing reaction
	\f[ A^a + B^0 -> A^{a+deltaA}-B^{deltaB} \f]
	where A and B may possess additional bonds not shown. An extra bond is
	added between A and B upon reaction.

	The reaction proceeds by testing for all possible (A,B) pairs and
	selects them only at a given rate. It works in parallel, by gathering
	first the successful pairs between neigboring CPUs and ensuring that
	each particle enters only in one new bond per reaction step.

    */

    class AssociationReaction : public Extension {

      public:

      AssociationReaction(shared_ptr<System> system, shared_ptr<VerletList> _verletList, shared_ptr<FixedPairList> _fixedPairList, shared_ptr<DomainDecomposition> _domdec);
      ~AssociationReaction();

      void setRate(real rate);
      real getRate();
      void setCutoff(real cutoff);
      real getCutoff();
      void setTypeA(size_t typeA);
      size_t getTypeA();
      void setTypeB(size_t typeB);
      size_t getTypeB();
      void setDeltaA(int deltaA);
      int getDeltaA();
      void setDeltaB(int deltaB);
      int getDeltaB();
      void setStateAMin(int stateAMin);
      int getStateAMin();
      void setInterval(int interval);
      int getInterval();

      void initialize();

      /** Actual reaction step */
      void react();

      void sendMultiMap(boost::unordered_multimap<longint, longint> &mm);
      void uniqueA(boost::unordered_multimap<longint, longint> &mm);
      void uniqueB(boost::unordered_multimap<longint, longint> &mm, boost::unordered_multimap<longint, longint> &nn);
      void applyAR();

      /** Register this class so it can be used from Python. */
      static void registerPython();

      private:

      boost::signals2::connection _initialize, _react;

      void reactPair(Particle& p1, Particle& p2);

      void connect();
      void disconnect();

      real rate; //!< reaction rate
      real cutoff; //!< reaction cutoff
      real cutoff_sqr; //!< reaction cutoff squared
      size_t typeA; //!< type of reactant A
      size_t typeB; //!< type of reactant B
      int deltaA; //!< state change for reactant A
      int deltaB; //!< state change for reactant B
      int stateAMin; //!< minimum state of reactant A
      int interval; //!< number of steps between reaction loops
      real dt; //!< timestep from the integrator
      shared_ptr<espressopp::interaction::Potential> potential;

      real current_cutoff;
      real current_cutoff_sqr;
      shared_ptr<VerletList> verletList;
      shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term
      shared_ptr<FixedPairList> fpl;
      shared_ptr<DomainDecomposition> domdec;

      /** container for (A,B) potential partners */
      boost::unordered_multimap<longint, longint> Alist;
      /** container for (A,B) effective partners */
      boost::unordered_multimap<longint, longint> Blist;
      static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
