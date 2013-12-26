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


namespace espresso {
  namespace integrator {

    const int AR_COMM_TAG = 0xac;
    using namespace storage;

    /** Reaction scheme for polymer growth and curing/crosslinking */

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
      shared_ptr<espresso::interaction::Potential> potential;

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
      static LOG4ESPP_DECL_LOGGER(thelogger);
    };
  }
}

#endif
